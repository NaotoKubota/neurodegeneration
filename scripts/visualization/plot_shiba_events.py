#!/usr/bin/env python3
"""Plot curated, species-aware Shiba splicing events across completed contrasts."""

from __future__ import annotations

import argparse
import json
import re
from math import ceil
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
EVENT_TYPES = {"SE", "FIVE", "THREE", "RI", "MXE", "AFE", "ALE", "MSE"}
SPECIES_ALIASES = {"Mus musculus": "mouse", "Homo sapiens": "human"}
REQUIRED_COLUMNS = {"pos_id", "ref_PSI", "alt_PSI", "dPSI", "q"}
NEURONAL_ANALYSIS_IDS = {
    "E-MTAB-6293",
    "GSE138527",
    "GSE139090",
    "GSE140205",
    "GSE234488",
    "GSE271392",
}
NEURONAL_CONTRAST_NAMES = {"APOE4_vs_APOE3_neurons"}
REFERENCE_COLOR = "#0072B2"
ALTERNATIVE_COLOR = "#D55E00"
POSITIVE_DPSI_COLOR = "#D55E00"
NEGATIVE_DPSI_COLOR = "#0072B2"


plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--events", type=Path, default=ROOT / "scripts/visualization/events_of_interest.json",
        help="Species-aware events JSON file.",
    )
    parser.add_argument(
        "--data-root", type=Path,
        default=Path("/rhome/naotok/bigdata/neurodegeneration/GEO_SRA"),
        help="Directory containing <analysis_id>/Shiba/<contrast_name> results.",
    )
    parser.add_argument(
        "--output-dir", type=Path,
        default=Path("~/bigdata/neurodegeneration/figure/shiba_events").expanduser(),
        help="Directory for figures and TSV manifests.",
    )
    parser.add_argument(
        "--contrast", action="append", dest="contrasts",
        help="Restrict to one contrast_name; repeat this option to select several.",
    )
    parser.add_argument(
        "--event", action="append", dest="event_names",
        help="Restrict to one event name; repeat this option to select several.",
    )
    return parser.parse_args()


def read_events(path: Path) -> list[dict[str, Any]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    events = payload.get("events")
    if not isinstance(events, list) or not events:
        raise ValueError(f"{path} must contain a non-empty 'events' list")

    seen_ids: dict[str, set[str]] = {"mouse": set(), "human": set()}
    for event in events:
        if not isinstance(event, dict) or not isinstance(event.get("name"), str) or not event["name"].strip():
            raise ValueError("Every event must have a non-empty string 'name'")
        if event.get("event_type") not in EVENT_TYPES:
            raise ValueError(f"{event['name']}: unsupported event_type {event.get('event_type')!r}")
        pos_ids = event.get("pos_id")
        if not isinstance(pos_ids, dict) or set(pos_ids) != {"mouse", "human"}:
            raise ValueError(f"{event['name']}: pos_id must have exactly mouse and human keys")
        for species, pos_id in pos_ids.items():
            if not isinstance(pos_id, str):
                raise ValueError(f"{event['name']}: {species} pos_id must be a string")
            if pos_id and pos_id in seen_ids[species]:
                raise ValueError(f"Duplicate {species} pos_id: {pos_id}")
            seen_ids[species].add(pos_id)
    return events


def species_for_analysis(analysis_id: str) -> str:
    metadata_path = ROOT / "metadata/analysis" / f"{analysis_id}.tsv"
    metadata = pd.read_csv(metadata_path, sep="\t", dtype=str)
    species_values = set(metadata["species"].dropna().str.strip())
    if len(species_values) != 1:
        raise ValueError(f"Expected one species in {metadata_path}, found {sorted(species_values)}")
    species = SPECIES_ALIASES.get(species_values.pop())
    if species is None:
        raise ValueError(f"Unsupported species in {metadata_path}")
    return species


def read_shiba_table(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t")
    if len(table.columns) == 1:
        table = pd.read_csv(path, sep=r"\s+", engine="python")
    missing_columns = REQUIRED_COLUMNS - set(table.columns)
    if missing_columns:
        raise ValueError(f"Missing required columns: {sorted(missing_columns)}")
    return table


def sample_groups(analysis_id: str, contrast: pd.Series) -> dict[str, str]:
    experiment_path = ROOT / "config/splicing" / f"experiment_Shiba_{analysis_id}.tsv"
    experiment = pd.read_csv(experiment_path, sep="\t", dtype=str)
    required = {"sample", "group"}
    if missing_columns := required - set(experiment.columns):
        raise ValueError(f"{experiment_path} missing columns: {sorted(missing_columns)}")
    selected = experiment.loc[experiment["group"].isin([contrast["reference_group"], contrast["alternative_group"]])]
    return dict(zip(selected["sample"], selected["group"]))


def neuronal_contrasts(contrasts: pd.DataFrame) -> pd.DataFrame:
    mask = contrasts["analysis_id"].isin(NEURONAL_ANALYSIS_IDS) | contrasts["contrast_name"].isin(
        NEURONAL_CONTRAST_NAMES
    )
    return contrasts.loc[mask].copy()


def slug(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._-") or "event"


def status_row(event: dict[str, Any], contrast: pd.Series, species: str, status: str, reason: str, **values: Any) -> dict[str, Any]:
    return {
        "event_name": event["name"],
        "gene_name": event.get("gene_name", ""),
        "event_type": event["event_type"],
        "species": species,
        "pos_id": event["pos_id"][species],
        "analysis_id": contrast["analysis_id"],
        "contrast_name": contrast["contrast_name"],
        "reference_group": contrast["reference_group"],
        "alternative_group": contrast["alternative_group"],
        "status": status,
        "reason": reason,
        **values,
    }


def collect_event_results(
    events: list[dict[str, Any]], contrasts: pd.DataFrame, data_root: Path
) -> tuple[pd.DataFrame, pd.DataFrame]:
    records: list[dict[str, Any]] = []
    sample_records: list[dict[str, Any]] = []
    species_cache: dict[str, str | Exception] = {}
    table_cache: dict[tuple[str, str], pd.DataFrame | Exception] = {}
    groups_cache: dict[tuple[str, str], dict[str, str] | Exception] = {}

    for _, contrast in contrasts.iterrows():
        analysis_id = contrast["analysis_id"]
        if analysis_id not in species_cache:
            try:
                species_cache[analysis_id] = species_for_analysis(analysis_id)
            except Exception as error:
                species_cache[analysis_id] = error
        species_result = species_cache[analysis_id]

        for event in events:
            if isinstance(species_result, Exception):
                records.append(status_row(event, contrast, "unknown", "skipped", "species_error", detail=str(species_result)))
                continue
            species = species_result
            pos_id = event["pos_id"][species]
            if not pos_id:
                records.append(status_row(event, contrast, species, "skipped", "missing_species_pos_id"))
                continue

            result_path = data_root / analysis_id / "Shiba" / contrast["contrast_name"] / "results/splicing" / f"PSI_{event['event_type']}.txt"
            cache_key = (str(result_path), event["event_type"])
            if cache_key not in table_cache:
                try:
                    table_cache[cache_key] = read_shiba_table(result_path)
                except Exception as error:
                    table_cache[cache_key] = error
            table_result = table_cache[cache_key]
            if isinstance(table_result, Exception):
                records.append(status_row(event, contrast, species, "skipped", "result_unavailable", detail=str(table_result)))
                continue

            matches = table_result.loc[table_result["pos_id"] == pos_id]
            if len(matches) != 1:
                reason = "event_not_found" if matches.empty else "duplicate_pos_id"
                records.append(status_row(event, contrast, species, "skipped", reason))
                continue
            row = matches.iloc[0]
            records.append(status_row(
                event, contrast, species, "found", "", ref_PSI=float(row["ref_PSI"]) * 100,
                alt_PSI=float(row["alt_PSI"]) * 100, dPSI=float(row["dPSI"]) * 100,
                q=row["q"], p_ttest=row.get("p_ttest", pd.NA), result_path=str(result_path),
            ))
            group_key = (analysis_id, contrast["contrast_name"])
            if group_key not in groups_cache:
                try:
                    groups_cache[group_key] = sample_groups(analysis_id, contrast)
                except Exception as error:
                    groups_cache[group_key] = error
            group_result = groups_cache[group_key]
            if isinstance(group_result, Exception):
                records[-1]["reason"] = "sample_metadata_unavailable"
                records[-1]["detail"] = str(group_result)
                continue
            for sample, group in group_result.items():
                psi_column = f"{sample}_PSI"
                if psi_column not in row.index:
                    continue
                sample_records.append({
                    "event_name": event["name"], "event_type": event["event_type"], "species": species,
                    "pos_id": pos_id, "analysis_id": analysis_id, "contrast_name": contrast["contrast_name"],
                    "sample": sample, "group": group, "group_role": (
                        "reference" if group == contrast["reference_group"] else "alternative"
                    ), "PSI": float(row[psi_column]) * 100,
                })
    return pd.DataFrame(records), pd.DataFrame(sample_records)


def plot_dpsi(event: dict[str, Any], found: pd.DataFrame, output_base: Path) -> None:
    ordered = found.sort_values(["species", "contrast_name"]).reset_index(drop=True)
    fig, axis = plt.subplots(figsize=(7.08, max(2.2, 0.38 * len(ordered) + 0.8)))
    positions = list(range(len(ordered)))
    for position, row in ordered.iterrows():
        axis.hlines(position, row["ref_PSI"], row["alt_PSI"], color="#777777", linewidth=0.8, zorder=1)
        axis.scatter(row["ref_PSI"], position, s=38, color=REFERENCE_COLOR, edgecolor="white", linewidth=0.5, zorder=3)
        axis.scatter(row["alt_PSI"], position, s=38, color=ALTERNATIVE_COLOR, edgecolor="white", linewidth=0.5, zorder=3)
    axis.set_yticks(positions, ordered["contrast_name"])
    axis.invert_yaxis()
    axis.set_xlim(-5, 105)
    axis.set_xlabel("PSI (%)")
    axis.set_title(f"{event['name']}: reference and alternative PSI", loc="left", fontweight="bold")
    for index, row in ordered.iterrows():
        q_label = f"q={float(row['q']):.2g}" if pd.notna(row["q"]) else "q=NA"
        label = f"dPSI={float(row['dPSI']):+.1f}; {q_label}"
        axis.text(103, index, label, va="center", ha="right", fontsize=7)
    axis.spines[["top", "right"]].set_visible(False)
    axis.legend(
        handles=[Patch(facecolor=REFERENCE_COLOR, label="Reference"), Patch(facecolor=ALTERNATIVE_COLOR, label="Alternative")],
        loc="lower right", frameon=False, ncol=2,
    )
    fig.tight_layout()
    for suffix in ("png", "pdf"):
        fig.savefig(output_base.with_suffix(f".{suffix}"), dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_sample_psi(event: dict[str, Any], samples: pd.DataFrame, output_base: Path) -> None:
    contrast_names = list(samples["contrast_name"].drop_duplicates())
    columns = 3
    rows = ceil(len(contrast_names) / columns)
    fig, axes = plt.subplots(
        rows, columns, figsize=(7.08, max(2.2, 2.15 * rows)), squeeze=False
    )
    for axis, contrast_name in zip(axes.flat, contrast_names):
        subset = samples.loc[samples["contrast_name"] == contrast_name]
        groups = [
            subset.loc[subset["group_role"] == role, "group"].iat[0]
            for role in ("reference", "alternative")
            if (subset["group_role"] == role).any()
        ]
        values = [subset.loc[subset["group"] == group, "PSI"].dropna().astype(float).tolist() for group in groups]
        positions = list(range(1, len(groups) + 1))
        roles = [subset.loc[subset["group"] == group, "group_role"].iat[0] for group in groups]
        colors = [REFERENCE_COLOR if role == "reference" else ALTERNATIVE_COLOR for role in roles]
        boxes = axis.boxplot(
            values, positions=positions, widths=0.09, showfliers=False, patch_artist=True,
            medianprops={"color": "#222222", "linewidth": 0.8},
            whiskerprops={"color": "#555555", "linewidth": 0.6},
            capprops={"color": "#555555", "linewidth": 0.6},
        )
        for box, color in zip(boxes["boxes"], colors):
            box.set(facecolor=color, edgecolor=color, alpha=0.3, linewidth=0.8)
        for position, color, group_values in zip(positions, colors, values):
            count = len(group_values)
            offsets = [position] if count == 1 else [position - 0.12 + 0.24 * index / (count - 1) for index in range(count)]
            axis.scatter(offsets, group_values, color=color, edgecolor="white", linewidth=0.4, s=28, zorder=3)
        labels = [f"{role.title()}\n(n={len(group_values)})" for role, group_values in zip(roles, values)]
        axis.set_xticks(positions, labels)
        axis.set_ylim(-5, 105)
        axis.set_ylabel("PSI (%)")
        axis.set_title(contrast_name, loc="left", fontweight="bold")
        axis.spines[["top", "right"]].set_visible(False)
    for axis in axes.flat[len(contrast_names):]:
        axis.set_visible(False)
    fig.suptitle(f"{event['name']}: sample-level PSI", x=0.12, ha="left", fontweight="bold", fontsize=10)
    fig.legend(
        handles=[Patch(facecolor=REFERENCE_COLOR, label="Reference"), Patch(facecolor=ALTERNATIVE_COLOR, label="Alternative")],
        loc="upper right", bbox_to_anchor=(0.98, 1.02), frameon=False, ncol=2,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    for suffix in ("png", "pdf"):
        fig.savefig(output_base.with_suffix(f".{suffix}"), dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    events = read_events(args.events)
    if args.event_names:
        requested = set(args.event_names)
        events = [event for event in events if event["name"] in requested]
        if not events:
            raise SystemExit("No requested event names were found in the JSON file")

    contrasts = pd.read_csv(ROOT / "config/contrasts.tsv", sep="\t", dtype=str)
    if args.contrasts:
        contrasts = contrasts.loc[contrasts["contrast_name"].isin(args.contrasts)]
        if contrasts.empty:
            raise SystemExit("No requested contrast names were found in config/contrasts.tsv")
    contrasts = neuronal_contrasts(contrasts)
    if contrasts.empty:
        raise SystemExit("No requested contrasts are neuronal cell or nervous-system tissue comparisons")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    results, sample_results = collect_event_results(events, contrasts, args.data_root)
    results.to_csv(args.output_dir / "event_contrast_manifest.tsv", sep="\t", index=False)
    sample_results.to_csv(args.output_dir / "event_sample_psi.tsv", sep="\t", index=False)
    results.groupby(["event_name", "status", "reason"], dropna=False).size().rename("count").reset_index().to_csv(
        args.output_dir / "event_contrast_summary.tsv", sep="\t", index=False
    )

    for event in events:
        found = results.loc[(results["event_name"] == event["name"]) & (results["status"] == "found")].copy()
        if found.empty:
            continue
        event_dir = args.output_dir / f"{slug(event['name'])}_{event['event_type']}"
        event_dir.mkdir(exist_ok=True)
        plot_dpsi(event, found, event_dir / "dpsi_summary")
        event_samples = sample_results.loc[sample_results["event_name"] == event["name"]]
        if not event_samples.empty:
            plot_sample_psi(event, event_samples, event_dir / "sample_psi")

    print(f"Processed {len(contrasts)} neuronal contrasts and wrote manifest and figures to {args.output_dir}")


if __name__ == "__main__":
    main()