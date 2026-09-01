#!/usr/bin/env python3
"""Generate SnakeNgs and Shiba inputs from run-level manifests."""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATA_ROOT = Path("/rhome/naotok/bigdata/neurodegeneration/GEO_SRA")
REFERENCES = {
    "Homo sapiens": {
        "star_index": "/rhome/naotok/bigdata/STAR/index_hg38_ensembl_v106",
        "gtf": "/rhome/naotok/bigdata/genome_annotation/Human/Homo_sapiens.GRCh38.106.gtf",
        "tag": "human",
    },
    "Mus musculus": {
        "star_index": "/rhome/naotok/bigdata/STAR/index_mm10_ensembl_v102",
        "gtf": "/rhome/naotok/bigdata/genome_annotation/Mouse/Mus_musculus.GRCm38.102.gtf",
        "tag": "mouse",
    },
}


def safe(value: str) -> str:
    value = value.replace("ε", "e")
    return re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_")


def sample_group(dataset: str, row: dict[str, str]) -> str:
    factors: list[str] = []
    title = row.get("title", "")
    keys_by_study = {
        "E-MTAB-6293": ("genotype", "age"),
        "GSE138527": ("genotype",),
        "GSE139090": ("genotype", "age"),
        "GSE271392": ("genotype", "tissue", "age"),
        "GSE102956": ("genotype", "cell_type"),
        "GSE163857": ("genotype", "background", "treatment", "sex"),
        "GSE140205": ("genotype", "sex"),
        "GSE234350": ("genotype", "cell_type"),
        "GSE234480": ("genotype",),
        "GSE234488": ("genotype",),
        "GSE254205": ("genotype", "treatment"),
    }
    for key in keys_by_study.get(dataset, ("genotype", "condition")):
        if row.get(key):
            factors.append(row[key])
    if dataset == "GSE254205" and "LD " in title:
        factors.append("LD_high" if "LD high" in title else "LD_low")
    if not factors:
        factors.append(title or row["sample"])
    return "__".join(safe(value) for value in factors)


def yaml_quote(value: str) -> str:
    return json.dumps(value, ensure_ascii=False)


def write_analysis(
    analysis_id: str,
    dataset: str,
    species: str,
    rows: list[dict[str, str]],
    data_root: Path,
) -> None:
    ref = REFERENCES[species]
    layouts = {row["layout"] for row in rows}
    if len(layouts) != 1:
        raise ValueError(f"{analysis_id}: mixed layouts {layouts}")
    by_sample: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_sample[row["sample"]].append(row)
    workdir = data_root / analysis_id

    config_dir = ROOT / "config" / "qc_mapping"
    config_dir.mkdir(parents=True, exist_ok=True)
    with (config_dir / f"config_{analysis_id}.yaml").open("w", encoding="utf-8") as handle:
        handle.write(f"workdir: {yaml_quote(str(workdir))}\n")
        handle.write("samples:\n")
        for sample, sample_rows in sorted(by_sample.items()):
            runs = ", ".join(yaml_quote(row["run"]) for row in sample_rows)
            handle.write(f"  {yaml_quote(sample)}: [{runs}]\n")
        handle.write(f"star_index: {yaml_quote(ref['star_index'])}\n")
        handle.write(f"gtf: {yaml_quote(ref['gtf'])}\n")
        handle.write(f"layout: {yaml_quote(next(iter(layouts)))}\n")

    sample_columns = [
        "sample", "group", "species", "title", "genotype", "age", "sex", "tissue",
        "cell_type", "disease", "treatment", "background", "batch", "runs",
    ]
    table_dir = ROOT / "metadata" / "analysis"
    table_dir.mkdir(parents=True, exist_ok=True)
    with (table_dir / f"{analysis_id}.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, sample_columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for sample, sample_rows in sorted(by_sample.items()):
            row = sample_rows[0]
            writer.writerow({
                **{key: row.get(key, "") for key in sample_columns},
                "sample": sample,
                "group": sample_group(dataset, row),
                "runs": ",".join(item["run"] for item in sample_rows),
            })

    splice_dir = ROOT / "config" / "splicing"
    splice_dir.mkdir(parents=True, exist_ok=True)
    experiment = splice_dir / f"experiment_Shiba_{analysis_id}.tsv"
    with experiment.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("sample", "bam", "group"))
        for sample, sample_rows in sorted(by_sample.items()):
            row = sample_rows[0]
            bam = workdir / "star" / sample / f"{sample}_Aligned.out.bam"
            writer.writerow((sample, bam, sample_group(dataset, row)))

    groups = sorted({sample_group(dataset, rows_[0]) for rows_ in by_sample.values()})
    with (splice_dir / f"groups_{analysis_id}.txt").open("w", encoding="utf-8") as handle:
        handle.write("\n".join(groups) + "\n")

    base_config = splice_dir / f"config_Shiba_{analysis_id}.yaml"
    with base_config.open("w", encoding="utf-8") as handle:
        handle.write(f"workdir: {yaml_quote(str(workdir / 'Shiba'))}\n")
        handle.write("container: docker://naotokubota/shiba:latest\n")
        handle.write(f"gtf: {yaml_quote(ref['gtf'])}\n")
        handle.write(f"experiment_table: {yaml_quote(str(experiment))}\n")
        handle.write("unannotated: true\nminimum_anchor_length: 6\nminimum_intron_length: 70\n")
        handle.write("maximum_intron_length: 500000\nstrand: XS\nonly_psi: false\nonly_psi_group: false\n")
        handle.write("fdr: 0.05\ndelta_psi: 0.1\nreference_group: EDIT_REFERENCE_GROUP\n")
        handle.write("alternative_group: EDIT_ALTERNATIVE_GROUP\nminimum_reads: 10\n")
        handle.write("individual_psi: true\nttest: true\nexcel: true\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", type=Path, default=ROOT / "metadata")
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    args = parser.parse_args()
    for manifest in sorted(args.metadata.glob("*.tsv")):
        with manifest.open(encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        if not rows or "dataset" not in rows[0]:
            continue
        dataset = rows[0]["dataset"]
        by_species: dict[str, list[dict[str, str]]] = defaultdict(list)
        for row in rows:
            by_species[row["species"]].append(row)
        for species, species_rows in by_species.items():
            by_layout: dict[str, list[dict[str, str]]] = defaultdict(list)
            for row in species_rows:
                by_layout[row["layout"]].append(row)
            for layout, analysis_rows in by_layout.items():
                suffixes: list[str] = []
                if len(by_species) > 1:
                    suffixes.append(REFERENCES[species]["tag"])
                if len(by_layout) > 1:
                    suffixes.append(layout)
                suffix = "_" + "_".join(suffixes) if suffixes else ""
                analysis_id = dataset + suffix
                write_analysis(analysis_id, dataset, species, analysis_rows, args.data_root)
                print(f"{analysis_id}: {len({row['sample'] for row in analysis_rows})} samples")


if __name__ == "__main__":
    main()
