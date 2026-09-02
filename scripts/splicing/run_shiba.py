#!/usr/bin/env python3
"""Validate a Shiba contrast, render its config, and invoke Snakemake."""

from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
import tempfile
from pathlib import Path
from typing import Optional


ROOT = Path(__file__).resolve().parents[2]


def safe_path_component(value: str) -> str:
    """Return a filesystem-safe contrast identifier."""
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._-")
    if not cleaned:
        raise ValueError("Contrast name does not contain a usable path component")
    return cleaned


def resolve_contrast_name(
    analysis_id: str,
    reference_group: str,
    alternative_group: str,
    explicit_name: Optional[str],
) -> str:
    if explicit_name:
        return safe_path_component(explicit_name)

    contrasts_file = ROOT / "config" / "contrasts.tsv"
    with contrasts_file.open(encoding="utf-8", newline="") as handle:
        matches = [
            row["contrast_name"]
            for row in csv.DictReader(handle, delimiter="\t")
            if row["analysis_id"] == analysis_id
            and row["reference_group"] == reference_group
            and row["alternative_group"] == alternative_group
        ]
    if len(matches) > 1:
        raise ValueError("Multiple contrast names match the requested groups")
    if matches:
        return safe_path_component(matches[0])
    return safe_path_component(f"{alternative_group}_vs_{reference_group}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("analysis_id")
    parser.add_argument("reference_group")
    parser.add_argument("alternative_group")
    parser.add_argument("--contrast-name")
    parser.add_argument("--cores", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "32")))
    parser.add_argument("--snakefile", default="/rhome/naotok/bigdata/Shiba/stable_release/Shiba/snakeshiba.smk")
    args = parser.parse_args()
    groups_file = ROOT / "config" / "splicing" / f"groups_{args.analysis_id}.txt"
    groups = set(groups_file.read_text(encoding="utf-8").splitlines())
    requested = {args.reference_group, args.alternative_group}
    if not requested <= groups:
        raise SystemExit(f"Unknown group(s): {sorted(requested - groups)}; see {groups_file}")
    template = ROOT / "config" / "splicing" / f"config_Shiba_{args.analysis_id}.yaml"
    contrast_name = resolve_contrast_name(
        args.analysis_id,
        args.reference_group,
        args.alternative_group,
        args.contrast_name,
    )
    config = template.read_text(encoding="utf-8")
    config = config.replace("EDIT_REFERENCE_GROUP", args.reference_group)
    config = config.replace("EDIT_ALTERNATIVE_GROUP", args.alternative_group)
    base_workdir_match = re.search(r'^workdir:\s*["\']?([^"\'\n]+)["\']?\s*$', config, re.MULTILINE)
    if not base_workdir_match:
        raise SystemExit(f"Missing workdir in {template}")
    contrast_workdir = Path(base_workdir_match.group(1)) / contrast_name
    config = re.sub(
        r'^workdir:\s*.*$',
        f'workdir: "{contrast_workdir}"',
        config,
        count=1,
        flags=re.MULTILINE,
    )
    # Cluster /tmp can be capacity constrained; keep the transient config beside
    # the versioned template and always remove it in the finally block.
    with tempfile.NamedTemporaryFile(
        "w", prefix=f".{args.analysis_id}.", suffix=".yaml", dir=template.parent, delete=False
    ) as handle:
        handle.write(config)
        rendered = handle.name
    try:
        snakemake_bin = os.environ.get(
            "SNAKEMAKE_BIN",
            "/opt/linux/rocky/8.x/x86_64/pkgs/snakemake/9.16.3/env/bin/snakemake",
        )
        subprocess.run([
            snakemake_bin, "-s", args.snakefile, "--configfile", rendered,
            "--cores", str(args.cores), "--use-singularity",
            "--singularity-args", "--bind /rhome/naotok:/rhome/naotok", "--rerun-incomplete",
        ], check=True)
    finally:
        Path(rendered).unlink(missing_ok=True)


if __name__ == "__main__":
    main()
