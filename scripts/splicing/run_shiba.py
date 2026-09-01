#!/usr/bin/env python3
"""Validate a Shiba contrast, render its config, and invoke Snakemake."""

from __future__ import annotations

import argparse
import os
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("analysis_id")
    parser.add_argument("reference_group")
    parser.add_argument("alternative_group")
    parser.add_argument("--cores", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "32")))
    parser.add_argument("--snakefile", default="/rhome/naotok/bigdata/Shiba/stable_release/Shiba/snakeshiba.smk")
    args = parser.parse_args()
    groups_file = ROOT / "config" / "splicing" / f"groups_{args.analysis_id}.txt"
    groups = set(groups_file.read_text(encoding="utf-8").splitlines())
    requested = {args.reference_group, args.alternative_group}
    if not requested <= groups:
        raise SystemExit(f"Unknown group(s): {sorted(requested - groups)}; see {groups_file}")
    template = ROOT / "config" / "splicing" / f"config_Shiba_{args.analysis_id}.yaml"
    config = template.read_text(encoding="utf-8")
    config = config.replace("EDIT_REFERENCE_GROUP", args.reference_group)
    config = config.replace("EDIT_ALTERNATIVE_GROUP", args.alternative_group)
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
