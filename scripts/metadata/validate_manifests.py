#!/usr/bin/env python3
"""Static integrity checks for generated bulk RNA-seq manifests/configs."""

from __future__ import annotations

import csv
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
FORBIDDEN = re.compile(r"single[- ]?(?:cell|nucle(?:us|i))|snrna|scrna|10x genomics|chromium", re.I)


def main() -> None:
    manifests = sorted((ROOT / "metadata").glob("*.tsv"))
    if len(manifests) != 11:
        raise SystemExit(f"Expected 11 included study manifests, found {len(manifests)}")
    total_samples = 0
    total_runs = 0
    for path in manifests:
        with path.open(encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        if not rows:
            raise SystemExit(f"Empty manifest: {path}")
        for row in rows:
            if row["library_strategy"].lower() != "rna-seq":
                raise SystemExit(f"Non-RNA-seq row: {path}/{row['run']}")
            if "single cell" in row["library_source"].lower():
                raise SystemExit(f"Single-cell source: {path}/{row['run']}")
            if FORBIDDEN.search(" ".join((row["title"], row["source"]))):
                raise SystemExit(f"Single-cell/nucleus metadata: {path}/{row['run']}")
            if row["layout"] not in {"single", "paired"}:
                raise SystemExit(f"Invalid layout: {path}/{row['run']}")
            if not row["sample"] or not row["run"] or not row["species"]:
                raise SystemExit(f"Missing identifier/species: {path}")
        samples = {row["sample"] for row in rows}
        runs = {row["run"] for row in rows}
        if len(runs) != len(rows):
            raise SystemExit(f"Duplicate run accession: {path}")
        total_samples += len(samples)
        total_runs += len(runs)
        print(f"OK {path.name}: {len(samples)} samples, {len(runs)} runs")
    print(f"OK total: {total_samples} study-samples, {total_runs} runs")


if __name__ == "__main__":
    main()
