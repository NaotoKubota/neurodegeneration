#!/usr/bin/env python3
"""Build run-level manifests for the selected bulk RNA-seq studies.

Metadata are fetched from GEO SOFT / NCBI SRA RunInfo, or from the
ArrayExpress SDRF for E-MTAB-6293.  Single-cell, single-nucleus and non-RNA-seq
samples are rejected even when they share a mixed GEO Series with bulk RNA-seq.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path


GEO_STUDIES = (
    "GSE138527", "GSE139090", "GSE271392", "GSE102956", "GSE163857",
    "GSE140205", "GSE234350", "GSE234480", "GSE234488", "GSE254205",
)

EXCLUDED_STUDIES = {
    "GSE164507": "single-nucleus RNA-seq",
    "GSE213446": "single-nucleus RNA-seq",
    "GSE223719": "single-cell RNA-seq",
    "GSE237718": "single-nucleus RNA-seq",
    "GSE242153": "single-nucleus RNA-seq",
    "GSE248020": "single-cell RNA-seq",
    "GSE295612": "single-cell RNA-seq",
}

BASE_COLUMNS = (
    "dataset", "sample", "run", "experiment", "species", "layout", "title",
    "source", "library_strategy", "library_source", "library_selection",
    "instrument", "fastq_url", "fastq_md5", "fastq_bytes", "characteristics",
)

SINGLE_CELL_RE = re.compile(
    r"single[- ]?(?:cell|nucle(?:us|i))|snrna|scrna|10x genomics|chromium",
    re.IGNORECASE,
)


def fetch(url: str, attempts: int = 5) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": "bulk-rnaseq-manifest/1.0"})
    for attempt in range(attempts):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                return response.read()
        except Exception:
            if attempt + 1 == attempts:
                raise
            time.sleep(2 ** attempt)
    raise RuntimeError("unreachable")


def clean_key(value: str) -> str:
    value = value.strip().lower().replace("ε", "e")
    value = re.sub(r"[^a-z0-9]+", "_", value).strip("_")
    aliases = {
        "sex": "sex", "gender": "sex", "age": "age", "genotype": "genotype",
        "genotype_variation": "genotype", "timepoint": "age",
        "tissue": "tissue", "organism_part": "tissue", "cell_type": "cell_type",
        "tissue_cell_type": "cell_type", "disease": "disease",
        "disease_status": "disease", "treatment": "treatment",
        "background": "background", "sequencing_batch": "batch",
    }
    return aliases.get(value, value)


def parse_soft(data: bytes) -> tuple[dict[str, list[str]], list[dict[str, object]]]:
    text = gzip.decompress(data).decode("utf-8", errors="replace")
    series: dict[str, list[str]] = defaultdict(list)
    samples: list[dict[str, object]] = []
    current: dict[str, object] | None = None
    for line in text.splitlines():
        if line.startswith("^SAMPLE = "):
            if current:
                samples.append(current)
            current = {"sample": line.split(" = ", 1)[1], "characteristics": {}}
            continue
        if line.startswith("^PLATFORM = "):
            if current:
                samples.append(current)
                current = None
            continue
        if line.startswith("!Series_"):
            key, _, value = line.partition(" = ")
            series[key[8:]].append(value)
            continue
        if current is None or not line.startswith("!Sample_"):
            continue
        key, _, value = line.partition(" = ")
        key = key[8:]
        if key == "characteristics_ch1":
            char_key, sep, char_value = value.partition(":")
            if sep:
                chars = current["characteristics"]
                assert isinstance(chars, dict)
                chars[clean_key(char_key)] = char_value.strip()
        elif key in {"extract_protocol_ch1", "growth_protocol_ch1", "treatment_protocol_ch1"}:
            current[key] = str(current.get(key, "")) + " " + value
        else:
            current[key] = value
    if current:
        samples.append(current)
    return dict(series), samples


def study_accession(series: dict[str, list[str]]) -> str:
    relations = series.get("relation", [])
    for pattern in (r"term=(SRP\d+)", r"bioproject/(PRJNA\d+)"):
        for relation in relations:
            match = re.search(pattern, relation, re.IGNORECASE)
            if match:
                return match.group(1).upper()
    raise ValueError(f"No SRA/BioProject relation in {relations}")


def runinfo(accession: str) -> list[dict[str, str]]:
    encoded = urllib.parse.quote(accession)
    url = f"https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc={encoded}"
    text = fetch(url).decode("utf-8-sig", errors="replace")
    rows = list(csv.DictReader(io.StringIO(text)))
    if not rows or "Run" not in rows[0]:
        raise ValueError(f"No SRA RunInfo returned for {accession}")
    return rows


def is_bulk_rnaseq(sample: dict[str, object], runs: list[dict[str, str]]) -> bool:
    strategy = str(sample.get("library_strategy", ""))
    if strategy and strategy.lower() != "rna-seq":
        return False
    run_strategies = {run.get("LibraryStrategy", "").lower() for run in runs}
    if run_strategies and run_strategies != {"rna-seq"}:
        return False
    sources = {run.get("LibrarySource", "").lower() for run in runs}
    if any("single cell" in source for source in sources):
        return False
    searchable = " ".join(
        str(sample.get(key, ""))
        for key in ("title", "source_name_ch1", "extract_protocol_ch1")
    )
    return SINGLE_CELL_RE.search(searchable) is None


def geo_rows(accession: str) -> list[dict[str, str]]:
    bucket = accession[:-3] + "nnn"
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{bucket}/{accession}/soft/{accession}_family.soft.gz"
    series, samples = parse_soft(fetch(url))
    runs = runinfo(study_accession(series))
    by_sample: dict[str, list[dict[str, str]]] = defaultdict(list)
    by_experiment: dict[str, list[dict[str, str]]] = defaultdict(list)
    for run in runs:
        by_sample[run.get("SampleName", "")].append(run)
        by_experiment[run.get("Experiment", "")].append(run)

    output: list[dict[str, str]] = []
    for sample in samples:
        gsm = str(sample["sample"])
        sample_runs = by_sample.get(gsm, [])
        if not sample_runs:
            relation = str(sample.get("relation", ""))
            match = re.search(r"term=(SRX\d+)", relation)
            if match:
                sample_runs = by_experiment.get(match.group(1), [])
        if not sample_runs:
            print(f"warning: {accession}/{gsm} has no mapped run", file=sys.stderr)
            continue
        sample["title"] = sample.get("title", "")
        if not is_bulk_rnaseq(sample, sample_runs):
            continue
        chars = sample.get("characteristics", {})
        assert isinstance(chars, dict)
        title = str(sample.get("title", ""))
        if "genotype" not in chars:
            for pattern, label in (
                (r"\bAPOE[ _]?3(?:\b|_)", "APOE3"), (r"\bAPOE[ _]?4(?:\b|_)", "APOE4"),
                (r"\bWTQ?\b", "WT"), (r"\bKihet\b", "SCA7_KI"),
                (r"\b92Q\b", "SCA7_92Q"),
            ):
                if re.search(pattern, title, re.IGNORECASE):
                    chars["genotype"] = label
                    break
        for run in sample_runs:
            fastq_url = run.get("download_path", "")
            output.append({
                "dataset": accession,
                "sample": gsm,
                "run": run.get("Run", ""),
                "experiment": run.get("Experiment", ""),
                "species": run.get("ScientificName", str(sample.get("organism_ch1", ""))),
                "layout": run.get("LibraryLayout", "").lower(),
                "title": str(sample.get("title", "")),
                "source": str(sample.get("source_name_ch1", "")),
                "library_strategy": run.get("LibraryStrategy", ""),
                "library_source": run.get("LibrarySource", ""),
                "library_selection": run.get("LibrarySelection", ""),
                "instrument": run.get("Model", ""),
                "fastq_url": fastq_url,
                "fastq_md5": "",
                "fastq_bytes": "",
                "characteristics": json.dumps(chars, ensure_ascii=False, sort_keys=True),
                **{key: str(value) for key, value in chars.items()},
            })
    return output


def arrayexpress_rows() -> list[dict[str, str]]:
    accession = "E-MTAB-6293"
    url = f"https://www.ebi.ac.uk/biostudies/files/{accession}/{accession}.sdrf.txt"
    text = fetch(url).decode("utf-8-sig", errors="replace")
    output: list[dict[str, str]] = []
    for row in csv.DictReader(io.StringIO(text), delimiter="\t"):
        chars: dict[str, str] = {}
        for key, value in row.items():
            match = re.match(r"(?:Characteristics|Factor Value)\[(.+)\]", key)
            if match and value:
                chars[clean_key(match.group(1))] = value
        age_unit = row.get("Unit[time unit]", "")
        if chars.get("age") and age_unit:
            chars["age"] = f"{chars['age']} {age_unit}"
        run = row["Comment[ENA_RUN]"]
        fastq_url = row["Comment[FASTQ_URI]"].replace("ftp://", "https://")
        output.append({
            "dataset": accession, "sample": row["Comment[ENA_SAMPLE]"], "run": run,
            "experiment": row["Comment[ENA_EXPERIMENT]"], "species": row["Characteristics[organism]"],
            "layout": row["Comment[LIBRARY_LAYOUT]"].lower(), "title": row["Source Name"],
            "source": row["Characteristics[organism part]"],
            "library_strategy": row["Comment[LIBRARY_STRATEGY]"],
            "library_source": row["Comment[LIBRARY_SOURCE]"],
            "library_selection": row["Comment[LIBRARY_SELECTION]"], "instrument": "Illumina HiSeq 2000",
            "fastq_url": fastq_url, "fastq_md5": "", "fastq_bytes": "",
            "characteristics": json.dumps(chars, ensure_ascii=False, sort_keys=True),
            **chars,
        })
    return output


def write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    extra = sorted({key for row in rows for key in row if key not in BASE_COLUMNS})
    columns = list(BASE_COLUMNS) + extra
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, columns, delimiter="\t", extrasaction="ignore", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda row: (row["sample"], row["run"])))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("metadata"))
    parser.add_argument("--study", action="append", choices=("E-MTAB-6293",) + GEO_STUDIES)
    args = parser.parse_args()
    studies = args.study or ["E-MTAB-6293", *GEO_STUDIES]
    for accession in studies:
        rows = arrayexpress_rows() if accession == "E-MTAB-6293" else geo_rows(accession)
        if not rows:
            raise RuntimeError(f"Bulk RNA-seq filtering removed every sample in {accession}")
        write_manifest(args.output / f"{accession}.tsv", rows)
        print(f"{accession}: {len({row['sample'] for row in rows})} samples, {len(rows)} runs")


if __name__ == "__main__":
    main()
