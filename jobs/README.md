# Slurm job records

The TSV files in this directory are immutable submission audit records. They
include completed jobs as well as failed or cancelled submissions that were
replaced by later jobs.

- `submitted_20260831.tsv` and `contrasts_20260831.tsv`: initial submissions.
- `resubmitted_20260901.tsv` and `contrasts_20260901.tsv`: submissions made
  after incomplete FASTQ files were removed.
- `reconfigured_20260901.tsv`: submissions using serialized STAR jobs and
  contrast-specific Shiba working directories.
- `mouse_lock_retry_20260901.tsv`: mouse analyses resubmitted after stale
  Snakemake locks were cleared; this supersedes the affected mouse rows in
  `reconfigured_20260901.tsv`.
- `incomplete_retry_20260902.tsv`: initial retry submissions for incomplete
  datasets. These were cancelled after detecting that a host installation
  shadowed the requested container version of ngsfetch.
- `incomplete_retry_v012_20260902.tsv`: replacement submissions that run
  ngsfetch v0.1.2 with a clean container environment.
- `GSE234488_shiba_retry_20260902.tsv`: final standalone Shiba retry for
  GSE234488.

Slurm logs, quarantined FASTQ files, and partial download files are generated
artifacts and are not tracked in Git.
