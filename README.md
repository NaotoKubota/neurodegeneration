# SCA2 / SCA7 / APOE4 bulk RNA-seq analysis

This repository provides a public bulk RNA-seq analysis pipeline modeled after
the directory structure and tool selection in
`/rhome/naotok/NeuronalMaturationSplicing`. It uses `ngsfetch` v0.1.2 for data
retrieval; `SnakeNgs` with fastp, STAR, samtools, Picard, deepTools, and MultiQC
for preprocessing; featureCounts and DESeq2 for differential gene expression;
and Shiba for differential splicing.

## Included datasets

The inclusion status and exclusion reason for every requested dataset are
recorded in `config/datasets.tsv`. All scRNA-seq and snRNA-seq datasets are
excluded. For the mixed-series dataset `GSE254205`, only the 17 bulk RNA-seq
samples from iPSC-derived microglia are retained. `GSE163857` is split into
separate human and mouse analysis units, while `GSE102956` is split by
single-end and paired-end library layout.

The current official SDRF for E-MTAB-6293 contains 64 runs: 16 at P1, 16 at
3 weeks, and 32 at 6 weeks. Therefore, the manifests include all 64 official
runs covering the three time points, rather than the 32 runs stated in the
initial dataset description.

## Directory structure

```text
config/
  datasets.tsv                 # Dataset inclusion and exclusion decisions
  contrasts.tsv                # Disease/control contrasts used downstream
  qc_mapping/                  # Auto-generated SnakeNgs configurations
  splicing/                    # Shiba configurations, experiment tables, groups
metadata/
  <accession>.tsv              # Official run-level metadata
  analysis/                    # Biological sample-level metadata
scripts/
  metadata/                    # Metadata and configuration generation
  download/                    # ngsfetch data retrieval
  qc_mapping/                  # fastp -> STAR -> alignment QC
  expression/                  # featureCounts -> DESeq2
  splicing/                    # Shiba
jobs/                          # Submitted Slurm job IDs and dependencies
logs/slurm/                    # Slurm output logs
```

## Regenerating metadata and configurations

```bash
cd /rhome/naotok/neurodegeneration
python3 scripts/metadata/build_manifests.py --output metadata
python3 scripts/metadata/make_configs.py
```

Reference genomes, the data root, and the SnakeNgs and Shiba locations default
to paths available in the existing HPCC environment. For another environment,
use `--data-root` when generating configurations or set `PROJECT_ROOT`,
`DATA_ROOT`, `SNAKENGS_ROOT`, and `SNAKEMAKE_BIN` when running the workflows.

## Example workflow

Available analysis-unit identifiers correspond to the filenames under
`metadata/analysis/`, without the `.tsv` extension.

```bash
# 1. Download only the selected bulk RNA-seq runs
sbatch scripts/download/Download_bulk_RNAseq.bash GSE139090

# 2. Run QC, STAR alignment, BAM/bigWig generation, and MultiQC
sbatch scripts/qc_mapping/qc_mapping_bulk_RNAseq.bash GSE139090

# 3. Generate gene-level counts
sbatch scripts/expression/expression_counts.bash GSE139090

# 4. Run DESeq2 using columns from metadata/analysis as factors or covariates
sbatch scripts/expression/DESeq2_bulk_RNAseq.bash \
  GSE139090 '~ genotype' 'genotype,fxSCA7 92Q,wild type' \
  'age=12 weeks' SCA7_vs_WT_12wk

# 5. Run Shiba using exact group names from the corresponding groups file
cat config/splicing/groups_GSE139090.txt
sbatch scripts/splicing/Shiba_bulk_RNAseq.bash \
  GSE139090 wild_type__12_weeks fxSCA7_92Q__12_weeks
```

Each Shiba contrast uses an independent output directory at
`<analysis_dir>/Shiba/<contrast_name>/`. STAR mapping is limited to one
concurrent Snakemake job per analysis to avoid memory overcommitment.

Within each submitted analysis chain, Slurm `afterok` dependencies enforce the
order `download -> QC/STAR -> featureCounts -> DESeq2`. Shiba depends directly
on successful QC and alignment. Downloads are additionally serialized across
datasets to avoid NCBI API rate limiting. Initial job IDs are recorded in
`jobs/submitted_20260831.tsv` and `jobs/contrasts_20260831.tsv`; replacement
jobs submitted after FASTQ cleanup are recorded in
`jobs/resubmitted_20260901.tsv` and `jobs/contrasts_20260901.tsv`. Subsequent
retry and resource-adjustment records are described in `jobs/README.md`.

## Analysis considerations

- Analyze GSE271392 separately by tissue and age. Do not combine its medulla or
  cervical spinal cord samples with cerebellar datasets.
- Do not combine the human and mouse subsets of GSE163857.
- The cervical spinal cord P8 comparison in GSE271392 is not included in the
  statistical contrasts because the wild-type group contains only one sample.
- DESeq2 designs can include relevant covariates, for example
  `~ sex + genotype` or `~ batch + genotype`.

## Visualizing selected Shiba events

Use `scripts/visualization/plot_shiba_events.py` to plot curated splicing
events from completed Shiba contrasts. It reads every contrast in
`config/contrasts.tsv` by default and searches for results at
`<data-root>/<analysis_id>/Shiba/<contrast_name>/results/splicing/`.

The input list is `scripts/visualization/events_of_interest.json`. Each event
has a display title and an exact Shiba `pos_id` for each species. Add a human
coordinate to `pos_id.human` only after it has been curated; an empty value is
recorded as `missing_species_pos_id` and is never matched against mouse output.

```json
{
  "name": "Example exon",
  "gene_name": "EXAMPLE",
  "event_type": "SE",
  "pos_id": {
    "mouse": "SE@chr1@...",
    "human": "SE@chr1@..."
  },
  "note": ""
}
```

The script determines the target species from the sole `species` value in
`metadata/analysis/<analysis_id>.tsv`, then matches only that species' `pos_id`
exactly. It needs Python packages `pandas` and `matplotlib`.

```bash
# Scan all registered completed contrasts.
python3 scripts/visualization/plot_shiba_events.py

# Recreate figures for selected contrasts and events.
python3 scripts/visualization/plot_shiba_events.py \
  --contrast SCA7_vs_WT_cervical_5wk \
  --event 'Bak1 E5'
```

Outputs are written to `~/bigdata/neurodegeneration/figure/shiba_events/` by
default. Use `--output-dir` to place a run elsewhere. Each found event receives
`sample_psi` and `dpsi_summary` figures in PNG and PDF. Figures use Arial with
7--10 pt text, embedded TrueType fonts in PDF, thin axes, and the
colorblind-accessible blue/orange Okabe-Ito palette for publication-ready
export.
`event_contrast_manifest.tsv` records one event per contrast, including skipped
results and their reason; `event_sample_psi.tsv` holds plotted sample PSI
values. All PSI values, including `ref_PSI`, `alt_PSI`, and `dPSI`, are written
as percentages; dPSI is expressed in percentage points. `event_contrast_summary.tsv`
provides a count summary.
