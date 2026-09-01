# SCA2 / SCA7 / APOE4 bulk RNA-seq analysis

This repository provides a public bulk RNA-seq analysis pipeline modeled after
the directory structure and tool selection in
`/rhome/naotok/NeuronalMaturationSplicing`. It uses `ngsfetch` for data
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

Within each submitted analysis chain, Slurm `afterok` dependencies enforce the
order `download -> QC/STAR -> featureCounts -> DESeq2`. Shiba depends directly
on successful QC and alignment. Downloads are additionally serialized across
datasets to avoid NCBI API rate limiting. Submitted job IDs are recorded in
`jobs/submitted_20260831.tsv` and `jobs/contrasts_20260831.tsv`.

## Analysis considerations

- Analyze GSE271392 separately by tissue and age. Do not combine its medulla or
  cervical spinal cord samples with cerebellar datasets.
- Do not combine the human and mouse subsets of GSE163857.
- The cervical spinal cord P8 comparison in GSE271392 is not included in the
  statistical contrasts because the wild-type group contains only one sample.
- DESeq2 designs can include relevant covariates, for example
  `~ sex + genotype` or `~ batch + genotype`.
