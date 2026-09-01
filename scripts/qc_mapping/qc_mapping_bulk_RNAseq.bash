#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=80G
#SBATCH --job-name="QC bulk RNAseq"
#SBATCH -p epyc

set -euo pipefail

analysis_id=${1:?Usage: sbatch qc_mapping_bulk_RNAseq.bash ANALYSIS_ID}
project_root=${PROJECT_ROOT:-/rhome/naotok/neurodegeneration}
snakengs_root=${SNAKENGS_ROOT:-/rhome/naotok/SnakeNgs}
snakemake_bin=${SNAKEMAKE_BIN:-/opt/linux/rocky/8.x/x86_64/pkgs/snakemake/9.16.3/env/bin/snakemake}
config="${project_root}/config/qc_mapping/config_${analysis_id}.yaml"

[[ -s "${config}" ]] || { echo "Missing config: ${config}" >&2; exit 1; }
[[ -x "${snakemake_bin}" ]] || { echo "Missing Snakemake: ${snakemake_bin}" >&2; exit 1; }
module load singularity
"${snakemake_bin}" -s "${snakengs_root}/snakefile/preprocessing_RNAseq.smk" \
    --configfile "${config}" \
    --cores "${SLURM_CPUS_PER_TASK:-64}" \
    --use-singularity \
    --singularity-args "--bind /rhome/naotok:/rhome/naotok" \
    --rerun-incomplete
