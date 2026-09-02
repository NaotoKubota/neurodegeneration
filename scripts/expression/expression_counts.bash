#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --job-name="featureCounts bulk"
#SBATCH -p epyc

set -euo pipefail
analysis_id=${1:?Usage: sbatch expression_counts.bash ANALYSIS_ID}
project_root=${PROJECT_ROOT:-/rhome/naotok/neurodegeneration}
data_root=${DATA_ROOT:-/rhome/naotok/bigdata/neurodegeneration/GEO_SRA}
snakemake_bin=${SNAKEMAKE_BIN:-/opt/linux/rocky/8.x/x86_64/pkgs/snakemake/9.16.3/env/bin/snakemake}
singularity_tmpdir=${SINGULARITY_TMPDIR:-${data_root}/.singularity/tmp}
[[ -x "${snakemake_bin}" ]] || { echo "Missing Snakemake: ${snakemake_bin}" >&2; exit 1; }
module load singularity
mkdir -p "${singularity_tmpdir}"
export SINGULARITY_TMPDIR="${singularity_tmpdir}"
"${snakemake_bin}" -s "${project_root}/scripts/expression/Snakefile_expression" \
    --configfile "${project_root}/config/qc_mapping/config_${analysis_id}.yaml" \
    --cores "${SLURM_CPUS_PER_TASK:-16}" --use-singularity \
    --singularity-args "--bind /rhome/naotok:/rhome/naotok" --rerun-incomplete
