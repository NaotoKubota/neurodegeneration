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
snakemake_bin=${SNAKEMAKE_BIN:-/opt/linux/rocky/8.x/x86_64/pkgs/snakemake/9.16.3/env/bin/snakemake}
[[ -x "${snakemake_bin}" ]] || { echo "Missing Snakemake: ${snakemake_bin}" >&2; exit 1; }
module load singularity
"${snakemake_bin}" -s "${project_root}/scripts/expression/Snakefile_expression" \
    --configfile "${project_root}/config/qc_mapping/config_${analysis_id}.yaml" \
    --cores "${SLURM_CPUS_PER_TASK:-16}" --use-singularity \
    --singularity-args "--bind /rhome/naotok:/rhome/naotok" --rerun-incomplete
