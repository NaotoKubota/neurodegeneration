#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --job-name="Shiba bulk RNAseq"
#SBATCH -p epyc

set -euo pipefail
module load singularity
project_root=${PROJECT_ROOT:-/rhome/naotok/neurodegeneration}
data_root=${DATA_ROOT:-/rhome/naotok/bigdata/neurodegeneration/GEO_SRA}
singularity_tmpdir=${SINGULARITY_TMPDIR:-${data_root}/.singularity/tmp}
job_tmpdir=${PIPELINE_TMPDIR:-${data_root}/.tmp/${SLURM_JOB_ID:-$$}}
mkdir -p "${singularity_tmpdir}" "${job_tmpdir}"
export SINGULARITY_TMPDIR="${singularity_tmpdir}"
export TMPDIR="${job_tmpdir}"
python3 "${project_root}/scripts/splicing/run_shiba.py" "$@"
