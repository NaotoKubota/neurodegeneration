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
python3 "${project_root}/scripts/splicing/run_shiba.py" "$@"
