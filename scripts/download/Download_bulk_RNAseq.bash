#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=20G
#SBATCH --job-name="Fetch bulk RNAseq"
#SBATCH -p epyc

set -euo pipefail

analysis_id=${1:?Usage: sbatch Download_bulk_RNAseq.bash ANALYSIS_ID}
project_root=${PROJECT_ROOT:-/rhome/naotok/neurodegeneration}
data_root=${DATA_ROOT:-/rhome/naotok/bigdata/neurodegeneration/GEO_SRA}
manifest="${project_root}/metadata/analysis/${analysis_id}.tsv"
singularity_tmpdir=${SINGULARITY_TMPDIR:-${data_root}/.singularity/tmp}

[[ -s "${manifest}" ]] || { echo "Missing manifest: ${manifest}" >&2; exit 1; }
module load singularity
mkdir -p "${data_root}/${analysis_id}" "${singularity_tmpdir}"
export SINGULARITY_TMPDIR="${singularity_tmpdir}"

awk -F '\t' 'NR == 1 { for (i = 1; i <= NF; i++) { gsub(/\r/, "", $i); if ($i == "runs") column = i }; next }
    { gsub(/\r/, "", $column); print $column }' "${manifest}" \
    | tr ',' '\n' | sort -u | while read -r run; do
    [[ -n "${run}" ]] || continue
    singularity exec --cleanenv --no-home --pwd /tmp \
        --bind "${data_root}:${data_root}" \
        docker://naotokubota/ngsfetch:v0.1.2 \
        ngsfetch -i "${run}" -o "${data_root}/${analysis_id}" -p "${SLURM_CPUS_PER_TASK:-16}" --attempts 10
done
