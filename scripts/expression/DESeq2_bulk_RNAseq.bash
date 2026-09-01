#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --job-name="DESeq2 bulk"
#SBATCH -p epyc

set -euo pipefail
analysis_id=${1:?Usage: sbatch DESeq2_bulk_RNAseq.bash ANALYSIS_ID DESIGN CONTRAST [SUBSET] [NAME]}
design=${2:?DESIGN is required, e.g. '~ sex + genotype'}
contrast=${3:?CONTRAST is required, e.g. 'genotype,APOE4,APOE3'}
subset=${4:-}
name=${5:-$(printf '%s' "${contrast}" | tr ', /' '___')}
project_root=${PROJECT_ROOT:-/rhome/naotok/neurodegeneration}
data_root=${DATA_ROOT:-/rhome/naotok/bigdata/neurodegeneration/GEO_SRA}
workdir="${data_root}/${analysis_id}"
module load singularity
singularity exec --bind /rhome/naotok:/rhome/naotok \
    docker://quay.io/biocontainers/bioconductor-deseq2:1.42.0--r43hf17093f_0 \
    Rscript "${project_root}/scripts/expression/differential_expression.R" \
    "${workdir}/expression/featureCounts.txt" \
    "${project_root}/metadata/analysis/${analysis_id}.tsv" \
    "${workdir}/expression/DESeq2/${name}" "${design}" "${contrast}" "${subset}"
