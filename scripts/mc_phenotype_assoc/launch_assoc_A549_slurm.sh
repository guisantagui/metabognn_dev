#!/bin/bash
#SBATCH --job-name=mc_pheno_A549
#SBATCH --array=0-1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs/mc_pheno_A549_%A_%a.out
# MetChange / phenotype association for A549.
#
#   sbatch scripts/mc_phenotype_assoc/launch_assoc_A549_slurm.sh
#
# Inputs must exist first (once per cell line):
#   cell=A549 bash scripts/mc_phenotype_assoc/prepare_assoc_inputs.sh
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
export cell=A549
exec bash scripts/mc_phenotype_assoc/launch_assoc_slurm_generic.sh
