#!/bin/bash
#SBATCH --job-name=mc_pheno_HEPG2
#SBATCH --array=0-1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs/mc_pheno_HEPG2_%A_%a.out
# MetChange / phenotype association for HEPG2.
#
#   sbatch scripts/mc_phenotype_assoc/launch_assoc_HEPG2_slurm.sh
#
# Inputs must exist first (once per cell line):
#   cell=HEPG2 bash scripts/mc_phenotype_assoc/prepare_assoc_inputs.sh
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
export cell=HEPG2
exec bash scripts/mc_phenotype_assoc/launch_assoc_slurm_generic.sh
