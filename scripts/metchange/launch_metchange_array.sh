#!/bin/bash --login
#SBATCH -N 1
#SBATCH --cpus-per-task=128
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100GB
#SBATCH --time=02-00:00:00
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH -p batch
#SBATCH --qos=normal

set -euo pipefail

module --force purge || true

module load env/development/2024a
module load math/Gurobi/12.0.1-GCCcore-13.3.0

export MAMBA_ROOT_PREFIX="/work/projects/car_t_methylome/software/micromamba"

eval "$(micromamba shell hook --shell bash)"

micromamba activate metabolic_drug_repurposing

# Arguments passed by coordinator
gem_f="$1"
medium_f="$2"
out_dir_expr="$3"
out_dir_metchange="$4"
expr_f="$5"
input_type="$6"

# Get name base
b_name=$(basename "$expr_f" .csv)

# Select the input corresponding to this array task
EXPR="${out_dir_expr}/${b_name}_batch_${SLURM_ARRAY_TASK_ID}.csv"


echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Expression: $EXPR"
echo "GEM: $gem_f"
echo "Medium: $medium_f"

python run_metchange_exe.py \
    --expr "$EXPR" \
    --input_type "$input_type" \
    --gem "$gem_f" \
    --medium "$medium_f" \
    --n_workers 124 \
    --out_dir "$out_dir_metchange"