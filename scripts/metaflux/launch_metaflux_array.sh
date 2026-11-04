#!/bin/bash --login
#SBATCH -N 1
#SBATCH --cpus-per-task=64
#SBATCH --ntasks-per-node=1
#SBATCH --mem=64GB
#SBATCH --time=00-04:00:00
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH -p batch
#SBATCH --qos=normal

set -euo pipefail

module --force purge || true

module load env/development/2024a

export MAMBA_ROOT_PREFIX="/work/projects/car_t_methylome/software/micromamba"

eval "$(micromamba shell hook --shell bash)"

micromamba activate metabolic_drug_repurposing

# Arguments passed by launch_metaflux_bigjob.sh
gem_f="$1"
medium_f="$2"
out_dir_expr="$3"
out_dir_metaflux="$4"
expr_f="$5"
mras_floor="$6"
status_dir="$7"

b_name=$(basename "$expr_f" .csv)
EXPR="${out_dir_expr}/${b_name}_batch_${SLURM_ARRAY_TASK_ID}.csv"

echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Expression: $EXPR"

# ~0.7 s per sample and core with clarabel
python run_metaflux.py \
    --expr "$EXPR" \
    --expr_scale linear \
    --gem "$gem_f" \
    --medium "$medium_f" \
    --solver clarabel \
    --mras_floor "$mras_floor" \
    --n_workers "$SLURM_CPUS_PER_TASK" \
    --status_dir "$status_dir" \
    --out_dir "$out_dir_metaflux"
