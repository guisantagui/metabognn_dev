#!/bin/bash --login
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=metpath_A549
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=64GB
#SBATCH -c 64
#SBATCH --time=00-12:00:00
#Define sdout path
#SBATCH --output=metpath_A549.out
#Define sderr path
#SBATCH --error=metpath_A549.err
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal

# MetPath for one cell line: (1) production / consumption pathways from one
# reference flux state, once; (2) scoring of the expression log2 fold changes
# (limma_DE.R on the expression) on them. The environment needs clarabel
# (pip install clarabel) for the reference flux state.

set -euo pipefail

module --force purge || true

module load env/development/2024a

export MAMBA_ROOT_PREFIX="/work/projects/car_t_methylome/software/micromamba"

eval "$(micromamba shell hook --shell bash)"

micromamba activate metabolic_drug_repurposing

project_dir="/home/users/gsantamaria/projects/metabognn_dev"
cell="A549"
gem_f="/home/users/gsantamaria/projects/metabo_pert/data/GEMs/Human2/Human-GEM/model/Human-GEM.xml"
medium_f="/home/users/gsantamaria/projects/metabo_pert/data/GEMs/gem_mediums/Human2-GEM/rpmi_medium.json"
hgnc_f="${project_dir}/data/utilities/hgnc_complete_set.txt"
baseline_f="/scratch/users/gsantamaria/test_large_files/data/perturbation/lincs/cell_matrices/${cell}/baseline_ctrl.tsv"
de_f="${project_dir}/results/labels/${cell}/expr/log2FC.csv"
# metaflux: METAFlux on the controls' mean expression (cell-line specific, ~2 CPU-hours of pathways);
# fba: maximal biomass with minimal norm, as in MetPath (not cell-specific, ~45 CPU-hours).
flux_state="metaflux"
mras_floor=0.25
distance=2
cutoff_fraction=0.05
fc_scale="ratio"     # ratio: MetPath's scores (1 = unchanged); log2: mean log2FC (0 = unchanged)
out_dir="${project_dir}/results/metpath/${cell}/${flux_state}_d${distance}/"
paths_dir="${out_dir}/pathways"
scores_dir="${out_dir}/scores_${fc_scale}"

if [ ! -f "${paths_dir}/metpath_pathways.tsv" ]; then
    python build_metpath_pathways.py \
        --gem "$gem_f" \
        --medium "$medium_f" \
        --flux_state "$flux_state" \
        --baseline "$baseline_f" \
        --hgnc "$hgnc_f" \
        --mras_floor "$mras_floor" \
        --distance "$distance" \
        --cutoff_fraction "$cutoff_fraction" \
        --n_workers "$SLURM_CPUS_PER_TASK" \
        --out_dir "$paths_dir"
else
    echo "${paths_dir}/metpath_pathways.tsv already exists."
fi

python score_metpath.py \
    --pathways "${paths_dir}/metpath_pathways.tsv" \
    --gem "$gem_f" \
    --de "$de_f" \
    --hgnc "$hgnc_f" \
    --fc_scale "$fc_scale" \
    --n_perms 1000 \
    --out_dir "$scores_dir"

micromamba deactivate
