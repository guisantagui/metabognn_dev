#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=200GB
#SBATCH -c 128
#SBATCH --time=02-00:00:00
#Define sdout path
#Define sderr path
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal

# Results are checkpointed in <out_dir>/checkpoints: if the job hits the time
# limit, resubmit it as is and it continues where it stopped.

# Script arguments
################################################################################

project_dir="$1"
cell="$2"
workers="$3"
mode="$4"
penalty="$5"
medium_f="$6"

# Validate arguments
if [[ $# -lt 4 ]]; then
    echo "Usage: sbatch $0 project_dir cell workers mode [penalty] [medium_f]"
    exit 1
fi

if [[ "$mode" == "mc" && ( -z "$penalty" || -z "$medium_f" ) ]]; then
    echo "ERROR: mode 'mc' requires penalty and medium_f."
    exit 1
fi

conda activate r-4.3.1
# One BLAS thread per worker: the workers already use all the cores
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

# Get paths
################################################################################

#mode="mc"
#penalty="rank"
#medium_f="/home/users/gsantamaria/projects/metabo_pert/data/GEMs/gem_mediums/Human2-GEM/rpmi_medium.json"


echo "Running limma of ${cell}."

mc_dir="${project_dir}/results/metChange/"
cell_mat_dir="${project_dir}/data/perturbation/lincs/cell_matrices/"
mdat="${cell_mat_dir}/${cell}/mdat.csv"
out_dir="${project_dir}/results/labels/"
mkdir -p "$out_dir"

if [ "$mode" == "mc" ]; then
    medium=$(basename "${medium_f%.json}")
    expr="${mc_dir}/${cell}/${penalty}_concat/expr_ensembl_${penalty}_metChange_${medium}_concat.csv"
    out_dir_cell="${out_dir}/${cell}/metchange/${penalty}/"
    reader="lines"   # fread segfaults on the MetChange matrix (~150k columns)
elif [ "$mode" == "expr" ]; then
    expr="${cell_mat_dir}/${cell}/expr_ensembl.csv"
    out_dir_cell="${out_dir}/${cell}/expr/"
    reader="fread"
else
    echo "ERROR: Unknown mode '$mode'. Use 'mc' or 'expr'."
    exit 1
fi

mkdir -p "$out_dir_cell"

# Run limma
################################################################################

if [[ -f "$mdat" && -f "$expr" ]]; then
    Rscript limma_DE.R --expr "$expr" \
        --mdat "$mdat" \
        --workers "$workers" \
        --reader "$reader" \
        --out_dir "$out_dir_cell" \
        || { echo "limma failed or stopped for ${cell}: resubmit to resume from checkpoints"; exit 1; }
else
    echo "ERROR: either '${mdat}' or '${expr}' don't exist."
    exit 1
fi
conda deactivate

export MAMBA_ROOT_PREFIX="/work/projects/car_t_methylome/software/micromamba"
eval "$(micromamba shell hook --shell bash)"


# Generate report file of number of DE genes per perturbation
################################################################################
micromamba activate metabolic_drug_repurposing

python check_limma_nsign_vs_nsamps.py --dir "$out_dir_cell"

micromamba deactivate