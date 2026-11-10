#!/bin/bash --login
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=get_lincs_master
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16GB
#SBATCH -c 1
#SBATCH --time=00-01:00:00
#Define sdout path
#SBATCH --output=./logs/get_lincs_master.out
#Define sderr path
#SBATCH --error=./logs/get_lincs_master.err
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal

mkdir -p "logs/"

cell_file="/home/users/gsantamaria/projects/metabognn_dev/data/LINCS/core_cell_lines.tsv"
project_dir="/home/users/gsantamaria/projects/metabognn_dev/"
project_dir_scratch="/scratch/users/gsantamaria/test_large_files/metabognn_dev/"
gem_f="${project_dir}/data/GEMs/Human-GEM/model/Human-GEM.xml"
out_dir="${project_dir_scratch}/results/metChange/"

max_fluxes_file="${out_dir}/lp1_fluxes.npy"

if [ ! -f "$max_fluxes_file" ]; then
    python run_metchange_exe_lp1.py --medium "$medium_f" \
        --gem "$gem_f" \
        --out_dir "$out_dir"
else
    echo "${max_fluxes_file} already exists."
fi

while IFS=$'\t' read -r cell; do
    echo "Launching ${cell} extraction..."
    sbatch \
        --job-name="mc_master_${cell}_rankAll" \
        --output="logs/mc_master_${cell}_rankAll.out" \
        --error="logs/mc_master_${cell}_rankAll.err" \
        launch_metchange_bigjob_generic.sh.sh \
        "$cell"
    fi
done < "$cell_file"