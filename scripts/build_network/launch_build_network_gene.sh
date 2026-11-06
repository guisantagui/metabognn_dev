#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=build_network_gene
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH -c 2
#SBATCH --mem=32GB
#SBATCH --time=0-02:00:00
#Define sdout path
#SBATCH --output=./build_network_gene.out
#Define sderr path
#SBATCH --error=./build_network_gene.err
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal

# Ablation: proteins only, transcriptional edges included: drug target -> ... -> TF -> any gene, no GEM layer.
# Same sources, evidence thresholds and curation in the three networks (launch_build_network_{metabolic,tf,gene}.sh);
# only the layers kept and the pruning differ. See build_network.py --help.

conda activate metabolic_drug_repurposing

project_dir="/home/users/gsantamaria/projects/metabognn_dev"
data_dir="${project_dir}/data"
out_dir="${data_dir}/networks/gene"
drug_targets="${data_dir}/drug_targets/combined/drug_targets.tsv"

python build_network.py --data-dir "$data_dir" \
        --out-dir "$out_dir" \
        --network gene \
        --drug-targets "$drug_targets"

conda deactivate
