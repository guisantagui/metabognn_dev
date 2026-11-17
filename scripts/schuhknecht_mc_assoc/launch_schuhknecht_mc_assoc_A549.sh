#!/bin/bash -l
#Usage # sbatch launch_schuhknecht_mc_assoc_A549.sh [compartments ...]   (from scripts/schuhknecht_mc_assoc)
#   MetChange vs Schuhknecht association on the limma labels, for each compartment setting given (default: c, the
#   earlier analysis, and intracellular = all compartments but e, copies averaged):
#     1. prepare_assoc_inputs.py   matched matrices (MetChange log2FC / q, Schuhknecht log2FC / p / q, expression
#                                  producers baseline from the expression log2FC, drug match)
#     2. overlap_enrichment.py     MetChange-significant vs Schuhknecht-significant metabolites per drug, against a
#                                  mismatched-drug null and a size-matched expression baseline
#     3. metchange_vs_expression_benchmark.R + sign_heterogeneity_analysis.R (needs R with data.table, tidyr, tibble,
#        dplyr, lme4, performance, ggplot2; skipped with a message if Rscript is not available: the folder from step
#        1 is small, rsync it and run step 3 locally with the command printed at the end)
#   Output: <scratch>/results/mc_schuhk_assoc/limma_<penalty>_<compartments>/ (overlap_enrichment/, benchmark/)
#SBATCH --job-name=mc_schuhk_assoc
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 4
#SBATCH --mem=32GB
#SBATCH --time=0-03:00:00
#SBATCH --output=./mc_schuhk_assoc_%j.out
#SBATCH --error=./mc_schuhk_assoc_%j.err
#SBATCH -p batch
#SBATCH --qos=normal

project_dir="/home/users/gsantamaria/projects/metabognn_dev"
project_dir_scratch="/scratch/users/gsantamaria/test_large_files/metabognn_dev"
cell="A549"
penalty="rank"
mc_dir="${project_dir_scratch}/results/labels/${cell}/metchange/${penalty}"    # corrected labels (scratch)
expr_dir="${project_dir_scratch}/results/labels/${cell}/expr"
schk_dir="${project_dir}/data/validation/schuhknecht"
net_dir="${project_dir}/data/networks/combined"
out_root="${project_dir_scratch}/results/mc_schuhk_assoc"
export LINCS_COMPOUNDINFO="${project_dir}/data/LINCS/compoundinfo_beta.txt"

settings=("$@")
[[ ${#settings[@]} -gt 0 ]] || settings=(c intracellular)

for f in "${mc_dir}/log2FC.csv" "${mc_dir}/qvals.csv" "${expr_dir}/log2FC.csv" "${net_dir}/nodes.tsv" \
         "${schk_dir}/41587_2024_2524_MOESM2_ESM.xlsx" "${schk_dir}/cpd_info_wSchuhknecht.csv" \
         "${schk_dir}/drug_match_table.csv"; do
    [[ -s "$f" ]] || { echo "missing input: $f"; exit 1; }
done

conda activate metabolic_drug_repurposing
have_r=0
command -v Rscript >/dev/null 2>&1 || module load lang/R >/dev/null 2>&1
command -v Rscript >/dev/null 2>&1 && Rscript -e 'for (p in c("data.table","tidyr","tibble","dplyr","lme4","performance","ggplot2")) if (!requireNamespace(p, quietly = TRUE)) quit(status = 1)' \
    && have_r=1
status=0
for comp in "${settings[@]}"; do
    out="${out_root}/limma_${penalty}_${comp//,/_}"
    echo "=== compartments ${comp} -> ${out}"
    python prepare_assoc_inputs.py --mc-dir "$mc_dir" --expr-dir "$expr_dir" --schuhknecht-dir "$schk_dir" \
        --network-dir "$net_dir" --compartments "$comp" --out-dir "$out" || { status=1; continue; }
    python overlap_enrichment.py --in-dir "$out" || status=1
    if [[ $have_r -eq 1 ]]; then
        Rscript -e 'source("metchange_vs_expression_benchmark.R"); source("sign_heterogeneity_analysis.R")' \
            "$out" "${out}/benchmark" > "${out}/benchmark.log" 2>&1 || { status=1; echo "R step failed: ${out}/benchmark.log"; }
    else
        echo "R (with the packages) not available here: rsync ${out} and run, from scripts/schuhknecht_mc_assoc:"
        echo "  Rscript -e 'source(\"metchange_vs_expression_benchmark.R\"); source(\"sign_heterogeneity_analysis.R\")' <dir> <dir>/benchmark > <dir>/benchmark.log 2>&1"
    fi
done
conda deactivate
exit $status
