#!/bin/bash -l
#Usage # sbatch launch_check_metabolic_coverage_A549.sh [runs_dir] [out_dir]
#   Full-data metabolic coverage check (report 03, 10.5): expression DE coverage and MetChange's own significant
#   features per active perturbation, and each run's spec_prob_timedose gap over its control split by them.
#   runs_dir: where the _min15 runs are (default: home results, where runs made before d11744b live)
#SBATCH --job-name=mc_coverage
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 1
#SBATCH --mem=32GB
#SBATCH --time=0-01:30:00
#SBATCH --output=./mc_coverage_%j.out
#SBATCH --error=./mc_coverage_%j.err
#SBATCH -p batch
#SBATCH --qos=normal

project_dir="/home/users/gsantamaria/projects/metabognn_dev"
project_dir_scratch="/scratch/users/gsantamaria/test_large_files/metabognn_dev"
cell="A549"

runs_dir="${1:-${project_dir}/results/predict_targets/${cell}}"
out_dir="${2:-${runs_dir}/metabolic_coverage_full}"

# labels as the runs read them: expr always from scratch; MetChange runs up to d11744b from home (q-values are
# unaffected by the sign correction, so either copy gives the same counts)
expr_qvals="${project_dir_scratch}/results/labels/${cell}/expr/qvals.csv"
mc_dir="${project_dir}/results/labels/${cell}/metchange/rank"
gem_genes="${project_dir}/data/GEMs/Human-GEM/model/genes.tsv"
mc_nodes="${project_dir}/data/networks/combined/nodes.tsv"

for f in "$expr_qvals" "${mc_dir}/qvals.csv" "${mc_dir}/log2FC.csv" "$gem_genes" "$mc_nodes"; do
    [[ -s "$f" ]] || { echo "missing input: $f"; exit 1; }
done

# real runs (name=dir); each needs <dir>_shuff. Missing pairs are skipped with a warning.
runs=()
for pair in \
    "mc=metchange_rank_noGlobals_lr1e-3_min15" \
    "expr=expr_noGlobals_lr1e-3_min15" \
    "mc_s1=metchange_rank_noGlobals_lr1e-3_min15_seed1" \
    "expr_s1=expr_noGlobals_lr1e-3_min15_seed1"; do
    name="${pair%%=*}"; d="${runs_dir}/${pair#*=}"
    if [[ -s "${d}/cv_drug_ranks.tsv" && -s "${d}_shuff/cv_drug_ranks.tsv" ]]; then
        runs+=(--run "${name}=${d}")
    else
        echo "WARNING: skipping ${name}: ${d} or ${d}_shuff has no cv_drug_ranks.tsv"
    fi
done
[[ ${#runs[@]} -gt 0 ]] || { echo "no runs found in ${runs_dir}"; exit 1; }

mkdir -p "$out_dir"
echo "runs: ${runs[*]}"
echo "out:  ${out_dir}"

conda activate metabolic_drug_repurposing
python check_metabolic_coverage.py \
    --qvals "$expr_qvals" \
    --gem-genes "$gem_genes" \
    --mc-qvals "${mc_dir}/qvals.csv" \
    --mc-log2fc "${mc_dir}/log2FC.csv" \
    --mc-nodes "$mc_nodes" \
    "${runs[@]}" \
    --out-dir "$out_dir"
status=$?
conda deactivate
exit $status
