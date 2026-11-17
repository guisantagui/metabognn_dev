#!/bin/bash -l
#Usage # sbatch launch_benchmark_schuhknecht_A549.sh [run ...] [-- extra benchmark_schuhknecht.py options]
#   Schuhknecht A549 benchmark (benchmark_schuhknecht.py; reports/2026-11-17_signature_selection_scale.md, section 6):
#   queries each fold's models with the measured intracellular metabolome change of each drug, for every run given
#   (default: corrected MetChange _min15, seeds 0 and 1) and its _shuff control. Each run first checks that it
#   reproduces its own held-out LINCS scores from the inputs here (reproduction_check.txt) and stops otherwise.
#   Output: <run dir>/schuhknecht_benchmark[_<tag>]/ (summary.tsv, per_drug.tsv, queries.tsv, scores.npz).
#   Needs data/validation/schuhknecht/ in project_dir (41587_2024_2524_MOESM2_ESM.xlsx, cpd_info_wSchuhknecht.csv,
#   drug_match_table.csv) and openpyxl in the env (first run caches the sheets as pickles in the out dir).
#   Example variants: sbatch launch_benchmark_schuhknecht_A549.sh metchange_rank_noGlobals_lr1e-3_min15 -- \
#       --restrict-nodes all --tag allnodes
#SBATCH --job-name=schuhknecht_bench
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 8
#SBATCH --mem=48GB
#SBATCH --time=0-03:00:00
#SBATCH --output=./schuhknecht_bench_%j.out
#SBATCH --error=./schuhknecht_bench_%j.err
#SBATCH -p batch
#SBATCH --qos=normal

project_dir="/home/users/gsantamaria/projects/metabognn_dev"
project_dir_scratch="/scratch/users/gsantamaria/test_large_files/metabognn_dev"
cell="A549"
runs_dir="${project_dir_scratch}/results/predict_targets/${cell}"
schk_dir="${project_dir}/data/validation/schuhknecht"

runs=()
extra=()
while [[ $# -gt 0 ]]; do
    if [[ "$1" == "--" ]]; then shift; extra=("$@"); break; fi
    runs+=("$1"); shift
done
[[ ${#runs[@]} -gt 0 ]] || runs=(metchange_rank_noGlobals_lr1e-3_min15 metchange_rank_noGlobals_lr1e-3_min15_seed1)

for f in 41587_2024_2524_MOESM2_ESM.xlsx cpd_info_wSchuhknecht.csv drug_match_table.csv; do
    [[ -s "${schk_dir}/${f}" ]] || { echo "missing input: ${schk_dir}/${f}"; exit 1; }
done

conda activate metabolic_drug_repurposing
status=0
for r in "${runs[@]}"; do
    for d in "${runs_dir}/${r}" "${runs_dir}/${r}_shuff"; do
        if [[ ! -d "$d" ]]; then echo "WARNING: ${d} not found, skipped"; status=1; continue; fi
        echo "=== ${d}"
        python benchmark_schuhknecht.py --run-dir "$d" --schuhknecht-dir "$schk_dir" --device cpu "${extra[@]}" \
            || status=1
    done
done
conda deactivate
exit $status
