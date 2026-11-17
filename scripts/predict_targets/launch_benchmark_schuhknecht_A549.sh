#!/bin/bash -l
#Usage # sbatch launch_benchmark_schuhknecht_A549.sh [run ...] [-- extra benchmark_schuhknecht.py options]
#   Schuhknecht A549 benchmark (benchmark_schuhknecht.py; reports/2026-11-17_signature_selection_scale.md, section 6):
#   queries each fold's models with the measured intracellular metabolome change of each drug, for every run given
#   (default: corrected MetChange _min15, seeds 0 and 1) and its _shuff control. The first real run and the first
#   control check that they reproduce their own held-out LINCS scores from the inputs here (reproduction_check.txt;
#   ~7 min each, mostly reading the labels) and stop otherwise; later runs share those inputs and skip it
#   (CHECK_ALL=1 sbatch ... to check every run). Add --only-shuff / --only-real before the runs to do one side only.
#   Output: <run dir>/schuhknecht_benchmark[_<tag>]/ (summary.tsv, per_drug.tsv, queries.tsv, scores.npz).
#   Needs data/validation/schuhknecht/ in project_dir (41587_2024_2524_MOESM2_ESM.xlsx, cpd_info_wSchuhknecht.csv,
#   drug_match_table.csv) and openpyxl in the env (the sheets are read once into <runs_dir>/schuhknecht_cache).
#   GPU job: ~1-2 min of scoring per run (CPU batch nodes took 20-45 min, thread oversubscription).
#   Example variants: sbatch launch_benchmark_schuhknecht_A549.sh metchange_rank_noGlobals_lr1e-3_min15 -- \
#       --restrict-nodes all --tag allnodes
#SBATCH --job-name=schuhknecht_bench
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 8
#SBATCH --mem=48GB
#SBATCH --time=0-02:00:00
#SBATCH --output=./schuhknecht_bench_%j.out
#SBATCH --error=./schuhknecht_bench_%j.err
#SBATCH -p gpu
#SBATCH -G 1
#SBATCH --qos=normal

project_dir="/home/users/gsantamaria/projects/metabognn_dev"
project_dir_scratch="/scratch/users/gsantamaria/test_large_files/metabognn_dev"
cell="A549"
runs_dir="${project_dir_scratch}/results/predict_targets/${cell}"
schk_dir="${project_dir}/data/validation/schuhknecht"

runs=()
extra=()
sides=(real shuff)
while [[ $# -gt 0 ]]; do
    case "$1" in
        --) shift; extra=("$@"); break ;;
        --only-shuff) sides=(shuff); shift ;;
        --only-real) sides=(real); shift ;;
        *) runs+=("$1"); shift ;;
    esac
done
[[ ${#runs[@]} -gt 0 ]] || runs=(metchange_rank_noGlobals_lr1e-3_min15 metchange_rank_noGlobals_lr1e-3_min15_seed1)

for f in 41587_2024_2524_MOESM2_ESM.xlsx cpd_info_wSchuhknecht.csv drug_match_table.csv; do
    [[ -s "${schk_dir}/${f}" ]] || { echo "missing input: ${schk_dir}/${f}"; exit 1; }
done

conda activate metabolic_drug_repurposing
# torch / BLAS threads = allocated cores (torch's default counts every core of the node)
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8} MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
device=$(python -c "import torch; print('cuda' if torch.cuda.is_available() else 'cpu')")
echo "device: ${device}"
cache="${runs_dir}/schuhknecht_cache"
status=0
checked_real=0
checked_shuff=0
for r in "${runs[@]}"; do
    for side in "${sides[@]}"; do
        if [[ $side == real ]]; then d="${runs_dir}/${r}"; else d="${runs_dir}/${r}_shuff"; fi
        if [[ ! -d "$d" ]]; then echo "WARNING: ${d} not found, skipped"; status=1; continue; fi
        check=()
        if [[ "${CHECK_ALL:-0}" != 1 ]]; then
            if [[ $side == real ]]; then
                [[ $checked_real == 1 ]] && check=(--check-n 0); checked_real=1
            else
                [[ $checked_shuff == 1 ]] && check=(--check-n 0); checked_shuff=1
            fi
        fi
        echo "=== ${d}"
        python benchmark_schuhknecht.py --run-dir "$d" --schuhknecht-dir "$schk_dir" --device "$device" \
            --cache-dir "$cache" "${check[@]}" "${extra[@]}" || status=1
    done
done
conda deactivate
exit $status
