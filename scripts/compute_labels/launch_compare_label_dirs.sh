#!/bin/bash --login
#Usage # sbatch launch_compare_label_dirs.sh [cell] [penalty] [max_perts]
# Compares the sign-corrected labels (home, fix_limma_reference_sign.py) with the rerun of the fixed limma_DE.R
# (scratch), for expression and MetChange. max_perts > 0 compares a random subset first (quick pass).
#SBATCH --job-name="cmp_labels"
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32GB
#SBATCH -c 1
#SBATCH --time=00-03:00:00
#SBATCH --output="logs/cmp_labels.out"
#SBATCH --error="logs/cmp_labels.err"
#SBATCH -p batch
#SBATCH --qos=normal

cell="${1:-A549}"
penalty="${2:-rank}"
max_perts="${3:-0}"

project_dir="/home/users/gsantamaria/projects/metabognn_dev"
project_dir_scratch="/scratch/users/gsantamaria/test_large_files/metabognn_dev"
lab_home="${project_dir}/results/labels/${cell}"
lab_scratch="${project_dir_scratch}/results/labels/${cell}"
out_root="${project_dir_scratch}/results/labels/${cell}/compare_corrected_vs_rerun"

export MAMBA_ROOT_PREFIX="/work/projects/car_t_methylome/software/micromamba"
eval "$(micromamba shell hook --shell bash)"
micromamba activate metabolic_drug_repurposing

# A = corrected (what the MetChange runs up to d11744b read), B = rerun (what launchers read now)
for pair in "expr:expr" "mc:metchange/${penalty}"; do
    tag="${pair%%:*}"; sub="${pair#*:}"
    echo "===== ${tag}: ${lab_home}/${sub}  vs  ${lab_scratch}/${sub}"
    ls -la "${lab_home}/${sub}" "${lab_scratch}/${sub}"
    # provenance: fix_limma_reference_sign.py copies say "sign-fixed"; a limma_DE.R rerun writes its own line
    ref_a="$(cat "${lab_home}/${sub}/reference_level.txt" 2>/dev/null || echo '<none>')"
    ref_b="$(cat "${lab_scratch}/${sub}/reference_level.txt" 2>/dev/null || echo '<none>')"
    echo "A reference_level.txt: ${ref_a}"
    echo "B reference_level.txt: ${ref_b}"
    if [[ "$ref_a" == "$ref_b" && "$ref_a" != "<none>" ]]; then
        echo "WARNING ${tag}: A and B have the same provenance; this is not corrected vs rerun. Skipped."
        continue
    fi
    python compare_label_dirs.py --a "${lab_home}/${sub}" --b "${lab_scratch}/${sub}" \
        --out "${out_root}/${tag}" --max-perts "$max_perts" \
        || echo "compare_label_dirs.py failed for ${tag}"
done

micromamba deactivate
