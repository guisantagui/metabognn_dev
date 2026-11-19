# Shared configuration for the MetChange / phenotype association launchers.
# Sourced, not executed.
#
# Labels live on scratch on the cluster and under the repo on the laptop; the
# directory structure below the base is identical, so only the base differs.
# Override anything by exporting it before calling a launcher.

repo_dir="${repo_dir:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"

# base holding results/labels/<cell>/...
if [ -d /scratch/users/gsantamaria/test_large_files/metabognn_dev ]; then
    labels_base="${labels_base:-/scratch/users/gsantamaria/test_large_files/metabognn_dev}"
else
    labels_base="${labels_base:-$repo_dir}"
fi

cell="${cell:-HEPG2}"
mc_labels="${mc_labels:-$labels_base/results/labels/$cell/metchange/rank}"
expr_labels="${expr_labels:-$labels_base/results/labels/$cell/expr}"

# outputs stay with the rest of the results tree (small: TSVs + one npz)
out_base="${out_base:-$repo_dir/results/mc_phenotype_assoc/$cell}"

sets_dir="${sets_dir:-$repo_dir/data/validation/metabolite_sets}"
met_sets="${met_sets:-$sets_dir/metabolite_sets.tsv}"
gene_sets="${gene_sets:-$sets_dir/gene_sets.tsv}"
pheno_dir="${pheno_dir:-$repo_dir/data/validation/phenotypes}"

py="${py:-python3}"
assoc="$repo_dir/scripts/mc_phenotype_assoc/mc_phenotype_assoc.py"
strat="$repo_dir/scripts/mc_phenotype_assoc/stratified_check.py"

n_perm="${n_perm:-5000}"
n_boot="${n_boot:-2000}"
n_random_sets="${n_random_sets:-200}"
qvalue="${qvalue:-0.05}"

say() { printf '\n=== %s\n' "$*"; }

# run_assoc <out-name> <labels-dir> <sets-file> <phenotype-file> <case> <control> [extra args...]
run_assoc() {
    local name="$1" labels="$2" sets="$3" pheno="$4" case_="$5" ctrl="$6"; shift 6
    local out="$out_base/$name"
    say "$name  ($case_ vs $ctrl)"
    if [ ! -d "$labels" ]; then
        echo "SKIP $name: no label dir $labels" >&2
        return 0
    fi
    mkdir -p "$out"
    "$py" -I "$assoc" \
        --labels "$labels" \
        --sets "$sets" \
        --phenotype "$pheno" \
        --case "$case_" --control "$ctrl" \
        --qvalue "$qvalue" \
        --n-perm "$n_perm" --n-boot "$n_boot" --n-random-sets "$n_random_sets" \
        --out-dir "$out" "$@" 2>&1 | grep -vE '^  [0-9]+ rows \('
}

# run_strat <out-name> <case> [sets spec]
run_strat() {
    local name="$1" case_="$2" spec="$3"
    local out="$out_base/$name"
    [ -f "$out/per_drug_scores.tsv" ] || { echo "SKIP strat $name: no scores" >&2; return 0; }
    for st in $(head -1 "$out/per_drug_scores.tsv" | tr '\t' '\n' | tail -n +8); do
        say "stratified: $name / $case_ / $st"
        "$py" -I "$strat" --scores "$out/per_drug_scores.tsv" \
            --case "$case_" --statistic "$st" \
            ${spec:+--sets "$spec"} \
            --out "$out/stratified_${case_}_${st}.tsv"
    done
}
