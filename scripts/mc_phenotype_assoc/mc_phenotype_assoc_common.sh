# Shared configuration for the MetChange / phenotype association launchers.
# Sourced, not executed.
#
# CELL LINE: set `cell` before sourcing (the slurm launchers export it). Every
# path below derives from it, so nothing else needs editing to change cell line.
#
# WHERE THINGS LIVE
#   labels      <scratch_proj_dir>/results/labels/<cell>/{metchange/rank,expr}
#   outputs     <scratch_proj_dir>/results/mc_phenotype_assoc/<cell>/
#   phenotypes  <scratch_proj_dir>/results/mc_phenotype_assoc/<cell>/phenotypes/
#   sets        <scratch_proj_dir>/results/mc_phenotype_assoc/metabolite_sets/
#
# Phenotype tables are PER CELL LINE because they are keyed on LINCS
# perturbation names and built from that cell line's compound mapping - the
# mappable compound set differs between cell lines, so sharing them would
# silently drop compounds. Metabolite / gene sets are per PROJECT: they are
# built from Human-GEM alone, and members absent from a given label matrix are
# dropped at scoring time (mc_phenotype_assoc.py reports how many it found).
#
# Everything is overridable by exporting it first.

repo_dir="${repo_dir:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"

# Project root on scratch when it exists (cluster), the repo otherwise (laptop).
# Keeps large intermediates off the home partition.
# Candidate project roots come from project_paths.conf (one source of truth; it
# used to be hardcoded here, in submit_assoc.sh and in
# check_label_availability.py, and the three drifted apart).
if [ -z "${scratch_proj_dir:-}" ]; then
    _conf="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/project_paths.conf"
    if [ -f "$_conf" ]; then
        while IFS='=' read -r _k _v; do
            case "$_k" in
                scratch_proj_dir)
                    [ -d "$_v" ] && [ -z "${scratch_proj_dir:-}" ] && scratch_proj_dir="$_v" ;;
            esac
        done < <(grep -E '^[a-z_]+=' "$_conf")
    fi
    : "${scratch_proj_dir:=$repo_dir}"
fi
labels_base="${labels_base:-$scratch_proj_dir}"

cell="${cell:?set cell before sourcing, e.g. cell=HEPG2}"

mc_labels="${mc_labels:-$labels_base/results/labels/$cell/metchange/rank}"
expr_labels="${expr_labels:-$labels_base/results/labels/$cell/expr}"

assoc_base="${assoc_base:-$scratch_proj_dir/results/mc_phenotype_assoc}"
out_base="${out_base:-$assoc_base/$cell}"
# Phenotypes: a per-cell directory wins when it exists (built from that cell
# line's own structure mapping); otherwise the shared union tables are used,
# built once with cell=ALL from lincs_unique_perts.tsv. The union tables are
# cell-independent on purpose - the association script keeps only the drugs that
# have labels in the cell line being run - so they are the cheap way to cover
# many cell lines. Which one was used is printed by `provenance`.
if [ -z "${pheno_dir:-}" ]; then
    if [ "$cell" = "ALL" ]; then
        pheno_dir="$assoc_base/phenotypes_ALL"
        pheno_scope="union build target"
    elif [ -f "$assoc_base/$cell/phenotypes/attene_ramos_mmp.tsv" ]; then
        pheno_dir="$assoc_base/$cell/phenotypes"
        pheno_scope="per-cell"
    elif [ -f "$assoc_base/phenotypes_ALL/attene_ramos_mmp.tsv" ]; then
        pheno_dir="$assoc_base/phenotypes_ALL"
        pheno_scope="union (cell=ALL)"
    else
        pheno_dir="$assoc_base/$cell/phenotypes"
        pheno_scope="per-cell (absent)"
    fi
else
    pheno_scope="explicit"
fi
sets_dir="${sets_dir:-$assoc_base/metabolite_sets}"

met_sets="${met_sets:-$sets_dir/metabolite_sets.tsv}"
gene_sets="${gene_sets:-$sets_dir/gene_sets.tsv}"
negctrl_sets="${negctrl_sets:-$sets_dir/metabolite_sets_negctrl.tsv}"

py="${py:-python3}"
here_dir="$repo_dir/scripts/mc_phenotype_assoc"
assoc="$here_dir/mc_phenotype_assoc.py"
strat="$here_dir/stratified_check.py"
matched="$here_dir/matched_check.py"

n_perm="${n_perm:-5000}"
n_boot="${n_boot:-2000}"
n_random_sets="${n_random_sets:-200}"
qvalue="${qvalue:-0.05}"

say() { printf '\n=== %s\n' "$*"; }

provenance() {
    cat <<PROV

--------------------------------------------------------------------
cell line        : $cell
project root     : $scratch_proj_dir
metchange labels : $mc_labels
expression labels: $expr_labels
phenotypes       : $pheno_dir  [$pheno_scope]
metabolite sets  : $sets_dir
outputs          : $out_base
resampling       : n_perm=$n_perm n_boot=$n_boot n_random_sets=$n_random_sets
--------------------------------------------------------------------
PROV
}

# Fail loudly rather than silently reusing another cell line's inputs: a missing
# phenotype table is the failure mode that yields a plausible-looking wrong
# answer instead of an error.
require_inputs() {
    local missing=0
    if [ ! -f "$met_sets" ]; then
        echo "MISSING metabolite sets: $met_sets" >&2; missing=1
    fi
    for f in "$pheno_dir/attene_ramos_mmp.tsv" "$pheno_dir/hallinger_ocr.tsv"; do
        if [ ! -f "$f" ]; then
            echo "MISSING phenotype table: $f" >&2; missing=1
        fi
    done
    if [ "$missing" -ne 0 ]; then
        cat <<MSG >&2

Build the inputs for cell=$cell first:

    cell=$cell bash $here_dir/prepare_assoc_inputs.sh

That writes the project-level metabolite / gene sets and the per-cell phenotype
tables. The phenotype tables need a structure mapping of the external datasets
onto THIS cell line's LINCS perturbations (map_smiles.R, needs R + rcdk); the
prepare script says so explicitly if the mapping is absent.
MSG
        return 1
    fi
    return 0
}

# run_assoc <out-name> <labels-dir> <sets-file> <phenotype-file> <case> <control> [extra args...]
run_assoc() {
    local name="$1" labels="$2" sets="$3" pheno="$4" case_="$5" ctrl="$6"; shift 6
    local out="$out_base/$name"
    say "$name  ($case_ vs $ctrl, cell=$cell)"
    if [ ! -d "$labels" ]; then
        echo "SKIP $name: no label dir $labels" >&2
        return 0
    fi
    if [ ! -f "$pheno" ]; then
        echo "SKIP $name: no phenotype table $pheno" >&2
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
    local name="$1" case_="$2" spec="${3:-}"
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

# run_matched <out-name> <case> <statistic>
run_matched() {
    local name="$1" case_="$2" st="$3"
    local out="$out_base/$name"
    [ -f "$out/per_drug_scores.tsv" ] || { echo "SKIP matched $name: no scores" >&2; return 0; }
    say "responsiveness-matched: $name / $case_ / $st"
    "$py" -I "$matched" --scores "$out/per_drug_scores.tsv" \
        --case "$case_" --statistic "$st" \
        --out "$out/matched_${case_}_${st}.tsv"
}
