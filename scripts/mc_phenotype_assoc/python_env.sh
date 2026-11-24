# Resolve a Python interpreter that actually has the project's dependencies.
# Sourced, not executed. Sets $py and defines require_python.
#
# WHY THIS EXISTS. On the cluster the dependencies live in the micromamba
# environment `metabolic_drug_repurposing`. A batch job does not inherit an
# activated environment, so plain `python3` is the system interpreter and
# mc_phenotype_assoc.py dies with "No module named numpy" - after Slurm has
# granted the allocation, and for the array tasks one at a time. On the laptop the
# system python3 already has numpy and scipy, so activation has to be OPTIONAL
# rather than assumed.
#
# RESOLUTION ORDER, first that works wins:
#   1. an exported $py, if it can import numpy and scipy (explicit override);
#   2. the current python3, if it can (laptop, or an already-activated shell);
#   3. the micromamba environment $mamba_env - its interpreter located by asking
#      micromamba rather than guessing a layout, else via the shell hook;
#   4. $mamba_env_prefix/bin/python3, or the same under common micromamba roots.
#
# Nothing here aborts on failure, because this file is also sourced by scripts
# that only print paths. `require_python` is the loud check: call it before doing
# real work, and it refuses with instructions rather than starting a job that will
# crash partway through.
#
# Override with: py=/path/to/python3, mamba_env=<name>, or mamba_env_prefix=<dir>
# (the last two can also be set in project_paths.conf).

mamba_env="${mamba_env:-metabolic_drug_repurposing}"

# Keep this list in step with what the scripts actually import. mc_phenotype_assoc.py
# needs numpy and scipy; the builders additionally need openpyxl or rdkit, which are
# checked by those scripts themselves because they fail immediately and clearly.
_py_modules="${_py_modules:-numpy scipy}"

_py_usable() {
    [ -n "${1:-}" ] || return 1
    command -v "$1" >/dev/null 2>&1 || [ -x "$1" ] || return 1
    "$1" -I -c "import $(echo "$_py_modules" | tr ' ' ',')" >/dev/null 2>&1
}

_resolve_python() {
    # 1. explicit override
    if [ -n "${py:-}" ] && _py_usable "$py"; then
        py_source="exported py"
        return 0
    fi
    # 2. whatever python3 is already on PATH
    if _py_usable python3; then
        py="$(command -v python3)"
        py_source="python3 on PATH"
        return 0
    fi
    # 3. micromamba: ask it for the env's interpreter rather than guessing a layout.
    #    `micromamba run` is used only to LOCATE the binary; the binary is then called
    #    directly, so no activation state is needed at run time.
    if command -v micromamba >/dev/null 2>&1; then
        local cand hook
        cand="$(micromamba run -n "$mamba_env" python3 -c 'import sys; print(sys.executable)' 2>/dev/null \
                | tail -1)"
        if _py_usable "$cand"; then
            py="$cand"
            py_source="micromamba env $mamba_env"
            return 0
        fi
        # activation, in case the env needs its LD_LIBRARY_PATH as well as its binary
        if hook="$(micromamba shell hook --shell bash 2>/dev/null)"; then
            eval "$hook" 2>/dev/null || true
            if micromamba activate "$mamba_env" 2>/dev/null && _py_usable python3; then
                py="$(command -v python3)"
                py_source="micromamba activate $mamba_env"
                return 0
            fi
        fi
    fi
    # 4. a prefix given explicitly, or the usual micromamba / conda roots
    local p
    for p in "${mamba_env_prefix:-}" \
             "${MAMBA_ROOT_PREFIX:-}/envs/$mamba_env" \
             "${CONDA_PREFIX:-}/../$mamba_env" \
             "$HOME/micromamba/envs/$mamba_env" \
             "$HOME/.micromamba/envs/$mamba_env" \
             "$HOME/miniforge3/envs/$mamba_env" \
             "$HOME/mambaforge/envs/$mamba_env"; do
        [ -n "$p" ] || continue
        if _py_usable "$p/bin/python3"; then
            py="$p/bin/python3"
            py_source="prefix $p"
            return 0
        fi
    done
    py="${py:-python3}"
    py_source="UNRESOLVED"
    return 1
}

# mamba_env / mamba_env_prefix may be supplied by project_paths.conf, like scratch_proj_dir.
_pe_conf="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/project_paths.conf"
if [ -f "$_pe_conf" ]; then
    while IFS='=' read -r _k _v; do
        case "$_k" in
            mamba_env_prefix) [ -z "${mamba_env_prefix:-}" ] && [ -n "$_v" ] && mamba_env_prefix="$_v" ;;
        esac
    done < <(grep -E '^[a-z_]+=' "$_pe_conf")
fi

_resolve_python || true
export py

# Loud check. Call this before any real work.
require_python() {
    if _py_usable "$py"; then
        return 0
    fi
    cat <<MSG >&2

CANNOT FIND A PYTHON WITH: $_py_modules
  tried            : exported \$py, python3 on PATH, micromamba env '$mamba_env',
                     \$mamba_env_prefix and the usual micromamba roots
  current \$py      : $py  [$py_source]
  python3 on PATH  : $(command -v python3 2>/dev/null || echo 'none')

On the cluster the dependencies are in the micromamba environment
'$mamba_env'. A Slurm job does not inherit an activated environment, so either
activate it before submitting, or let these scripts find it:

    micromamba activate $mamba_env          # then sbatch / bash as usual

or point them at it explicitly, once:

    export mamba_env_prefix=\$(micromamba info --base)/envs/$mamba_env
    # or add a line to scripts/mc_phenotype_assoc/project_paths.conf:
    #     mamba_env_prefix=/path/to/envs/$mamba_env

or override the interpreter outright:

    export py=/path/to/envs/$mamba_env/bin/python3

Refusing to start work that would fail partway through.
MSG
    return 1
}
