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
#   3. the micromamba environment $mamba_env, reached the way it is reached by
#      hand on this cluster - MAMBA_ROOT_PREFIX exported first, then the shell
#      hook, then activate:
#          export MAMBA_ROOT_PREFIX=/work/projects/.../software/micromamba
#          eval "$(micromamba shell hook --shell bash)"
#          micromamba activate metabolic_drug_repurposing
#      The root prefix matters: the environment lives under a shared project
#      directory, NOT the default ~/micromamba, so without it micromamba looks in
#      the wrong place and reports no such environment. Roots come from
#      project_paths.conf (mamba_root_prefix), so the path is written once.
#   4. failing that, $mamba_env_prefix/bin/python3 or <root>/envs/$mamba_env/bin/python3
#      used directly, which needs neither the hook nor micromamba on PATH.
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
    # 3. If micromamba is already on PATH and configured - a shell rc, a .mambarc,
    #    or MAMBA_ROOT_PREFIX inherited from the submitting shell - then it knows
    #    its own root and which environments exist. Ask it rather than guessing: a
    #    root that is not in the candidate list below would otherwise be missed
    #    entirely, which is how the first version of this file failed its own test.
    local root mm hook cand mm_base
    mm="$(command -v micromamba 2>/dev/null || true)"
    if [ -n "$mm" ]; then
        cand="$("$mm" run -n "$mamba_env" python3 -c 'import sys; print(sys.executable)' \
                2>/dev/null | tail -1)"
        if _py_usable "$cand"; then
            py="$cand"
            py_source="micromamba run -n $mamba_env"
            return 0
        fi
        # its base directory becomes the first candidate root for the loop below
        mm_base="$("$mm" info --base 2>/dev/null | tr -d '[:space:]')"
        if [ -z "$mm_base" ]; then
            mm_base="$("$mm" info 2>/dev/null \
                       | sed -n 's/.*base environment *: *\([^ ]*\).*/\1/p' | tail -1)"
        fi
        [ -n "$mm_base" ] && _mamba_roots="$mm_base ${_mamba_roots:-}"
    fi

    # 4. Each candidate root in turn. MAMBA_ROOT_PREFIX must be exported BEFORE
    #    micromamba is asked anything, because the environment lives under a shared
    #    project directory rather than ~/micromamba.
    for root in "${MAMBA_ROOT_PREFIX:-}" ${_mamba_roots:-} \
                "$HOME/micromamba" "$HOME/.micromamba" \
                "$HOME/miniforge3" "$HOME/mambaforge"; do
        [ -n "$root" ] || continue
        [ -d "$root" ] || continue
        export MAMBA_ROOT_PREFIX="$root"

        # (a) the hook-then-activate route, which is the procedure known to work
        #     here. Activation also exports the environment's own PATH and
        #     LD_LIBRARY_PATH, which a bare interpreter path would not.
        if [ -z "$mm" ]; then
            for cand in "$root/bin/micromamba" "$root/micromamba" \
                        "$root/condabin/micromamba"; do
                [ -x "$cand" ] && mm="$cand" && break
            done
        fi
        if [ -n "$mm" ] && hook="$("$mm" shell hook --shell bash 2>/dev/null)"; then
            eval "$hook" 2>/dev/null || true
            if micromamba activate "$mamba_env" 2>/dev/null && _py_usable python3; then
                py="$(command -v python3)"
                py_source="micromamba activate $mamba_env (root $root)"
                return 0
            fi
        fi

        # (b) the environment's interpreter called directly - no hook, and
        #     micromamba need not be on PATH at all.
        if _py_usable "$root/envs/$mamba_env/bin/python3"; then
            py="$root/envs/$mamba_env/bin/python3"
            py_source="prefix $root/envs/$mamba_env"
            return 0
        fi
    done

    # 5. an environment prefix given outright, wherever it lives
    if _py_usable "${mamba_env_prefix:-}/bin/python3"; then
        py="$mamba_env_prefix/bin/python3"
        py_source="mamba_env_prefix $mamba_env_prefix"
        return 0
    fi
    if _py_usable "${CONDA_PREFIX:-}/../$mamba_env/bin/python3"; then
        py="$CONDA_PREFIX/../$mamba_env/bin/python3"
        py_source="sibling of CONDA_PREFIX"
        return 0
    fi
    py="${py:-python3}"
    py_source="UNRESOLVED"
    return 1
}

# mamba_env / mamba_env_prefix may be supplied by project_paths.conf, like scratch_proj_dir.
_pe_conf="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/project_paths.conf"
_mamba_roots=""
if [ -f "$_pe_conf" ]; then
    while IFS='=' read -r _k _v; do
        case "$_k" in
            # repeatable: every mamba_root_prefix line becomes a candidate, so a
            # path for another machine can sit here harmlessly (only existing
            # directories are tried), exactly as for scratch_proj_dir.
            mamba_root_prefix) [ -n "$_v" ] && _mamba_roots="$_mamba_roots $_v" ;;
            mamba_env_prefix)  [ -z "${mamba_env_prefix:-}" ] && [ -n "$_v" ] && mamba_env_prefix="$_v" ;;
            mamba_env)         [ -z "${mamba_env_explicit:-}" ] && [ -n "$_v" ] && mamba_env="$_v" ;;
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

The usual procedure here, which these scripts now reproduce by themselves:

    export MAMBA_ROOT_PREFIX=/work/projects/car_t_methylome/software/micromamba
    eval \"\$(micromamba shell hook --shell bash)\"
    micromamba activate $mamba_env

If that root is wrong on this machine, correct the mamba_root_prefix line in
scripts/mc_phenotype_assoc/project_paths.conf (several lines are allowed; only
directories that exist are tried):

    mamba_root_prefix=/your/micromamba/root

roots tried: ${MAMBA_ROOT_PREFIX:-unset}${_mamba_roots:+ $_mamba_roots}

or override the interpreter outright:

    export py=/path/to/envs/$mamba_env/bin/python3

Refusing to start work that would fail partway through.
MSG
    return 1
}
