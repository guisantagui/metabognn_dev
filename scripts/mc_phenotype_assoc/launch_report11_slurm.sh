#!/bin/bash
# Batch wrapper for reproduce_report11.sh, so the console is not blocked.
#
#   cd scripts/mc_phenotype_assoc
#   sbatch launch_report11_slurm.sh                  # HEPG2, sections 8 9 10 11
#   sbatch launch_report11_slurm.sh 9 10 11          # selected sections
#   sbatch --export=ALL,cell=A549 launch_report11_slurm.sh
#
# Watch it:  squeue -u $USER      then      tail -f logs/report11_<jobid>.out
#
# WHY --output IS WHAT IT IS. Slurm expands --output BEFORE the job script runs, so it
# cannot contain a shell variable and is resolved relative to the SUBMIT directory. Hence
# the relative `logs/` path and the `cd` into this script's directory below: submit from
# anywhere else and the log lands somewhere else. submit_assoc.sh documents the same trap.
#
# RESOURCES. Everything is single-process NumPy on local files: no MPI, no GPU, no Gurobi.
# The cost is the size-matched random-set null, 200 draws x ~45 sizes per arm, which is
# CPU-bound and single-threaded. Memory is dominated by holding one label matrix
# (~1.4 GB for HEPG2 metchange, ~2.7 GB expr) plus a float copy, so 48G is comfortable
# and 24G would probably do. Four arms x ~12 min, plus the derived-matrix builds.
#SBATCH --job-name=report11
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=logs/report11_%j.out
set -uo pipefail

# Locate this directory. NOT from BASH_SOURCE / $0: Slurm copies the batch script to
# /var/spool/slurmd/job<id>/ before running it, so dirname of the script is the spool
# directory (job 6428066 died exactly that way: "reproduce_report11.sh: No such file").
# Walk up from the submit directory instead, as launch_assoc_HEPG2_slurm.sh does.
marker="scripts/mc_phenotype_assoc/reproduce_report11.sh"
here=""
for start in "${SLURM_SUBMIT_DIR:-}" "$PWD" "$(dirname "${BASH_SOURCE[0]}")"; do
    [ -n "$start" ] || continue
    d="$(cd "$start" 2>/dev/null && pwd)" || continue
    while [ -n "$d" ] && [ "$d" != "/" ]; do
        if [ -f "$d/$marker" ]; then here="$d/scripts/mc_phenotype_assoc"; break 2; fi
        d="$(dirname "$d")"
    done
done
if [ -z "$here" ]; then
    echo "cannot locate $marker upward from SLURM_SUBMIT_DIR=${SLURM_SUBMIT_DIR:-unset} or $PWD." >&2
    echo "Submit from inside the repository (e.g. scripts/mc_phenotype_assoc)." >&2
    exit 1
fi
mkdir -p "$here/logs"
cd "$here"

echo "host        : $(hostname)"
echo "job         : ${SLURM_JOB_ID:-<interactive>}"
echo "started     : $(date '+%F %T')"
echo "cell        : ${cell:-HEPG2}"
echo "sections    : ${*:-all}"
echo

# Keep the BLAS libraries single-threaded. The work here is one big NumPy loop per arm,
# not linear algebra, so letting OpenBLAS spawn threads buys nothing and oversubscribes
# the allocation.
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

bash "$here/reproduce_report11.sh" "$@"
rc=$?

echo
echo "finished    : $(date '+%F %T')   exit=$rc"
if [ "$rc" -ne 0 ]; then
    echo "FAILED. The two lines to check first are the 'labels base' and 'DESIGN' lines"
    echo "near the top of this log: a labels base pointing at the repo rather than scratch"
    echo "means mc_phenotype_assoc_common.sh did not resolve, and you can force it with"
    echo "  sbatch --export=ALL,labels_base=/scratch/users/\$USER/test_large_files/metabognn_dev ..."
fi
exit $rc
