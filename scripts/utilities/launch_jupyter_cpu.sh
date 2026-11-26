#!/usr/bin/bash --login
#SBATCH --job-name=Jupyter
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=batch
#SBATCH --time=0-08:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

# Clean environment (important on HPC)
module --force purge || true

module load env/development/2024a
module load math/Gurobi/12.0.1-GCCcore-13.3.0

export MAMBA_ROOT_PREFIX="/work/projects/car_t_methylome/software/micromamba"

# Activate micromamba
if command -v micromamba >/dev/null 2>&1; then
  eval "$(micromamba shell hook --shell bash)"
  micromamba activate jupyter_launch   # ONLY env that needs Jupyter
else
  export PATH="/work/projects/car_t_methylome/software/micromamba/envs/jupyter_launch/bin:$PATH"
fi

# --- Sanity checks ---
echo "----- Sanity checks -----"
echo "which python: $(which python || true)"
echo "python -V: $(python -V 2>&1 || true)"
echo "which jupyter: $(which jupyter || true)"
echo "which jupyter-lab: $(which jupyter-lab || true)"
echo "jupyter --version: $(jupyter --version 2>&1 || true)"
echo "------------------------"

# Ensure jupyter-lab exists
if ! command -v jupyter-lab >/dev/null 2>&1; then
  echo "***ERROR*** jupyter-lab not found in jupyter_launch."
  echo "Install it with: micromamba install -n jupyter_launch -c conda-forge jupyterlab"
  exit 1
fi

# --- Jupyter config ---
loopback_device="127.0.0.1"
port="18888"
connection_instructions="connection_instructions.log"

# Start Jupyter Lab
jupyter lab --ip=${loopback_device} --port=${port} --no-browser &
lab_pid=$!

# Save connection instructions
echo "# Connection instructions" > "${connection_instructions}"
echo "" >> "${connection_instructions}"
echo "Run this on your local machine:" >> "${connection_instructions}"
echo "ssh -N -J ${USER}@access-iris.uni.lu:8022 -L ${port}:${loopback_device}:${port} ${USER}@$(hostname -i)" >> "${connection_instructions}"
echo "" >> "${connection_instructions}"
echo "Then open in browser:" >> "${connection_instructions}"

# Wait and capture token
sleep 3

URL=$(jupyter lab list 2>&1 | grep -Eo 'http://[^ ]+\?token=[a-z0-9]+' | head -n 1 | sed -r 's/([0-9]{1,3}\.){3}[0-9]{1,3}/127.0.0.1/g' || true)

if [ -n "${URL}" ]; then
  echo "${URL}" >> "${connection_instructions}"
fi

# Debug info
echo -e '\n=== AVAILABLE KERNELS ==='
jupyter kernelspec list || true

echo -e '\n=== JUPYTER PATHS ==='
jupyter --paths || true

# Wait for user to stop Jupyter
wait ${lab_pid}