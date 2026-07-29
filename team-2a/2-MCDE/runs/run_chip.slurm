#!/bin/bash
#SBATCH --cluster=chip-gpu
#SBATCH --job-name=MCDE_abbytest            # Job name
#SBATCH --mem=48G                      # Job memory request
#SBATCH --gres=gpu:1                   # Number of requested GPU(s)
#SBATCH --time=10:30:00                  # Time limit days-hrs:min:sec
#SBATCH --constraint=rtx_8000           # Specific hardware constraint
#SBATCH --error=slurm/error               # Error file name
#SBATCH --output=slurm/output

echo "Job started"

module load Anaconda3/2024.02-1

# Proper conda initialization for batch jobs
source /usr/ebuild/installs/software/Anaconda3/2024.02-1/etc/profile.d/conda.sh

echo "Activating conda environment..."
conda activate /umbc/rs/cybertrn/users/anam2/conda_envs/bride
echo "Conda activated"

# sanity checks (VERY important for debugging)
echo "Python location:"
which python
python --version
python -c "import sys; print(sys.executable)"

echo "Running script..."

srun python ../../scripts/MCDE_Main_Script_ErgInit_Peter_3e9.py

echo "Script finished"

conda deactivate
echo "Job done"
