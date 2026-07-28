#!/bin/bash
#SBATCH --cluster=chip-gpu
#SBATCH --job-name=MCDE_Process            # Job name
#SBATCH --mem=48G                      # Job memory request
#SBATCH --gres=gpu:1                   # Number of requested GPU(s)
#SBATCH --time=10:30:00                  # Time limit days-hrs:min:sec
#SBATCH --constraint=rtx_8000           # Specific hardware constraint
#SBATCH --error=slurm/error               # Error file name
#SBATCH --output=slurm/output

if [ -f "model-final" ] || [ -d "model-final" ]
then
    scancel $SLURM_ARRAY_JOB_ID
else

    module load Anaconda3/2024.02-1
    source /usr/ebuild/software/emerald/software/Anaconda3/2024.02-1/bin/activate
    echo -e "Output\n---\n"
    echo "activating conda environment"
    eval "$(conda shell.bash hook)"

    unset PYTHONPATH
    source /usr/ebuild/installs/software/Anaconda3/2024.02-1/etc/profile.d/conda.sh

    ### choose which environment carefully
    conda activate mcde
    echo "conda environment activated"

    # debugging flags
    #export NCCL_DEBUG=INFO
    #export PYTHONFAULTHANDLER=1

    ### choose file name and then args 
    srun python3 ../../scripts/MCDE_Main_Script_ErgInit_Peter_3e9.py
    conda deactivate
    echo "conda environment deactivated"
fi

