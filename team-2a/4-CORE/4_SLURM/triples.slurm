#!/bin/bash
#SBATCH --job-name=CORE                # Job name
#SBATCH --output=6_OUTPUT/1-CORE_out/triples/slurm.op     
#SBATCH --error=6_OUTPUT/1-CORE_out/triples/slurm.err
#SBATCH --ntasks=1                     # Number of tasks (1 in this case)
#SBATCH --cpus-per-task=4  			   # Number of CPU cores per task
#SBATCH --time=00:20:00                # Maximum run time (5 minutes, adjust as needed)
#SBATCH --mem=2000                     # Memory
#SBATCH --gres=gpu:1                   # Request 1 GPU
#SBATCH --cluster=chip-gpu             # Specify the cluster (adjust based on your cluster's name)
#SBATCH --constraint='L40S|RTX_8000|RTX_6000'              # Request L40S GPUs
#SBATCH --partition=gpu                # Use the GPU partition (adjust based on your cluster's partition names)

module load CUDA/12.9.0
module load GCCcore/11.2.0

PROGRAM=./coreGPU
srun ./1_BUILD/coreGPU 3_CFG/triples.cfg
# For use with MPI
# mpirun "${PROGRAM}" FLAG FLAG FLAG
