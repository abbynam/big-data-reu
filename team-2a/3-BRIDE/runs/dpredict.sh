#!/bin/bash
#SBATCH --job-name=demo_run                
#SBATCH --output=slurm/output
#SBATCH --error=slurm/error
#SBATCH --account=pi_gobbert
#SBATCH --ntasks=1                     
#SBATCH --cpus-per-task=1              
#SBATCH --time=10:20:00                
#SBATCH --mem=48G                     
#SBATCH --cluster=chip-gpu
#SBATCH --gres=gpu:1                   
#SBATCH --constraint='L40S|RTX_8000|RTX_6000'  # Request L40S, RTX_8000, RTX_6000 GPUs
#SBATCH --partition=gpu 

# variables
run_id='demo'  # CHANGE THIS
# shouldn't change variables below
PROGRAM_BASE=/umbc/rs/cybertrn/reu2026/team2/research/3-BRIDE
config_path='../../config/'
config_path+=${run_id}
config_path+='.yaml'

# DON'T CHANGE ANYTHING BELOW
# activate conda env
module load Anaconda3/2024.02-1
source /usr/ebuild/software/emerald/software/Anaconda3/2024.02-1/bin/activate
echo "activating conda environment..."
eval "$(conda shell.bash hook)"
conda activate /umbc/rs/cybertrn/reu2024/team2/envs/ada_main  # choose which environment carefully  
echo "conda environment activated."

# debugging flags
export NCCL_DEBUG=INFO
export PYTHONFAULTHANDLER=1

# # write this script to the outputsß
# scontrol write bash_script $SLURM_JOB_ID slurm_output/script.txt

# run
srun python3 ../../predict.py -c $config_path
echo "finished running predict.py, now running eval.py ..."
# get the largest version number
cd ${PROGRAM_BASE}/logs/csv_logs/${run_id}/lightning_logs/
version=$(ls | grep '^version_[0-9]\+$' | cut -c2- | sort -n | tail -n1 | sed 's|^|v|')
echo "using version: "${version}
cd - 
# do rest of run
#srun python3 ../../eval/eval.py -o b -p eval/${run_id}/ -t logs/csv_logs/${run_id}/lightning_logs/${version}/
echo "eval.py skipped"
#echo "finished running eval.py"
conda deactivate
echo "conda environment deactivated."
