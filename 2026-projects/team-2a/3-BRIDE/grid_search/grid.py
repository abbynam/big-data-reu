import os
import yaml
import subprocess
import getpass
from sklearn.model_selection import ParameterGrid


#Edit the parameters as needed to grid search, add any additional numbers to the rows as necessary
#Edit the Slurm script time as needed 
#Under GRID COMBOS section, certain values are hardcoded and held constant for each run. change as needed.


# ==PARAMETERS=====
# add what numbers to test for hyperparameters
param_grid = {
    'batch_size': [1024, 2048],  # change as needed 
    'val_split': [0.2],  # change as needed 
    'num_layers': [2, 3],  # change as needed
    'neurons': [64, 128],  # change as needed 
    'lr': [0.01, 0.001],  # change as needed
    'lr_step':[10] #change as needed
    'lr_gamma':[0.1] #change as needed 
    'dropout': [0.2, 0.5],  # change numbers and add to rows as necessary 
    'data_name': ['barajas']  # change dataset names as needed
}

# === NEEDED PATHS ===
Template_path = '../../config/template.yaml'
config_output = './config'  # place output in the config folder 
slurm_out = '../../runs/'

# ===LOAD IN TEMPLATE ====
with open(Template_path) as f:
    template = yaml.safe_load(f)

# ===CREATING RUNID =========
username = 'USERNAME'  # ADD YOUR USERNAME HERE
date_and_time = 'month_date'  # ex: 0503 for May third 
run_name = 'whatever_your_run_name_is'  # name whatever this run name is

# ===CHECKING DIRECTORY EXISTENCE ====
os.makedirs(config_output, exist_ok=True)

# ===GRID COMBOS ===== 
grid = ParameterGrid(param_grid)

for i, params in enumerate(grid):
    tc = template.copy()
    run_id = f"{username}{date_and_time}{run_name}_{i+1}"
    tc['run_id'] = run_id
    tc['pred_ckpt'] = ''
    tc['resume_ckpt'] = ''
    tc['mdl_key'] = 'deep_impr_fcn'  # CHANGE MODEL AS NEEDED
    tc['data']['train_data_path'] = f"/umbc/rs/cybertrn/reu2025/team2/research/base/pp2/data/{params['data_name']}/train/"
    tc['data']['test_data_path'] = f"/umbc/rs/cybertrn/reu2025/team2/research/base/pp2/data/{params['data_name']}/test/"
    tc['data']['batch_size'] = params['batch_size']
    tc['data']['val_split'] = params['val_split']
    tc['fit']['max_epochs'] = params['max_epochs']
    tc['fit']['n_nodes'] = 1  # Usually 1
    tc['fit']['n_devices'] = 4  # num of GPUs
    # tc['fit']['patience'] = params['']  # if needed
    tc['fit']['ckpt_freq'] = 200
    tc['fit']['max_epochs'] = 60 #change as needed
    tc['model']['num_layers'] = params['num_layers']
    tc['model']['neurons'] = params['neurons']
    tc['model']['lr'] = params['lr']
    tc['model']['lr_step'] = params['lr_step'] 
    tc['model']['lr_gam'] = params['lr_gam']
    tc['model']['dropout'] = params['dropout']
    tc['model']['activation'] = 'relu'  # change if needed
    tc['model']['optimizer'] = 'adam'  # change if needed

    # write config YAML
    yaml_path = f"{config_output}/{run_id}.yaml"
    with open(yaml_path, 'w') as f_out:
        yaml.dump(tc, f_out)

    # build SLURM script
    slurm_script = f"""#!/bin/bash
#SBATCH --job-name={run_id}
#SBATCH --mem=48G
#SBATCH --nodes=1                # num nodes: MUST match .yaml file
#SBATCH --gres=gpu:4             # num gpus per node: MUST match .yaml file AND ntasks-per-node=
#SBATCH --ntasks-per-node=4      # num gpus per node: MUST match .yaml file AND gres=gpu:
#SBATCH --time=23:00:00          # Time limit days-hrs:min:sec
#SBATCH --constraint=rtx_6000    # see hpcf website
#SBATCH --error=slurm_output/slurm.err
#SBATCH --output=slurm_output/slurm.out

# variables
run_id={run_id}  # CHANGE THIS
# shouldn't change variables below
config_path='../../config/'${{run_id}}'.yaml'

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

# run
srun python3 ../../train.py -c $config_path
conda deactivate
echo "conda environment deactivated."
"""

    # check writing out slurm file 
    slurm_directory = slurm_out + run_id
    os.makedirs(slurm_directory, exist_ok=True)
    slurm_path = f"{slurm_directory}/{run_id}.sh"
    with open(slurm_path, 'w') as the_slurm:
        the_slurm.write(slurm_script)

    # sbatch Job
    subprocess.run(['sbatch', slurm_path], check=True)
