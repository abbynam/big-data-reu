import os
import sys
import yaml
import subprocess
import getpass
import argparse
from sklearn.model_selection import ParameterGrid
from datetime import datetime
from pathlib import Path

def print_usage():
    print("Usage: python grid_search.py <run_name> <batch_sizes> <lrs> <dropouts>")
    print("Usage: python grid_search.py --lr=0.1,0.2 --dropout=0.1,0.2 --run_name='abc' --batch_size=1,2")
    print("Optional: --username --model --template_path --config_dir --slurm_dir")


def grid_search(run_name, batch_size, lr, dropout,
                username='student', model='cnn_1d_2025',
                template_path='template.yaml',
                config_dir='../config',
                slurm_dir='../runs'):

    #Edit the parameters as needed to grid search, add any additional numbers to the rows as necessary
    #Edit the Slurm script time as needed 
    #Under GRID COMBOS section, certain values are hardcoded and held constant for each run. change as needed.


    # ==PARAMETERS=====
    #CHANGE THESE 

    # python running dir: grid_search/
    base_dir = Path(__file__).resolve().parent

    # add what numbers to test for hyperparameters
    param_grid = {
        'batch_size': [batch_size],  # change as needed
        'patience' : [3000], #change as needed
        'input_size':[15],
        'hidden_layers': [[2048]], 
        'lr': [lr],  # change as needed
        'lr_step':[10], #change as needed
        'lr_gam':[0.1], #change as needed 
        'dropout': [dropout],  # change numbers and add to rows as necessary 
        'data_name': ['barajas'],  # change dataset names as needed
        'l2':[0.001]
    }

    # ===CHECKING DIRECTORY EXISTENCE ====
    os.makedirs(base_dir/config_dir, exist_ok=True)

    # ===LOAD IN TEMPLATE ====
    with open(base_dir/config_dir/template_path) as f: #grid_search folder is stored in base
        template = yaml.safe_load(f)

    # ===GRID COMBOS =====
    date_and_time = datetime.now().strftime("%y%m%d_%H%M%S") 
    grid = ParameterGrid(param_grid)

    for i, params in enumerate(grid):
        tc = template.copy()
        run_id = f"{username}{date_and_time}{run_name}_{i+1}"
        run_dir = (base_dir/"../runs"/run_id)
        tc['run_id'] = run_id
        tc['pred_ckpt'] = ''
        tc['resume_ckpt'] = ''
        tc['mdl_key'] = model  # CHANGE MODEL AS NEEDED
        tc['data']['train_data_path'] = f"/umbc/rs/cybertrn/reu2025/team2/research/base/pp2/data/{params['data_name']}/train/"
        tc['data']['test_data_path'] = f"/umbc/rs/cybertrn/reu2025/team2/research/base/pp2/data/{params['data_name']}/test/"
        tc['data']['batch_size'] = params['batch_size']
        tc['data']['val_split'] =  0.2 #change if needed
        tc['fit']['max_epochs'] = 200  #change if needed
        tc['fit']['n_nodes'] = 1  # Usually 1
        tc['fit']['n_devices'] = 4  # num of GPUs
        tc['fit']['patience'] = params['patience']  # if needed
        tc['fit']['ckpt_freq'] = 200 #CHANGE IF NEEDED
        tc['model']['input_size'] = params['input_size']
        tc['model']['num_classes'] = 13
        tc['model']['hidden_layers'] = params['hidden_layers']
        tc['model']['activation'] = 'relu'
        tc['model']['lr'] = params['lr']
        tc['model']['lr_step'] = params['lr_step'] 
        tc['model']['lr_gam'] = params['lr_gam']
        tc['model']['penalty'] = 0
        tc['model']['dropout'] = params['dropout']
        tc['model']['l2'] = params['l2'] 
        # write config YAML
        yaml_path = f"{config_dir}/{run_id}.yaml"
        with open(yaml_path, 'w') as f_out:
            yaml.dump(tc, f_out)

        # build SLURM script
        slurm_script = f"""#!/bin/bash
#SBATCH --job-name={run_id}
#SBATCH --cluster=chip-gpu
#SBATCH --mem=48G
#SBATCH --nodes=1                # num nodes: MUST match .yaml file
#SBATCH --gres=gpu:4             # num gpus per node: MUST match .yaml file AND ntasks-per-node=
#SBATCH --ntasks-per-node=4      # num gpus per node: MUST match .yaml file AND gres=gpu:
#SBATCH --time=23:00:00          # Time limit days-hrs:min:sec
#SBATCH --constraint=rtx_6000    # see hpcf website
#SBATCH --error={run_dir}/slurm_output/slurm.err
#SBATCH --output={run_dir}/slurm_output/slurm.out

# variables
run_id={run_id}  
# shouldn't change variables below
config_path='{base_dir}/../config/'${{run_id}}'.yaml'

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
export MKL_THREADING_LAYER=GNU

# run
srun python3 {base_dir}/../train.py -c $config_path
conda deactivate
echo "conda environment deactivated."
"""

        # check writing out slurm file 
        slurm_directory = slurm_dir + run_id
        os.makedirs(slurm_directory, exist_ok=True)
        # Give user used hyperparameter combination for the run
        param_path = f"{run_dir}/used_hyperparams.yaml"

        # ===CHECKING DIRECTORY EXISTENCE ====
        os.makedirs(os.path.dirname(param_path), exist_ok=True)

        with open(param_path, 'w') as param_file:
            yaml.dump(params, param_file)

        slurm_path = f"{run_dir}/{run_id}.sh"
        with open(slurm_path, 'w') as the_slurm:
            the_slurm.write(slurm_script)
        
        # sbatch Job
        subprocess.run(['sbatch', slurm_path], check=True)


def main():
    parser = argparse.ArgumentParser(description="Grid search launcher")

    parser.add_argument('--run_name', required=True, help='Name for the run')
    parser.add_argument('--batch_size', default='1024', help='Comma-separated list of batch sizes')
    parser.add_argument('--lr', default='0.01', help='Comma-separated list of learning rates')
    parser.add_argument('--dropout', default='0.2', help='Comma-separated list of dropout values')
    parser.add_argument('--username', default=os.getenv('USER', 'student'), help='Username for run_id')

    args = parser.parse_args()

    batch_sizes = args.batch_size.split(',')
    lrs = args.lr.split(',')
    dropouts = args.dropout.split(',')

    for bs in batch_sizes:
        for lr in lrs:
            for do in dropouts:
                grid_search(
                    run_name=args.run_name,
                    batch_size=bs,
                    lr=lr,
                    dropout=do,
                    username=args.username
                )


if __name__ == '__main__':
    main()

