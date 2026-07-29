#!/bin/bash

#SBATCH --job-name=harrison_initial		#
#SBATCH --output=root/slurm_op/slurm%a.out	#
#SBATCH --error=root/slurm_err/slurm%a.err	#
#SBATCH --cluster=chip-cpu			# chip-cpu
#SBATCH --account=pi_gobbert			# or: cybertrn
#SBATCH --partition=2024			# 
#SBATCH --qos=shared				# shared
#SBATCH --time=05:00:00				# up as needed
#SBATCH --cpus-per-task=1			# 1
#SBATCH --ntasks=1				# 1
#SBATCH --nodes=1				# 1
#SBATCH --mem=250G				# up to 256G on any node, leaving a small reservation for the OS
#SBATCH --array=1-360%360                       # number of nodes to use and max to load at once

# load root
module load ROOT/6.24.06-foss-2021b

# run the program to make root data
srun ../../config/POLARIS3 "$(pwd)/parameters.mac" root/Team2_Data_150

# ------subsequent run, dependent job

# combine all files into one
#hadd root/compiled.root root/Team2*.root

# fix permissions so the group can collaborate
chmod g+w root/*.root

# convert the compiled file to csv
#python ../../scripts/converter_core.py "$(pwd)"

