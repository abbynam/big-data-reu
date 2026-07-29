This file documents step by step instructions on how to create an MCDE run. 

Start an interactive job if you have not already:  

```
srun --cluster=chip-cpu --account=pi_gobbert --partition=2024 --qos=shared --mem=40G --pty /bin/bash
```

or

```
srun --cluster=chip-cpu --account=cybertrn --partition=2024 --qos=shared --mem=40G --pty /bin/bash
```

Activate the appropriate conda environment.  


```
unset $PYTHONPATH
module purge 
echo $PYTHONPATH
conda activate /umbc/rs/cybertrn/users/anam2/conda_envs/bride
```
Copy an existing example run and rename it or create a new directory in the runs folder, and make sure to copy the ```run_chip.slurm``` file into it.  

Then, ```sbatch run_chip.slurm``` to get output in ```slurm```.

```oldslurm``` used
	```
	/umbc/rs/cybertrn/reu2025/team2/PG_DICOM_cat/1-item1_PG_cat.csv
	/umbc/rs/cybertrn/reu2025/team2/PG_DICOM_cat/1-item1_511_cat.csv
	```  
	
```newslurm``` used the water phantom contour in ```../../csv_in```

