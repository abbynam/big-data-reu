This file details step-by-step instructions on how to use BRIDE to test different machine learning models.  

Copy the predict.sh file into a new directory. Change the job name and **CHANGE THE RUN ID**.  

Create a new config file by copying the ```template.yaml``` file, name it the same as your run. Change the model by editing the ```mdl_key``` (reference ```utils.py```). 
Make sure the run_id matches the run_id in your .sh file). Modify the train_data_path and test_data_path.

Then: ```sbatch predict.sh```  

You can interactively see the training of the model in progress by using: 
```
tail -f slurm.out (see the progress of each epoch)
tail -f slurm.err (see the new best score of the validation loss per epoch)
```

