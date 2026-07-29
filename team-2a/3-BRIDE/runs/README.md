This file details how to use BRIDE to apply a model to experimental triples data, using ```farshad``` as the original directory with the raw triple scatter data stored in CSVs and ```farshad-converted``` and ```farshad-reoriented``` as the new directories that later slurm files will reference.

**Creating training datasets**  
Navigate to ```3-BRIDE/raw_data/farshad-converted```.  

Create a new directory using mkdir with the same name as the CSV (from ```raw_data/farshad```) that you want to use. 
- *FOR REORIENTATION ONLY:* Make this new directory in ```farshad-reoriented```. 
- *FOR REORIENTATION ONLY:* In the new directory, run:
  ```python ../../../reorientData.py ../../farshad/<folder name>/AllEventsCombined_M400-scatters-UnPhy_3x.csv <90 or 270 gantry angle>```
In the new directory, run:
```
python ../../../csvToNpy.py ../../farshad/<folder name>/AllEventsCombined_M400-scatters-UnPhy_3x.csv```
to convert the CSV to numpy arrays.
```

- *FOR REORIENTATION ONLY:* (example) ```python ../../../csvToNpy.py AllEventsCombined_M400-scatters-UnPhy_3x_reoriented90.csv```
  - You should see ```X.npy``` and ```Y.npy``` in your subdirectory in ```farshad-converted```.

**Predictions**  
This is essentially the same workflow as BRIDE (See Testing_README.md) Again, the ```farshad``` datasets are used as an example for the naming convention.

Create a new directory using mkdir in 3-BRIDE/runs following the convention: ```1_farshad_pred```. 
- *FOR REORIENTATION ONLY:* Follow the convention ```##_farshad_reorient_pred```.

In runs, make a copy of the ```demo/predict.sh``` file: 
```
cp -p demo/predict.sh #_farshad_pred/#_farshad_pred.sh
```
This predict.sh file is different from the typical BRIDE slurm script because it skips the evaluation step.  

Edit the predict.sh file so the run_id follows the convention: ```1_farshad_pred```.  

Create a copy of the ```demo.yaml``` file in config and make sure you change the ```run_id```:  
```
cp demo.yaml #_farshad_pred.yaml
```

Change the test data paths to your directory with the ```X.npy``` and ```Y.npy``` arrays in ```3-BRIDE/raw_data/farshad-converted```. The train data path does not matter; you can either comment it out or just leave it as is.  
- Make sure the REORIENTED data path is to: ```farshad-reoriented```
  
Then: ```sbatch predict.sh```.  

If you get an error in the slurm.err regarding missing lightning logs: 
```
/umbc/rs/cybertrn/reu2026/team2/research/3-BRIDE/logs/csv_logs/##_farshad_pred/lightning_logs: No such file or directory
```
 
you need to create that directory:  
```
mkdir -p /umbc/rs/cybertrn/reu2026/team2/research/3-BRIDE/logs/csv_logs/##_farshad_pred/lightning_logs
```

You can find the results in:  

```
/umbc/rs/cybertrn/reu2026/team2/research/3_BRIDE/eval/##_farshad_pred
```

**Repaired**  

Navigate to your runs folder: ```3_BRIDE/runs/#_farshad_pred```  

Then, run: 
```
python ../../predToCsv.py ../../raw_data/farshad-converted/<original data folder name>/X.npy ../../eval/#_farshad_pred/y_pred.npy
```
- *FOR REORIENTATION ONLY:* Reverse the reorientation on the new CSV before reconstruction. (Gantry angle argument takes -90 or -270). For example:
  ```
  python reorientData.py runs/11_farshad_reorient_pred/repaired.csv -90
  ```

