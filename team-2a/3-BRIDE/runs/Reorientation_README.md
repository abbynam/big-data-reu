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

- **FOR REORIENTATION ONLY:** (example) python ../../../csvToNpy.py AllEventsCombined_M400-scatters-UnPhy_3x_reoriented90.csv
  - You should see X.npy and Y.npy in your subdirectory in farshad-converted. 

