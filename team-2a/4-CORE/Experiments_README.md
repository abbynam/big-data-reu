This file explains how to create data files for CORE when working with experimental triples data. 

1. Navigate to ```4-CORE/2_DATA```.
2. Make a symbolic link: ```ln -s ../../3-BRIDE/runs/#_farshad_pred/repaired.csv #_farshad_repaired.csv```
   - This represents the repaired data. 
3. Make another symbolic link: ln -s ../../3-BRIDE/raw_data/farshad/5-251206-SOBP-7.5Gray-G270-trial2(replace)/AllEventsCombined_M400-scatters-UnPhy_3x.csv ##_farshad_control.csv
   - This represents the control data. 
