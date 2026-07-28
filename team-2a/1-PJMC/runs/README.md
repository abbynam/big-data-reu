This file details step-by-step instructions on how to create a PJMC run.

**Editing the bash file**  
Before you do anything in PJMC, make sure your ```~/bash.rc``` file includes the following: 

        # from Vijay
        module load ROOT/6.24.06-foss-2021b
        source /umbc/rs/pi_gobbert/common/PromptGamma/geant4-install-seq/bin/geant4.sh
        ###DATA FILES
        export G4ENSDFSTATEDATA=/umbc/rs/pi_gobbert/common/PromptGamma/geant4-install-seq/share/Geant4-10.7.3/data/G4ENSDFSTATE2.3
        export G4PARTICLEXSDATA=/umbc/rs/pi_gobbert/common/PromptGamma/geant4-install-seq/share/Geant4-10.7.3/data/G4PARTICLEXS3.1.1
        export G4LEDATA=/umbc/rs/pi_gobbert/common/PromptGamma/geant4-install-seq/share/Geant4-10.7.3/data/G4EMLOW7.13
        export G4LEVELGAMMADATA=/umbc/rs/pi_gobbert/common/PromptGamma/geant4-install-seq/share/Geant4-10.7.3/data/PhotonEvaporation5.7
        export G4NEUTRONHPDATA=/umbc/rs/pi_gobbert/common/PromptGamma/geant4-install-seq/share/Geant4-10.7.3/data/G4NDL4.6
        export G4ABLADATA=/umbc/rs/pi_gobbert/common/PromptGamma/geant4-install-seq/share/Geant4-10.7.3/data/G4ABLA3.1
        export G4INCLDATA=/umbc/rs/pi_gobbert/common/PromptGamma/geant4-install-seq/share/Geant4-10.7.3/data/G4INCL1.0
        export G4PIIDATA=/umbc/rs/pi_gobbert/common/PromptGamma/geant4-install-seq/share/Geant4-10.7.3/data/G4PII1.3
        export G4RADIOACTIVEDATA=/umbc/rs/pi_gobbert/common/PromptGamma/geant4-install-seq/share/Geant4-10.7.3/data/RadioactiveDecay5.6
        export G4REALSURFACEDATA=/umbc/rs/pi_gobbert/common/PromptGamma/geant4-install-seq/share/Geant4-10.7.3/data/RealSurface2.2
        export G4SAIDXSDATA=/umbc/rs/pi_gobbert/common/PromptGamma/geant4-install-seq/share/Geant4-10.7.3/data/G4SAIDDATA2.0

Then run the command: source ~/bashrc to update. If it worked, you should see no output.

**Creating a run**  
Start an interactive job using the following command: 

```
srun --cluster=chip-cpu --account=cybertrn --qos=shared --time=02:50:00 --mem=16G --pty /bin/bash
```

Then, navigate to ```1-PJMC/runs``` and copy the example run. Rename it appropriately. 

```
cp -r example abby_test1
```

Edit the ```run.slurm``` file. You should change the following fields:  

- job-name 
- array (1-3%3 is sufficient for test runs) 
- account (to cybertrn)

If necessary, also edit and change the appropriate fields in the ```parameters.mac``` file:  

- /beam/position
- /beam/energy
- /run/beamOn → represents the number of events

Submit the ```run.slurm``` file as a job. You should see 3 jobs in the queue (if using array=1-3%3)

```
sbatch run.slurm
```

Navigate into the root directory. If you see that you do not have a ```compiled.root``` file, go back to the previous folder level and run the following: 

```
hadd root/compiled.root root/Team2*.root
```

(You can also rewrite the compiled.root by running: ```hadd -f root/compiled.root root/Team2*.root```)

Turn ```compiled.root``` into a CSV by running the following. 
Make sure you run this first: 

```
pip install uproot
```

Then:
```
python ../../scripts/converter_core.py "$(pwd)"
```

The csv directory of your run folder should contain a 511 (annihilation data, high energy collisions) and PG (lower energy collisions) CSV. 
