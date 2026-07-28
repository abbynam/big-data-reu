1. Copy the template slurm file

	cp 4_SLURM/template.slurm 4_SLURM/job_name.slurm

2. Rename Job to something informative

	#SBATCH --job-name=triples_run25_sbp

2. Copy triples config, edit filename and contents to match job name

	cp 3_CFG/triples.cfg 3_CFG/triples_run25_sbp.cfg

	OUTPUT_FOLDER_PATH = 6_OUTPUT/1-CORE_out/triples_run25_sbp/

3. Choose parameters for the run, mainly

	Specify the input file:
	EVENT_FILE_PATH = ./2_DATA/something.csv

	Choose one of the four reconstruction algorithms:
	IMAGE_ALGORITHM = KEM

3. Use sbatch to submit slurm job

	CORE$ sbatch 4_SLURM/triples_run25_sbp.slurm

	The reconstruction output will end up in your OUTPUT_FOLDER_PATH

4. Run the plotter to build images from reconstruction data

	Make sure the output folder exists:
	mkdir -p 6_OUTPUT/2-plotter/triples_run25_sbp/

	python3 5_SCRIPTS/recon_plotter_ummc.py 6_OUTPUT/1-CORE_out/triples_run25_sbp/output.dat 6_OUTPUT/2-plotter/triples_run25_sbp

