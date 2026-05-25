#!/bin/bash
#SBATCH --job-name=openmm
#SBATCH --error=%j.err.txt
#SBATCH --output=%j.out.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=gpu
#SBATCH --time=4:15:00
#SBATCH --signal=B:TERM@120

cd /home/10034103@uvic.local/EVB-LEH/calculs/scripts

conda activate md-openmm
RUN_NAME="${RUN_NAME:-run_001}"
python 03_openmm_md_simulation_holo.py --run-name "$RUN_NAME"
