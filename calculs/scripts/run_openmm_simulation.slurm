#!/bin/bash
#SBATCH --job-name=openmm
#SBATCH --error=%j.err.txt
#SBATCH --output=%j.out.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
##SBATCH --partition=normal2

source /home/10034103/miniforge3/etc/profile.d/conda.sh
conda init
conda activate md-openmm

python run_openmm_simulations.py
