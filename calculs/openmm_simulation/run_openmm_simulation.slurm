#!/bin/bash
#SBATCH --job-name=openmm
#SBATCH --error=%j.err.txt
#SBATCH --output=%j.out.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
##SBATCH --partition=normal2

conda activate md-openmm
python run_openmm_simulationsqueue.py
