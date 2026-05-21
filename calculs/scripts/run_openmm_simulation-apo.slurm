#!/bin/bash
#SBATCH --job-name=openmm
#SBATCH --error=%j.err.txt
#SBATCH --output=%j.out.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=gpu
#SBATCH --time=4:00:00

conda init
conda activate md-openmm
python run_openmm_simulation-apo.py
