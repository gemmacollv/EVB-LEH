#!/bin/bash
#SBATCH --job-name=openmm
#SBATCH --error=%j.err.txt
#SBATCH --output=%j.out.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=gpu
#SBATCH --time=4:00:00

conda activate md-openmm

python openmm_md_simulation-ligand.py
