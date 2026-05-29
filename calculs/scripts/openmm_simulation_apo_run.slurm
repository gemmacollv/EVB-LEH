#!/bin/bash
#SBATCH --job-name=openmm-apo-run
#SBATCH --error=%j.err.txt
#SBATCH --output=%j.out.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=normal2
#SBATCH --time=4:15:00
#SBATCH --signal=B:TERM@120

cd /home/10034103@uvic.local/EVB-LEH/calculs/scripts

source /home/10034103/miniforge3/etc/profile.d/conda.sh
conda activate md-openmm
RUN_NAME="${RUN_NAME:-run_001}"
CHECKPOINT_IN="${CHECKPOINT_IN:-}"

if [ -n "$CHECKPOINT_IN" ]; then
    python 03_openmm_md_simulation_apo_run.py --run-name "$RUN_NAME" --checkpoint-in "$CHECKPOINT_IN"
else
    python 03_openmm_md_simulation_apo_run.py --run-name "$RUN_NAME"
fi
