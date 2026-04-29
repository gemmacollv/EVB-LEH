#!/usr/bin/env bash
set -euo pipefail

python scripts/openmm_md.py \
  --steps 5000 \
  --equilibration-steps 1000 \
  --temperature 300.0 \
  --friction 1.0 \
  --timestep-fs 2.0 \
  --padding-nm 1.2 \
  --report-interval 500 \
  --output-dir outputs \
  "$@"
