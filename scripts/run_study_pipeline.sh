#!/usr/bin/env bash
set -euo pipefail

INPUT_PDB="inputs/trp_cage_1l2y_model1.pdb"
STUDY_DIR="study_runs"
PH="7.0"
TEMPERATURE="300.0"
PADDING_NM="1.2"
IONIC_STRENGTH="0.15"
FRICTION="1.0"
TIMESTEP_FS="2.0"
REPORT_INTERVAL="500"
MINIMIZE_ITERATIONS="5000"
NVT_STEPS="25000"
NPT_STEPS="100000"
PRODUCTION_STEPS="250000"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-pdb)
      INPUT_PDB="$2"
      shift 2
      ;;
    --study-dir)
      STUDY_DIR="$2"
      shift 2
      ;;
    --ph)
      PH="$2"
      shift 2
      ;;
    --temperature)
      TEMPERATURE="$2"
      shift 2
      ;;
    --padding-nm)
      PADDING_NM="$2"
      shift 2
      ;;
    --ionic-strength)
      IONIC_STRENGTH="$2"
      shift 2
      ;;
    --friction)
      FRICTION="$2"
      shift 2
      ;;
    --timestep-fs)
      TIMESTEP_FS="$2"
      shift 2
      ;;
    --report-interval)
      REPORT_INTERVAL="$2"
      shift 2
      ;;
    --minimize-iterations)
      MINIMIZE_ITERATIONS="$2"
      shift 2
      ;;
    --nvt-steps)
      NVT_STEPS="$2"
      shift 2
      ;;
    --npt-steps)
      NPT_STEPS="$2"
      shift 2
      ;;
    --production-steps)
      PRODUCTION_STEPS="$2"
      shift 2
      ;;
    *)
      echo "Argument no reconegut: $1" >&2
      exit 2
      ;;
  esac
done

echo "1/5 Netejant PDB..."
python scripts/01_clean_pdb.py \
  --input-pdb "$INPUT_PDB" \
  --study-dir "$STUDY_DIR"

echo "2/5 Afegint estat de protonacio i preparant carregues..."
python scripts/02_add_protonation.py \
  --study-dir "$STUDY_DIR" \
  --ph "$PH" \
  --temperature "$TEMPERATURE"

echo "3/5 Preparant variants: proteina + lligands i nomes proteina..."
python scripts/03_prepare_variants.py \
  --study-dir "$STUDY_DIR"

echo "4/5 Executant dinamica molecular..."
python scripts/04_run_md.py \
  --study-dir "$STUDY_DIR" \
  --ph "$PH" \
  --padding-nm "$PADDING_NM" \
  --ionic-strength "$IONIC_STRENGTH" \
  --temperature "$TEMPERATURE" \
  --friction "$FRICTION" \
  --timestep-fs "$TIMESTEP_FS" \
  --report-interval "$REPORT_INTERVAL" \
  --minimize-iterations "$MINIMIZE_ITERATIONS" \
  --nvt-steps "$NVT_STEPS" \
  --npt-steps "$NPT_STEPS" \
  --production-steps "$PRODUCTION_STEPS"

echo "5/5 Analitzant RMSD, RMSF i radi de gir..."
python scripts/05_analyze_basic.py \
  --study-dir "$STUDY_DIR"

echo "Pipeline completat. Resultats a: $STUDY_DIR"
