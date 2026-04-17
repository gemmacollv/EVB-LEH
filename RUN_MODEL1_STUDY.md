# Estudi Pas A Pas del Trp-cage 1L2Y Model 1

S'ha triat una estructura de l'ensamblat NMR del Trp-cage:

- PDB: `1L2Y`
- model seleccionat: `model 1`
- fitxer utilitzat: [trp_cage_1l2y_model1.pdb](/Users/gemmacollvila/github/EVB-LEH/inputs/trp_cage_1l2y_model1.pdb)

## Pas 1. Preparació del sistema

```bash
bash scripts/run_01_prepare.sh
```

Sortida esperada:

- `runs/01_prepare/prepared_solvated.pdb`
- `runs/01_prepare/summary.txt`

## Pas 2. Minimització

```bash
bash scripts/run_02_minimize.sh
```

Sortida esperada:

- `runs/02_minimize/minimized.pdb`
- `runs/02_minimize/summary.txt`

## Pas 3. Equilibratge NVT

```bash
bash scripts/run_03_nvt.sh
```

Sortida esperada:

- `runs/03_nvt/nvt_final.pdb`
- `runs/03_nvt/nvt_state.csv`
- `runs/03_nvt/nvt_trajectory.dcd`

## Pas 4. Equilibratge NPT

```bash
bash scripts/run_04_npt.sh
```

Sortida esperada:

- `runs/04_npt/npt_final.pdb`
- `runs/04_npt/npt_state.csv`
- `runs/04_npt/npt_trajectory.dcd`

## Pas 5. Producció

```bash
bash scripts/run_05_production.sh
```

Sortida esperada:

- `runs/05_production/production_final.pdb`
- `runs/05_production/production_state.csv`
- `runs/05_production/production_trajectory.dcd`

## Notes

- Aquest estudi està configurat sobre una proteïna real molt petita per fer proves ràpides.
- El model triat és explícitament el `model 1` de l'entrada NMR `1L2Y`.
- Si es vol comparar amb altres estructures de l'ensamblat NMR, caldria generar nous fitxers d'entrada per a cada model.
