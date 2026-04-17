# Estudi Pas A Pas del Trp-cage 1L2Y Model 1

S'ha triat una estructura de l'ensamblat NMR del Trp-cage:

- PDB: `1L2Y`
- model seleccionat: `model 1`
- fitxer utilitzat: [trp_cage_1l2y_model1.pdb](/Users/gemmacollvila/github/EVB-LEH/inputs/trp_cage_1l2y_model1.pdb)

## Flux recomanat

Els passos que segueix ara el repositori són:

1. netejar PDB
2. afegir estat de protonació i preparar carregues
3. preparar una variant `prot + lligands` i una `nomes prot`
4. executar la dinàmica molecular
5. analitzar RMSD, RMSF i radi de gir

## Pas 1. Netejar PDB

```bash
bash scripts/run_clean_pdb.sh
```

Sortida esperada:

- `study_runs/01_clean_pdb/cleaned.pdb`
- `study_runs/01_clean_pdb/summary.txt`

## Pas 2. Afegir protonació i carregues

```bash
bash scripts/run_add_protonation.sh
```

Sortida esperada:

- `study_runs/02_protonation/protonated.pdb`
- `study_runs/02_protonation/summary.txt`

## Pas 3. Preparar variants

```bash
bash scripts/run_prepare_variants.sh
```

Sortida esperada:

- `study_runs/03_variants/protein_only.pdb`
- `study_runs/03_variants/protein_plus_ligands.pdb`
- `study_runs/03_variants/summary.txt`

Nota:

- en el cas del Trp-cage `1L2Y` no hi ha lligands, així que les dues variants coincideixen

## Pas 4. Dinàmica molecular

```bash
bash scripts/run_md_study.sh
```

Sortida esperada:

- `study_runs/04_md/solvated_initial.pdb`
- `study_runs/04_md/minimized.pdb`
- `study_runs/04_md/after_nvt.pdb`
- `study_runs/04_md/after_npt.pdb`
- `study_runs/04_md/production_final.pdb`
- trajectòries i fitxers CSV d'estat

## Pas 5. Analisi

```bash
bash scripts/run_basic_analysis.sh
```

Sortida esperada:

- `study_runs/05_analysis/rmsd.csv`
- `study_runs/05_analysis/rmsf_ca.csv`
- `study_runs/05_analysis/radius_of_gyration.csv`
- `study_runs/05_analysis/summary.txt`

## Notes

- Aquest estudi està configurat sobre una proteïna real molt petita per fer proves ràpides.
- El model triat és explícitament el `model 1` de l'entrada NMR `1L2Y`.
- Si es vol comparar amb altres estructures de l'ensamblat NMR, caldria generar nous fitxers d'entrada per a cada model.
- L'analisi necessita `mdtraj`, que no està instal·lat actualment en aquest entorn.
