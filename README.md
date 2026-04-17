# EVB-LEH

Exemple mínim de dinàmica molecular amb OpenMM sobre una proteïna real petita.

## Ús

Instal·la OpenMM i executa:

```bash
python openmm_md.py --steps 5000 --output-dir outputs
```

El script carrega el Trp-cage miniprotein (`1L2Y`, model 1), afegeix hidrògens,
solvata el sistema amb aigua TIP3P, minimitza l'energia i executa una simulació
curta de dinàmica molecular. Els fitxers de sortida es guardaran a `outputs/`.

## Protocol Trp-cage

Hi ha una guia pas a pas específica per a la proteïna petita d'exemple a:

- [PROTOCOL_TRP_CAGE_MD.md](/Users/gemmacollvila/github/EVB-LEH/PROTOCOL_TRP_CAGE_MD.md)
- [RUN_MODEL1_STUDY.md](/Users/gemmacollvila/github/EVB-LEH/RUN_MODEL1_STUDY.md)

## Scripts per passos

També hi ha una versió separada del protocol en petits scripts, seguint aquest esquema:

1. netejar PDB
2. afegir protonació i carregues
3. preparar variants amb i sense lligands
4. executar MD
5. analitzar RMSD, RMSF i radi de gir

```bash
bash scripts/run_clean_pdb.sh
bash scripts/run_add_protonation.sh
bash scripts/run_prepare_variants.sh
bash scripts/run_md_study.sh
bash scripts/run_basic_analysis.sh
```

L'entrada inicial està a [trp_cage_1l2y_model1.pdb](/Users/gemmacollvila/github/EVB-LEH/inputs/trp_cage_1l2y_model1.pdb) i cada pas escriu la seva sortida a `study_runs/`.
