# EVB-LEH

Exemple mínim de dinàmica molecular amb OpenMM sobre una proteïna real petita.

## Ús

Instal·la OpenMM i executa:

```bash
python scripts/openmm_md.py --input-pdb inputs/trp_cage_1l2y_model1.pdb --steps 5000 --output-dir outputs
```

Per executar la simulació amb un altre PDB, només cal canviar `--input-pdb`:

```bash
python scripts/openmm_md.py --input-pdb inputs/el_teu_fitxer.pdb --steps 5000 --output-dir outputs/el_teu_sistema
```

El script carrega el PDB indicat, afegeix hidrògens, solvata el sistema amb aigua TIP3P, minimitza l'energia i executa una simulació curta de dinàmica molecular. Els fitxers de sortida es guardaran al directori indicat amb
`--output-dir`.

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
STUDY_DIR=study_runs/el_teu_sistema

bash scripts/run_clean_pdb.sh --input-pdb inputs/el_teu_fitxer.pdb --study-dir "$STUDY_DIR"
bash scripts/run_add_protonation.sh --study-dir "$STUDY_DIR"
bash scripts/run_prepare_variants.sh --study-dir "$STUDY_DIR"
bash scripts/run_md_study.sh --study-dir "$STUDY_DIR"
bash scripts/run_basic_analysis.sh --study-dir "$STUDY_DIR"
```

Si no passes cap `--input-pdb`, el projecte continua utilitzant el Trp-cage
d'exemple. Si no passes cap `--study-dir`, cada pas escriu la seva sortida a
`study_runs/`. Per tenir resultats separats per cada PDB, usa un `--study-dir`
diferent per sistema.

També pots executar tot el flux recomanat d'una sola vegada:

```bash
bash scripts/run_study_pipeline.sh \
  --input-pdb inputs/el_teu_fitxer.pdb \
  --study-dir "$STUDY_DIR"
```
