# Protocol Pas A Pas per a la Dinàmica Molecular del Trp-cage 1L2Y

Aquest document descriu un protocol curt i reproduïble per estudiar amb OpenMM una proteïna real molt petita: el Trp-cage miniprotein (`PDB: 1L2Y`, model 1).

És un sistema molt útil per:

- provar pipelines de dinàmica molecular
- validar preparació, minimització i equilibratge
- fer assajos ràpids abans de passar a proteïnes més grans

## 1. Sistema de partida

La proteïna d'entrada és:

- Trp-cage miniprotein
- `1L2Y`, model 1
- 20 residus

El fitxer inicial del repositori és:

- [trp_cage_1l2y_model1.pdb](/Users/gemmacollvila/github/EVB-LEH/inputs/trp_cage_1l2y_model1.pdb)

## 2. Flux de treball

Per seguir un esquema més semblant a un estudi real, el flux recomanat és:

1. netejar PDB
2. afegir estat de protonació i preparar carregues
3. preparar una variant `prot + lligands` i una `nomes prot`
4. executar la dinàmica molecular
5. analitzar RMSD, RMSF i radi de gir

## 3. Neteja del PDB

Script:

- [01_clean_pdb.py](/Users/gemmacollvila/github/EVB-LEH/scripts/01_clean_pdb.py)

Llançament:

```bash
bash scripts/run_clean_pdb.sh
```

Línies de codi:

```bash
python scripts/01_clean_pdb.py \
  --input-pdb inputs/trp_cage_1l2y_model1.pdb \
  --output-dir study_runs/01_clean_pdb
```

## 4. Estat de protonació i carregues

En aquest pas es fa:

- addició d'hidrògens a `pH 7.0`
- definició de l'estat de protonació
- preparació del sistema perquè el force field assigni carregues parcials

Script:

- [02_add_protonation.py](/Users/gemmacollvila/github/EVB-LEH/scripts/02_add_protonation.py)

Llançament:

```bash
bash scripts/run_add_protonation.sh
```

Línies de codi:

```bash
python scripts/02_add_protonation.py \
  --input-pdb study_runs/01_clean_pdb/cleaned.pdb \
  --output-dir study_runs/02_protonation \
  --ph 7.0 \
  --temperature 300.0
```

## 5. Variants: proteina + lligands i nomes proteina

Script:

- [03_prepare_variants.py](/Users/gemmacollvila/github/EVB-LEH/scripts/03_prepare_variants.py)

Llançament:

```bash
bash scripts/run_prepare_variants.sh
```

Línies de codi:

```bash
python scripts/03_prepare_variants.py \
  --input-pdb study_runs/02_protonation/protonated.pdb \
  --output-dir study_runs/03_variants
```

En el Trp-cage `1L2Y` no hi ha lligands, així que les dues variants acaben sent equivalents.

## 6. Dinàmica molecular

Objectiu:

- solvatació
- minimització
- equilibratge NVT
- equilibratge NPT
- producció

Script:

- [04_run_md.py](/Users/gemmacollvila/github/EVB-LEH/scripts/04_run_md.py)

Llançament:

```bash
bash scripts/run_md_study.sh
```

Línies de codi:

```bash
python scripts/04_run_md.py \
  --input-pdb study_runs/03_variants/protein_only.pdb \
  --output-dir study_runs/04_md \
  --ph 7.0 \
  --padding-nm 1.2 \
  --ionic-strength 0.15 \
  --temperature 300.0 \
  --friction 1.0 \
  --timestep-fs 2.0 \
  --report-interval 500 \
  --minimize-iterations 5000 \
  --nvt-steps 25000 \
  --npt-steps 100000 \
  --production-steps 250000
```

## 7. Analisi

Objectiu:

- RMSD
- RMSF
- radi de gir

Script:

- [05_analyze_basic.py](/Users/gemmacollvila/github/EVB-LEH/scripts/05_analyze_basic.py)

Llançament:

```bash
bash scripts/run_basic_analysis.sh
```

Línies de codi:

```bash
python scripts/05_analyze_basic.py \
  --topology study_runs/04_md/production_final.pdb \
  --trajectory study_runs/04_md/production_production_trajectory.dcd \
  --output-dir study_runs/05_analysis
```

Nota:

- aquest pas requereix `mdtraj`
- en aquest entorn no està instal·lat actualment

## 8. Ordre recomanat

Executa els passos en aquest ordre:

1. `bash scripts/run_clean_pdb.sh`
2. `bash scripts/run_add_protonation.sh`
3. `bash scripts/run_prepare_variants.sh`
4. `bash scripts/run_md_study.sh`
5. `bash scripts/run_basic_analysis.sh`

Versió completa en línia de comandes:

```bash
python scripts/01_clean_pdb.py --input-pdb inputs/trp_cage_1l2y_model1.pdb --output-dir study_runs/01_clean_pdb
python scripts/02_add_protonation.py --input-pdb study_runs/01_clean_pdb/cleaned.pdb --output-dir study_runs/02_protonation --ph 7.0 --temperature 300.0
python scripts/03_prepare_variants.py --input-pdb study_runs/02_protonation/protonated.pdb --output-dir study_runs/03_variants
python scripts/04_run_md.py --input-pdb study_runs/03_variants/protein_only.pdb --output-dir study_runs/04_md --ph 7.0 --padding-nm 1.2 --ionic-strength 0.15 --temperature 300.0 --friction 1.0 --timestep-fs 2.0 --report-interval 500 --minimize-iterations 5000 --nvt-steps 25000 --npt-steps 100000 --production-steps 250000
python scripts/05_analyze_basic.py --topology study_runs/04_md/production_final.pdb --trajectory study_runs/04_md/production_production_trajectory.dcd --output-dir study_runs/05_analysis
```

## 9. Sortides

Cada etapa escriu els seus resultats dins `study_runs/`:

- PDBs intermedis
- fitxers `csv` amb energies i temperatura
- trajectòries `dcd`
- resum curt de cada pas

## 10. Què analitzar després

Per aquest sistema petit, les anàlisis més útils són:

- RMSD global
- RMSF per residu
- radi de gir
- estabilitat estructural general

## 11. Utilitat del sistema

El Trp-cage és especialment útil perquè:

- és petit i ràpid de simular
- és una proteïna real
- permet verificar que el pipeline funciona de punta a punta
- és ideal com a exemple docent i de desenvolupament
