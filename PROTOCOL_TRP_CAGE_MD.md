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

## 2. Preparació del sistema

En aquest projecte, la preparació fa:

- càrrega del PDB
- addició d'hidrògens a `pH 7.0`
- solvatació en una caixa periòdica
- aigua `TIP3P`
- neutralització i força iònica per defecte

Script:

- [01_prepare_system.py](/Users/gemmacollvila/github/EVB-LEH/scripts/01_prepare_system.py)

Llançament:

```bash
bash scripts/run_01_prepare.sh
```

## 3. Minimització d'energia

Objectiu:

- eliminar contactes estèrics
- relaxar aigua i ions
- obtenir una geometria inicial estable

Script:

- [02_minimize.py](/Users/gemmacollvila/github/EVB-LEH/scripts/02_minimize.py)

Llançament:

```bash
bash scripts/run_02_minimize.sh
```

## 4. Equilibratge NVT

Objectiu:

- estabilitzar temperatura
- començar relaxació dinàmica sense canvi de volum

Script:

- [03_equilibrate_nvt.py](/Users/gemmacollvila/github/EVB-LEH/scripts/03_equilibrate_nvt.py)

Llançament:

```bash
bash scripts/run_03_nvt.sh
```

## 5. Equilibratge NPT

Objectiu:

- estabilitzar densitat i volum
- adaptar el solvent al sistema

Script:

- [04_equilibrate_npt.py](/Users/gemmacollvila/github/EVB-LEH/scripts/04_equilibrate_npt.py)

Llançament:

```bash
bash scripts/run_04_npt.sh
```

## 6. Producció

Objectiu:

- generar una trajectòria útil per analitzar estabilitat i flexibilitat

Script:

- [05_production.py](/Users/gemmacollvila/github/EVB-LEH/scripts/05_production.py)

Llançament:

```bash
bash scripts/run_05_production.sh
```

## 7. Ordre recomanat

Executa els passos en aquest ordre:

1. `bash scripts/run_01_prepare.sh`
2. `bash scripts/run_02_minimize.sh`
3. `bash scripts/run_03_nvt.sh`
4. `bash scripts/run_04_npt.sh`
5. `bash scripts/run_05_production.sh`

## 8. Sortides

Cada etapa escriu els seus resultats dins `runs/`:

- PDBs intermedis
- fitxers `csv` amb energies i temperatura
- trajectòries `dcd`
- resum curt de cada pas

## 9. Què analitzar després

Per aquest sistema petit, les anàlisis més útils són:

- RMSD global
- RMSF per residu
- estabilitat estructural
- evolució de l'energia potencial
- compactació de la proteïna

## 10. Utilitat del sistema

El Trp-cage és especialment útil perquè:

- és petit i ràpid de simular
- és una proteïna real
- permet verificar que el pipeline funciona de punta a punta
- és ideal com a exemple docent i de desenvolupament
