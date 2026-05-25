# Dinàmica molecular amb OpenMM sobre LEH

Aquest fitxer descriu el protocol sobre l’estudi de la dinàmica molecular de la _Limonene-1,2-epoxide hydrolase_ (LEH) (PDB: 1NWW)

Primer de tot s'ha d'importar el PDB de la proteïna que estudiarem, en el nostre cas la LEH (1NWW): https://www.rcsb.org/structure/1NWW
Un cop importat el PDB, els passos que seguirem per tal de realitzar la dinàmica molecular són els següents:

1. Natejar PDB
2. Afegir protonació i càrregues
3. Preparar prtoeïna amb i sense holo
4. Executar _Molecular Dynamics_ (MD)
5. Analitzar RMSD, RMSF i radi de gir

## 1. Nateja del PDB

El següent codi (01_clean_pdb.py) agafa un fitxer .pddb i elimina HOH, MES i HPN, on nommés deixa la part proteica per poder protonar, solvatar i simular.

HOH = aigües cristal·logràfiques
MES = buffer
HPN = holo/inhibidor o molècula no estàndard


```bash
python calculs/run_clean_pdb.sh --input-pdb inputs/leh_1nww.pdb"
```

## 2. Adició de protons (protonació)

Addició d'hidrògens a pH 7.0, definició de l'estat de protonació i preparació del sistema perquè el force field (amberxl...) assigni càrregues parcials

```bash
bash scripts/run_add_protonation.sh"
```

## 3. Preparar prtoeïna amb i sense holo
```bash
bash scripts/run_prepare_variants.sh"
```

## 4. Executar Molecular Dynamics (MD)
```bash
bash scripts/run_md_study.sh"
```

## 5. analitzar RMSD, RMSF i radi de gir
```bash
bash scripts/run_basic_analysis.sh"
```


