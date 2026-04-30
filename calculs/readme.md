# Dinàmica molecular amb OpenMM sobre LEH

ProtocolAquest fitxer sobre l’estudi de la dinàmica molecular de la LEH

Primer de tot s'ha d'importar el PDB de la proteïna que estudiarem, en el nostre cas la LEH (1NWW): https://www.rcsb.org/structure/1NWW

1. Natejar PDB
2. Afegir protonació i carregues
3. Preparar prtoeïna amb i sense lligands
4. Executar Molecular Dynamics (MD)
5. analitzar RMSD, RMSF i radi de gir

# 1. Nateja del PDB

El següent codi (01_clean_pdb.py) agafa un fitxer .pddb i elimina HOH, MES i HPN, on nommés deixa la part proteica per poder protonar, solvatar i simular.

HOH = aigües cristal·logràfiques
MES = buffer
HPN = lligand/inhibidor o molècula no estàndard


```bash
python calculs/run_clean_pdb.sh --input-pdb inputs/leh_1nww.pdb --study-dir "$STUDY_DIR"
```

# 2. Adició de protons (protonació)
```bash
bash scripts/run_add_protonation.sh --study-dir "$STUDY_DIR"
```

# 3. Preparar prtoeïna amb i sense lligands
```bash
bash scripts/run_prepare_variants.sh --study-dir "$STUDY_DIR"
```

# 4. Executar Molecular Dynamics (MD)
```bash
bash scripts/run_md_study.sh --study-dir "$STUDY_DIR"
```

# 5. analitzar RMSD, RMSF i radi de gir
```bash
bash scripts/run_basic_analysis.sh --study-dir "$STUDY_DIR"```





