# Dinàmica molecular amb OpenMM sobre LEH

ProtocolAquest fitxer sobre l’estudi de la dinàmica molecular de la LEH

Primer de tot s'ha d'importar el PDB de la proteïna que estudiarem, en el nostre cas la LEH (1NWW): https://www.rcsb.org/structure/1NWW

1. Natejar PDB
2. Afegir protonació i carregues
3. Preparar prtoeïna amb i sense lligands
4. Executar Molecular Dynamics (MD)
5. analitzar RMSD, RMSF i radi de gir

# Nateja del PDB

El següent codi (01_clean_pdb.py) agafa un fitxer .pddb i elimina HOH, MES i HPN, on nommés deixa la part proteica per poder protonar, solvatar i simular.

HOH = aigües cristal·logràfiques
MES = buffer
HPN = lligand/inhibidor o molècula no estàndard


Instal·la OpenMM i executa:

```bash
python calculs/01_claean_pdb.py --input-pdb inputs/leh_1nww.pdb --steps 5000 --output-dir outputs
```

Fes-ho en petits scripts de python

dia 2 17/04/2026
Passos seguits

1. netejar 
2. protonar
3. construir els inpcrd i prm top

fes la protonació amb de el pdb protonated

scripts molt fàcil en python de cada pas




bash scripts/run_clean_pdb.sh --input-pdb inputs/el_teu_fitxer.pdb --study-dir "$STUDY_DIR"

bash scripts/run_add_protonation.sh --study-dir "$STUDY_DIR"

bash scripts/run_prepare_variants.sh --study-dir "$STUDY_DIR"

bash scripts/run_md_study.sh --study-dir "$STUDY_DIR"

bash scripts/run_basic_analysis.sh --study-dir "$STUDY_DIR"