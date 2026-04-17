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

## Protocol LEH

Hi ha una guia pas a pas per plantejar un estudi de dinàmica molecular de la LEH a:

- [PROTOCOL_LEH_MD.md](/Users/gemmacollvila/github/EVB-LEH/PROTOCOL_LEH_MD.md)
