# EVB-LEH

Exemple mínim de dinàmica molecular amb OpenMM.

## Ús

Instal·la OpenMM i executa:

```bash
python openmm_md.py --steps 5000 --output-dir outputs
```

El script crearà una caixa d'aigua, farà minimització d'energia i una simulació curta
de dinàmica molecular. Els fitxers de sortida es guardaran a `outputs/`.
