"""Genera fitxers AMBER a partir d'un PDB extern.

Aquest script:
1. Llegeix un PDB des de disc.
2. Construeix el sistema amb force field AMBER.
3. Exporta els fitxers `prmtop` i `inpcrd`.

No incrusta cap PDB dins del codi.
"""

from __future__ import annotations

from pathlib import Path

import parmed as pmd
from openmm.app import ForceField, HBonds, NoCutoff, PDBFile, PME


# Fitxer d'entrada i directori de sortida
INPUT_PDB = Path("study_runs/02_protonation/protonated.pdb")
OUTPUT_DIR = Path("study_runs/amber_from_pdb")

# Force field AMBER
PROTEIN_FORCEFIELD = "amber14-all.xml"
WATER_FORCEFIELD = "amber14/tip3p.xml"

# Parametres del sistema
NONBONDED_CUTOFF_NM = 1.0
USE_CONSTRAINTS = True


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    pdb = PDBFile(str(INPUT_PDB))
    forcefield = ForceField(PROTEIN_FORCEFIELD, WATER_FORCEFIELD)

    has_periodic_box = pdb.topology.getPeriodicBoxVectors() is not None
    create_system_kwargs = {
        "topology": pdb.topology,
        "nonbondedMethod": PME if has_periodic_box else NoCutoff,
        "constraints": HBonds if USE_CONSTRAINTS else None,
        "rigidWater": USE_CONSTRAINTS,
    }
    if has_periodic_box:
        create_system_kwargs["nonbondedCutoff"] = NONBONDED_CUTOFF_NM

    system = forcefield.createSystem(**create_system_kwargs)

    structure = pmd.openmm.load_topology(pdb.topology, system, xyz=pdb.positions)

    prmtop_path = OUTPUT_DIR / f"{INPUT_PDB.stem}.prmtop"
    inpcrd_path = OUTPUT_DIR / f"{INPUT_PDB.stem}.inpcrd"
    summary_path = OUTPUT_DIR / "summary.txt"

    structure.save(str(prmtop_path), overwrite=True)
    structure.save(str(inpcrd_path), overwrite=True)

    summary_path.write_text(
        "\n".join(
            [
                f"Input PDB: {INPUT_PDB}",
                f"PRMTOP: {prmtop_path}",
                f"INPCRD: {inpcrd_path}",
                f"Protein force field: {PROTEIN_FORCEFIELD}",
                f"Water force field: {WATER_FORCEFIELD}",
                f"Periodic box detected: {has_periodic_box}",
                f"Nonbonded method: {'PME' if has_periodic_box else 'NoCutoff'}",
                f"Constraints: {'HBonds' if USE_CONSTRAINTS else 'None'}",
                f"Atoms totals: {pdb.topology.getNumAtoms()}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(f"PDB d'entrada: {INPUT_PDB}")
    print(f"Fitxer AMBER topology: {prmtop_path}")
    print(f"Fitxer AMBER coordinates: {inpcrd_path}")


if __name__ == "__main__":
    main()
