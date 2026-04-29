"""Genera fitxers AMBER a partir d'un PDB extern.

Aquest script:
1. Llegeix un PDB des de disc.
2. Construeix el sistema amb force field AMBER.
3. Exporta els fitxers `prmtop` i `inpcrd`.

No incrusta cap PDB dins del codi.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import parmed as pmd
from openmm import unit
from openmm.app import ForceField, HBonds, NoCutoff, PDBFile, PME


# Force field AMBER
PROTEIN_FORCEFIELD = "amber14-all.xml"
WATER_FORCEFIELD = "amber14/tip3p.xml"

# Parametres del sistema
NONBONDED_CUTOFF_NM = 1.0
USE_CONSTRAINTS = True


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Genera fitxers AMBER a partir d'un PDB extern."
    )
    parser.add_argument(
        "--study-dir",
        type=Path,
        default=Path("study_runs"),
        help="Directori base de l'estudi.",
    )
    parser.add_argument(
        "--input-pdb",
        type=Path,
        default=None,
        help="PDB d'entrada. Per defecte: STUDY_DIR/02_protonation/protonated.pdb.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directori de sortida. Per defecte: STUDY_DIR/amber_from_pdb.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.input_pdb is None:
        args.input_pdb = args.study_dir / "02_protonation" / "protonated.pdb"
    if args.output_dir is None:
        args.output_dir = args.study_dir / "amber_from_pdb"
    args.output_dir.mkdir(parents=True, exist_ok=True)

    pdb = PDBFile(str(args.input_pdb))
    forcefield = ForceField(PROTEIN_FORCEFIELD, WATER_FORCEFIELD)

    has_periodic_box = pdb.topology.getPeriodicBoxVectors() is not None
    create_system_kwargs = {
        "topology": pdb.topology,
        "nonbondedMethod": PME if has_periodic_box else NoCutoff,
        "constraints": HBonds if USE_CONSTRAINTS else None,
        "rigidWater": USE_CONSTRAINTS,
    }
    if has_periodic_box:
        create_system_kwargs["nonbondedCutoff"] = NONBONDED_CUTOFF_NM * unit.nanometer

    system = forcefield.createSystem(**create_system_kwargs)

    structure = pmd.openmm.load_topology(pdb.topology, system, xyz=pdb.positions)

    prmtop_path = args.output_dir / f"{args.input_pdb.stem}.prmtop"
    inpcrd_path = args.output_dir / f"{args.input_pdb.stem}.inpcrd"
    summary_path = args.output_dir / "summary.txt"

    structure.save(str(prmtop_path), overwrite=True)
    structure.save(str(inpcrd_path), overwrite=True)

    summary_path.write_text(
        "\n".join(
            [
                f"Input PDB: {args.input_pdb}",
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

    print(f"PDB d'entrada: {args.input_pdb}")
    print(f"Fitxer AMBER topology: {prmtop_path}")
    print(f"Fitxer AMBER coordinates: {inpcrd_path}")


if __name__ == "__main__":
    main()
