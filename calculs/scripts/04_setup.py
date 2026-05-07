"""Prepara un sistema OpenMM i exporta fitxers AMBER.

Aquest script fa nomes el setup:
1. Llegeix el PDB preparat de l'estudi.
2. Afegeix hidrogens, solvent i ions.
3. Exporta `system.prmtop` i `system.inpcrd`.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import parmed as pmd
from openmm import unit
from openmm.app import PME

from md_common import create_forcefield, prepare_modeller, write_pdb


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepara un sistema solvatar i exporta fitxers AMBER."
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
        help="PDB d'entrada. Per defecte: STUDY_DIR/03_variants/protein_only.pdb.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directori de sortida. Per defecte: STUDY_DIR/amber_from_pdb.",
    )
    parser.add_argument("--prefix", default="system", help="Prefix dels fitxers AMBER.")
    parser.add_argument("--ph", type=float, default=7.0)
    parser.add_argument("--padding-nm", type=float, default=1.2)
    parser.add_argument("--ionic-strength", type=float, default=0.15, help="M")
    parser.add_argument("--temperature", type=float, default=300.0)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.input_pdb is None:
        args.input_pdb = args.study_dir / "03_variants" / "protein_only.pdb"
    if args.output_dir is None:
        args.output_dir = args.study_dir / "amber_from_pdb"
    args.input_pdb = args.input_pdb.resolve()
    args.output_dir = args.output_dir.resolve()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    modeller = prepare_modeller(
        input_pdb=args.input_pdb,
        padding_nm=args.padding_nm,
        ionic_strength_molar=args.ionic_strength,
        ph=args.ph,
    )

    # Export without constraints so ParmEd can write all AMBER bond parameters.
    forcefield = create_forcefield()
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=None,
        rigidWater=False,
    )
    structure = pmd.openmm.load_topology(
        modeller.topology,
        system,
        xyz=modeller.positions,
    )

    prmtop_path = args.output_dir / f"{args.prefix}.prmtop"
    inpcrd_path = args.output_dir / f"{args.prefix}.inpcrd"
    solvated_pdb_path = args.output_dir / f"{args.prefix}_solvated.pdb"
    summary_path = args.output_dir / "summary.txt"

    structure.save(str(prmtop_path), overwrite=True)
    structure.save(str(inpcrd_path), overwrite=True)
    write_pdb(modeller.topology, modeller.positions, solvated_pdb_path)

    box_vectors = modeller.topology.getPeriodicBoxVectors()
    box_vectors_nm = [vector.value_in_unit(unit.nanometer) for vector in box_vectors]
    summary_path.write_text(
        "\n".join(
            [
                f"Input PDB: {args.input_pdb}",
                f"PRMTOP: {prmtop_path}",
                f"INPCRD: {inpcrd_path}",
                f"Solvated PDB: {solvated_pdb_path}",
                "Protein force field: amber14-all.xml",
                "Water force field: amber14/tip3p.xml",
                "Nonbonded method: PME",
                "Constraints: None during AMBER export; simulation can apply HBonds.",
                f"pH: {args.ph}",
                f"Padding (nm): {args.padding_nm}",
                f"Ionic strength (M): {args.ionic_strength}",
                f"Temperature (K): {args.temperature}",
                f"Atoms totals: {modeller.topology.getNumAtoms()}",
                f"Box vectors (nm): {box_vectors_nm}",
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
