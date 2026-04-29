from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

from openmm.app import ForceField

from md_common import load_modeller, nonprotein_residues, protein_residue_label, write_pdb


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Pas 2: afegir estat de protonacio i preparar carregues via force field."
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
        help="PDB netejat. Per defecte: STUDY_DIR/01_clean_pdb/cleaned.pdb.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directori de sortida. Per defecte: STUDY_DIR/02_protonation.",
    )
    parser.add_argument("--ph", type=float, default=7.0)
    parser.add_argument("--temperature", type=float, default=300.0)
    parser.add_argument(
        "--keep-nonprotein",
        action="store_true",
        help="Conserva residus no proteics. Per defecte s'eliminen abans de protonar.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.input_pdb is None:
        args.input_pdb = args.study_dir / "01_clean_pdb" / "cleaned.pdb"
    if args.output_dir is None:
        args.output_dir = args.study_dir / "02_protonation"
    args.output_dir.mkdir(parents=True, exist_ok=True)

    modeller = load_modeller(args.input_pdb)
    removed = []
    if not args.keep_nonprotein:
        removed = nonprotein_residues(modeller.topology)
        modeller.delete(removed)

    if modeller.topology.getNumAtoms() == 0:
        raise SystemExit(
            "No queden atoms de proteina per protonar. Revisa l'entrada o usa --keep-nonprotein."
        )

    forcefield = ForceField("amber14-all.xml", "amber14/tip3p.xml")
    modeller.addHydrogens(forcefield, pH=args.ph)

    output_pdb = args.output_dir / "protonated.pdb"
    write_pdb(modeller.topology, modeller.positions, output_pdb)

    removed_counts = Counter(residue.name for residue in removed)
    removed_summary = (
        ", ".join(f"{name}={count}" for name, count in sorted(removed_counts.items()))
        if removed_counts
        else "cap"
    )
    removed_examples = ", ".join(protein_residue_label(residue) for residue in removed[:20])
    (args.output_dir / "summary.txt").write_text(
        "\n".join(
            [
                f"Input PDB: {args.input_pdb}",
                f"Output PDB: {output_pdb}",
                f"pH: {args.ph}",
                "Force field: amber14-all.xml + amber14/tip3p.xml",
                f"Atoms totals: {modeller.topology.getNumAtoms()}",
                f"Residus no proteics eliminats: {removed_summary}",
                f"Exemples eliminats: {removed_examples or 'cap'}",
                "Nota: les carregues parcials s'assignen en crear el sistema amb el force field.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"Proteina protonada desada a {output_pdb}")


if __name__ == "__main__":
    main()
