from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

from openmm.app import Modeller, PDBFile

from md_common import nonprotein_residues, protein_residue_label

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_PREPARED_PDB_DIR = SCRIPT_DIR.parent / "prepared pdbs"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Pas 1: netejar i normalitzar el PDB.")
    parser.add_argument(
        "--study-dir",
        type=Path,
        default=DEFAULT_PREPARED_PDB_DIR,
        help="Directori base dels PDB preparats.",
    )
    parser.add_argument(
        "--input-pdb",
        type=Path,
        default=Path("input/leh_1nww.pdb"),
        help="Fitxer PDB d'entrada. Es pot passar qualsevol .pdb.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directori de sortida. Per defecte: STUDY_DIR/01_cleaned_pdb.",
    )
    parser.add_argument(
        "--keep-nonprotein",
        action="store_true",
        help="Conserva HETATM/aigues/holo. Per defecte es deixa nomes la proteina.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.output_dir is None:
        args.output_dir = args.study_dir / "01_cleaned_pdb"
    args.output_dir.mkdir(parents=True, exist_ok=True)

    pdb = PDBFile(str(args.input_pdb))
    modeller = Modeller(pdb.topology, pdb.positions)
    removed = []
    if not args.keep_nonprotein:
        removed = nonprotein_residues(modeller.topology)
        modeller.delete(removed)

    if modeller.topology.getNumAtoms() == 0:
        raise SystemExit(
            "No queden atoms despres de la neteja. Revisa el PDB o prova --keep-nonprotein."
        )

    output_pdb = args.output_dir / "protein_only_cleaned.pdb"
    with output_pdb.open("w", encoding="utf-8") as handle:
        PDBFile.writeFile(modeller.topology, modeller.positions, handle)

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
                f"Atoms entrada: {pdb.topology.getNumAtoms()}",
                f"Atoms sortida: {modeller.topology.getNumAtoms()}",
                f"Residus no proteics eliminats: {removed_summary}",
                f"Exemples eliminats: {removed_examples or 'cap'}",
                "Accio: normalitzacio del PDB i filtratge a proteina estandard.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"PDB netejat i desat a {output_pdb}")


if __name__ == "__main__":
    main()
