from __future__ import annotations

import argparse
from pathlib import Path

from openmm.app import PDBFile


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Pas 1: netejar i normalitzar el PDB.")
    parser.add_argument(
        "--input-pdb",
        type=Path,
        default=Path("inputs/trp_cage_1l2y_model1.pdb"),
        help="Fitxer PDB d'entrada. Es pot passar qualsevol .pdb.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("study_runs/01_clean_pdb"),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    pdb = PDBFile(str(args.input_pdb))
    output_pdb = args.output_dir / "cleaned.pdb"
    with output_pdb.open("w", encoding="utf-8") as handle:
        PDBFile.writeFile(pdb.topology, pdb.positions, handle)

    (args.output_dir / "summary.txt").write_text(
        "\n".join(
            [
                f"Input PDB: {args.input_pdb}",
                f"Output PDB: {output_pdb}",
                f"Atoms totals: {pdb.topology.getNumAtoms()}",
                "Accio: normalitzacio del PDB d'entrada.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"PDB netejat i desat a {output_pdb}")


if __name__ == "__main__":
    main()
