from __future__ import annotations

import argparse
from pathlib import Path

from md_common import DEFAULT_INPUT_PDB, prepare_modeller, write_pdb


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Pas 1: preparar, protonar i solvatar el sistema.")
    parser.add_argument("--input-pdb", type=Path, default=DEFAULT_INPUT_PDB)
    parser.add_argument("--output-dir", type=Path, default=Path("runs/01_prepare"))
    parser.add_argument("--padding-nm", type=float, default=1.2)
    parser.add_argument("--ionic-strength", type=float, default=0.15, help="M")
    parser.add_argument("--ph", type=float, default=7.0)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    modeller = prepare_modeller(
        input_pdb=args.input_pdb,
        padding_nm=args.padding_nm,
        ionic_strength_molar=args.ionic_strength,
        ph=args.ph,
    )
    out_pdb = args.output_dir / "prepared_solvated.pdb"
    write_pdb(modeller.topology, modeller.positions, out_pdb)

    summary = args.output_dir / "summary.txt"
    summary.write_text(
        "\n".join(
            [
                f"Input PDB: {args.input_pdb}",
                f"Output PDB: {out_pdb}",
                f"pH: {args.ph}",
                f"Padding (nm): {args.padding_nm}",
                f"Ionic strength (M): {args.ionic_strength}",
                f"Atoms totals: {modeller.topology.getNumAtoms()}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"Sistema preparat a {out_pdb}")


if __name__ == "__main__":
    main()
