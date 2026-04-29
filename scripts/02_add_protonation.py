from __future__ import annotations

import argparse
from pathlib import Path

from openmm.app import ForceField

from md_common import build_system, load_modeller, write_pdb


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
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.input_pdb is None:
        args.input_pdb = args.study_dir / "01_clean_pdb" / "cleaned.pdb"
    if args.output_dir is None:
        args.output_dir = args.study_dir / "02_protonation"
    args.output_dir.mkdir(parents=True, exist_ok=True)

    modeller = load_modeller(args.input_pdb)
    forcefield = ForceField("amber14-all.xml", "amber14/tip3p.xml")
    modeller.addHydrogens(forcefield, pH=args.ph)

    output_pdb = args.output_dir / "protonated.pdb"
    write_pdb(modeller.topology, modeller.positions, output_pdb)

    # system = build_system(modeller.topology, args.temperature, with_barostat=False)
    # (args.output_dir / "summary.txt").write_text(
    #     "\n".join(
    #         [
    #             f"Input PDB: {args.input_pdb}",
    #             f"Output PDB: {output_pdb}",
    #             f"pH: {args.ph}",
    #             "Force field: amber14-all.xml + amber14/tip3p.xml",
    #             f"Atoms totals: {modeller.topology.getNumAtoms()}",
    #             f"Forces al sistema: {system.getNumForces()}",
    #             "Nota: les carregues parcials s'assignen en crear el sistema amb el force field.",
    #         ]
    #     )
    #     + "\n",
    #     encoding="utf-8",
    # )
    print(f"Proteina protonada desada a {output_pdb}")


if __name__ == "__main__":
    main()
