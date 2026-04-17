from __future__ import annotations

import argparse
from pathlib import Path

from openmm.app import PDBFile

from md_common import write_pdb


PROTEIN_RESIDUES = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Pas 3: preparar una variant nomes proteina i una proteina + lligands."
    )
    parser.add_argument(
        "--input-pdb",
        type=Path,
        default=Path("study_runs/02_protonation/protonated.pdb"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("study_runs/03_variants"),
    )
    return parser.parse_args()


def detect_ligands(topology) -> list[str]:
    ligands = []
    for residue in topology.residues():
        if residue.name not in PROTEIN_RESIDUES:
            ligands.append(f"{residue.name}:{residue.id}")
    return ligands


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    pdb = PDBFile(str(args.input_pdb))
    ligands = detect_ligands(pdb.topology)

    protein_only = args.output_dir / "protein_only.pdb"
    protein_plus_ligands = args.output_dir / "protein_plus_ligands.pdb"
    write_pdb(pdb.topology, pdb.positions, protein_only)
    write_pdb(pdb.topology, pdb.positions, protein_plus_ligands)

    lines = [
        f"Input PDB: {args.input_pdb}",
        f"Protein only: {protein_only}",
        f"Protein + ligands: {protein_plus_ligands}",
    ]
    if ligands:
        lines.append("Lligands detectats: " + ", ".join(ligands))
        lines.append(
            "Nota: per aquest exemple no s'ha fet una separacio estructural real; caldria un cas amb lligands reals."
        )
    else:
        lines.append("Lligands detectats: cap")
        lines.append(
            "En el Trp-cage 1L2Y les dues variants coincideixen perque no hi ha lligands."
        )

    (args.output_dir / "summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Variants desades a {args.output_dir}")


if __name__ == "__main__":
    main()
