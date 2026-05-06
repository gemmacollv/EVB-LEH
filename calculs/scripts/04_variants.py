from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

from openmm.app import Modeller, PDBFile

from md_common import nonprotein_residues, protein_residue_label, write_pdb


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Pas 3: preparar una variant nomes proteina i una proteina + lligands."
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
        help="PDB protonat. Per defecte: STUDY_DIR/02_protonation/protonated.pdb.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directori de sortida. Per defecte: STUDY_DIR/03_variants.",
    )
    return parser.parse_args()


def detect_ligands(topology) -> list[str]:
    return [protein_residue_label(residue) for residue in nonprotein_residues(topology)]


def main() -> None:
    args = parse_args()
    if args.input_pdb is None:
        args.input_pdb = args.study_dir / "02_protonation" / "protonated.pdb"
    if args.output_dir is None:
        args.output_dir = args.study_dir / "03_variants"
    args.output_dir.mkdir(parents=True, exist_ok=True)

    pdb = PDBFile(str(args.input_pdb))
    modeller = Modeller(pdb.topology, pdb.positions)
    ligands = detect_ligands(pdb.topology)

    protein_only = args.output_dir / "protein_only.pdb"
    protein_plus_ligands = args.output_dir / "protein_plus_ligands.pdb"
    modeller.delete(nonprotein_residues(modeller.topology))
    if modeller.topology.getNumAtoms() == 0:
        raise SystemExit(
            "No s'ha pogut crear la variant protein_only: no hi ha atoms de proteina."
        )
    write_pdb(modeller.topology, modeller.positions, protein_only)
    write_pdb(pdb.topology, pdb.positions, protein_plus_ligands)

    ligand_counts = Counter(label.split(":", maxsplit=1)[0] for label in ligands)
    ligand_summary = (
        ", ".join(f"{name}={count}" for name, count in sorted(ligand_counts.items()))
        if ligand_counts
        else "cap"
    )

    lines = [
        f"Input PDB: {args.input_pdb}",
        f"Protein only: {protein_only}",
        f"Protein + ligands: {protein_plus_ligands}",
        f"Atoms protein_only: {modeller.topology.getNumAtoms()}",
    ]
    if ligands:
        lines.append("Residus no proteics detectats: " + ligand_summary)
        lines.append("Exemples: " + ", ".join(ligands[:20]))
        lines.append(
            "Nota: la MD per defecte usa protein_only; "
            "protein_plus_ligands requereix parametres dels lligands."
        )
    else:
        lines.append("Residus no proteics detectats: cap")
        lines.append("Les dues variants coincideixen perque no s'han detectat lligands.")

    (args.output_dir / "summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Variants desades a {args.output_dir}")


if __name__ == "__main__":
    main()