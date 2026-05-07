from __future__ import annotations

import argparse
import csv
from pathlib import Path



def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Pas 5: analisi basica de RMSD, RMSF i radi de gir."
    )
    parser.add_argument(
        "--study-dir",
        type=Path,
        default=Path("study_runs"),
        help="Directori base de l'estudi.",
    )
    parser.add_argument(
        "--topology",
        type=Path,
        default=None,
        help="Topologia/PDB final. Per defecte: STUDY_DIR/04_md/production_final.pdb.",
    )
    parser.add_argument(
        "--trajectory",
        type=Path,
        default=None,
        help="Trajectoria DCD. Per defecte: STUDY_DIR/04_md/production_trajectory.dcd.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directori de sortida. Per defecte: STUDY_DIR/05_analysis.",
    )
    return parser.parse_args()


def write_timeseries_csv(path: Path, header: tuple[str, str], values) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        for idx, value in enumerate(values):
            writer.writerow([idx, value])


def write_residue_csv(path: Path, residues, values) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["residue", "rmsf_nm"])
        for residue, value in zip(residues, values, strict=True):
            writer.writerow([residue, value])


def main() -> None:
    args = parse_args()
    if args.topology is None:
        args.topology = args.study_dir / "04_md" / "production_final.pdb"
    if args.trajectory is None:
        args.trajectory = args.study_dir / "04_md" / "production_trajectory.dcd"
    if args.output_dir is None:
        args.output_dir = args.study_dir / "05_analysis"
    args.output_dir.mkdir(parents=True, exist_ok=True)

    try:
        import mdtraj as md
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "Per executar l'analisi cal instal·lar mdtraj. Exemple: `pip install mdtraj`"
        ) from exc

    traj = md.load_dcd(str(args.trajectory), top=str(args.topology))
    protein_atoms = traj.topology.select("protein")
    ca_atoms = traj.topology.select("protein and name CA")

    aligned = traj.atom_slice(protein_atoms)
    aligned.superpose(aligned, 0)

    rmsd = md.rmsd(aligned, aligned, 0)
    rg = md.compute_rg(aligned)

    ca_traj = traj.atom_slice(ca_atoms)
    ca_traj.superpose(ca_traj, 0)
    rmsf = md.rmsf(ca_traj, ca_traj, 0)
    residues = [str(atom.residue) for atom in ca_traj.topology.atoms]

    write_timeseries_csv(args.output_dir / "rmsd.csv", ("frame", "rmsd_nm"), rmsd)
    write_timeseries_csv(args.output_dir / "radius_of_gyration.csv", ("frame", "rg_nm"), rg)
    write_residue_csv(args.output_dir / "rmsf_ca.csv", residues, rmsf)

    mean_rmsd = float(rmsd.mean()) if len(rmsd) else 0.0
    mean_rg = float(rg.mean()) if len(rg) else 0.0
    max_rmsf = float(rmsf.max()) if len(rmsf) else 0.0
    (args.output_dir / "summary.txt").write_text(
        "\n".join(
            [
                f"Topology: {args.topology}",
                f"Trajectory: {args.trajectory}",
                f"Frames: {traj.n_frames}",
                f"RMSD mitja (nm): {mean_rmsd:.6f}",
                f"Radi de gir mitja (nm): {mean_rg:.6f}",
                f"RMSF CA maxima (nm): {max_rmsf:.6f}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"Analisi desada a {args.output_dir}")


if __name__ == "__main__":
    main()
