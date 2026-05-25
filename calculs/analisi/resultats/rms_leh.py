#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[2]
DEFAULT_DATA_DIR = PROJECT_ROOT / "calculs" / "analisi" / "data"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "calculs" / "analisi" / "resultats" / "grafics"
DEFAULT_SETUP_DIR = PROJECT_ROOT / "calculs" / "prepared pdbs" / "02_openmm_md_setup"
DEFAULT_REPORT_INTERVAL = 10000
DEFAULT_TIMESTEP_FS = 4.0


@dataclass(frozen=True)
class SimulationKind:
    name: str
    topology: Path
    trajectory_name: str


SIMULATION_KINDS = {
    "apo": SimulationKind(
        name="apo",
        topology=DEFAULT_SETUP_DIR / "protein_apo" / "protein_apo.prmtop",
        trajectory_name="trajectory-apo.dcd",
    ),
    "holo": SimulationKind(
        name="holo",
        topology=DEFAULT_SETUP_DIR / "protein_holo" / "protein_holo.prmtop",
        trajectory_name="trajectory-holo.dcd",
    ),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Calcula RMSD, RMSF CA i radi de gir amb els fitxers OpenMM "
            "de calculs/analisi/data."
        )
    )
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--kind", choices=["all", "apo", "holo"], default="all")
    parser.add_argument(
        "--run",
        action="append",
        default=None,
        help="Run concret a analitzar, per exemple run_001. Es pot repetir.",
    )
    parser.add_argument("--report-interval", type=int, default=DEFAULT_REPORT_INTERVAL)
    parser.add_argument("--timestep-fs", type=float, default=DEFAULT_TIMESTEP_FS)
    return parser.parse_args()


def frame_times_ns(n_frames: int, report_interval: int, timestep_fs: float) -> np.ndarray:
    ns_per_frame = report_interval * timestep_fs / 1_000_000.0
    return np.arange(n_frames, dtype=float) * ns_per_frame


def select_atoms(topology, selector: str, label: str) -> np.ndarray:
    indices = topology.select(selector)
    if len(indices) == 0:
        raise ValueError(f"No s'han trobat atoms per a {label}: {selector}")
    return indices


def write_rows(path: Path, header: list[str], rows: list[list[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def save_line_plot(output_path: Path, x_values, y_values, title: str, xlabel: str, ylabel: str) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        print("Matplotlib no esta instal-lat; s'escriuen nomes els CSV.")
        return

    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(8, 5), dpi=300)
    plt.plot(x_values, y_values, linewidth=1.4)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def analyze_run(md, kind: SimulationKind, run_dir: Path, output_dir: Path, report_interval: int, timestep_fs: float) -> bool:
    trajectory_path = run_dir / kind.trajectory_name
    if not trajectory_path.exists():
        print(f"S'omet {kind.name}/{run_dir.name}: no existeix {trajectory_path.name}")
        return False
    if not kind.topology.exists():
        raise SystemExit(f"No existeix la topologia per {kind.name}: {kind.topology}")

    run_output_dir = output_dir / kind.name / run_dir.name
    run_output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Carregant {kind.name}/{run_dir.name}")
    print(f"  Topologia: {kind.topology}")
    print(f"  Trajectoria: {trajectory_path}")

    traj = md.load_dcd(str(trajectory_path), top=str(kind.topology))
    times_ns = frame_times_ns(traj.n_frames, report_interval, timestep_fs)

    protein_atoms = select_atoms(traj.topology, "protein", "proteina")
    ca_atoms = select_atoms(traj.topology, "protein and name CA", "C-alpha")

    protein_traj = traj.atom_slice(protein_atoms)
    protein_traj.superpose(protein_traj, frame=0)
    rmsd = md.rmsd(protein_traj, protein_traj, frame=0)
    rg = md.compute_rg(protein_traj)

    ca_traj = traj.atom_slice(ca_atoms)
    ca_traj.superpose(ca_traj, frame=0)
    rmsf = md.rmsf(ca_traj, ca_traj, frame=0)
    residues = [str(atom.residue) for atom in ca_traj.topology.atoms]

    write_rows(
        run_output_dir / "rmsd.csv",
        ["frame", "time_ns", "rmsd_nm"],
        [[idx, time_ns, value] for idx, (time_ns, value) in enumerate(zip(times_ns, rmsd, strict=True))],
    )
    write_rows(
        run_output_dir / "radius_of_gyration.csv",
        ["frame", "time_ns", "rg_nm"],
        [[idx, time_ns, value] for idx, (time_ns, value) in enumerate(zip(times_ns, rg, strict=True))],
    )
    write_rows(
        run_output_dir / "rmsf_ca.csv",
        ["residue", "rmsf_nm"],
        [[residue, value] for residue, value in zip(residues, rmsf, strict=True)],
    )

    save_line_plot(run_output_dir / "rmsd.png", times_ns, rmsd, f"{kind.name} {run_dir.name}: RMSD", "Temps (ns)", "RMSD (nm)")
    save_line_plot(run_output_dir / "radius_of_gyration.png", times_ns, rg, f"{kind.name} {run_dir.name}: radi de gir", "Temps (ns)", "Rg (nm)")
    save_line_plot(run_output_dir / "rmsf_ca.png", np.arange(len(rmsf)), rmsf, f"{kind.name} {run_dir.name}: RMSF CA", "Residus CA", "RMSF (nm)")

    summary = [
        f"System: {kind.name}",
        f"Run: {run_dir.name}",
        f"Topology: {kind.topology}",
        f"Trajectory: {trajectory_path}",
        f"Frames trajectoria: {traj.n_frames}",
        f"Temps final trajectoria (ns): {times_ns[-1]:.6f}" if len(times_ns) else "Temps final trajectoria (ns): 0.000000",
        f"RMSD mitja proteina (nm): {float(np.mean(rmsd)):.6f}" if len(rmsd) else "RMSD mitja proteina (nm): 0.000000",
        f"Radi de gir mitja (nm): {float(np.mean(rg)):.6f}" if len(rg) else "Radi de gir mitja (nm): 0.000000",
        f"RMSF CA maxim (nm): {float(np.max(rmsf)):.6f}" if len(rmsf) else "RMSF CA maxim (nm): 0.000000",
    ]
    (run_output_dir / "summary.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")
    print(f"Resultats desats a: {run_output_dir}")
    return True


def iter_run_dirs(data_dir: Path, kind_name: str, requested_runs: list[str] | None) -> list[Path]:
    kind_dir = data_dir / kind_name
    if requested_runs:
        return [kind_dir / run_name for run_name in requested_runs]
    return sorted(path for path in kind_dir.glob("run_*") if path.is_dir())


def main() -> None:
    args = parse_args()
    try:
        import mdtraj as md
    except ModuleNotFoundError as exc:
        raise SystemExit("Per executar aquest script cal instal-lar mdtraj. Exemple: conda install -c conda-forge mdtraj") from exc

    args.data_dir = args.data_dir.resolve()
    args.output_dir = args.output_dir.resolve()

    kind_names = ["apo", "holo"] if args.kind == "all" else [args.kind]
    analyzed_any = False
    for kind_name in kind_names:
        kind = SIMULATION_KINDS[kind_name]
        run_dirs = iter_run_dirs(args.data_dir, kind_name, args.run)
        if not run_dirs:
            print(f"No s'han trobat runs per {kind_name} dins {args.data_dir / kind_name}")
            continue
        for run_dir in run_dirs:
            if not run_dir.exists():
                print(f"S'omet {kind.name}/{run_dir.name}: no existeix")
                continue
            analyzed_any = analyze_run(md, kind, run_dir, args.output_dir, args.report_interval, args.timestep_fs) or analyzed_any

    if not analyzed_any:
        raise SystemExit("No s'ha analitzat cap trajectoria.")


if __name__ == "__main__":
    main()
