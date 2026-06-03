from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
DEFAULT_DATA_DIR = PROJECT_ROOT / "calculs" / "analisi" / "data"
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
        topology=DEFAULT_SETUP_DIR / "protein_only" / "protein_only.prmtop",
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
            "Analitza les trajectòries OpenMM de calculs/analisi/data: "
            "RMSD, RMSF, radi de gir i evolució dels ponts d'hidrogen."
        )
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=DEFAULT_DATA_DIR,
        help=f"Directori amb apo/run_* i holo/run_*. Per defecte: {DEFAULT_DATA_DIR}.",
    )
    parser.add_argument(
        "--kind",
        choices=["all", "apo", "holo"],
        default="all",
        help="Quin sistema analitzar. Per defecte: all.",
    )
    parser.add_argument(
        "--run",
        action="append",
        default=None,
        help="Run concret a analitzar, per exemple run_001. Es pot repetir.",
    )
    parser.add_argument(
        "--output-subdir",
        default=Path("/home/10034103@uvic.local/EVB-LEH/calculs/scripts/resultats"),
        help="Subcarpeta on escriure els resultats dins de cada run.",
    )
    parser.add_argument(
        "--report-interval",
        type=int,
        default=DEFAULT_REPORT_INTERVAL,
        help=f"Passos entre frames de la trajectòria. Per defecte: {DEFAULT_REPORT_INTERVAL}.",
    )
    parser.add_argument(
        "--timestep-fs",
        type=float,
        default=DEFAULT_TIMESTEP_FS,
        help=f"Timestep de la simulacio en fs. Per defecte: {DEFAULT_TIMESTEP_FS}.",
    )
    return parser.parse_args()


def write_rows(path: Path, header: list[str], rows: list[list[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def frame_times_ns(n_frames: int, report_interval: int, timestep_fs: float) -> np.ndarray:
    ns_per_frame = report_interval * timestep_fs / 1_000_000.0
    return np.arange(n_frames) * ns_per_frame


def select_atoms(topology, selector: str, label: str) -> np.ndarray:
    indices = topology.select(selector)
    if len(indices) == 0:
        raise ValueError(f"No s'han trobat àtoms per a la selecció {label}: {selector}")
    return indices


def compute_hbond_counts(md, traj) -> np.ndarray:
    try:
        hbonds_by_frame = md.wernet_nilsson(traj, periodic=True)
    except AttributeError as exc:
        raise RuntimeError(
            "Aquesta versió de mdtraj no té wernet_nilsson; cal actualitzar mdtraj "
            "o adaptar el mètode de ponts d'hidrogen."
        ) from exc
    return np.array([len(frame_hbonds) for frame_hbonds in hbonds_by_frame], dtype=int)


def save_plots(output_dir: Path, times_ns: np.ndarray, series: dict[str, object]) -> None:
    try:
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        print("Matplotlib no està instal·lat; s'escriuen només els CSV.")
        return

    plot_specs = {
        "rmsd": ("RMSD", "Distància (nm)", "rmsd.png"),
        "radius_of_gyration": ("Radi de gir", "Distància (nm)", "radius_of_gyration.png"),
        "hbonds": ("Ponts d'hidrogen totals", "Nombre de ponts d'hidrogen", "hydrogen_bonds.png"),
    }
    for key, (title, ylabel, filename) in plot_specs.items():
        if key not in series:
            continue
        plt.figure(figsize=(8, 5), dpi=300)
        plt.plot(times_ns, series[key], linewidth=1.4)
        plt.title(title)
        plt.xlabel("Temps (ns)")
        plt.ylabel(ylabel)
        plt.grid(True, linestyle="--", alpha=0.4)
        plt.tight_layout()
        plt.savefig(output_dir / filename)
        plt.close()

    if "rmsf" in series and "rmsf_residues" in series:
        residues = series["rmsf_residues"]
        rmsf = series["rmsf"]
        plt.figure(figsize=(9, 5), dpi=300)
        plt.plot(np.arange(1, len(rmsf) + 1), rmsf, linewidth=1.4)
        plt.title("RMSF")
        plt.xlabel("Residus C-α")
        plt.ylabel("Distància (nm)")
        if len(residues) <= 30:
            plt.xticks(np.arange(1, len(residues) + 1), residues, rotation=90)
        plt.grid(True, linestyle="--", alpha=0.4)
        plt.tight_layout()
        plt.savefig(output_dir / "rmsf_ca.png")
        plt.close()


def analyze_run(md, kind: SimulationKind, run_dir: Path, output_subdir: str, report_interval: int, timestep_fs: float) -> bool:
    trajectory_path = run_dir / kind.trajectory_name
    if not trajectory_path.exists():
        print(f"S'omet {run_dir}: no existeix {trajectory_path.name}")
        return False
    if not kind.topology.exists():
        raise SystemExit(f"No existeix la topologia per {kind.name}: {kind.topology}")

    output_dir = run_dir / output_subdir
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Carregant {kind.name}/{run_dir.name}: {trajectory_path}")
    traj = md.load_dcd(str(trajectory_path), top=str(kind.topology))
    times_ns = frame_times_ns(traj.n_frames, report_interval, timestep_fs)

    protein_atoms = select_atoms(traj.topology, "protein", "proteïna")
    ca_atoms = select_atoms(traj.topology, "protein and name CA", "C-α")

    protein_traj = traj.atom_slice(protein_atoms)
    protein_traj.superpose(protein_traj, frame=0)
    rmsd = md.rmsd(protein_traj, protein_traj, frame=0)
    rg = md.compute_rg(protein_traj)

    ca_traj = traj.atom_slice(ca_atoms)
    ca_traj.superpose(ca_traj, frame=0)
    rmsf = md.rmsf(ca_traj, ca_traj, frame=0)
    residues = [str(atom.residue) for atom in ca_traj.topology.atoms]

    hbond_counts = compute_hbond_counts(md, traj)

    write_rows(
        output_dir / "rmsd.csv",
        ["frame", "time_ns", "rmsd_nm"],
        [[idx, time_ns, value] for idx, (time_ns, value) in enumerate(zip(times_ns, rmsd, strict=True))],
    )
    write_rows(
        output_dir / "radius_of_gyration.csv",
        ["frame", "time_ns", "rg_nm"],
        [[idx, time_ns, value] for idx, (time_ns, value) in enumerate(zip(times_ns, rg, strict=True))],
    )
    write_rows(
        output_dir / "rmsf_ca.csv",
        ["residue", "rmsf_nm"],
        [[residue, value] for residue, value in zip(residues, rmsf, strict=True)],
    )
    write_rows(
        output_dir / "hydrogen_bonds.csv",
        ["frame", "time_ns", "n_hydrogen_bonds"],
        [[idx, time_ns, int(value)] for idx, (time_ns, value) in enumerate(zip(times_ns, hbond_counts, strict=True))],
    )

    save_plots(
        output_dir,
        times_ns,
        {
            "rmsd": rmsd,
            "radius_of_gyration": rg,
            "hbonds": hbond_counts,
            "rmsf": rmsf,
            "rmsf_residues": residues,
        },
    )

    summary = [
        f"System: {kind.name}",
        f"Run: {run_dir.name}",
        f"Topology: {kind.topology}",
        f"Trajectory: {trajectory_path}",
        f"Frames: {traj.n_frames}",
        f"Temps final (ns): {times_ns[-1]:.6f}" if len(times_ns) else "Temps final (ns): 0.000000",
        f"RMSD mitjà de la proteïna (nm): {float(np.mean(rmsd)):.6f}" if len(rmsd) else "RMSD mitjà de la proteïna (nm): 0.000000",
        f"Radi de gir mitjà (nm): {float(np.mean(rg)):.6f}" if len(rg) else "Radi de gir mitjà (nm): 0.000000",
        f"RMSF C-α màxim (nm): {float(np.max(rmsf)):.6f}" if len(rmsf) else "RMSF C-α màxim (nm): 0.000000",
        f"Ponts d'hidrogen mitjans: {float(np.mean(hbond_counts)):.6f}" if len(hbond_counts) else "Ponts d'hidrogen mitjans: 0.000000",
    ]
    (output_dir / "summary.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")
    print(f"Analisi desada a: {output_dir}")
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
        raise SystemExit(
            "Per executar l'anàlisi cal instal·lar mdtraj. Exemple: `pip install mdtraj`"
        ) from exc

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
                print(f"S'omet {run_dir}: no existeix")
                continue
            analyzed_any = analyze_run(
                md=md,
                kind=kind,
                run_dir=run_dir,
                output_subdir=args.output_subdir,
                report_interval=args.report_interval,
                timestep_fs=args.timestep_fs,
            ) or analyzed_any

    if not analyzed_any:
        raise SystemExit("No s'ha analitzat cap trajectòria.")


if __name__ == "__main__":
    main()
