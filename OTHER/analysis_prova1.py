from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[2]
DEFAULT_DATA_DIR = PROJECT_ROOT / "calculs" / "analisi" / "data"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "calculs" / "analisi" / "resultats"
DEFAULT_SETUP_DIR = PROJECT_ROOT / "calculs" / "prepared pdbs" / "02_openmm_md_setup"
DEFAULT_REPORT_INTERVAL = 10000
DEFAULT_TIMESTEP_FS = 4.0


@dataclass(frozen=True)
class SimulationKind:
    name: str
    topology: Path
    trajectory_name: str
    log_name: str


SIMULATION_KINDS = {
    "apo": SimulationKind(
        name="apo",
        topology=DEFAULT_SETUP_DIR / "protein_only" / "protein_only.prmtop",
        trajectory_name="trajectory-apo.dcd",
        log_name="log-apo.txt",
    ),
    "holo": SimulationKind(
        name="holo",
        topology=DEFAULT_SETUP_DIR / "protein_holo" / "protein_holo.prmtop",
        trajectory_name="trajectory-holo.dcd",
        log_name="log-holo.txt",
    ),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Analitza els resultats OpenMM de calculs/analisi/data. "
            "Calcula RMSD, RMSF CA, radi de gir, ponts d'hidrogen i resumeix el log."
        )
    )
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--kind", choices=["all", "apo", "holo"], default="all")
    parser.add_argument("--run", action="append", default=None)
    parser.add_argument("--report-interval", type=int, default=DEFAULT_REPORT_INTERVAL)
    parser.add_argument("--timestep-fs", type=float, default=DEFAULT_TIMESTEP_FS)
    parser.add_argument("--skip-hbonds", action="store_true")
    return parser.parse_args()


def frame_times_ns(n_frames: int, report_interval: int, timestep_fs: float) -> np.ndarray:
    ns_per_frame = report_interval * timestep_fs / 1_000_000.0
    return np.arange(n_frames) * ns_per_frame


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


def read_log_table(log_path: Path) -> list[dict[str, float]]:
    if not log_path.exists():
        return []
    rows: list[dict[str, float]] = []
    with log_path.open("r", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if not header:
            return []
        columns = [column.strip().lstrip("#").strip("\"") for column in header]
        for raw_row in reader:
            if raw_row:
                rows.append(dict(zip(columns, [float(value) for value in raw_row], strict=True)))
    return rows


def write_log_outputs(log_path: Path, output_dir: Path, report_interval: int, timestep_fs: float) -> list[str]:
    rows = read_log_table(log_path)
    if not rows:
        return [f"Log: no trobat o buit ({log_path})"]
    output_rows = []
    for row in rows:
        step = row.get("Step", 0.0)
        time_ns = step * timestep_fs / 1_000_000.0
        output_rows.append([
            int(step),
            time_ns,
            row.get("Potential Energy (kJ/mole)", ""),
            row.get("Temperature (K)", ""),
        ])
    write_rows(output_dir / "thermo.csv", ["step", "time_ns", "potential_energy_kj_mol", "temperature_k"], output_rows)
    steps = [row.get("Step", 0.0) for row in rows]
    temperatures = [row.get("Temperature (K)", np.nan) for row in rows]
    potentials = [row.get("Potential Energy (kJ/mole)", np.nan) for row in rows]
    return [
        f"Log: {log_path}",
        f"Log frames: {len(rows)}",
        f"Ultim pas log: {int(max(steps)) if steps else 0}",
        f"Interval report assumit: {report_interval} passos",
        f"Temperatura mitjana (K): {float(np.nanmean(temperatures)):.6f}",
        f"Energia potencial mitjana (kJ/mol): {float(np.nanmean(potentials)):.6f}",
    ]


def compute_hbond_counts(md, traj) -> np.ndarray:
    try:
        hbonds_by_frame = md.wernet_nilsson(traj, periodic=True)
    except Exception:
        hbonds_by_frame = md.wernet_nilsson(traj, periodic=False)
    return np.array([len(frame_hbonds) for frame_hbonds in hbonds_by_frame], dtype=int)


def save_plots(output_dir: Path, times_ns: np.ndarray, series: dict[str, object]) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        print("Matplotlib no esta instal-lat; s'escriuen nomes els CSV.")
        return
    plot_specs = {
        "rmsd": ("RMSD proteina", "Temps (ns)", "RMSD (nm)", "rmsd.png"),
        "radius_of_gyration": ("Radi de gir", "Temps (ns)", "Rg (nm)", "radius_of_gyration.png"),
        "hbonds": ("Ponts d'hidrogen", "Temps (ns)", "Nombre de ponts", "hydrogen_bonds.png"),
    }
    for key, (title, xlabel, ylabel, filename) in plot_specs.items():
        if key not in series:
            continue
        plt.figure(figsize=(8, 5), dpi=300)
        plt.plot(times_ns, series[key], linewidth=1.4)
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.grid(True, linestyle="--", alpha=0.4)
        plt.tight_layout()
        plt.savefig(output_dir / filename)
        plt.close()
    if "rmsf" in series and "rmsf_residues" in series:
        residues = series["rmsf_residues"]
        rmsf = series["rmsf"]
        x = np.arange(len(rmsf))
        plt.figure(figsize=(9, 5), dpi=300)
        plt.plot(x, rmsf, linewidth=1.4)
        plt.title("RMSF C-alpha")
        plt.xlabel("Residu")
        plt.ylabel("RMSF (nm)")
        if len(residues) <= 40:
            plt.xticks(x, residues, rotation=90)
        plt.grid(True, linestyle="--", alpha=0.4)
        plt.tight_layout()
        plt.savefig(output_dir / "rmsf_ca.png")
        plt.close()


def analyze_run(md, kind: SimulationKind, run_dir: Path, output_dir: Path, report_interval: int, timestep_fs: float, skip_hbonds: bool) -> bool:
    trajectory_path = run_dir / kind.trajectory_name
    log_path = run_dir / kind.log_name
    if not trajectory_path.exists():
        print(f"S'omet {kind.name}/{run_dir.name}: no existeix {trajectory_path.name}")
        return False
    if not kind.topology.exists():
        raise SystemExit(f"No existeix la topologia per {kind.name}: {kind.topology}")
    run_output_dir = output_dir / kind.name / run_dir.name
    run_output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Carregant {kind.name}/{run_dir.name}: {trajectory_path}")
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
    write_rows(run_output_dir / "rmsd.csv", ["frame", "time_ns", "rmsd_nm"], [[idx, time_ns, value] for idx, (time_ns, value) in enumerate(zip(times_ns, rmsd, strict=True))])
    write_rows(run_output_dir / "radius_of_gyration.csv", ["frame", "time_ns", "rg_nm"], [[idx, time_ns, value] for idx, (time_ns, value) in enumerate(zip(times_ns, rg, strict=True))])
    write_rows(run_output_dir / "rmsf_ca.csv", ["residue", "rmsf_nm"], [[residue, value] for residue, value in zip(residues, rmsf, strict=True)])
    plot_series: dict[str, object] = {"rmsd": rmsd, "radius_of_gyration": rg, "rmsf": rmsf, "rmsf_residues": residues}
    hbond_summary: list[str] = ["Ponts d'hidrogen: omesos"]
    if not skip_hbonds:
        hbond_counts = compute_hbond_counts(md, traj)
        write_rows(run_output_dir / "hydrogen_bonds.csv", ["frame", "time_ns", "n_hydrogen_bonds"], [[idx, time_ns, int(value)] for idx, (time_ns, value) in enumerate(zip(times_ns, hbond_counts, strict=True))])
        plot_series["hbonds"] = hbond_counts
        hbond_summary = [f"Ponts d'hidrogen mitjans: {float(np.mean(hbond_counts)):.6f}"]
    save_plots(run_output_dir, times_ns, plot_series)
    log_summary = write_log_outputs(log_path, run_output_dir, report_interval, timestep_fs)
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
        *hbond_summary,
        *log_summary,
    ]
    (run_output_dir / "summary.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")
    print(f"Analisi desada a: {run_output_dir}")
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
        raise SystemExit("Per executar l'analisi cal instal-lar mdtraj. Exemple: pip install mdtraj") from exc
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
            analyzed_any = analyze_run(md, kind, run_dir, args.output_dir, args.report_interval, args.timestep_fs, args.skip_hbonds) or analyzed_any
    if not analyzed_any:
        raise SystemExit("No s'ha analitzat cap trajectoria.")


if __name__ == "__main__":
    main()
