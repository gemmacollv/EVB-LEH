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
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "calculs" / "analisi" / "resultats" / "joined"
DEFAULT_SETUP_DIR = PROJECT_ROOT / "calculs" / "prepared pdbs" / "02_openmm_md_setup"
DEFAULT_REPORT_INTERVAL = 10000
DEFAULT_TIMESTEP_FS = 4.0
ACTIVE_SITE_RESIDUES = {
    "LEU53", "MET73", "ARG94", "ASP96", "LEU98", "ASP127", "PHE129",
    "LEU199", "MET219", "ARG240", "ASP242", "LEU244", "ASP273", "PHE275",
}
ACTIVE_SITE_CONTACT_CUTOFF_NM = 0.45
DEFAULT_CATALYTIC_DISTANCE_SPECS = [
    (
        "HPN_C1_ASP_OD",
        "resname HPN and name C1",
        "protein and resname ASP and (name OD1 or name OD2)",
    ),
    (
        "HPN_O1_ARG_NH_NE",
        "resname HPN and name O1",
        "protein and resname ARG and (name NH1 or name NH2 or name NE)",
    ),
    (
        "HPN_O1_ASP_OD",
        "resname HPN and name O1",
        "protein and resname ASP and (name OD1 or name OD2)",
    ),
]
DEFAULT_ATTACK_ANGLE_SPECS = [
    ("HPN_N1_C1_O1", "resname HPN and name N1", "resname HPN and name C1", "resname HPN and name O1"),
]


@dataclass(frozen=True)
class SimulationKind:
    name: str
    topology: Path
    trajectory_name: str
    log_name: str


SIMULATION_KINDS = {
    "apo": SimulationKind(
        name="apo",
        topology=DEFAULT_SETUP_DIR / "protein_apo" / "protein_apo.prmtop",
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
            "Uneix tots els runs apo/holo en memoria i calcula RMSD, RMSF C-α, "
            "radi de gir i resum termodinamic per cada sistema."
        )
    )
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--kind", choices=["all", "apo", "holo"], default="all")
    parser.add_argument(
        "--run",
        action="append",
        default=None,
        help="Run concret a incloure, per exemple run_001. Es pot repetir.",
    )
    parser.add_argument("--report-interval", type=int, default=DEFAULT_REPORT_INTERVAL)
    parser.add_argument("--timestep-fs", type=float, default=DEFAULT_TIMESTEP_FS)
    parser.add_argument(
        "--skip-hbonds",
        action="store_true",
        help="Omet el calcul de ponts d'hidrogen, que pot ser mes lent.",
    )
    parser.add_argument(
        "--ligand-resname",
        default="HPN",
        help="Nom del residu del lligand per filtrar ponts d'hidrogen lligand-entorn.",
    )
    parser.add_argument(
        "--only-ligand-hbonds",
        action="store_true",
        help="Calcula nomes els ponts d'hidrogen proteina-lligand per accelerar aquesta analisi.",
    )
    parser.add_argument(
        "--save-joined-dcd",
        action="store_true",
        help="Desa la trajectoria concatenada com joined-apo.dcd/joined-holo.dcd.",
    )
    parser.add_argument(
        "--skip-catalytic-figure",
        action="store_true",
        help="Omet la figura de distancies catalitiques i contactes lligand-centre actiu.",
    )
    parser.add_argument(
        "--catalytic-distance",
        action="append",
        default=None,
        metavar="NOM::SELECTOR1::SELECTOR2",
        help="Distancia addicional per a la Figura 8, com distancia minima entre dos selectors MDTraj.",
    )
    parser.add_argument(
        "--attack-angle",
        action="append",
        default=None,
        metavar="NOM::SELECTOR1::SELECTOR2::SELECTOR3",
        help="Angle addicional per a la Figura 8. El segon selector es pren com a vertex de l angle.",
    )
    parser.add_argument(
        "--no-default-catalytic-metrics",
        action="store_true",
        help="No calcula les distancies/angles catalitics per defecte basats en HPN.",
    )
    parser.add_argument(
        "--active-site-contact-cutoff-nm",
        type=float,
        default=ACTIVE_SITE_CONTACT_CUTOFF_NM,
        help="Tall en nm per considerar contacte entre lligand i residus del centre actiu.",
    )
    return parser.parse_args()


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


def save_rmsf_plot(output_path: Path, residues: list[str], rmsf_nm: np.ndarray, title: str) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        print("Matplotlib no esta instal-lat; s'escriuen nomes els CSV.")
        return

    x_values = np.arange(1, len(rmsf_nm) + 1)
    active_positions = [idx + 1 for idx, residue in enumerate(residues) if residue in ACTIVE_SITE_RESIDUES]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(9, 5), dpi=300)
    plt.plot(x_values, rmsf_nm, linewidth=1.4, label="C-α")
    for position in active_positions:
        plt.axvline(position, color="tab:red", linewidth=0.7, alpha=0.35)
    if active_positions:
        plt.scatter(active_positions, rmsf_nm[np.array(active_positions) - 1], color="tab:red", s=14, label="Centre actiu")
    plt.title(title)
    plt.xlabel("Residus C-α")
    plt.ylabel("Distancia (nm)")
    plt.grid(True, linestyle="--", alpha=0.4)
    if active_positions:
        plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def save_catalytic_preorganization_plot(
    output_path: Path,
    times_ns: np.ndarray,
    distance_series_nm: dict[str, np.ndarray],
    angle_series_deg: dict[str, np.ndarray],
) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        print("Matplotlib no esta instal-lat; s\x27escriuen nomes els CSV.")
        return

    if not distance_series_nm and not angle_series_deg:
        return

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax_distance = plt.subplots(figsize=(9, 5), dpi=300)
    for label, values in distance_series_nm.items():
        ax_distance.plot(times_ns, values, linewidth=1.2, label=label)
    ax_distance.set_title("Figura 8. Distancies catalitiques")
    ax_distance.set_xlabel("Temps (ns)")
    ax_distance.set_ylabel("Distancia (nm)")
    ax_distance.grid(True, linestyle="--", alpha=0.4)

    handles, labels = ax_distance.get_legend_handles_labels()
    if angle_series_deg:
        ax_angle = ax_distance.twinx()
        for label, values in angle_series_deg.items():
            line = ax_angle.plot(times_ns, values, linewidth=1.0, linestyle=":", label=f"{label} angle")[0]
            handles.append(line)
            labels.append(line.get_label())
        ax_angle.set_ylabel("Angle (graus)")

    if handles:
        ax_distance.legend(handles, labels, fontsize=7, loc="best")
    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)


def save_contact_bar_plot(output_path: Path, contact_rows: list[list[object]], title: str) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        print("Matplotlib no esta instal-lat; s\x27escriuen nomes els CSV.")
        return

    if not contact_rows:
        return

    labels = [str(row[0]) for row in contact_rows]
    occupancies = [float(row[3]) for row in contact_rows]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(9, 5), dpi=300)
    plt.bar(labels, occupancies, color="tab:green", alpha=0.8)
    plt.title(title)
    plt.xlabel("Residus del centre actiu")
    plt.ylabel("Ocupacio de contacte (%)")
    plt.xticks(rotation=45, ha="right")
    plt.grid(True, axis="y", linestyle="--", alpha=0.4)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

def iter_run_dirs(data_dir: Path, kind_name: str, requested_runs: list[str] | None) -> list[Path]:
    kind_dir = data_dir / kind_name
    if requested_runs:
        return [kind_dir / run_name for run_name in requested_runs]
    return sorted(path for path in kind_dir.glob("run_*") if path.is_dir())


def frame_times_ns(n_frames: int, report_interval: int, timestep_fs: float) -> np.ndarray:
    ns_per_frame = report_interval * timestep_fs / 1_000_000.0
    return np.arange(n_frames, dtype=float) * ns_per_frame


def select_atoms(topology, selector: str, label: str) -> np.ndarray:
    indices = topology.select(selector)
    if len(indices) == 0:
        raise ValueError(f"No s'han trobat atoms per a {label}: {selector}")
    return indices


def load_and_join_trajectories(md, kind: SimulationKind, run_dirs: list[Path]):
    if not kind.topology.exists():
        raise SystemExit(f"No existeix la topologia per {kind.name}: {kind.topology}")

    trajectories = []
    frame_map = []
    frame_offset = 0
    for run_dir in run_dirs:
        trajectory_path = run_dir / kind.trajectory_name
        if not trajectory_path.exists():
            print(f"S'omet {kind.name}/{run_dir.name}: no existeix {trajectory_path.name}")
            continue
        print(f"Carregant {kind.name}/{run_dir.name}: {trajectory_path}")
        traj = md.load_dcd(str(trajectory_path), top=str(kind.topology))
        trajectories.append(traj)
        for frame in range(traj.n_frames):
            frame_map.append([frame_offset + frame, run_dir.name, frame])
        frame_offset += traj.n_frames

    if not trajectories:
        return None, []

    joined = trajectories[0]
    for traj in trajectories[1:]:
        joined = joined.join(traj, check_topology=True, discard_overlapping_frames=False)
    return joined, frame_map


def read_log_rows(log_path: Path, run_name: str, timestep_fs: float) -> list[list[object]]:
    if not log_path.exists():
        return []
    rows = []
    with log_path.open("r", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if not header:
            return []
        columns = [column.strip().lstrip("#").strip('"') for column in header]
        for raw_row in reader:
            if not raw_row:
                continue
            values = dict(zip(columns, raw_row, strict=True))
            step = int(float(values.get("Step", 0)))
            rows.append(
                [
                    run_name,
                    step,
                    step * timestep_fs / 1_000_000.0,
                    values.get("Potential Energy (kJ/mole)", ""),
                    values.get("Temperature (K)", ""),
                ]
            )
    return rows


def write_thermo(kind: SimulationKind, run_dirs: list[Path], output_dir: Path, timestep_fs: float) -> list[str]:
    rows = []
    for run_dir in run_dirs:
        rows.extend(read_log_rows(run_dir / kind.log_name, run_dir.name, timestep_fs))
    if not rows:
        return ["Thermo: no hi ha logs disponibles"]

    write_rows(
        output_dir / "thermo_by_run.csv",
        ["run", "step", "time_ns_within_run", "potential_energy_kj_mol", "temperature_k"],
        rows,
    )
    temperatures = [float(row[4]) for row in rows if row[4] != ""]
    potentials = [float(row[3]) for row in rows if row[3] != ""]
    return [
        f"Thermo frames: {len(rows)}",
        f"Temperatura mitjana (K): {float(np.mean(temperatures)):.6f}" if temperatures else "Temperatura mitjana (K): n/a",
        f"Energia potencial mitjana (kJ/mol): {float(np.mean(potentials)):.6f}" if potentials else "Energia potencial mitjana (kJ/mol): n/a",
    ]


def compute_hbond_counts(md, traj) -> np.ndarray:
    try:
        hbonds_by_frame = md.wernet_nilsson(traj, periodic=True)
    except Exception:
        hbonds_by_frame = md.wernet_nilsson(traj, periodic=False)
    return np.array([len(frame_hbonds) for frame_hbonds in hbonds_by_frame], dtype=int)


def compute_ligand_hbond_counts(md, traj, ligand_resname: str) -> np.ndarray:
    ligand_atoms_original = set(traj.topology.select(f"resname {ligand_resname}"))
    if not ligand_atoms_original:
        raise ValueError(f"No s'han trobat atoms del lligand amb resname {ligand_resname}.")

    protein_atoms = traj.topology.select("protein")
    if len(protein_atoms) == 0:
        raise ValueError("No s'han trobat atoms de proteina per calcular ponts proteina-lligand.")

    try:
        nearby_by_frame = md.compute_neighbors(
            traj,
            0.45,
            query_indices=np.array(sorted(ligand_atoms_original), dtype=int),
            haystack_indices=protein_atoms,
            periodic=True,
        )
    except Exception:
        nearby_by_frame = md.compute_neighbors(
            traj,
            0.45,
            query_indices=np.array(sorted(ligand_atoms_original), dtype=int),
            haystack_indices=protein_atoms,
            periodic=False,
        )

    selected_atoms = set(ligand_atoms_original)
    for nearby_atoms in nearby_by_frame:
        for atom_index in nearby_atoms:
            residue = traj.topology.atom(int(atom_index)).residue
            selected_atoms.update(atom.index for atom in residue.atoms)

    if selected_atoms == ligand_atoms_original:
        return np.zeros(traj.n_frames, dtype=int)

    local_traj = traj.atom_slice(np.array(sorted(selected_atoms), dtype=int))
    ligand_atoms = set(local_traj.topology.select(f"resname {ligand_resname}"))

    try:
        hbonds = md.baker_hubbard(local_traj, freq=0.0, periodic=True)
    except Exception:
        hbonds = md.baker_hubbard(local_traj, freq=0.0, periodic=False)

    ligand_hbonds = []
    for donor, hydrogen, acceptor in hbonds:
        hbond_atoms = {int(donor), int(hydrogen), int(acceptor)}
        touches_ligand = bool(hbond_atoms & ligand_atoms)
        touches_environment = bool(hbond_atoms - ligand_atoms)
        if touches_ligand and touches_environment:
            ligand_hbonds.append([int(donor), int(hydrogen), int(acceptor)])

    if not ligand_hbonds:
        return np.zeros(local_traj.n_frames, dtype=int)

    ligand_hbonds = np.array(ligand_hbonds, dtype=int)
    h_acceptor_pairs = ligand_hbonds[:, [1, 2]]
    try:
        distances = md.compute_distances(local_traj, h_acceptor_pairs, periodic=True)
        angles = md.compute_angles(local_traj, ligand_hbonds, periodic=True)
    except Exception:
        distances = md.compute_distances(local_traj, h_acceptor_pairs, periodic=False)
        angles = md.compute_angles(local_traj, ligand_hbonds, periodic=False)
    present = (distances < 0.25) & (np.degrees(angles) > 120.0)
    return np.sum(present, axis=1).astype(int)


def analyze_joined(md, kind: SimulationKind, run_dirs: list[Path], output_dir: Path, args: argparse.Namespace) -> bool:
    joined, frame_map = load_and_join_trajectories(md, kind, run_dirs)
    if joined is None:
        print(f"No s'ha pogut carregar cap trajectoria per {kind.name}.")
        return False

    kind_output_dir = output_dir / kind.name
    kind_output_dir.mkdir(parents=True, exist_ok=True)
    if not args.only_ligand_hbonds:
        try:
            joined = joined.image_molecules(inplace=False)
        except Exception as exc:
            print(f"Avís: no s ha pogut recentrar/aplicar PBC a {kind.name}: {exc}")
    times_ns = frame_times_ns(joined.n_frames, args.report_interval, args.timestep_fs)
    system_label = kind.name.upper()

    write_rows(kind_output_dir / "frame_map.csv", ["joined_frame", "run", "frame_in_run"], frame_map)

    if args.only_ligand_hbonds:
        if kind.name != "holo":
            print("--only-ligand-hbonds nomes aplica al sistema holo; s'omet aquest sistema.")
            return False
        ligand_hbond_counts = compute_ligand_hbond_counts(md, joined, args.ligand_resname)
        write_rows(
            kind_output_dir / "ligand_hydrogen_bonds.csv",
            ["frame", "time_ns", "n_ligand_hydrogen_bonds"],
            [
                [idx, time_ns, int(value)]
                for idx, (time_ns, value) in enumerate(zip(times_ns, ligand_hbond_counts, strict=True))
            ],
        )
        save_line_plot(
            kind_output_dir / "ligand_hydrogen_bonds.png",
            times_ns,
            ligand_hbond_counts,
            f"{system_label}: ponts lligand-proteina",
            "Temps (ns)",
            "Nombre de ponts",
        )
        summary = [
            f"System: {kind.name}",
            f"Runs units: {', '.join(run_dir.name for run_dir in run_dirs)}",
            f"Frames totals: {joined.n_frames}",
            f"Temps final concatenat (ns): {times_ns[-1]:.6f}" if len(times_ns) else "Temps final concatenat (ns): 0.000000",
            f"Ponts d'hidrogen proteina-lligand mitjans ({args.ligand_resname}): {float(np.mean(ligand_hbond_counts)):.6f}",
        ]
        (kind_output_dir / "ligand_hydrogen_bonds_summary.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")
        print(f"Analisi de ponts proteina-lligand desada a: {kind_output_dir}")
        return True

    if args.save_joined_dcd:
        joined.save_dcd(str(kind_output_dir / f"joined-{kind.name}.dcd"))

    protein_atoms = select_atoms(joined.topology, "protein", "proteina")
    ca_atoms = select_atoms(joined.topology, "protein and name CA", "C-α")

    protein_traj = joined.atom_slice(protein_atoms)
    protein_traj.superpose(protein_traj, frame=0)
    rmsd = md.rmsd(protein_traj, protein_traj, frame=0)
    rg = md.compute_rg(protein_traj)

    ca_traj = joined.atom_slice(ca_atoms)
    ca_traj.superpose(ca_traj, frame=0)
    rmsf = md.rmsf(ca_traj, ca_traj, frame=0)
    residues = [str(atom.residue) for atom in ca_traj.topology.atoms]

    write_rows(
        kind_output_dir / "rmsd.csv",
        ["frame", "time_ns", "rmsd_nm"],
        [[idx, time_ns, value] for idx, (time_ns, value) in enumerate(zip(times_ns, rmsd, strict=True))],
    )
    write_rows(
        kind_output_dir / "radius_of_gyration.csv",
        ["frame", "time_ns", "rg_nm"],
        [[idx, time_ns, value] for idx, (time_ns, value) in enumerate(zip(times_ns, rg, strict=True))],
    )
    write_rows(
        kind_output_dir / "rmsf_ca.csv",
        ["residue", "rmsf_nm"],
        [[residue, value] for residue, value in zip(residues, rmsf, strict=True)],
    )

    save_line_plot(
        kind_output_dir / "rmsd.png",
        times_ns,
        rmsd,
        f"{system_label}: RMSD de la proteina",
        "Temps (ns)",
        "Distancia (nm)",
    )
    save_line_plot(
        kind_output_dir / "radius_of_gyration.png",
        times_ns,
        rg,
        f"{system_label}: radi de gir",
        "Temps (ns)",
        "Distancia (nm)",
    )
    save_rmsf_plot(
        kind_output_dir / "rmsf_ca.png",
        residues,
        rmsf,
        f"{system_label}: flexibilitat C-α",
    )

    hbond_summary = ["Ponts d'hidrogen: omesos"]
    if not args.skip_hbonds:
        hbond_counts = compute_hbond_counts(md, joined)
        write_rows(
            kind_output_dir / "hydrogen_bonds.csv",
            ["frame", "time_ns", "n_hydrogen_bonds"],
            [[idx, time_ns, int(value)] for idx, (time_ns, value) in enumerate(zip(times_ns, hbond_counts, strict=True))],
        )
        save_line_plot(
            kind_output_dir / "hydrogen_bonds.png",
            times_ns,
            hbond_counts,
            f"{system_label}: ponts d'hidrogen",
            "Temps (ns)",
            "Nombre de ponts",
        )
        hbond_summary = [f"Ponts d'hidrogen mitjans: {float(np.mean(hbond_counts)):.6f}"]

    if kind.name == "holo" and (not args.skip_hbonds or args.only_ligand_hbonds):
        ligand_hbond_counts = compute_ligand_hbond_counts(md, joined, args.ligand_resname)
        write_rows(
            kind_output_dir / "ligand_hydrogen_bonds.csv",
            ["frame", "time_ns", "n_ligand_hydrogen_bonds"],
            [
                [idx, time_ns, int(value)]
                for idx, (time_ns, value) in enumerate(zip(times_ns, ligand_hbond_counts, strict=True))
            ],
        )
        save_line_plot(
            kind_output_dir / "ligand_hydrogen_bonds.png",
            times_ns,
            ligand_hbond_counts,
            f"{system_label}: ponts lligand-proteina",
            "Temps (ns)",
            "Nombre de ponts",
        )
        hbond_summary.append(
            f"Ponts d'hidrogen proteina-lligand mitjans ({args.ligand_resname}): "
            f"{float(np.mean(ligand_hbond_counts)):.6f}"
        )

    thermo_summary = write_thermo(kind, run_dirs, kind_output_dir, args.timestep_fs)
    summary = [
        f"System: {kind.name}",
        f"Runs: {', '.join(run_dir.name for run_dir in run_dirs)}",
        f"Topology: {kind.topology}",
        f"Frames totals: {joined.n_frames}",
        f"Temps final (ns): {times_ns[-1]:.6f}" if len(times_ns) else "Temps final (ns): 0.000000",
        f"RMSD mitja proteina (nm): {float(np.mean(rmsd)):.6f}" if len(rmsd) else "RMSD mitjà proteina (nm): 0.000000",
        f"Radi de gir mitja (nm): {float(np.mean(rg)):.6f}" if len(rg) else "Radi de gir mitjà (nm): 0.000000",
        f"RMSF C-α maxim (nm): {float(np.max(rmsf)):.6f}" if len(rmsf) else "RMSF C-α màxim (nm): 0.000000",
        *hbond_summary,
        *thermo_summary,
    ]
    (kind_output_dir / "summary.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")
    print(f"Analisi concatenada desada a: {kind_output_dir}")
    return True


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
        existing_run_dirs = [run_dir for run_dir in run_dirs if run_dir.exists()]
        if not existing_run_dirs:
            print(f"No s'han trobat runs per {kind_name} dins {args.data_dir / kind_name}")
            continue
        analyzed_any = analyze_joined(md, kind, existing_run_dirs, args.output_dir, args) or analyzed_any

    if not analyzed_any:
        raise SystemExit("No s'ha analitzat cap trajectoria.")


if __name__ == "__main__":
    main()
