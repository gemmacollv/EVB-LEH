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
    "TYR48", "ASN50", "ARG94", "ASP96", "ASP127",
    "TYR194", "ASN196", "ARG240", "ASP242", "ASP273",
}
ACTIVE_SITE_CONTACT_CUTOFF_NM = 0.45
ACTIVE_SITE_DISPLAY_LABELS = {
    "TYR48": "Tyr53 A",
    "ASN50": "Asn55 A",
    "ARG94": "Arg99 A",
    "ASP96": "Asp101 A",
    "ASP127": "Asp132 A",
    "TYR194": "Tyr53 B",
    "ASN196": "Asn55 B",
    "ARG240": "Arg99 B",
    "ASP242": "Asp101 B",
    "ASP273": "Asp132 B",
}
GENERAL_BASE_ASP_SELECTOR = "protein and resname ASP and (resSeq 127 or resSeq 273)"
CATALYTIC_TYR_SELECTOR = "protein and resname TYR and (resSeq 48 or resSeq 194)"
CATALYTIC_ASN_SELECTOR = "protein and resname ASN and (resSeq 50 or resSeq 196)"
DEFAULT_CATALYTIC_DISTANCE_SPECS = [
    (
        "WAT_O_HPN_C1",
        "water and name O",
        "resname HPN and (name C1 or name C1x)",
    ),
    (
        "HPN_O1_TYR53_OH",
        "resname HPN and (name O1 or name O1x)",
        f"{CATALYTIC_TYR_SELECTOR} and name OH",
    ),
    (
        "HPN_O1_ASN55_ND2",
        "resname HPN and (name O1 or name O1x)",
        f"{CATALYTIC_ASN_SELECTOR} and name ND2",
    ),
]
CATALYTIC_DISTANCE_DISPLAY_LABELS = {
    "WAT_O_HPN_C1": "Aigua nucleòfila - HPN C1",
    "HPN_O1_TYR53_OH": "HPN O1 - Tyr53",
    "HPN_O1_ASN55_ND2": "HPN O1 - Asn55",
}
NUCLEOPHILIC_ATTACK_ANGLE_LABEL = "ASP132_OD_WAT_O_HPN_C1"
NUCLEOPHILIC_WATER_C1_CUTOFF_NM = 0.65
NUCLEOPHILIC_WATER_ASP_CUTOFF_NM = 0.45
_NUCLEOPHILIC_WATER_GEOMETRY_CACHE = {}

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
            "Uneix tots els runs apo/holo en memòria i calcula RMSD, RMSF C-α, "
            "radi de gir i resum termodinàmic per cada sistema."
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
        help="Omet el càlcul de ponts d'hidrogen, que pot ser més lent.",
    )
    parser.add_argument(
        "--ligand-resname",
        default="HPN",
        help="Nom del residu del lligand per filtrar ponts d'hidrogen lligand-entorn.",
    )
    parser.add_argument(
        "--only-ligand-hbonds",
        action="store_true",
        help="Calcula només els ponts d'hidrogen proteïna-lligand per accelerar aquesta anàlisi.",
    )
    parser.add_argument(
        "--save-joined-dcd",
        action="store_true",
        help="Desa la trajectòria concatenada com joined-apo.dcd/joined-holo.dcd.",
    )
    parser.add_argument(
        "--skip-catalytic-figure",
        action="store_true",
        help="Omet la figura de distàncies catalítiques i contactes lligand-centre actiu.",
    )
    parser.add_argument(
        "--only-catalytic-metrics",
        action="store_true",
        help="Calcula només les distàncies i l'angle catalítics del sistema holo.",
    )
    parser.add_argument(
        "--catalytic-distance",
        action="append",
        default=None,
        metavar="NOM::SELECTOR1::SELECTOR2",
        help="Distància addicional per a la Figura 8, com distància mínima entre dos selectors MDTraj.",
    )
    parser.add_argument(
        "--no-default-catalytic-metrics",
        action="store_true",
        help="No calcula les distàncies catalítiques per defecte basades en HPN.",
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
        print("Matplotlib no està instal·lat; s'escriuen només els CSV.")
        return

    output_path.parent.mkdir(parents=True, exist_ok=True)
    x_values = np.asarray(x_values, dtype=float)
    y_values = np.asarray(y_values, dtype=float)
    finite = y_values[np.isfinite(y_values)]
    y_min = float(np.min(finite)) if len(finite) else 0.0
    y_max = float(np.max(finite)) if len(finite) else 1.0
    y_span = max(y_max - y_min, 0.01)
    y_floor = max(0.0, y_min - 0.08 * y_span)

    fig, ax = plt.subplots(figsize=(9, 5.2), dpi=300)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.plot(x_values, y_values, linewidth=2.0, color="#1f77b4")
    ax.set_ylim(y_floor, y_max + 0.12 * y_span)
    ax.set_title(title, fontsize=14, weight="bold", pad=10)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, axis="y", linestyle="-", alpha=0.18)
    ax.grid(True, axis="x", linestyle=":", alpha=0.08)
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)
    ax.spines["left"].set_color("#b0b0b0")
    ax.spines["bottom"].set_color("#b0b0b0")
    ax.spines["top"].set_color("#b0b0b0")
    ax.spines["right"].set_color("#b0b0b0")
    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)


def save_rmsf_plot(output_path: Path, residues: list[str], rmsf_nm: np.ndarray, title: str) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        print("Matplotlib no està instal·lat; s'escriuen només els CSV.")
        return

    x_values = np.arange(1, len(rmsf_nm) + 1)

    active_points = [
        (idx + 1, ACTIVE_SITE_DISPLAY_LABELS[residue])
        for idx, residue in enumerate(residues)
        if residue in ACTIVE_SITE_DISPLAY_LABELS
    ]
    site_colors = ["#1f77b4", "#ff7f0e"]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(12, 6.2), dpi=300)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    y_min = float(np.min(rmsf_nm)) if len(rmsf_nm) else 0.0
    y_max = float(np.max(rmsf_nm)) if len(rmsf_nm) else 1.0
    y_span = max(y_max - y_min, 0.01)
    y_floor = max(0.0, y_min - 0.08 * y_span)
    y_ceiling = y_max + 0.18 * y_span
    ax.set_ylim(y_floor, y_ceiling)
    ax.plot(x_values, rmsf_nm, linewidth=2.0, color="#1f77b4")
    ax.plot(x_values, rmsf_nm, linewidth=0.8, color="#8ec7f0", alpha=0.9)
    label_lanes = [0.74, 0.56, 0.38]
    for color_idx, (position, label) in enumerate(active_points):
        color = site_colors[color_idx % len(site_colors)]
        point_y = float(rmsf_nm[position - 1])
        label_y = y_min + y_span * label_lanes[color_idx % len(label_lanes)]
        label_x = min(position + 1.2, len(rmsf_nm))
        ax.axvspan(position - 0.55, position + 0.55, color=color, alpha=0.18, linewidth=0)
        ax.axvline(position, color=color, linewidth=0.9, alpha=0.5)
        ax.scatter([position], [point_y], color=color, edgecolor="white", linewidth=0.6, s=32, zorder=4)
        ax.annotate(
            label,
            xy=(position, point_y),
            xytext=(label_x, label_y),
            textcoords="data",
            ha="left",
            va="center",
            fontsize=6.7,
            color="#263238",
            arrowprops={"arrowstyle": "-", "color": color, "alpha": 0.7, "linewidth": 0.8},
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": color, "alpha": 0.92, "linewidth": 0.75},
        )
    ax.set_title(title, fontsize=14, weight="bold", pad=10)
    ax.set_xlabel("Residus C-α")
    ax.set_ylabel("Distància (nm)")
    ax.grid(True, axis="y", linestyle="-", alpha=0.18)
    ax.grid(True, axis="x", linestyle=":", alpha=0.08)
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)
    ax.spines["left"].set_color("#b0b0b0")
    ax.spines["bottom"].set_color("#b0b0b0")
    ax.spines["top"].set_color("#b0b0b0")
    ax.spines["right"].set_color("#b0b0b0")
    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)


def save_catalytic_preorganization_plot(
    output_path: Path,
    times_ns: np.ndarray,
    distance_series_nm: dict[str, np.ndarray],
) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        print("Matplotlib no està instal·lat; s\x27escriuen només els CSV.")
        return

    if not distance_series_nm:
        return

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax_distance = plt.subplots(figsize=(9.5, 5.4), dpi=300)
    fig.patch.set_facecolor("white")
    ax_distance.set_facecolor("white")
    palette = ["#1f77b4", "#ff7f0e", "#2ca02c"]
    linestyles = ["-", "-", "-", "--"]
    for color_idx, (label, values) in enumerate(distance_series_nm.items()):
        ax_distance.plot(
            times_ns,
            values,
            linewidth=2.0,
            color=palette[color_idx % len(palette)],
            linestyle=linestyles[color_idx % len(linestyles)],
            label=CATALYTIC_DISTANCE_DISPLAY_LABELS.get(label, label),
        )
    ax_distance.set_title("Distàncies catalítiques", fontsize=14, weight="bold", pad=10)
    ax_distance.set_xlabel("Temps (ns)")
    ax_distance.set_ylabel("Distància (nm)")
    ax_distance.grid(True, axis="y", linestyle="-", alpha=0.18)
    ax_distance.grid(True, axis="x", linestyle=":", alpha=0.08)
    ax_distance.spines["top"].set_visible(True)
    ax_distance.spines["right"].set_visible(True)
    ax_distance.spines["left"].set_color("#b0b0b0")
    ax_distance.spines["bottom"].set_color("#b0b0b0")
    ax_distance.spines["top"].set_color("#b0b0b0")
    ax_distance.spines["right"].set_color("#b0b0b0")

    if distance_series_nm:
        ax_distance.legend(fontsize=6.7, loc="best", frameon=True, facecolor="white", edgecolor="#d6d6d6")
    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)


def save_contact_bar_plot(output_path: Path, contact_rows: list[list[object]], title: str) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        print("Matplotlib no està instal·lat; s\x27escriuen només els CSV.")
        return

    if not contact_rows:
        return

    labels = [str(row[0]) for row in contact_rows]
    occupancies = [float(row[3]) for row in contact_rows]
    bar_colors = ["#1f77b4", "#ff7f0e"]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(9.5, 5.4), dpi=300)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.bar(labels, occupancies, color=[bar_colors[idx % len(bar_colors)] for idx in range(len(labels))], alpha=0.84, edgecolor="white", linewidth=0.8)
    ax.set_title(title, fontsize=14, weight="bold", pad=10)
    ax.set_xlabel("Residus del centre actiu")
    ax.set_ylabel("Ocupació de contacte (%)")
    ax.tick_params(axis="x", rotation=45)
    for tick in ax.get_xticklabels():
        tick.set_ha("right")
    ax.grid(True, axis="y", linestyle="-", alpha=0.18)
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)
    ax.spines["left"].set_color("#b0b0b0")
    ax.spines["bottom"].set_color("#b0b0b0")
    ax.spines["top"].set_color("#b0b0b0")
    ax.spines["right"].set_color("#b0b0b0")
    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)

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
        raise ValueError(f"No s'han trobat àtoms per a {label}: {selector}")
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
        raise ValueError(f"No s'han trobat àtoms del lligand amb resname {ligand_resname}.")

    protein_atoms = traj.topology.select("protein")
    if len(protein_atoms) == 0:
        raise ValueError("No s'han trobat àtoms de proteïna per calcular ponts proteïna-lligand.")

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



def parse_metric_spec(raw_spec: str, expected_parts: int, option_name: str) -> tuple[str, ...]:
    parts = tuple(part.strip() for part in raw_spec.split("::"))
    if len(parts) != expected_parts or any(not part for part in parts):
        raise ValueError(f"Format invalid per {option_name}: {raw_spec}")
    return parts


def catalytic_distance_specs(args: argparse.Namespace) -> list[tuple[str, str, str]]:
    specs: list[tuple[str, str, str]] = []
    if not args.no_default_catalytic_metrics:
        specs.extend(DEFAULT_CATALYTIC_DISTANCE_SPECS)
    for raw_spec in args.catalytic_distance or []:
        specs.append(parse_metric_spec(raw_spec, 3, "--catalytic-distance"))
    return specs



def select_nucleophilic_water_geometry(md, traj) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, list[object]]:
    cache_key = id(traj)
    if cache_key in _NUCLEOPHILIC_WATER_GEOMETRY_CACHE:
        return _NUCLEOPHILIC_WATER_GEOMETRY_CACHE[cache_key]

    water_selector = "water and name O"
    c1_selector = "resname HPN and (name C1 or name C1x)"
    asp_selector = f"{GENERAL_BASE_ASP_SELECTOR} and (name OD1 or name OD2)"
    c1_cutoff_nm = NUCLEOPHILIC_WATER_C1_CUTOFF_NM
    asp_cutoff_nm = NUCLEOPHILIC_WATER_ASP_CUTOFF_NM
    water_atoms = traj.topology.select(water_selector)
    c1_atoms = traj.topology.select(c1_selector)
    asp_od_atoms = traj.topology.select(asp_selector)
    if len(water_atoms) == 0 or len(c1_atoms) == 0 or len(asp_od_atoms) == 0:
        raise ValueError("Selecció buida per identificar l'aigua nucleòfila.")

    selected_distances = np.full(traj.n_frames, np.nan, dtype=float)
    selected_waters = np.full(traj.n_frames, -1, dtype=int)
    selected_c1_atoms = np.full(traj.n_frames, -1, dtype=int)
    selected_asp_atoms = np.full(traj.n_frames, -1, dtype=int)
    chunk_size = 100

    for start in range(0, traj.n_frames, chunk_size):
        stop = min(start + chunk_size, traj.n_frames)
        chunk = traj[start:stop]
        try:
            nearby_waters_by_frame = md.compute_neighbors(
                chunk,
                c1_cutoff_nm,
                query_indices=c1_atoms,
                haystack_indices=water_atoms,
                periodic=True,
            )
        except Exception:
            nearby_waters_by_frame = md.compute_neighbors(
                chunk,
                c1_cutoff_nm,
                query_indices=c1_atoms,
                haystack_indices=water_atoms,
                periodic=False,
            )

        for local_frame, candidate_waters in enumerate(nearby_waters_by_frame):
            frame_index = start + local_frame
            candidates = np.asarray(candidate_waters, dtype=int)
            if len(candidates) == 0:
                continue

            frame_xyz = chunk.xyz[local_frame]
            candidate_xyz = frame_xyz[candidates]
            c1_xyz = frame_xyz[c1_atoms]
            asp_xyz = frame_xyz[asp_od_atoms]

            water_c1 = np.linalg.norm(candidate_xyz[:, np.newaxis, :] - c1_xyz[np.newaxis, :, :], axis=2)
            water_asp = np.linalg.norm(candidate_xyz[:, np.newaxis, :] - asp_xyz[np.newaxis, :, :], axis=2)
            nearest_c1_positions = np.argmin(water_c1, axis=1)
            nearest_asp_positions = np.argmin(water_asp, axis=1)
            nearest_c1_distances = np.min(water_c1, axis=1)
            nearest_asp_distances = np.min(water_asp, axis=1)
            nucleophilic_mask = (nearest_c1_distances <= c1_cutoff_nm) & (nearest_asp_distances <= asp_cutoff_nm)
            if not np.any(nucleophilic_mask):
                continue

            candidate_scores = nearest_c1_distances + nearest_asp_distances
            candidate_scores = np.where(nucleophilic_mask, candidate_scores, np.inf)
            best_candidate_position = int(np.argmin(candidate_scores))

            selected_distances[frame_index] = nearest_c1_distances[best_candidate_position]
            selected_waters[frame_index] = candidates[best_candidate_position]
            selected_c1_atoms[frame_index] = c1_atoms[nearest_c1_positions[best_candidate_position]]
            selected_asp_atoms[frame_index] = asp_od_atoms[nearest_asp_positions[best_candidate_position]]

    valid_frames = int(np.sum(selected_waters >= 0))
    metadata_row = [
        "WAT_O_HPN_C1",
        (
            f"{water_selector}; només aigües a <= {c1_cutoff_nm:.2f} nm de C1 "
            f"i <= {asp_cutoff_nm:.2f} nm d'Asp132 Oδ"
        ),
        c1_selector,
        valid_frames,
        len(c1_atoms),
        f"frames_valids={valid_frames}/{traj.n_frames}",
    ]
    result = (selected_distances, selected_waters, selected_c1_atoms, selected_asp_atoms, metadata_row)
    _NUCLEOPHILIC_WATER_GEOMETRY_CACHE[cache_key] = result
    return result


def finite_summary(values: np.ndarray) -> tuple[float, float, float] | None:
    finite_values = values[np.isfinite(values)]
    if len(finite_values) == 0:
        return None
    return float(np.mean(finite_values)), float(np.min(finite_values)), float(np.max(finite_values))


def compute_min_distances(md, traj, specs: list[tuple[str, str, str]]) -> tuple[dict[str, np.ndarray], list[list[object]]]:
    distance_series: dict[str, np.ndarray] = {}
    metadata_rows: list[list[object]] = []
    for label, selector_a, selector_b in specs:
        if label == "WAT_O_HPN_C1":
            try:
                distances, _waters, _c1_atoms, _asp_atoms, metadata_row = select_nucleophilic_water_geometry(md, traj)
            except ValueError as exc:
                print(f"S'omet distància catalítica {label}: {exc}")
                continue
            distance_series[label] = distances
            metadata_rows.append(metadata_row)
            continue

        atoms_a = traj.topology.select(selector_a)
        atoms_b = traj.topology.select(selector_b)
        if len(atoms_a) == 0 or len(atoms_b) == 0:
            print(f"S'omet distància catalítica {label}: selecció buida.")
            continue
        pairs = np.array([(int(atom_a), int(atom_b)) for atom_a in atoms_a for atom_b in atoms_b if atom_a != atom_b], dtype=int)
        if len(pairs) == 0:
            print(f"S'omet distància catalítica {label}: no hi ha parelles d'àtoms vàlides.")
            continue
        try:
            distances = md.compute_distances(traj, pairs, periodic=True)
        except Exception:
            distances = md.compute_distances(traj, pairs, periodic=False)
        min_distances = np.min(distances, axis=1)
        distance_series[label] = min_distances
        metadata_rows.append([label, selector_a, selector_b, len(atoms_a), len(atoms_b), len(pairs)])
    return distance_series, metadata_rows



def compute_nucleophilic_attack_angles(md, traj) -> tuple[dict[str, np.ndarray], list[list[object]]]:
    try:
        _distances, selected_waters, selected_c1_atoms, selected_asp_atoms, _metadata = select_nucleophilic_water_geometry(md, traj)
    except ValueError as exc:
        print(f"S'omet angle catalític {NUCLEOPHILIC_ATTACK_ANGLE_LABEL}: {exc}")
        return {}, []

    angles_deg = np.full(traj.n_frames, np.nan, dtype=float)
    valid_frames = (selected_waters >= 0) & (selected_c1_atoms >= 0) & (selected_asp_atoms >= 0)
    frame_indices = np.flatnonzero(valid_frames)
    if len(frame_indices) > 0:
        water_xyz = traj.xyz[frame_indices, selected_waters[valid_frames], :]
        c1_xyz = traj.xyz[frame_indices, selected_c1_atoms[valid_frames], :]
        asp_xyz = traj.xyz[frame_indices, selected_asp_atoms[valid_frames], :]

        asp_to_water = asp_xyz - water_xyz
        c1_to_water = c1_xyz - water_xyz
        dot_products = np.sum(asp_to_water * c1_to_water, axis=1)
        norms = np.linalg.norm(asp_to_water, axis=1) * np.linalg.norm(c1_to_water, axis=1)
        cosines = np.divide(dot_products, norms, out=np.full_like(dot_products, np.nan), where=norms > 0.0)
        angles_deg[frame_indices] = np.degrees(np.arccos(np.clip(cosines, -1.0, 1.0)))

    metadata_rows = [[
        NUCLEOPHILIC_ATTACK_ANGLE_LABEL,
        f"{GENERAL_BASE_ASP_SELECTOR} and (name OD1 or name OD2)",
        "water and name O; triada per proximitat conjunta a C1 i Asp132 Oδ",
        "resname HPN and (name C1 or name C1x)",
        len(np.unique(selected_asp_atoms[selected_asp_atoms >= 0])),
        len(np.unique(selected_waters[selected_waters >= 0])),
        len(np.unique(selected_c1_atoms[selected_c1_atoms >= 0])),
    ]]
    return {NUCLEOPHILIC_ATTACK_ANGLE_LABEL: angles_deg}, metadata_rows


def residue_label(residue) -> str:
    chain = getattr(residue.chain, "chain_id", None) or getattr(residue.chain, "id", "") or "?"
    return f"{residue.name}{residue.resSeq}:{chain}"


def active_site_residue_atoms(topology) -> dict[str, list[int]]:
    residues: dict[str, list[int]] = {}
    active_labels = {label.upper() for label in ACTIVE_SITE_RESIDUES}
    for residue in topology.residues:
        residue_key = f"{residue.name}{residue.resSeq}".upper()
        if residue_key in active_labels:
            residues[residue_label(residue)] = [atom.index for atom in residue.atoms]
    return residues


def compute_active_site_contacts(md, traj, ligand_resname: str, cutoff_nm: float) -> list[list[object]]:
    ligand_atoms = traj.topology.select(f"resname {ligand_resname}")
    if len(ligand_atoms) == 0:
        print(f"S'ometen contactes del centre actiu: no s'han trobat àtoms {ligand_resname}.")
        return []

    rows = []
    for label, residue_atoms in active_site_residue_atoms(traj.topology).items():
        pairs = np.array([(int(ligand_atom), int(residue_atom)) for ligand_atom in ligand_atoms for residue_atom in residue_atoms], dtype=int)
        if len(pairs) == 0:
            continue
        try:
            distances = md.compute_distances(traj, pairs, periodic=True)
        except Exception:
            distances = md.compute_distances(traj, pairs, periodic=False)
        min_distances = np.min(distances, axis=1)
        rows.append([
            label,
            float(np.mean(min_distances)),
            float(np.min(min_distances)),
            float(np.mean(min_distances < cutoff_nm) * 100.0),
        ])
    return sorted(rows, key=lambda row: float(row[3]), reverse=True)


def write_catalytic_metrics(md, joined, kind: SimulationKind, kind_output_dir: Path, times_ns: np.ndarray, args: argparse.Namespace) -> list[str]:
    if kind.name != "holo" or args.skip_catalytic_figure:
        return []

    distance_series, metadata_rows = compute_min_distances(md, joined, catalytic_distance_specs(args))
    summary: list[str] = []

    if distance_series:
        labels = list(distance_series)
        write_rows(
            kind_output_dir / "catalytic_atom_distances.csv",
            ["frame", "time_ns", *[f"{label}_nm" for label in labels]],
            [[idx, time_ns, *[distance_series[label][idx] for label in labels]] for idx, time_ns in enumerate(times_ns)],
        )
        distance_summary_rows = []
        for label, values in distance_series.items():
            stats = finite_summary(values)
            if stats is None:
                distance_summary_rows.append([label, np.nan, np.nan, np.nan])
                summary.append(f"Distància catalítica mitjana {label} (nm): n/a")
            else:
                mean_value, min_value, max_value = stats
                distance_summary_rows.append([label, mean_value, min_value, max_value])
                summary.append(f"Distància catalítica mitjana {label} (nm): {mean_value:.6f}")
        write_rows(
            kind_output_dir / "catalytic_atom_distance_summary.csv",
            ["metric", "mean_nm", "min_nm", "max_nm"],
            distance_summary_rows,
        )
    if metadata_rows:
        write_rows(
            kind_output_dir / "catalytic_atom_distance_selections.csv",
            ["metric", "selector_1", "selector_2", "n_atoms_1", "n_atoms_2", "n_pairs"],
            metadata_rows,
        )

    angle_series, angle_metadata_rows = compute_nucleophilic_attack_angles(md, joined)
    if angle_series:
        labels = list(angle_series)
        write_rows(
            kind_output_dir / "catalytic_attack_angles.csv",
            ["frame", "time_ns", *[f"{label}_deg" for label in labels]],
            [[idx, time_ns, *[angle_series[label][idx] for label in labels]] for idx, time_ns in enumerate(times_ns)],
        )
        angle_summary_rows = []
        for label, values in angle_series.items():
            stats = finite_summary(values)
            if stats is None:
                angle_summary_rows.append([label, np.nan, np.nan, np.nan])
                summary.append(f"Angle d'atac nucleòfil mitjà {label} (graus): n/a")
            else:
                mean_value, min_value, max_value = stats
                angle_summary_rows.append([label, mean_value, min_value, max_value])
                summary.append(f"Angle d'atac nucleòfil mitjà {label} (graus): {mean_value:.6f}")
        write_rows(
            kind_output_dir / "catalytic_attack_angle_summary.csv",
            ["metric", "mean_deg", "min_deg", "max_deg"],
            angle_summary_rows,
        )
    if angle_metadata_rows:
        write_rows(
            kind_output_dir / "catalytic_attack_angle_selections.csv",
            ["metric", "selector_1", "selector_2", "selector_3", "n_atoms_1", "n_atoms_2", "n_atoms_3"],
            angle_metadata_rows,
        )

    save_catalytic_preorganization_plot(kind_output_dir / "catalytic_preorganization.png", times_ns, distance_series)

    contact_rows = compute_active_site_contacts(md, joined, args.ligand_resname, args.active_site_contact_cutoff_nm)
    if contact_rows:
        write_rows(
            kind_output_dir / "active_site_ligand_contacts.csv",
            ["residue", "mean_min_distance_nm", "min_distance_nm", "contact_occupancy_percent"],
            contact_rows,
        )
        save_contact_bar_plot(kind_output_dir / "active_site_ligand_contacts.png", contact_rows, "Contactes lligand-centre actiu")
        summary.append(f"Residus del centre actiu en contacte: {len(contact_rows)}")
    return summary

def analyze_joined(md, kind: SimulationKind, run_dirs: list[Path], output_dir: Path, args: argparse.Namespace) -> bool:
    joined, frame_map = load_and_join_trajectories(md, kind, run_dirs)
    if joined is None:
        print(f"No s'ha pogut carregar cap trajectòria per {kind.name}.")
        return False

    kind_output_dir = output_dir / kind.name
    kind_output_dir.mkdir(parents=True, exist_ok=True)
    if not args.only_ligand_hbonds and not args.only_catalytic_metrics:
        try:
            joined = joined.image_molecules(inplace=False)
        except Exception as exc:
            print(f"Avís: no s'ha pogut recentrar/aplicar PBC a {kind.name}: {exc}")
    times_ns = frame_times_ns(joined.n_frames, args.report_interval, args.timestep_fs)
    system_label = kind.name.upper()

    write_rows(kind_output_dir / "frame_map.csv", ["joined_frame", "run", "frame_in_run"], frame_map)

    if args.only_catalytic_metrics:
        if kind.name != "holo":
            print("--only-catalytic-metrics només aplica al sistema holo; s'omet aquest sistema.")
            return False
        catalytic_summary = write_catalytic_metrics(md, joined, kind, kind_output_dir, times_ns, args)
        summary = [
            f"System: {kind.name}",
            f"Runs units: {', '.join(run_dir.name for run_dir in run_dirs)}",
            f"Frames totals: {joined.n_frames}",
            f"Temps final concatenat (ns): {times_ns[-1]:.6f}" if len(times_ns) else "Temps final concatenat (ns): 0.000000",
            *catalytic_summary,
        ]
        (kind_output_dir / "summary.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")
        print(f"Mètriques catalítiques desades a: {kind_output_dir}")
        return True

    if args.only_ligand_hbonds:
        if kind.name != "holo":
            print("--only-ligand-hbonds només aplica al sistema holo; s'omet aquest sistema.")
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
            "Ponts d'hidrogen lligand - proteïna",
            "Temps (ns)",
            "Ponts d'hidrogen",
        )
        summary = [
            f"System: {kind.name}",
            f"Runs units: {', '.join(run_dir.name for run_dir in run_dirs)}",
            f"Frames totals: {joined.n_frames}",
            f"Temps final concatenat (ns): {times_ns[-1]:.6f}" if len(times_ns) else "Temps final concatenat (ns): 0.000000",
            f"Ponts d'hidrogen lligand - proteïna mitjans ({args.ligand_resname}): {float(np.mean(ligand_hbond_counts)):.6f}",
        ]
        (kind_output_dir / "ligand_hydrogen_bonds_summary.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")
        print(f"Anàlisi de ponts proteïna-lligand desada a: {kind_output_dir}")
        return True

    if args.save_joined_dcd:
        joined.save_dcd(str(kind_output_dir / f"joined-{kind.name}.dcd"))

    protein_atoms = select_atoms(joined.topology, "protein", "proteïna")
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
        "RMSD",
        "Temps (ns)",
        "Distància (nm)",
    )
    save_line_plot(
        kind_output_dir / "radius_of_gyration.png",
        times_ns,
        rg,
        "Radi de gir",
        "Temps (ns)",
        "Distància (nm)",
    )
    save_rmsf_plot(
        kind_output_dir / "rmsf_ca.png",
        residues,
        rmsf,
        "RMSF",
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
            "Ponts d'hidrogen totals",
            "Temps (ns)",
            "Ponts d'hidrogen",
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
            "Ponts d'hidrogen lligand - proteïna",
            "Temps (ns)",
            "Ponts d'hidrogen",
        )
        hbond_summary.append(
            f"Ponts d'hidrogen lligand - proteïna mitjans ({args.ligand_resname}): "
            f"{float(np.mean(ligand_hbond_counts)):.6f}"
        )

    catalytic_summary = write_catalytic_metrics(md, joined, kind, kind_output_dir, times_ns, args)
    thermo_summary = write_thermo(kind, run_dirs, kind_output_dir, args.timestep_fs)
    summary = [
        f"System: {kind.name}",
        f"Runs: {', '.join(run_dir.name for run_dir in run_dirs)}",
        f"Topology: {kind.topology}",
        f"Frames totals: {joined.n_frames}",
        f"Temps final (ns): {times_ns[-1]:.6f}" if len(times_ns) else "Temps final (ns): 0.000000",
        f"RMSD mitjà de la proteïna (nm): {float(np.mean(rmsd)):.6f}" if len(rmsd) else "RMSD mitjà de la proteïna (nm): 0.000000",
        f"Radi de gir mitjà (nm): {float(np.mean(rg)):.6f}" if len(rg) else "Radi de gir mitjà (nm): 0.000000",
        f"RMSF C-α màxim (nm): {float(np.max(rmsf)):.6f}" if len(rmsf) else "RMSF C-α màxim (nm): 0.000000",
        *hbond_summary,
        *catalytic_summary,
        *thermo_summary,
    ]
    (kind_output_dir / "summary.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")
    print(f"Anàlisi concatenada desada a: {kind_output_dir}")
    return True


def main() -> None:
    args = parse_args()
    try:
        import mdtraj as md
    except ModuleNotFoundError as exc:
        raise SystemExit("Per executar aquest script cal instal·lar mdtraj. Exemple: conda install -c conda-forge mdtraj") from exc

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
        raise SystemExit("No s'ha analitzat cap trajectòria.")


if __name__ == "__main__":
    main()
