#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_INPUT_DIR = SCRIPT_DIR / "joined"
DEFAULT_OUTPUT_DIR = SCRIPT_DIR / "comparacio_apo_holo"
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compara els resultats APO i HOLO i genera gràfics "
            "superposats per veure les diferències principals."
        )
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=DEFAULT_INPUT_DIR,
        help=f"Directori amb joined/apo i joined/holo. Per defecte: {DEFAULT_INPUT_DIR}.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Directori on desar la comparació. Per defecte: {DEFAULT_OUTPUT_DIR}.",
    )
    return parser.parse_args()


def read_numeric_csv(path: Path, x_column: str, y_column: str) -> tuple[np.ndarray, np.ndarray]:
    if not path.exists():
        raise FileNotFoundError(f"No existeix: {path}")

    x_values = []
    y_values = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            x_values.append(float(row[x_column]))
            y_values.append(float(row[y_column]))
    return np.array(x_values, dtype=float), np.array(y_values, dtype=float)


def read_rmsf(path: Path) -> tuple[list[str], np.ndarray]:
    if not path.exists():
        raise FileNotFoundError(f"No existeix: {path}")

    residues = []
    values = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            residues.append(row["residue"])
            values.append(float(row["rmsf_nm"]))
    return residues, np.array(values, dtype=float)


def read_summary(path: Path) -> dict[str, str]:
    values = {}
    if not path.exists():
        return values
    for line in path.read_text(encoding="utf-8").splitlines():
        if ":" not in line:
            continue
        key, value = line.split(":", 1)
        values[key.strip()] = value.strip()
    return values


def write_rows(path: Path, header: list[str], rows: list[list[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def save_overlay_plot(
    output_path: Path,
    apo_x: np.ndarray,
    apo_y: np.ndarray,
    holo_x: np.ndarray,
    holo_y: np.ndarray,
    title: str,
    xlabel: str,
    ylabel: str,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    output_path.parent.mkdir(parents=True, exist_ok=True)
    combined = np.concatenate([apo_y, holo_y]) if len(apo_y) or len(holo_y) else np.array([0.0, 1.0])
    finite = combined[np.isfinite(combined)]
    y_min = float(np.min(finite)) if len(finite) else 0.0
    y_max = float(np.max(finite)) if len(finite) else 1.0
    y_span = max(y_max - y_min, 0.01)
    y_floor = max(0.0, y_min - 0.08 * y_span)

    fig, ax = plt.subplots(figsize=(9, 5.2), dpi=300)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.plot(apo_x, apo_y, label="APO", linewidth=2.0, color="#1f77b4")
    ax.plot(holo_x, holo_y, label="HOLO", linewidth=2.0, color="#ff7f0e")
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
    ax.legend(frameon=True, facecolor="white", edgecolor="#d6d6d6")
    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)


def save_rmsf_plot(
    output_path: Path,
    apo_residues: list[str],
    apo: np.ndarray,
    holo_residues: list[str],
    holo: np.ndarray,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n_points = min(len(apo), len(holo))
    x_values = np.arange(1, n_points + 1)
    active_points = [
        (idx + 1, residue, ACTIVE_SITE_DISPLAY_LABELS[residue])
        for idx, residue in enumerate(holo_residues[:n_points])
        if residue in ACTIVE_SITE_DISPLAY_LABELS
    ]
    if not active_points:
        active_points = [
            (idx + 1, residue, ACTIVE_SITE_DISPLAY_LABELS[residue])
            for idx, residue in enumerate(apo_residues[:n_points])
            if residue in ACTIVE_SITE_DISPLAY_LABELS
        ]
    site_colors = ["#1f77b4", "#ff7f0e"]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(12, 6.2), dpi=300)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    combined = np.concatenate([apo[:n_points], holo[:n_points]]) if n_points else np.array([0.0, 1.0])
    y_min = float(np.min(combined))
    y_max = float(np.max(combined))
    y_span = max(y_max - y_min, 0.01)
    y_floor = max(0.0, y_min - 0.08 * y_span)
    y_ceiling = y_max + 0.18 * y_span
    ax.set_ylim(y_floor, y_ceiling)
    ax.plot(x_values, apo[:n_points], label="APO", linewidth=2.0, color="#1f77b4")
    ax.plot(x_values, holo[:n_points], label="HOLO", linewidth=2.0, color="#ff7f0e")
    label_lanes = [0.74, 0.56, 0.38]
    for color_idx, (position, _residue, label) in enumerate(active_points):
        color = site_colors[color_idx % len(site_colors)]
        point_y = max(float(apo[position - 1]), float(holo[position - 1]))
        label_y = y_min + y_span * label_lanes[color_idx % len(label_lanes)]
        label_x = min(position + 1.2, n_points)
        ax.axvspan(position - 0.55, position + 0.55, color=color, alpha=0.18, linewidth=0)
        ax.axvline(position, color=color, linewidth=0.9, alpha=0.5)
        ax.scatter(
            [position, position],
            [apo[position - 1], holo[position - 1]],
            color=[color],
            edgecolor="white",
            linewidth=0.6,
            s=32,
            zorder=4,
        )
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
    ax.set_title("RMSF", fontsize=14, weight="bold", pad=10)
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
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)


def stats(values: np.ndarray) -> tuple[float, float, float]:
    return float(np.mean(values)), float(np.min(values)), float(np.max(values))


def main() -> None:
    args = parse_args()
    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()
    apo_dir = input_dir / "apo"
    holo_dir = input_dir / "holo"

    try:
        import matplotlib  # noqa: F401
    except ModuleNotFoundError as exc:
        raise SystemExit("Per executar aquest script cal instal·lar matplotlib.") from exc

    comparisons = [
        (
            "rmsd.csv",
            "time_ns",
            "rmsd_nm",
            "rmsd_apo_holo.png",
            "RMSD",
            "Temps (ns)",
            "Distància (nm)",
        ),
        (
            "radius_of_gyration.csv",
            "time_ns",
            "rg_nm",
            "radi_gir_apo_holo.png",
            "Radi de gir",
            "Temps (ns)",
            "Distància (nm)",
        ),
        (
            "hydrogen_bonds.csv",
            "time_ns",
            "n_hydrogen_bonds",
            "ponts_hidrogen_apo_holo.png",
            "Ponts d'hidrogen totals",
            "Temps (ns)",
            "Ponts d'hidrogen",
        ),
    ]

    summary_rows = []
    for filename, x_col, y_col, plot_name, title, xlabel, ylabel in comparisons:
        apo_x, apo_y = read_numeric_csv(apo_dir / filename, x_col, y_col)
        holo_x, holo_y = read_numeric_csv(holo_dir / filename, x_col, y_col)
        save_overlay_plot(output_dir / plot_name, apo_x, apo_y, holo_x, holo_y, title, xlabel, ylabel)

        apo_mean, apo_min, apo_max = stats(apo_y)
        holo_mean, holo_min, holo_max = stats(holo_y)
        summary_rows.append([filename, "APO", apo_mean, apo_min, apo_max])
        summary_rows.append([filename, "HOLO", holo_mean, holo_min, holo_max])
        summary_rows.append([filename, "HOLO-APO mitjana", holo_mean - apo_mean, "", ""])

    apo_residues, apo_rmsf = read_rmsf(apo_dir / "rmsf_ca.csv")
    holo_residues, holo_rmsf = read_rmsf(holo_dir / "rmsf_ca.csv")
    save_rmsf_plot(output_dir / "rmsf_apo_holo.png", apo_residues, apo_rmsf, holo_residues, holo_rmsf)

    n_residues = min(len(apo_rmsf), len(holo_rmsf))
    rmsf_diff = holo_rmsf[:n_residues] - apo_rmsf[:n_residues]
    write_rows(
        output_dir / "rmsf_diferencies.csv",
        ["residue_apo", "residue_holo", "apo_rmsf_nm", "holo_rmsf_nm", "holo_minus_apo_nm"],
        [
            [apo_residues[idx], holo_residues[idx], apo_rmsf[idx], holo_rmsf[idx], rmsf_diff[idx]]
            for idx in range(n_residues)
        ],
    )
    summary_rows.append(["rmsf_ca.csv", "APO", *stats(apo_rmsf)])
    summary_rows.append(["rmsf_ca.csv", "HOLO", *stats(holo_rmsf)])
    summary_rows.append(["rmsf_ca.csv", "HOLO-APO mitjana", float(np.mean(rmsf_diff)), "", ""])

    write_rows(
        output_dir / "resum_comparatiu.csv",
        ["metric", "system", "mean", "min", "max"],
        summary_rows,
    )

    apo_summary = read_summary(apo_dir / "summary.txt")
    holo_summary = read_summary(holo_dir / "summary.txt")
    report = [
        "Comparació APO vs HOLO",
        f"Input: {input_dir}",
        f"Output: {output_dir}",
        "",
        "Gràfics generats:",
        "- rmsd_apo_holo.png",
        "- radi_gir_apo_holo.png",
        "- rmsf_apo_holo.png",
        "- ponts_hidrogen_apo_holo.png",
        "",
        "Valors del resum original:",
    ]
    for key in sorted(set(apo_summary) | set(holo_summary)):
        report.append(f"{key}: APO={apo_summary.get(key, 'n/a')} | HOLO={holo_summary.get(key, 'n/a')}")

    (output_dir / "resum_comparatiu.txt").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(f"Comparació desada a: {output_dir}")


if __name__ == "__main__":
    main()
