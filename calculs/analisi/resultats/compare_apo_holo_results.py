#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_INPUT_DIR = SCRIPT_DIR / "joined"
DEFAULT_OUTPUT_DIR = SCRIPT_DIR / "comparacio_apo_holo"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compara els resultats APO i HOLO i genera grafics "
            "superposats per veure les diferencies principals."
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
        help=f"Directori on desar la comparacio. Per defecte: {DEFAULT_OUTPUT_DIR}.",
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
    plt.figure(figsize=(8, 5), dpi=300)
    plt.plot(apo_x, apo_y, label="APO", linewidth=1.4)
    plt.plot(holo_x, holo_y, label="HOLO", linewidth=1.4)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def save_rmsf_plot(output_path: Path, apo: np.ndarray, holo: np.ndarray) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n_points = min(len(apo), len(holo))
    x_values = np.arange(1, n_points + 1)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(9, 5), dpi=300)
    plt.plot(x_values, apo[:n_points], label="APO", linewidth=1.4)
    plt.plot(x_values, holo[:n_points], label="HOLO", linewidth=1.4)
    plt.title("Comparació RMSF")
    plt.xlabel("Residus C-α")
    plt.ylabel("RMSF (nm)")
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


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
        raise SystemExit("Per executar aquest script cal instal-lar matplotlib.") from exc

    comparisons = [
        (
            "rmsd.csv",
            "time_ns",
            "rmsd_nm",
            "rmsd_apo_holo.png",
            "Comparacio RMSD",
            "Temps (ns)",
            "RMSD (nm)",
        ),
        (
            "radius_of_gyration.csv",
            "time_ns",
            "rg_nm",
            "radi_gir_apo_holo.png",
            "Comparacio radi de gir",
            "Temps (ns)",
            "Radi de gir (nm)",
        ),
        (
            "hydrogen_bonds.csv",
            "time_ns",
            "n_hydrogen_bonds",
            "ponts_hidrogen_apo_holo.png",
            "Comparacio ponts d'hidrogen",
            "Temps (ns)",
            "Nombre de ponts",
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
    save_rmsf_plot(output_dir / "rmsf_apo_holo.png", apo_rmsf, holo_rmsf)

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
        "Comparacio APO vs HOLO",
        f"Input: {input_dir}",
        f"Output: {output_dir}",
        "",
        "Grafics generats:",
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
    print(f"Comparacio desada a: {output_dir}")


if __name__ == "__main__":
    main()
