from __future__ import annotations

import argparse
from pathlib import Path

from md_common import add_standard_args, build_simulation, build_system, load_modeller, summarize_energy, write_pdb


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Pas 2: minimització d'energia.")
    parser.add_argument("--input-pdb", type=Path, default=Path("runs/01_prepare/prepared_solvated.pdb"))
    parser.add_argument("--output-dir", type=Path, default=Path("runs/02_minimize"))
    parser.add_argument("--max-iterations", type=int, default=5000)
    add_standard_args(parser)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    modeller = load_modeller(args.input_pdb)
    system = build_system(modeller.topology, args.temperature, with_barostat=False)
    simulation = build_simulation(
        modeller.topology, system, args.temperature, args.friction, args.timestep_fs
    )
    simulation.context.setPositions(modeller.positions)
    simulation.minimizeEnergy(maxIterations=args.max_iterations)

    positions = simulation.context.getState(getPositions=True).getPositions()
    out_pdb = args.output_dir / "minimized.pdb"
    write_pdb(simulation.topology, positions, out_pdb)
    (args.output_dir / "summary.txt").write_text(summarize_energy(simulation) + "\n", encoding="utf-8")
    print(f"Sistema minimitzat a {out_pdb}")


if __name__ == "__main__":
    main()
