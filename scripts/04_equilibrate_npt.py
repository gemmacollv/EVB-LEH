from __future__ import annotations

import argparse
from pathlib import Path

from openmm import unit

from md_common import (
    add_standard_args,
    attach_reporters,
    build_simulation,
    build_system,
    load_modeller,
    summarize_energy,
    write_pdb,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Pas 4: equilibratge NPT.")
    parser.add_argument("--input-pdb", type=Path, default=Path("runs/03_nvt/nvt_final.pdb"))
    parser.add_argument("--output-dir", type=Path, default=Path("runs/04_npt"))
    parser.add_argument("--steps", type=int, default=100000)
    add_standard_args(parser)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    modeller = load_modeller(args.input_pdb)
    system = build_system(modeller.topology, args.temperature, with_barostat=True)
    simulation = build_simulation(
        modeller.topology, system, args.temperature, args.friction, args.timestep_fs
    )
    simulation.context.setPositions(modeller.positions)
    simulation.context.setVelocitiesToTemperature(args.temperature * unit.kelvin)
    attach_reporters(simulation, args.output_dir, args.steps, args.report_interval, "npt")
    simulation.step(args.steps)

    positions = simulation.context.getState(getPositions=True).getPositions()
    out_pdb = args.output_dir / "npt_final.pdb"
    write_pdb(simulation.topology, positions, out_pdb)
    (args.output_dir / "summary.txt").write_text(summarize_energy(simulation) + "\n", encoding="utf-8")
    print(f"Equilibratge NPT completat a {out_pdb}")


if __name__ == "__main__":
    main()
