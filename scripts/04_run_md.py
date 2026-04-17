from __future__ import annotations

import argparse
from pathlib import Path

from openmm import unit
from openmm.app import PDBFile

from md_common import attach_reporters, build_simulation, build_system, prepare_modeller, summarize_energy, write_pdb


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Pas 4: executar la dinamica molecular.")
    parser.add_argument(
        "--input-pdb",
        type=Path,
        default=Path("study_runs/03_variants/protein_only.pdb"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("study_runs/04_md"),
    )
    parser.add_argument("--ph", type=float, default=7.0)
    parser.add_argument("--padding-nm", type=float, default=1.2)
    parser.add_argument("--ionic-strength", type=float, default=0.15, help="M")
    parser.add_argument("--temperature", type=float, default=300.0)
    parser.add_argument("--friction", type=float, default=1.0, help="1/ps")
    parser.add_argument("--timestep-fs", type=float, default=2.0)
    parser.add_argument("--report-interval", type=int, default=500)
    parser.add_argument("--minimize-iterations", type=int, default=5000)
    parser.add_argument("--nvt-steps", type=int, default=25000)
    parser.add_argument("--npt-steps", type=int, default=100000)
    parser.add_argument("--production-steps", type=int, default=250000)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    modeller = prepare_modeller(
        input_pdb=args.input_pdb,
        padding_nm=args.padding_nm,
        ionic_strength_molar=args.ionic_strength,
        ph=args.ph,
    )
    solvated_pdb = args.output_dir / "solvated_initial.pdb"
    write_pdb(modeller.topology, modeller.positions, solvated_pdb)

    total_steps = args.nvt_steps + args.npt_steps + args.production_steps
    current_positions = modeller.positions

    # Minimitzacio
    min_system = build_system(modeller.topology, args.temperature, with_barostat=False)
    min_sim = build_simulation(
        modeller.topology, min_system, args.temperature, args.friction, args.timestep_fs
    )
    min_sim.context.setPositions(current_positions)
    min_sim.minimizeEnergy(maxIterations=args.minimize_iterations)
    current_positions = min_sim.context.getState(getPositions=True).getPositions()
    write_pdb(modeller.topology, current_positions, args.output_dir / "minimized.pdb")

    # NVT
    nvt_system = build_system(modeller.topology, args.temperature, with_barostat=False)
    nvt_sim = build_simulation(
        modeller.topology, nvt_system, args.temperature, args.friction, args.timestep_fs
    )
    nvt_sim.context.setPositions(current_positions)
    nvt_sim.context.setVelocitiesToTemperature(args.temperature * unit.kelvin)
    attach_reporters(nvt_sim, args.output_dir, total_steps, args.report_interval, "md")
    if args.nvt_steps > 0:
        nvt_sim.step(args.nvt_steps)
    current_positions = nvt_sim.context.getState(getPositions=True).getPositions()
    write_pdb(modeller.topology, current_positions, args.output_dir / "after_nvt.pdb")

    # NPT + produccio en una sola simulacio per mantenir la continuitat
    npt_system = build_system(modeller.topology, args.temperature, with_barostat=True)
    md_sim = build_simulation(
        modeller.topology, npt_system, args.temperature, args.friction, args.timestep_fs
    )
    md_sim.context.setPositions(current_positions)
    md_sim.context.setVelocitiesToTemperature(args.temperature * unit.kelvin)
    attach_reporters(md_sim, args.output_dir, total_steps, args.report_interval, "production")
    if args.npt_steps > 0:
        md_sim.step(args.npt_steps)
    write_pdb(
        modeller.topology,
        md_sim.context.getState(getPositions=True).getPositions(),
        args.output_dir / "after_npt.pdb",
    )
    if args.production_steps > 0:
        md_sim.step(args.production_steps)

    final_state = md_sim.context.getState(getPositions=True, getEnergy=True)
    final_pdb = args.output_dir / "production_final.pdb"
    write_pdb(modeller.topology, final_state.getPositions(), final_pdb)

    (args.output_dir / "summary.txt").write_text(
        "\n".join(
            [
                f"Input PDB: {args.input_pdb}",
                f"Solvated PDB: {solvated_pdb}",
                f"Final PDB: {final_pdb}",
                f"NVT steps: {args.nvt_steps}",
                f"NPT steps: {args.npt_steps}",
                f"Production steps: {args.production_steps}",
                summarize_energy(md_sim),
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"MD completada. Resultats a {args.output_dir}")


if __name__ == "__main__":
    main()
