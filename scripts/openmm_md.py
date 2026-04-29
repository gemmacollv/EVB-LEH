"""Dinamica molecular amb OpenMM a partir d'un PDB extern.

Per defecte utilitza el PDB d'exemple del Trp-cage, pero es pot executar amb
qualsevol altre PDB passant `--input-pdb`.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from sys import stdout

from openmm import LangevinMiddleIntegrator, MonteCarloBarostat, unit
from openmm.app import (
    DCDReporter,
    ForceField,
    HBonds,
    Modeller,
    PDBFile,
    PME,
    Simulation,
    StateDataReporter,
)


DEFAULT_INPUT_PDB = Path("inputs/trp_cage_1l2y_model1.pdb")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Executa una simulacio de dinamica molecular amb OpenMM sobre un PDB extern."
    )
    parser.add_argument(
        "--input-pdb",
        type=Path,
        default=DEFAULT_INPUT_PDB,
        help="Fitxer PDB d'entrada. Per defecte usa el Trp-cage d'exemple.",
    )
    parser.add_argument("--steps", type=int, default=5000, help="Passos de produccio.")
    parser.add_argument(
        "--equilibration-steps",
        type=int,
        default=1000,
        help="Passos curts d'equilibratge abans de la produccio.",
    )
    parser.add_argument("--temperature", type=float, default=300.0, help="Temperatura en Kelvin.")
    parser.add_argument(
        "--friction",
        type=float,
        default=1.0,
        help="Coeficient de friccio del termostat en 1/ps.",
    )
    parser.add_argument("--timestep-fs", type=float, default=2.0, help="Pas temporal en fs.")
    parser.add_argument(
        "--padding-nm",
        type=float,
        default=1.2,
        help="Gruix de solvent al voltant de la proteina en nm.",
    )
    parser.add_argument(
        "--ionic-strength",
        type=float,
        default=0.0,
        help="Forca ionica en molar.",
    )
    parser.add_argument("--ph", type=float, default=7.0, help="pH per afegir hidrogens.")
    parser.add_argument("--report-interval", type=int, default=500, help="Interval de report.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs"),
        help="Directori de sortida.",
    )
    return parser.parse_args()


def build_simulation(
    input_pdb: Path,
    temperature_kelvin: float,
    friction_per_ps: float,
    timestep_fs: float,
    padding_nm: float,
    ionic_strength_molar: float,
    ph: float,
) -> tuple[Simulation, Modeller]:
    pdb = PDBFile(str(input_pdb))
    modeller = Modeller(pdb.topology, pdb.positions)

    forcefield = ForceField("amber14-all.xml", "amber14/tip3p.xml")
    modeller.addHydrogens(forcefield, pH=ph)
    modeller.addSolvent(
        forcefield,
        model="tip3p",
        padding=padding_nm * unit.nanometer,
        neutralize=True,
        ionicStrength=ionic_strength_molar * unit.molar,
    )

    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=HBonds,
        rigidWater=True,
    )
    system.addForce(MonteCarloBarostat(1.0 * unit.bar, temperature_kelvin * unit.kelvin))

    integrator = LangevinMiddleIntegrator(
        temperature_kelvin * unit.kelvin,
        friction_per_ps / unit.picosecond,
        timestep_fs * unit.femtoseconds,
    )
    simulation = Simulation(modeller.topology, system, integrator)
    return simulation, modeller


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    simulation, modeller = build_simulation(
        input_pdb=args.input_pdb,
        temperature_kelvin=args.temperature,
        friction_per_ps=args.friction,
        timestep_fs=args.timestep_fs,
        padding_nm=args.padding_nm,
        ionic_strength_molar=args.ionic_strength,
        ph=args.ph,
    )

    simulation.context.setPositions(modeller.positions)
    simulation.minimizeEnergy()
    simulation.context.setVelocitiesToTemperature(args.temperature * unit.kelvin)

    solvated_pdb_path = args.output_dir / "solvated_initial.pdb"
    with solvated_pdb_path.open("w", encoding="utf-8") as handle:
        PDBFile.writeFile(modeller.topology, modeller.positions, handle)

    state_report_path = args.output_dir / "state.csv"
    traj_report_path = args.output_dir / "trajectory.dcd"
    final_pdb_path = args.output_dir / "final.pdb"
    total_steps = args.steps + args.equilibration_steps

    simulation.reporters.append(
        StateDataReporter(
            str(state_report_path),
            args.report_interval,
            step=True,
            time=True,
            potentialEnergy=True,
            kineticEnergy=True,
            totalEnergy=True,
            temperature=True,
            density=True,
            volume=True,
            progress=True,
            remainingTime=True,
            speed=True,
            totalSteps=total_steps,
            separator=",",
        )
    )
    simulation.reporters.append(DCDReporter(str(traj_report_path), args.report_interval))
    simulation.reporters.append(
        StateDataReporter(
            stdout,
            args.report_interval,
            step=True,
            potentialEnergy=True,
            temperature=True,
            speed=True,
            progress=True,
            totalSteps=total_steps,
        )
    )

    if args.equilibration_steps > 0:
        simulation.step(args.equilibration_steps)
    simulation.step(args.steps)

    final_state = simulation.context.getState(getPositions=True, getEnergy=True)
    with final_pdb_path.open("w", encoding="utf-8") as handle:
        PDBFile.writeFile(simulation.topology, final_state.getPositions(), handle)

    potential_energy = final_state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    print(f"\nSimulacio completada. Energia potencial final: {potential_energy:.3f} kJ/mol")
    print(f"PDB d'entrada: {args.input_pdb}")
    print(f"Resultats desats a: {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
