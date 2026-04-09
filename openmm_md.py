"""Exemple senzill de dinàmica molecular amb OpenMM.

El script:
1. Crea una molècula d'aigua a partir d'un PDB incrustat.
2. Genera una caixa periòdica de solvent TIP3P-FB.
3. Minimitza l'energia.
4. Executa una simulació curta de dinàmica molecular.
5. Desa trajectòria, estructura final i estadístiques bàsiques.
"""

from __future__ import annotations

import argparse
from io import StringIO
from pathlib import Path

from openmm import LangevinMiddleIntegrator, MonteCarloBarostat, unit
from openmm.app import (
    DCDReporter,
    ForceField,
    Modeller,
    PDBFile,
    PDBReporter,
    PME,
    Simulation,
    StateDataReporter,
)


WATER_PDB = """\
ATOM      1  O   HOH A   1       0.000   0.000   0.000  1.00  0.00           O
ATOM      2  H1  HOH A   1       0.095   0.000   0.000  1.00  0.00           H
ATOM      3  H2  HOH A   1      -0.024   0.092   0.000  1.00  0.00           H
TER
END
"""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Executa una simulació de dinàmica molecular bàsica amb OpenMM."
    )
    parser.add_argument(
        "--steps",
        type=int,
        default=5000,
        help="Nombre de passos d'integració de producció.",
    )
    parser.add_argument(
        "--temperature",
        type=float,
        default=300.0,
        help="Temperatura en Kelvin.",
    )
    parser.add_argument(
        "--friction",
        type=float,
        default=1.0,
        help="Coeficient de fricció en 1/ps per al termòstat de Langevin.",
    )
    parser.add_argument(
        "--timestep-fs",
        type=float,
        default=2.0,
        help="Pas temporal en femtosegons.",
    )
    parser.add_argument(
        "--padding-nm",
        type=float,
        default=1.0,
        help="Gruix addicional de solvent al voltant de la molècula en nm.",
    )
    parser.add_argument(
        "--report-interval",
        type=int,
        default=500,
        help="Interval de report en passos.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs"),
        help="Directori on es guardaran els resultats.",
    )
    return parser.parse_args()


def build_solvated_system(
    temperature_kelvin: float, padding_nm: float
) -> tuple[Simulation, Modeller]:
    pdb = PDBFile(StringIO(WATER_PDB))
    modeller = Modeller(pdb.topology, pdb.positions)

    forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
    modeller.addSolvent(
        forcefield,
        model="tip3pfb",
        padding=padding_nm * unit.nanometer,
        ionicStrength=0.0 * unit.molar,
        neutralize=False,
    )

    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=None,
        rigidWater=True,
    )
    system.addForce(MonteCarloBarostat(1.0 * unit.bar, temperature_kelvin * unit.kelvin))

    integrator = LangevinMiddleIntegrator(
        temperature_kelvin * unit.kelvin,
        1.0 / unit.picosecond,
        2.0 * unit.femtoseconds,
    )
    simulation = Simulation(modeller.topology, system, integrator)
    return simulation, modeller


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    simulation, modeller = build_solvated_system(args.temperature, args.padding_nm)
    integrator = simulation.integrator
    integrator.setFriction(args.friction / unit.picosecond)
    integrator.setStepSize(args.timestep_fs * unit.femtoseconds)

    simulation.context.setPositions(modeller.positions)
    simulation.minimizeEnergy()
    simulation.context.setVelocitiesToTemperature(args.temperature * unit.kelvin)

    state_report_path = args.output_dir / "state.csv"
    traj_report_path = args.output_dir / "trajectory.dcd"
    final_pdb_path = args.output_dir / "final.pdb"

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
            totalSteps=args.steps,
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
            totalSteps=args.steps,
        )
    )

    simulation.step(args.steps)

    with final_pdb_path.open("w", encoding="utf-8") as handle:
        PDBFile.writeFile(
            simulation.topology,
            simulation.context.getState(getPositions=True).getPositions(),
            handle,
        )

    final_state = simulation.context.getState(getEnergy=True)
    potential_energy = final_state.getPotentialEnergy().value_in_unit(
        unit.kilojoule_per_mole
    )
    print(f"\nSimulació completada. Energia potencial final: {potential_energy:.3f} kJ/mol")
    print(f"Resultats desats a: {args.output_dir.resolve()}")


if __name__ == "__main__":
    from sys import stdout

    main()
