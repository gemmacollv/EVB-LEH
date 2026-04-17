from __future__ import annotations

import argparse
from pathlib import Path

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


ROOT = Path(__file__).resolve().parent.parent
DEFAULT_INPUT_PDB = ROOT / "inputs" / "trp_cage_1l2y_model1.pdb"


def add_standard_args(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    parser.add_argument("--temperature", type=float, default=300.0)
    parser.add_argument("--friction", type=float, default=1.0, help="1/ps")
    parser.add_argument("--timestep-fs", type=float, default=2.0)
    parser.add_argument("--report-interval", type=int, default=500)
    return parser


def create_forcefield() -> ForceField:
    return ForceField("amber14-all.xml", "amber14/tip3p.xml")


def load_modeller(pdb_path: Path) -> Modeller:
    pdb = PDBFile(str(pdb_path))
    return Modeller(pdb.topology, pdb.positions)


def prepare_modeller(
    input_pdb: Path,
    padding_nm: float,
    ionic_strength_molar: float,
    ph: float,
) -> Modeller:
    modeller = load_modeller(input_pdb)
    forcefield = create_forcefield()
    modeller.addHydrogens(forcefield, pH=ph)
    modeller.addSolvent(
        forcefield,
        model="tip3p",
        padding=padding_nm * unit.nanometer,
        neutralize=True,
        ionicStrength=ionic_strength_molar * unit.molar,
    )
    return modeller


def build_system(topology, temperature_kelvin: float, with_barostat: bool):
    forcefield = create_forcefield()
    system = forcefield.createSystem(
        topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=HBonds,
        rigidWater=True,
    )
    if with_barostat:
        system.addForce(MonteCarloBarostat(1.0 * unit.bar, temperature_kelvin * unit.kelvin))
    return system


def build_simulation(topology, system, temperature_kelvin: float, friction_per_ps: float, timestep_fs: float) -> Simulation:
    integrator = LangevinMiddleIntegrator(
        temperature_kelvin * unit.kelvin,
        friction_per_ps / unit.picosecond,
        timestep_fs * unit.femtoseconds,
    )
    return Simulation(topology, system, integrator)


def write_pdb(topology, positions, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        PDBFile.writeFile(topology, positions, handle)


def attach_reporters(
    simulation: Simulation,
    output_dir: Path,
    total_steps: int,
    report_interval: int,
    prefix: str,
) -> None:
    simulation.reporters.append(
        StateDataReporter(
            str(output_dir / f"{prefix}_state.csv"),
            report_interval,
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
    simulation.reporters.append(
        DCDReporter(str(output_dir / f"{prefix}_trajectory.dcd"), report_interval)
    )


def summarize_energy(simulation: Simulation) -> str:
    state = simulation.context.getState(getEnergy=True)
    energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    return f"Energia potencial final: {energy:.3f} kJ/mol"
