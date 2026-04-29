from __future__ import annotations

import argparse
from pathlib import Path

from openmm import LangevinMiddleIntegrator, unit
from openmm.app import (
    AmberInpcrdFile,
    AmberPrmtopFile,
    DCDReporter,
    HBonds,
    NoCutoff,
    PDBFile,
    PME,
    Simulation,
    StateDataReporter,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Executa una simulacio OpenMM a partir de fitxers AMBER."
    )
    parser.add_argument(
        "--study-dir",
        type=Path,
        default=Path("study_runs"),
        help="Directori base de l'estudi.",
    )
    parser.add_argument(
        "--prmtop",
        type=Path,
        default=None,
        help="Fitxer PRMTOP. Per defecte: STUDY_DIR/amber_from_pdb/protonated.prmtop.",
    )
    parser.add_argument(
        "--inpcrd",
        type=Path,
        default=None,
        help="Fitxer INPCRD. Per defecte: STUDY_DIR/amber_from_pdb/protonated.inpcrd.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directori de sortida. Per defecte: STUDY_DIR/amber_simulation.",
    )
    parser.add_argument("--steps", type=int, default=5000)
    parser.add_argument("--temperature", type=float, default=300.0)
    parser.add_argument("--friction", type=float, default=1.0, help="1/ps")
    parser.add_argument("--timestep-fs", type=float, default=2.0)
    parser.add_argument("--report-interval", type=int, default=500)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.prmtop is None:
        args.prmtop = args.study_dir / "amber_from_pdb" / "protonated.prmtop"
    if args.inpcrd is None:
        args.inpcrd = args.study_dir / "amber_from_pdb" / "protonated.inpcrd"
    if args.output_dir is None:
        args.output_dir = args.study_dir / "amber_simulation"
    args.output_dir.mkdir(parents=True, exist_ok=True)

    prmtop = AmberPrmtopFile(str(args.prmtop))
    inpcrd = AmberInpcrdFile(str(args.inpcrd))

    has_periodic_box = inpcrd.boxVectors is not None
    create_system_kwargs = {
        "nonbondedMethod": PME if has_periodic_box else NoCutoff,
        "constraints": HBonds,
    }
    if has_periodic_box:
        create_system_kwargs["nonbondedCutoff"] = 1.0 * unit.nanometer
    system = prmtop.createSystem(**create_system_kwargs)

    integrator = LangevinMiddleIntegrator(
        args.temperature * unit.kelvin,
        args.friction / unit.picosecond,
        args.timestep_fs * unit.femtoseconds,
    )
    simulation = Simulation(prmtop.topology, system, integrator)
    simulation.context.setPositions(inpcrd.positions)
    if inpcrd.boxVectors is not None:
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
    simulation.minimizeEnergy()
    simulation.context.setVelocitiesToTemperature(args.temperature * unit.kelvin)

    simulation.reporters.append(
        DCDReporter(str(args.output_dir / "trajectory.dcd"), args.report_interval)
    )
    simulation.reporters.append(
        StateDataReporter(
            str(args.output_dir / "state.csv"),
            args.report_interval,
            step=True,
            potentialEnergy=True,
            temperature=True,
            separator=",",
        )
    )

    simulation.step(args.steps)

    state = simulation.context.getState(getPositions=True, getEnergy=True)
    with (args.output_dir / "final.pdb").open("w", encoding="utf-8") as handle:
        PDBFile.writeFile(prmtop.topology, state.getPositions(), handle)

    energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    print(f"Energia final: {energy:.3f} kJ/mol")
    print(f"Resultats a: {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
