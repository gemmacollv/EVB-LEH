from __future__ import annotations

import argparse
from pathlib import Path

from openmm import LangevinMiddleIntegrator, MonteCarloBarostat, unit
from openmm.app import (
    AmberInpcrdFile,
    AmberPrmtopFile,
    CheckpointReporter,
    DCDReporter,
    HBonds,
    NoCutoff,
    PDBFile,
    PME,
    Simulation,
    StateDataReporter,
)

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_STUDY_DIR = SCRIPT_DIR.parent / "study_runs"
DEFAULT_PRMTOP = DEFAULT_STUDY_DIR / "solvated_pdb" / "protein_only.prmtop"
DEFAULT_INPCRD = DEFAULT_STUDY_DIR / "solvated_pdb" / "protein_only.inpcrd"
DEFAULT_OUTPUT_DIR = DEFAULT_STUDY_DIR / "amber_simulation_apo"

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Executa una simulacio OpenMM a partir de fitxers AMBER."
    )
    parser.add_argument(
        "--study-dir",
        type=Path,
        default=DEFAULT_STUDY_DIR,
        help=f"Directori base de l'estudi. Per defecte: {DEFAULT_STUDY_DIR}.",
    )
    parser.add_argument(
        "--prmtop",
        type=Path,
        default=None,
        help=f"Fitxer PRMTOP. Per defecte: {DEFAULT_PRMTOP}.",
    )
    parser.add_argument(
        "--inpcrd",
        type=Path,
        default=None,
        help=f"Fitxer INPCRD. Per defecte: {DEFAULT_INPCRD}.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help=f"Directori de sortida. Per defecte: {DEFAULT_OUTPUT_DIR}.",
    )
    parser.add_argument("--steps", type=int, default=1000000)
    parser.add_argument("--equilibration-steps", type=int, default=1000)
    parser.add_argument("--temperature", type=float, default=300.0)
    parser.add_argument("--friction", type=float, default=1.0, help="1/ps")
    parser.add_argument("--timestep-fs", type=float, default=4.0)
    parser.add_argument("--hydrogen-mass", type=float, default=1.5, help="amu; posa 0 per desactivar HMR")
    parser.add_argument("--pressure", type=float, default=1.0, help="atm; posa 0 per desactivar el barostat")
    parser.add_argument("--barostat-interval", type=int, default=25)
    parser.add_argument("--ewald-error-tolerance", type=float, default=0.0005)
    parser.add_argument("--report-interval", type=int, default=10000)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.prmtop is None:
        args.prmtop = DEFAULT_PRMTOP
    if args.inpcrd is None:
        args.inpcrd = DEFAULT_INPCRD
    if args.output_dir is None:
        args.output_dir = DEFAULT_OUTPUT_DIR
    args.prmtop = args.prmtop.resolve()
    args.inpcrd = args.inpcrd.resolve()
    args.output_dir = args.output_dir.resolve()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    dcd_path = args.output_dir / "trajectory-apo.dcd"
    log_path = args.output_dir / "log-apo.txt"
    checkpoint_path = args.output_dir / "checkpoint-apo.chk"

    prmtop = AmberPrmtopFile(str(args.prmtop))
    inpcrd = AmberInpcrdFile(str(args.inpcrd))

    has_periodic_box = inpcrd.boxVectors is not None
    create_system_kwargs = {
        "nonbondedMethod": PME if has_periodic_box else NoCutoff,
        "constraints": HBonds,
        "rigidWater": True,
        "ewaldErrorTolerance": args.ewald_error_tolerance,
    }
    if has_periodic_box:
        create_system_kwargs["nonbondedCutoff"] = 1.0 * unit.nanometer
    if args.hydrogen_mass > 0:
        create_system_kwargs["hydrogenMass"] = args.hydrogen_mass * unit.amu

    system = prmtop.createSystem(**create_system_kwargs)
    if has_periodic_box and args.pressure > 0:
        system.addForce(
            MonteCarloBarostat(
                args.pressure * unit.atmospheres,
                args.temperature * unit.kelvin,
                args.barostat_interval,
            )
        )

    integrator = LangevinMiddleIntegrator(
        args.temperature * unit.kelvin,
        args.friction / unit.picosecond,
        args.timestep_fs * unit.femtoseconds,
    )
    simulation = Simulation(prmtop.topology, system, integrator)
    restarting = checkpoint_path.exists()
    if restarting:
        print(f"Carregant checkpoint: {checkpoint_path}")
        simulation.context.loadCheckpoint(checkpoint_path.read_bytes())
    else:
        print("Comencant simulacio nova.")
        simulation.context.setPositions(inpcrd.positions)
        if inpcrd.boxVectors is not None:
            simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
        print("Minimitzant energia...")
        simulation.minimizeEnergy()
        print("Equilibrant...")
        simulation.context.setVelocitiesToTemperature(args.temperature * unit.kelvin)
        simulation.step(args.equilibration_steps)

    simulation.reporters.append(
        DCDReporter(str(dcd_path), args.report_interval, append=restarting and dcd_path.exists())
    )
    simulation.reporters.append(
        StateDataReporter(
            str(log_path),
            args.report_interval,
            step=True,
            potentialEnergy=True,
            temperature=True,
            separator="\t",
            append=restarting and log_path.exists(),
        )
    )

    simulation.reporters.append(
        CheckpointReporter(str(checkpoint_path), args.report_interval)
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