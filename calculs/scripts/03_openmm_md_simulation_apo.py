from __future__ import annotations

import argparse
import signal
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
DEFAULT_STUDY_DIR = SCRIPT_DIR.parent / "prepared pdbs"
DEFAULT_ANALYSIS_DAT_DIR = SCRIPT_DIR.parent / "analisi" / "data"
DEFAULT_PRMTOP = DEFAULT_STUDY_DIR / "02_openmm_md_setup" / "protein_apo" / "protein_apo.prmtop"
DEFAULT_INPCRD = DEFAULT_STUDY_DIR / "02_openmm_md_setup" / "protein_apo" / "protein_apo.inpcrd"
DEFAULT_RUN_NAME = "run_001"
DEFAULT_OUTPUT_BASE_DIR = DEFAULT_ANALYSIS_DAT_DIR / "apo"
DEFAULT_OUTPUT_DIR = DEFAULT_OUTPUT_BASE_DIR / DEFAULT_RUN_NAME
DEFAULT_STEPS = 6000000
DEFAULT_STEP_CHUNK = 10000

STOP_REQUESTED = False


def request_stop(signum, _frame) -> None:
    global STOP_REQUESTED
    STOP_REQUESTED = True
    print(f"Rebuda senyal {signum}; s'acabara el bloc actual i es guardara final.pdb.", flush=True)


def install_signal_handlers() -> None:
    signal.signal(signal.SIGTERM, request_stop)
    signal.signal(signal.SIGINT, request_stop)



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
    parser.add_argument(
        "--run-name",
        default=DEFAULT_RUN_NAME,
        help=(
            "Nom del grup de resultats si no es passa --output-dir. "
            f"Per defecte: {DEFAULT_RUN_NAME}."
        ),
    )
    parser.add_argument(
        "--steps",
        type=int,
        default=DEFAULT_STEPS,
        help=f"Nombre màxim de passos de producció. Per defecte: {DEFAULT_STEPS}.",
    )
    parser.add_argument(
        "--step-chunk",
        type=int,
        default=DEFAULT_STEP_CHUNK,
        help=f"Passos per bloc abans de comprovar el temps. Per defecte: {DEFAULT_STEP_CHUNK}.",
    )
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


def write_final_pdb_image(topology, positions, output_path: Path) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        print("Matplotlib no esta instal-lat; no es pot guardar final.png.")
        return

    coords = positions.value_in_unit(unit.nanometer)
    protein_ca = []
    holo_atoms = []
    amino_acids = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
        "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
        "TYR", "VAL", "HID", "HIE", "HIP", "CYX", "ASH", "GLH", "LYN",
    }
    solvent_and_ions = {
        "HOH", "WAT", "TIP3", "SOL", "NA", "CL", "K", "MG", "CA", "ZN",
    }

    for atom in topology.atoms():
        residue_name = atom.residue.name.upper()
        coord = coords[atom.index]
        point = (coord.x, coord.y, coord.z)
        if residue_name in amino_acids and atom.name == "CA":
            protein_ca.append(point)
        elif residue_name not in amino_acids and residue_name not in solvent_and_ions:
            holo_atoms.append(point)

    if not protein_ca and not holo_atoms:
        print("No hi ha àtoms per dibuixar final.png.")
        return

    fig = plt.figure(figsize=(7, 6), dpi=220)
    ax = fig.add_subplot(111, projection="3d")
    all_points = []
    if protein_ca:
        xs, ys, zs = zip(*protein_ca)
        ax.plot(xs, ys, zs, color="#1f77b4", linewidth=1.2, label="Proteina C-α")
        ax.scatter(xs, ys, zs, color="#1f77b4", s=5)
        all_points.extend(protein_ca)
    if holo_atoms:
        xs, ys, zs = zip(*holo_atoms)
        ax.scatter(xs, ys, zs, color="#d62728", s=18, label="Holo/no proteic")
        all_points.extend(holo_atoms)

    xs, ys, zs = zip(*all_points)
    x_range = max(xs) - min(xs) or 1.0
    y_range = max(ys) - min(ys) or 1.0
    z_range = max(zs) - min(zs) or 1.0
    ax.set_xlim(min(xs), max(xs))
    ax.set_ylim(min(ys), max(ys))
    ax.set_zlim(min(zs), max(zs))
    try:
        ax.set_box_aspect((x_range, y_range, z_range))
    except Exception:
        pass
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("y (nm)")
    ax.set_zlabel("z (nm)")
    ax.set_title("Estructura final")
    ax.legend(loc="upper right")
    plt.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)
    print(f"Imatge final desada a: {output_path}")


def save_final_outputs(simulation, prmtop, checkpoint_path: Path, output_dir: Path, steps_done: int) -> None:
    checkpoint_path.write_bytes(simulation.context.createCheckpoint())
    print(f"Passos de producció executats en aquest job: {steps_done}")

    state = simulation.context.getState(getPositions=True, getEnergy=True)
    final_pdb_path = output_dir / "final.pdb"
    final_image_path = output_dir / "final.png"
    with final_pdb_path.open("w", encoding="utf-8") as handle:
        PDBFile.writeFile(prmtop.topology, state.getPositions(), handle)
    write_final_pdb_image(prmtop.topology, state.getPositions(), final_image_path)

    energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    print(f"Checkpoint final desat a: {checkpoint_path}")
    print(f"PDB final desat a: {final_pdb_path}")
    print(f"Energia final: {energy:.3f} kJ/mol")
    print(f"Resultats a: {output_dir.resolve()}")


def main() -> None:
    args = parse_args()
    if args.prmtop is None:
        args.prmtop = DEFAULT_PRMTOP
    if args.inpcrd is None:
        args.inpcrd = DEFAULT_INPCRD
    if args.output_dir is None:
        args.output_dir = DEFAULT_OUTPUT_BASE_DIR / args.run_name
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

    if args.step_chunk <= 0:
        raise SystemExit("--step-chunk ha de ser mes gran que 0.")

    install_signal_handlers()

    steps_done = 0
    try:
        while (args.steps is None or steps_done < args.steps) and not STOP_REQUESTED:
            if args.steps is None:
                chunk_steps = args.step_chunk
            else:
                chunk_steps = min(args.step_chunk, args.steps - steps_done)
            if chunk_steps <= 0:
                break

            simulation.step(chunk_steps)
            steps_done += chunk_steps
    finally:
        save_final_outputs(simulation, prmtop, checkpoint_path, args.output_dir, steps_done)


if __name__ == "__main__":
    main()