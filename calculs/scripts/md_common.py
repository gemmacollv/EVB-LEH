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

# Defineix la carpeta principal del projecte i el PDB d’entrada per defecte.
# En el aquest cas leh_1nww.pdb

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_INPUT_PDB = ROOT / "input" / "leh_1nww.pdb"

# Llista dels 20 aminoàcids estàndard que es permet saber si un residu és proteïna o no.

PROTEIN_RESIDUES = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
}

# Crea etiquetes per poder distingir els diferets residus
def protein_residue_label(residue) -> str:
    chain_id = residue.chain.id or "?"
    return f"{residue.name}:{chain_id}:{residue.id}"

# Recorre tots els residus d’una topologia OpenMM i retorna els que no són aminoàcids estàndard.
def nonprotein_residues(topology) -> list:
    return [
        residue
        for residue in topology.residues()
        if residue.name not in PROTEIN_RESIDUES
    ]


def add_standard_args(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    parser.add_argument("--temperature", type=float, default=300.0)
    parser.add_argument("--friction", type=float, default=1.0, help="1/ps")
    parser.add_argument("--timestep-fs", type=float, default=2.0)
    parser.add_argument("--report-interval", type=int, default=500)
    return parser

# Crea el camp de forces, on amber14-all.xml defineix els paràmetres de la proteïna i amber14/tip3p.xml el model de l'aigua
def create_forcefield() -> ForceField:
    return ForceField("amber14-all.xml", "amber14/tip3p.xml")


def load_modeller(pdb_path: Path) -> Modeller:
    pdb = PDBFile(str(pdb_path))
    return Modeller(pdb.topology, pdb.positions)

# Càrrega un fitxer .pdb i el transforma en un objecte modeller d'OpenMM. 
# Prepara el sistema per abans de la simulació.
    # Carrega el PDB.
    # Crea el force field.
    # Afegeix hidrògens segons el pH.
    # Afegeix aigua TIP3P.
    # Neutralitza la càrrega.
    # Afegeix ions segons la força iònica.

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
