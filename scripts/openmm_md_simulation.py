from pathlib import Path

from openmm import LangevinMiddleIntegrator, unit
from openmm.app import (
    AmberInpcrdFile,
    AmberPrmtopFile,
    DCDReporter,
    HBonds,
    NoCutoff,
    PDBFile,
    Simulation,
    StateDataReporter,
)


prmtop = AmberPrmtopFile("study_runs/amber_from_pdb/protonated.prmtop")
inpcrd = AmberInpcrdFile("study_runs/amber_from_pdb/protonated.inpcrd")
output_dir = Path("study_runs/amber_simulation")
output_dir.mkdir(parents=True, exist_ok=True)

system = prmtop.createSystem(nonbondedMethod=NoCutoff, constraints=HBonds)
integrator = LangevinMiddleIntegrator(
    300 * unit.kelvin,
    1.0 / unit.picosecond,
    2.0 * unit.femtoseconds,
)
simulation = Simulation(prmtop.topology, system, integrator)
simulation.context.setPositions(inpcrd.positions)
simulation.minimizeEnergy()
simulation.context.setVelocitiesToTemperature(300 * unit.kelvin)

simulation.reporters.append(DCDReporter(str(output_dir / "trajectory.dcd"), 500))
simulation.reporters.append(
    StateDataReporter(
        str(output_dir / "state.csv"),
        500,
        step=True,
        potentialEnergy=True,
        temperature=True,
        separator=",",
    )
)

simulation.step(5000)

state = simulation.context.getState(getPositions=True, getEnergy=True)
with (output_dir / "final.pdb").open("w", encoding="utf-8") as handle:
    PDBFile.writeFile(prmtop.topology, state.getPositions(), handle)

energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
print(f"Energia final: {energy:.3f} kJ/mol")
print(f"Resultats a: {output_dir.resolve()}")
