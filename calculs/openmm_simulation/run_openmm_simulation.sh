#!/bin/bash
#SBATCH --job-name=openmm
#SBATCH --error=%j.err.txt
#SBATCH --output=%j.out.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=normal2

python3 <<EOF
from openmm import *
from openmm.app import *
from openmm.unit import *
# Input files
pdb = PDBFile('protonated-processed.pdb')
forcefield = ForceField('amber19-all.xml', 'amber19/tip3pfb.xml')
# System configuration
nonbondedMethod = PME
nonbondedCutoff = 1.0 * nanometers
constraints = HBonds
rigidWater = True
constraintTolerance = 1e-6
hydrogenMass = 1.5 * amu
ewaldErrorTolerance = 0.0005
# Integration options
dt = 0.004 * picoseconds
temperature = 300 * kelvin
friction = 1.0 / picosecond
pressure = 1.0 * atmospheres
barostatInterval = 25
# Simulation options
steps = 10000
equilibrationSteps = 1000
dcdReporter = DCDReporter('trajectory.dcd', 1000)
dataReporter = StateDataReporter(
    'log.txt',
    1000,
    totalSteps=steps,
    step=True,
    speed=True,
    progress=True,
    potentialEnergy=True,
    temperature=True,
    separator='\t'
)
checkpointReporter = CheckpointReporter('checkpoint.chk', 10000)
# Prepare simulation
print('Building system...')
topology = pdb.topology
positions = pdb.positions
system = forcefield.createSystem(
    topology,
    nonbondedMethod=nonbondedMethod,
    nonbondedCutoff=nonbondedCutoff,
    constraints=constraints,
    rigidWater=rigidWater,
    ewaldErrorTolerance=ewaldErrorTolerance,
    hydrogenMass=hydrogenMass
)
system.addForce(
    MonteCarloBarostat(
        pressure,
        temperature,
        barostatInterval
    )
)
integrator = LangevinMiddleIntegrator(
    temperature,
    friction,
    dt
)
integrator.setConstraintTolerance(constraintTolerance)
simulation = Simulation(
    topology,
    system,
    integrator
)
simulation.context.setPositions(positions)
# Energy minimization
print('Performing energy minimization...')
simulation.minimizeEnergy()
# Equilibration
print('Equilibrating...')
simulation.context.setVelocitiesToTemperature(temperature)
simulation.step(equilibrationSteps)
# Production run
print('Simulating...')
simulation.reporters.append(dcdReporter)
simulation.reporters.append(dataReporter)
simulation.reporters.append(checkpointReporter)
simulation.currentStep = 0
simulation.step(steps)
print('Simulation complete.')
EOF