"""Dinàmica molecular amb OpenMM sobre una proteïna real molt petita.

Aquest exemple utilitza el Trp-cage miniprotein (PDB: 1L2Y, model 1),
una proteïna de 20 residus que és prou petita per a proves ràpides.
"""

from __future__ import annotations

import argparse
from io import StringIO
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


TRP_CAGE_PDB = """\
ATOM      1  N   ASN A   1      -8.901   4.127  -0.555  1.00  0.00           N
ATOM      2  CA  ASN A   1      -8.608   3.135  -1.618  1.00  0.00           C
ATOM      3  C   ASN A   1      -7.117   2.964  -1.897  1.00  0.00           C
ATOM      4  O   ASN A   1      -6.634   1.849  -1.758  1.00  0.00           O
ATOM      5  CB  ASN A   1      -9.437   3.396  -2.889  1.00  0.00           C
ATOM      6  CG  ASN A   1     -10.915   3.130  -2.611  1.00  0.00           C
ATOM      7  OD1 ASN A   1     -11.269   2.700  -1.524  1.00  0.00           O
ATOM      8  ND2 ASN A   1     -11.806   3.406  -3.543  1.00  0.00           N
ATOM     17  N   LEU A   2      -6.379   4.031  -2.228  1.00  0.00           N
ATOM     18  CA  LEU A   2      -4.923   4.002  -2.452  1.00  0.00           C
ATOM     19  C   LEU A   2      -4.136   3.187  -1.404  1.00  0.00           C
ATOM     20  O   LEU A   2      -3.391   2.274  -1.760  1.00  0.00           O
ATOM     21  CB  LEU A   2      -4.411   5.450  -2.619  1.00  0.00           C
ATOM     22  CG  LEU A   2      -4.795   6.450  -1.495  1.00  0.00           C
ATOM     23  CD1 LEU A   2      -3.612   6.803  -0.599  1.00  0.00           C
ATOM     24  CD2 LEU A   2      -5.351   7.748  -2.084  1.00  0.00           C
ATOM     36  N   TYR A   3      -4.354   3.455  -0.111  1.00  0.00           N
ATOM     37  CA  TYR A   3      -3.690   2.738   0.981  1.00  0.00           C
ATOM     38  C   TYR A   3      -4.102   1.256   1.074  1.00  0.00           C
ATOM     39  O   TYR A   3      -3.291   0.409   1.442  1.00  0.00           O
ATOM     40  CB  TYR A   3      -3.964   3.472   2.302  1.00  0.00           C
ATOM     41  CG  TYR A   3      -2.824   3.339   3.290  1.00  0.00           C
ATOM     42  CD1 TYR A   3      -2.746   2.217   4.138  1.00  0.00           C
ATOM     43  CD2 TYR A   3      -1.820   4.326   3.332  1.00  0.00           C
ATOM     44  CE1 TYR A   3      -1.657   2.076   5.018  1.00  0.00           C
ATOM     45  CE2 TYR A   3      -0.725   4.185   4.205  1.00  0.00           C
ATOM     46  CZ  TYR A   3      -0.639   3.053   5.043  1.00  0.00           C
ATOM     47  OH  TYR A   3       0.433   2.881   5.861  1.00  0.00           O
ATOM     57  N   ILE A   4      -5.342   0.925   0.689  1.00  0.00           N
ATOM     58  CA  ILE A   4      -5.857  -0.449   0.613  1.00  0.00           C
ATOM     59  C   ILE A   4      -5.089  -1.221  -0.470  1.00  0.00           C
ATOM     60  O   ILE A   4      -4.621  -2.334  -0.226  1.00  0.00           O
ATOM     61  CB  ILE A   4      -7.386  -0.466   0.343  1.00  0.00           C
ATOM     62  CG1 ILE A   4      -8.197   0.540   1.197  1.00  0.00           C
ATOM     63  CG2 ILE A   4      -7.959  -1.884   0.501  1.00  0.00           C
ATOM     64  CD1 ILE A   4      -8.019   0.412   2.715  1.00  0.00           C
ATOM     76  N   GLN A   5      -4.907  -0.601  -1.645  1.00  0.00           N
ATOM     77  CA  GLN A   5      -4.122  -1.167  -2.743  1.00  0.00           C
ATOM     78  C   GLN A   5      -2.629  -1.321  -2.390  1.00  0.00           C
ATOM     79  O   GLN A   5      -1.986  -2.240  -2.884  1.00  0.00           O
ATOM     80  CB  GLN A   5      -4.292  -0.313  -4.013  1.00  0.00           C
ATOM     81  CG  GLN A   5      -4.244  -1.171  -5.290  1.00  0.00           C
ATOM     82  CD  GLN A   5      -5.576  -1.860  -5.585  1.00  0.00           C
ATOM     83  OE1 GLN A   5      -5.769  -3.044  -5.335  1.00  0.00           O
ATOM     84  NE2 GLN A   5      -6.532  -1.146  -6.152  1.00  0.00           N
ATOM     93  N   TRP A   6      -2.074  -0.459  -1.528  1.00  0.00           N
ATOM     94  CA  TRP A   6      -0.716  -0.631  -0.993  1.00  0.00           C
ATOM     95  C   TRP A   6      -0.631  -1.766   0.044  1.00  0.00           C
ATOM     96  O   TRP A   6       0.295  -2.579  -0.004  1.00  0.00           O
ATOM     97  CB  TRP A   6      -0.221   0.703  -0.417  1.00  0.00           C
ATOM     98  CG  TRP A   6       1.148   0.652   0.194  1.00  0.00           C
ATOM     99  CD1 TRP A   6       2.319   0.664  -0.482  1.00  0.00           C
ATOM    100  CD2 TRP A   6       1.508   0.564   1.606  1.00  0.00           C
ATOM    101  NE1 TRP A   6       3.371   0.560   0.411  1.00  0.00           N
ATOM    102  CE2 TRP A   6       2.928   0.515   1.710  1.00  0.00           C
ATOM    103  CE3 TRP A   6       0.779   0.524   2.812  1.00  0.00           C
ATOM    104  CZ2 TRP A   6       3.599   0.445   2.938  1.00  0.00           C
ATOM    105  CZ3 TRP A   6       1.439   0.433   4.053  1.00  0.00           C
ATOM    106  CH2 TRP A   6       2.842   0.407   4.120  1.00  0.00           C
ATOM    117  N   LEU A   7      -1.600  -1.860   0.967  1.00  0.00           N
ATOM    118  CA  LEU A   7      -1.641  -2.932   1.963  1.00  0.00           C
ATOM    119  C   LEU A   7      -1.847  -4.319   1.342  1.00  0.00           C
ATOM    120  O   LEU A   7      -1.144  -5.248   1.742  1.00  0.00           O
ATOM    121  CB  LEU A   7      -2.710  -2.645   3.033  1.00  0.00           C
ATOM    122  CG  LEU A   7      -2.301  -1.579   4.069  1.00  0.00           C
ATOM    123  CD1 LEU A   7      -3.475  -1.323   5.018  1.00  0.00           C
ATOM    124  CD2 LEU A   7      -1.093  -2.007   4.914  1.00  0.00           C
ATOM    136  N   LYS A   8      -2.753  -4.481   0.360  1.00  0.00           N
ATOM    137  CA  LYS A   8      -3.024  -5.791  -0.269  1.00  0.00           C
ATOM    138  C   LYS A   8      -1.796  -6.427  -0.937  1.00  0.00           C
ATOM    139  O   LYS A   8      -1.719  -7.648  -1.030  1.00  0.00           O
ATOM    140  CB  LYS A   8      -4.224  -5.697  -1.232  1.00  0.00           C
ATOM    141  CG  LYS A   8      -3.930  -5.009  -2.577  1.00  0.00           C
ATOM    142  CD  LYS A   8      -3.682  -5.986  -3.736  1.00  0.00           C
ATOM    143  CE  LYS A   8      -3.494  -5.199  -5.039  1.00  0.00           C
ATOM    144  NZ  LYS A   8      -4.563  -5.483  -6.023  1.00  0.00           N
ATOM    158  N   ASP A   9      -0.828  -5.607  -1.355  1.00  0.00           N
ATOM    159  CA  ASP A   9       0.466  -6.016  -1.905  1.00  0.00           C
ATOM    160  C   ASP A   9       1.481  -6.464  -0.832  1.00  0.00           C
ATOM    161  O   ASP A   9       2.545  -6.971  -1.194  1.00  0.00           O
ATOM    162  CB  ASP A   9       1.033  -4.839  -2.724  1.00  0.00           C
ATOM    163  CG  ASP A   9       0.672  -4.906  -4.210  1.00  0.00           C
ATOM    164  OD1 ASP A   9      -0.532  -5.051  -4.522  1.00  0.00           O
ATOM    165  OD2 ASP A   9       1.627  -4.815  -5.017  1.00  0.00           O
ATOM    170  N   GLY A  10       1.185  -6.278   0.464  1.00  0.00           N
ATOM    171  CA  GLY A  10       2.060  -6.618   1.593  1.00  0.00           C
ATOM    172  C   GLY A  10       2.628  -5.412   2.353  1.00  0.00           C
ATOM    173  O   GLY A  10       3.496  -5.594   3.208  1.00  0.00           O
ATOM    177  N   GLY A  11       2.172  -4.187   2.055  1.00  0.00           N
ATOM    178  CA  GLY A  11       2.626  -2.967   2.723  1.00  0.00           C
ATOM    179  C   GLY A  11       4.157  -2.802   2.654  1.00  0.00           C
ATOM    180  O   GLY A  11       4.710  -2.829   1.551  1.00  0.00           O
ATOM    184  N   PRO A  12       4.871  -2.651   3.794  1.00  0.00           N
ATOM    185  CA  PRO A  12       6.333  -2.533   3.806  1.00  0.00           C
ATOM    186  C   PRO A  12       7.058  -3.729   3.165  1.00  0.00           C
ATOM    187  O   PRO A  12       8.139  -3.562   2.601  1.00  0.00           O
ATOM    188  CB  PRO A  12       6.740  -2.387   5.279  1.00  0.00           C
ATOM    189  CG  PRO A  12       5.460  -1.952   5.987  1.00  0.00           C
ATOM    190  CD  PRO A  12       4.362  -2.615   5.160  1.00  0.00           C
ATOM    198  N   SER A  13       6.463  -4.929   3.205  1.00  0.00           N
ATOM    199  CA  SER A  13       7.049  -6.179   2.704  1.00  0.00           C
ATOM    200  C   SER A  13       6.897  -6.369   1.185  1.00  0.00           C
ATOM    201  O   SER A  13       7.025  -7.488   0.697  1.00  0.00           O
ATOM    202  CB  SER A  13       6.458  -7.371   3.472  1.00  0.00           C
ATOM    203  OG  SER A  13       6.763  -7.264   4.850  1.00  0.00           O
ATOM    209  N   SER A  14       6.637  -5.290   0.434  1.00  0.00           N
ATOM    210  CA  SER A  14       6.389  -5.315  -1.015  1.00  0.00           C
ATOM    211  C   SER A  14       7.332  -4.405  -1.823  1.00  0.00           C
ATOM    212  O   SER A  14       7.082  -4.123  -2.993  1.00  0.00           O
ATOM    213  CB  SER A  14       4.914  -4.993  -1.265  1.00  0.00           C
ATOM    214  OG  SER A  14       4.431  -5.743  -2.358  1.00  0.00           O
ATOM    220  N   GLY A  15       8.419  -3.920  -1.202  1.00  0.00           N
ATOM    221  CA  GLY A  15       9.451  -3.116  -1.870  1.00  0.00           C
ATOM    222  C   GLY A  15       8.984  -1.725  -2.316  1.00  0.00           C
ATOM    223  O   GLY A  15       9.539  -1.177  -3.267  1.00  0.00           O
ATOM    227  N   ARG A  16       7.956  -1.164  -1.660  1.00  0.00           N
ATOM    228  CA  ARG A  16       7.289   0.084  -2.054  1.00  0.00           C
ATOM    229  C   ARG A  16       6.855   0.916  -0.829  1.00  0.00           C
ATOM    230  O   ARG A  16       6.222   0.366   0.076  1.00  0.00           O
ATOM    231  CB  ARG A  16       6.110  -0.243  -2.994  1.00  0.00           C
ATOM    232  CG  ARG A  16       5.046  -1.171  -2.378  1.00  0.00           C
ATOM    233  CD  ARG A  16       3.923  -1.592  -3.338  1.00  0.00           C
ATOM    234  NE  ARG A  16       4.251  -2.811  -4.100  1.00  0.00           N
ATOM    235  CZ  ARG A  16       4.859  -2.914  -5.274  1.00  0.00           C
ATOM    236  NH1 ARG A  16       5.289  -1.864  -5.937  1.00  0.00           N
ATOM    237  NH2 ARG A  16       5.035  -4.095  -5.809  1.00  0.00           N
ATOM    251  N   PRO A  17       7.156   2.230  -0.780  1.00  0.00           N
ATOM    252  CA  PRO A  17       6.782   3.088   0.345  1.00  0.00           C
ATOM    253  C   PRO A  17       5.261   3.331   0.395  1.00  0.00           C
ATOM    254  O   PRO A  17       4.586   3.165  -0.624  1.00  0.00           O
ATOM    255  CB  PRO A  17       7.554   4.394   0.119  1.00  0.00           C
ATOM    256  CG  PRO A  17       7.677   4.474  -1.401  1.00  0.00           C
ATOM    257  CD  PRO A  17       7.820   3.010  -1.816  1.00  0.00           C
ATOM    265  N   PRO A  18       4.710   3.739   1.555  1.00  0.00           N
ATOM    266  CA  PRO A  18       3.287   4.031   1.686  1.00  0.00           C
ATOM    267  C   PRO A  18       2.901   5.305   0.913  1.00  0.00           C
ATOM    268  O   PRO A  18       3.684   6.256   0.871  1.00  0.00           O
ATOM    269  CB  PRO A  18       3.035   4.190   3.187  1.00  0.00           C
ATOM    270  CG  PRO A  18       4.385   4.655   3.729  1.00  0.00           C
ATOM    271  CD  PRO A  18       5.393   3.949   2.823  1.00  0.00           C
ATOM    279  N   PRO A  19       1.688   5.360   0.336  1.00  0.00           N
ATOM    280  CA  PRO A  19       1.185   6.543  -0.353  1.00  0.00           C
ATOM    281  C   PRO A  19       0.715   7.607   0.655  1.00  0.00           C
ATOM    282  O   PRO A  19      -0.124   7.324   1.513  1.00  0.00           O
ATOM    283  CB  PRO A  19       0.048   6.014  -1.229  1.00  0.00           C
ATOM    284  CG  PRO A  19      -0.519   4.852  -0.412  1.00  0.00           C
ATOM    285  CD  PRO A  19       0.716   4.275   0.272  1.00  0.00           C
ATOM    293  N   SER A  20       1.271   8.822   0.549  1.00  0.00           N
ATOM    294  CA  SER A  20       0.852  10.027   1.285  1.00  0.00           C
ATOM    295  C   SER A  20      -0.406  10.657   0.683  1.00  0.00           C
ATOM    296  O   SER A  20      -0.387  10.916  -0.540  1.00  0.00           O
ATOM    297  CB  SER A  20       1.972  11.071   1.284  1.00  0.00           C
ATOM    298  OG  SER A  20       3.120  10.541   1.911  1.00  0.00           O
ATOM    299  OXT SER A  20      -1.341  10.903   1.473  1.00  0.00           O
TER     305      SER A  20
END
"""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Executa una simulació de dinàmica molecular amb OpenMM sobre el Trp-cage 1L2Y."
    )
    parser.add_argument("--steps", type=int, default=5000, help="Passos de producció.")
    parser.add_argument(
        "--equilibration-steps",
        type=int,
        default=1000,
        help="Passos curts d'equilibratge abans de la producció.",
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
        help="Coeficient de fricció del termòstat en 1/ps.",
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
        default=1.2,
        help="Gruix de solvent al voltant de la proteïna en nm.",
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
        help="Directori de sortida.",
    )
    return parser.parse_args()


def build_simulation(
    temperature_kelvin: float,
    friction_per_ps: float,
    timestep_fs: float,
    padding_nm: float,
) -> tuple[Simulation, Modeller]:
    pdb = PDBFile(StringIO(TRP_CAGE_PDB))
    modeller = Modeller(pdb.topology, pdb.positions)

    forcefield = ForceField("amber14-all.xml", "amber14/tip3p.xml")
    modeller.addHydrogens(forcefield, pH=7.0)
    modeller.addSolvent(
        forcefield,
        model="tip3p",
        padding=padding_nm * unit.nanometer,
        neutralize=True,
        ionicStrength=0.0 * unit.molar,
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
        temperature_kelvin=args.temperature,
        friction_per_ps=args.friction,
        timestep_fs=args.timestep_fs,
        padding_nm=args.padding_nm,
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
            totalSteps=args.steps + args.equilibration_steps,
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
            totalSteps=args.steps + args.equilibration_steps,
        )
    )

    if args.equilibration_steps > 0:
        simulation.step(args.equilibration_steps)
    simulation.step(args.steps)

    final_state = simulation.context.getState(getPositions=True, getEnergy=True)
    with final_pdb_path.open("w", encoding="utf-8") as handle:
        PDBFile.writeFile(simulation.topology, final_state.getPositions(), handle)

    potential_energy = final_state.getPotentialEnergy().value_in_unit(
        unit.kilojoule_per_mole
    )
    print(f"\nSimulació completada. Energia potencial final: {potential_energy:.3f} kJ/mol")
    print("Proteïna d'entrada: Trp-cage miniprotein (PDB 1L2Y, model 1)")
    print(f"Resultats desats a: {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
