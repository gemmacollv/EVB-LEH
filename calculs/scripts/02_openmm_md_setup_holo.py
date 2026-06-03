from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

CONDA_EXE = Path("/home/10034103/miniforge3/bin/conda")
REEXEC_ENV_VAR = "EVB_LEH_MD_OPENMM_REEXEC"


def reexec_with_md_openmm() -> None:
    if os.environ.get(REEXEC_ENV_VAR) == "1" or not CONDA_EXE.exists():
        return
    os.environ[REEXEC_ENV_VAR] = "1"
    command = [
        str(CONDA_EXE),
        "run",
        "-n",
        "md-openmm",
        "python",
        str(Path(__file__).resolve()),
        *sys.argv[1:],
    ]
    print("Falten dependencies en aquest Python; reexecutant amb md-openmm...", flush=True)
    os.execv(str(CONDA_EXE), command)


if os.environ.get(REEXEC_ENV_VAR) != "1":
    try:
        import parmed
        import openff
        import openmm
        import openmmforcefields
        import rdkit
    except ImportError:
        reexec_with_md_openmm()

import parmed as pmd
from openmm import unit
from openmm.app import ForceField, Modeller, PME, PDBFile

from md_common import write_pdb

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_PREPARED_PDB_DIR = SCRIPT_DIR.parent / "prepared pdbs"
DEFAULT_STUDY_DIR = DEFAULT_PREPARED_PDB_DIR
DEFAULT_INPUT_PDB = DEFAULT_PREPARED_PDB_DIR / "01_cleaned_pdb_holo" / "cleaned_holo.pdb"
DEFAULT_OUTPUT_DIR = DEFAULT_PREPARED_PDB_DIR / "02_openmm_md_setup" / "protein_holo"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepara la proteïna amb HPN, solvata, ionitza i exporta PRMTOP/INPCRD."
    )
    parser.add_argument("--study-dir", type=Path, default=DEFAULT_STUDY_DIR)
    parser.add_argument("--input-pdb", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--prefix", default="protein_holo")
    parser.add_argument("--holo-resname", default="HPN")
    parser.add_argument("--holo-smiles", default="CCCCCCC(=O)N")
    parser.add_argument("--holo-forcefield", default="openff-2.2.1")
    parser.add_argument("--ph", type=float, default=7.0)
    parser.add_argument("--padding-nm", type=float, default=1.2)
    parser.add_argument("--ionic-strength", type=float, default=0.15, help="M")
    parser.add_argument("--temperature", type=float, default=300.0)
    return parser.parse_args()


def import_holo_tools():
    try:
        from openff.toolkit import Molecule
        from openmmforcefields.generators import SMIRNOFFTemplateGenerator
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError as exc:
        reexec_with_md_openmm()
        raise SystemExit(
            "Falten dependencies per parametritzar HPN dins aquest Python.\n"
            f"Python actual: {sys.executable}\n"
            "Executa dins md-openmm: conda run -n md-openmm python "
            "calculs/scripts/02_openmm_md_setup_holo.py"
        ) from exc
    return Molecule, SMIRNOFFTemplateGenerator, Chem, AllChem


def atom_serial(line: str) -> int:
    return int(line[6:11])


def conect_serials(line: str) -> list[int]:
    serials = []
    for start in range(6, len(line), 5):
        token = line[start:start + 5].strip()
        if token:
            serials.append(int(token))
    return serials


def holo_pdb_blocks(pdb_path: Path, holo_resname: str) -> list[str]:
    text = pdb_path.read_text(encoding="utf-8").splitlines(keepends=True)
    residues: dict[tuple[str, str, str], list[str]] = {}

    for line in text:
        if not line.startswith("HETATM") or line[17:20].strip() != holo_resname:
            continue
        key = (line[21].strip(), line[22:26].strip(), line[26].strip())
        residues.setdefault(key, []).append(line)

    blocks = []
    for atom_lines in residues.values():
        serials = {atom_serial(line) for line in atom_lines}
        conect_lines = []
        for line in text:
            if not line.startswith("CONECT"):
                continue
            ids = conect_serials(line)
            if ids and all(serial in serials for serial in ids):
                conect_lines.append(line)
        blocks.append("".join(atom_lines + conect_lines) + "END\n")
    return blocks


def create_holo_molecules(pdb_path: Path, holo_resname: str, holo_smiles: str):
    Molecule, _, Chem, AllChem = import_holo_tools()
    template = Chem.MolFromSmiles(holo_smiles)
    if template is None:
        raise SystemExit(f"SMILES invalid per al holo {holo_resname}: {holo_smiles}")

    molecules = []
    for index, block in enumerate(holo_pdb_blocks(pdb_path, holo_resname), start=1):
        pdb_mol = Chem.MolFromPDBBlock(block, sanitize=False, removeHs=False)
        if pdb_mol is None:
            raise SystemExit(f"No s ha pogut llegir el bloc PDB del holo {holo_resname} #{index}.")
        try:
            rdkit_mol = AllChem.AssignBondOrdersFromTemplate(template, pdb_mol)
            rdkit_mol = Chem.AddHs(rdkit_mol, addCoords=True)
            Chem.SanitizeMol(rdkit_mol)
        except Exception as exc:
            raise SystemExit(
                f"No s ha pogut assignar la quimica del SMILES al holo {holo_resname} #{index}: {exc}"
            ) from exc

        molecule = Molecule.from_rdkit(
            rdkit_mol,
            allow_undefined_stereo=True,
            hydrogens_are_explicit=True,
        )
        molecule.name = holo_resname
        molecules.append(molecule)

    if not molecules:
        raise SystemExit(f"No s ha trobat cap residu {holo_resname} a {pdb_path}.")
    return molecules


def create_holo_forcefield(holo_molecules, holo_forcefield: str):
    _, SMIRNOFFTemplateGenerator, _, _ = import_holo_tools()
    forcefield = ForceField("amber14-all.xml", "amber14/tip3p.xml")
    smirnoff = SMIRNOFFTemplateGenerator(
        molecules=holo_molecules,
        forcefield=holo_forcefield,
    )
    forcefield.registerTemplateGenerator(smirnoff.generator)
    return forcefield


def add_hydrogenated_holos(modeller: Modeller, holo_molecules) -> None:
    for molecule in holo_molecules:
        modeller.add(
            molecule.to_topology().to_openmm(),
            molecule.conformers[0].to_openmm(),
        )


def normalize_residue_ids(topology) -> None:
    for residue in topology.residues():
        residue.id = str(residue.id)


def main() -> None:
    args = parse_args()
    if args.input_pdb is None:
        args.input_pdb = DEFAULT_INPUT_PDB
    if args.output_dir is None:
        args.output_dir = args.study_dir / "02_openmm_md_setup" / "protein_holo"

    args.input_pdb = args.input_pdb.resolve()
    args.output_dir = args.output_dir.resolve()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    if not args.input_pdb.exists():
        raise SystemExit(f"No existeix el PDB d entrada: {args.input_pdb}")

    pdb = PDBFile(str(args.input_pdb))
    holo_molecules = create_holo_molecules(
        args.input_pdb,
        args.holo_resname,
        args.holo_smiles,
    )
    forcefield = create_holo_forcefield(holo_molecules, args.holo_forcefield)

    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.delete([
        residue
        for residue in modeller.topology.residues()
        if residue.name == args.holo_resname
    ])
    modeller.addHydrogens(forcefield, pH=args.ph)
    add_hydrogenated_holos(modeller, holo_molecules)
    modeller.addSolvent(
        forcefield,
        model="tip3p",
        padding=args.padding_nm * unit.nanometer,
        neutralize=True,
        ionicStrength=args.ionic_strength * unit.molar,
    )
    normalize_residue_ids(modeller.topology)

    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=None,
        rigidWater=False,
    )
    structure = pmd.openmm.load_topology(
        modeller.topology,
        system,
        xyz=modeller.positions,
    )

    prmtop_path = args.output_dir / f"{args.prefix}.prmtop"
    inpcrd_path = args.output_dir / f"{args.prefix}.inpcrd"
    solvated_pdb_path = args.output_dir / f"{args.prefix}_solvated.pdb"
    summary_path = args.output_dir / "summary_holo.txt"

    structure.save(str(prmtop_path), overwrite=True)
    structure.save(str(inpcrd_path), overwrite=True)
    write_pdb(modeller.topology, modeller.positions, solvated_pdb_path)

    box_vectors = modeller.topology.getPeriodicBoxVectors()
    box_vectors_nm = [vector.value_in_unit(unit.nanometer) for vector in box_vectors]
    summary = [
        "System: holo",
        f"Input PDB: {args.input_pdb}",
        f"Holo residue: {args.holo_resname}",
        f"Holo molecules: {len(holo_molecules)}",
        f"Holo SMILES: {args.holo_smiles}",
        f"Holo force field: {args.holo_forcefield}",
        f"PRMTOP: {prmtop_path}",
        f"INPCRD: {inpcrd_path}",
        f"Solvated PDB: {solvated_pdb_path}",
        "Protein force field: amber14-all.xml",
        "Water force field: amber14/tip3p.xml",
        "Nonbonded method: PME",
        "Constraints: None during AMBER export; simulation can apply HBonds.",
        f"pH: {args.ph}",
        f"Padding (nm): {args.padding_nm}",
        f"Ionic strength (M): {args.ionic_strength}",
        f"Temperature (K): {args.temperature}",
        f"Atoms totals: {modeller.topology.getNumAtoms()}",
        f"Box vectors (nm): {box_vectors_nm}",
    ]
    summary_path.write_text("\n".join(summary) + "\n", encoding="utf-8")

    print("Sistema protein_holo preparat")
    print(f"PDB entrada: {args.input_pdb}")
    print(f"Fitxer AMBER topology: {prmtop_path}")
    print(f"Fitxer AMBER coordinates: {inpcrd_path}")


if __name__ == "__main__":
    main()
