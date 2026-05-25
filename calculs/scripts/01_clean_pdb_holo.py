from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

PROTEIN_RECORDS = {"ATOM"}
KEEP_RESIDUES_DEFAULT = ("HPN",)
DROP_RESIDUES_DEFAULT = ("HOH", "WAT", "MES")
SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_PREPARED_PDB_DIR = SCRIPT_DIR.parent / "prepared pdbs"

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Neteja el PDB per conservar proteina i un o mes holo, eliminant aigua i buffer."
    )
    parser.add_argument("--study-dir", type=Path, default=DEFAULT_PREPARED_PDB_DIR)
    parser.add_argument("--input-pdb", type=Path, default=Path("input/leh_1nww.pdb"))
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--keep-residue", action="append", default=None)
    parser.add_argument("--drop-residue", action="append", default=None)
    return parser.parse_args()

def residue_name(line: str) -> str:
    return line[17:20].strip()


def atom_serial(line: str) -> int | None:
    try:
        return int(line[6:11])
    except ValueError:
        return None

def conect_serials(line: str) -> list[int]:
    serials = []
    for start in range(6, len(line), 5):
        token = line[start:start + 5].strip()
        if token:
            try:
                serials.append(int(token))
            except ValueError:
                pass
    return serials


def main() -> None:
    args = parse_args()
    if args.output_dir is None:
        args.output_dir = args.study_dir / "01_cleaned_pdb_holo"
    args.output_dir.mkdir(parents=True, exist_ok=True)

    keep_residues = set(args.keep_residue or KEEP_RESIDUES_DEFAULT)
    drop_residues = set(args.drop_residue or DROP_RESIDUES_DEFAULT)

    kept_lines: list[str] = []
    kept_serials: set[int] = set()
    deferred_conect: list[str] = []
    kept_hetero = Counter()
    removed_hetero = Counter()

    for line in args.input_pdb.read_text(encoding="utf-8").splitlines(keepends=True):
        record = line[:6].strip()

        if record in PROTEIN_RECORDS:
            kept_lines.append(line)
            serial = atom_serial(line)
            if serial is not None:
                kept_serials.add(serial)
            continue

        if record == "HETATM":
            name = residue_name(line)
            if name in keep_residues:
                kept_lines.append(line)
                kept_hetero[name] += 1
                serial = atom_serial(line)
                if serial is not None:
                    kept_serials.add(serial)
            elif name in drop_residues:
                removed_hetero[name] += 1
            else:
                removed_hetero[name] += 1
            continue

        if record == "ANISOU":
            serial = atom_serial(line)
            if serial is not None and serial in kept_serials:
                kept_lines.append(line)
            continue

        if record == "CONECT":
            deferred_conect.append(line)
            continue

        if record in {"TER", "END"}:
            kept_lines.append(line)
            continue

        if record in {"HET", "HETNAM", "FORMUL"}:
            if any(name in line for name in keep_residues):
                kept_lines.append(line)
            continue

        if record.startswith("REMARK") or record in {
            "HEADER", "TITLE", "COMPND", "SOURCE", "KEYWDS", "EXPDTA",
            "AUTHOR", "JRNL", "DBREF", "SEQRES", "HELIX", "SHEET", "SITE",
            "CRYST1", "ORIGX1", "ORIGX2", "ORIGX3", "SCALE1", "SCALE2", "SCALE3"
        }:
            kept_lines.append(line)
            continue

        kept_lines.append(line)

    conect_lines = []
    for line in deferred_conect:
        serials = conect_serials(line)
        if serials and all(s in kept_serials for s in serials):
            conect_lines.append(line)

    output_pdb = args.output_dir / "cleaned_holo.pdb"
    if kept_lines and kept_lines[-1].startswith("END"):
        final_lines = kept_lines[:-1] + conect_lines + [kept_lines[-1]]
    else:
        final_lines = kept_lines + conect_lines + ["END\n"]
    output_pdb.write_text("".join(final_lines), encoding="utf-8")

    summary = [
        f"Input PDB: {args.input_pdb}",
        f"Output PDB: {output_pdb}",
        "Residus HETATM conservats: " + (", ".join(f"{k}={v}" for k, v in sorted(kept_hetero.items())) or "cap"),
        "Residus HETATM eliminats: " + (", ".join(f"{k}={v}" for k, v in sorted(removed_hetero.items()) if v) or "cap"),
        "Residus conservats per defecte: " + ", ".join(sorted(keep_residues)),
        "Residus eliminats per defecte: " + ", ".join(sorted(drop_residues)),
    ]
    (args.output_dir / "summary.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")
    print(f"PDB netejat i desat a {output_pdb}")


if __name__ == "__main__":
    main()
