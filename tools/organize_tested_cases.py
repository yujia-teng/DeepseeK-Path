#!/usr/bin/env python3
"""
Organize PNG-backed structure test cases by HPKOT lattice class.

The script scans the seven lattice-family folders in the sibling ``structure``
directory, keeps only cases that have the expected PNG products, classifies the
matching structure with spglib/seekpath, writes a CSV manifest, and optionally
copies the files into:

    output_root / family / hpkot_key / SGnumber-SGsymbol

Spin-Laue labels are intentionally not inferred. The manifest contains an empty
``manual_spin_laue`` column for later annotation.
"""

from __future__ import annotations

import argparse
import csv
import shutil
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import seekpath
import spglib
from pymatgen.core import Structure


FAMILIES = (
    "triclinic",
    "monoclinic",
    "orthorhombic",
    "tetragonal",
    "trigonal",
    "hexagnoal",
    "cubic",
)

STRUCTURE_SUFFIXES = (".vasp", ".cif", ".mcif")
PNG_MARKERS = ("_spinbz_top_", "_spinflip_", "_spinbz_", "_ibz_")
CLASSIFY_SUFFIX_PRIORITY = {".vasp": 0, ".cif": 1, ".mcif": 2}

warnings.filterwarnings("ignore", message="We strongly encourage explicit.*encoding")
warnings.filterwarnings("ignore", message="dict interface is deprecated")
warnings.filterwarnings("ignore", message="Issues encountered while parsing CIF:.*")
warnings.filterwarnings("ignore", message="Cannot determine chemical composition from CIF!.*")


@dataclass(frozen=True)
class ClassifiedCase:
    family: str
    case_stem: str
    hpkot_key: str
    sg_number: int
    sg_symbol: str
    laue_group: str
    target_folder: Path
    files_to_copy: tuple[Path, ...]


@dataclass(frozen=True)
class Problem:
    family: str
    case_stem: str
    reason: str


def default_structure_root() -> Path:
    repo_root = Path(__file__).resolve().parents[1]
    return repo_root.parent / "structure"


def default_output_root(source_root: Path) -> Path:
    return source_root / "tested_by_lattice_class"


def png_case_stem(path: Path) -> str | None:
    name = path.name
    for marker in PNG_MARKERS:
        if marker in name:
            return name.split(marker, 1)[0]
    return None


def safe_sg_symbol(symbol: str) -> str:
    return (
        symbol.strip()
        .replace("/", "_")
        .replace("\\", "_")
        .replace(" ", "")
        .replace(":", "_")
    )


def laue_group_from_point_group(point_group: str) -> str:
    pg = str(point_group).strip().replace(" ", "")
    mapping = {
        "1": "-1",
        "-1": "-1",
        "2": "2/m",
        "m": "2/m",
        "2/m": "2/m",
        "222": "mmm",
        "mm2": "mmm",
        "mmm": "mmm",
        "4": "4/m",
        "-4": "4/m",
        "4/m": "4/m",
        "422": "4/mmm",
        "4mm": "4/mmm",
        "-42m": "4/mmm",
        "-4m2": "4/mmm",
        "4/mmm": "4/mmm",
        "3": "-3",
        "-3": "-3",
        "32": "-3m",
        "3m": "-3m",
        "-3m": "-3m",
        "6": "6/m",
        "-6": "6/m",
        "6/m": "6/m",
        "622": "6/mmm",
        "6mm": "6/mmm",
        "-6m2": "6/mmm",
        "-62m": "6/mmm",
        "6/mmm": "6/mmm",
        "23": "m-3",
        "m-3": "m-3",
        "432": "m-3m",
        "-43m": "m-3m",
        "m-3m": "m-3m",
    }
    return mapping.get(pg, "")


def structure_to_seekpath_cell(structure: Structure):
    lattice = structure.lattice.matrix
    positions = structure.frac_coords
    numbers = [site.specie.Z for site in structure]
    return lattice, positions, numbers


def classify_structure(path: Path, symprec: float, angle_tolerance: float):
    structure = Structure.from_file(str(path))
    cell = structure_to_seekpath_cell(structure)
    dataset = spglib.get_symmetry_dataset(
        cell, symprec=symprec, angle_tolerance=angle_tolerance
    )
    if dataset is None:
        raise ValueError("spglib could not determine symmetry")

    sp_result = seekpath.get_path(
        cell,
        with_time_reversal=True,
        recipe="hpkot",
        threshold=symprec,
        symprec=symprec,
        angle_tolerance=angle_tolerance,
    )
    return {
        "hpkot_key": sp_result.get(
            "bravais_lattice_extended", sp_result.get("bravais_lattice", "")
        ),
        "sg_number": int(dataset.number),
        "sg_symbol": str(dataset.international),
        "point_group": str(dataset.pointgroup),
    }


def find_structure_files(family_dir: Path, stem: str) -> list[Path]:
    matches: list[Path] = []
    for path in family_dir.iterdir():
        if not path.is_file():
            continue
        if path.suffix.lower() in STRUCTURE_SUFFIXES and path.stem == stem:
            matches.append(path)
        elif stem.upper().startswith("POSCAR") and path.name == stem:
            matches.append(path)
        elif path.name.upper().startswith("POSCAR") and path.stem == stem:
            matches.append(path)
    return sorted(matches, key=lambda p: p.name.lower())


def choose_classification_file(paths: list[Path]) -> Path:
    return sorted(
        paths,
        key=lambda p: (CLASSIFY_SUFFIX_PRIORITY.get(p.suffix.lower(), 99), p.name.lower()),
    )[0]


def matching_pngs(family_dir: Path, stem: str) -> list[Path]:
    matches = []
    for path in family_dir.glob("*.png"):
        if png_case_stem(path) == stem:
            matches.append(path)
    return sorted(matches, key=lambda p: p.name.lower())


def discover_cases(
    source_root: Path,
    output_root: Path,
    symprec: float,
    angle_tolerance: float,
) -> tuple[list[ClassifiedCase], list[Problem]]:
    cases: list[ClassifiedCase] = []
    problems: list[Problem] = []

    for family in FAMILIES:
        family_dir = source_root / family
        if not family_dir.is_dir():
            problems.append(Problem(family, "", f"missing folder: {family_dir}"))
            continue

        stems = sorted(
            {
                stem
                for png in family_dir.glob("*.png")
                for stem in [png_case_stem(png)]
                if stem
            }
        )

        for stem in stems:
            pngs = matching_pngs(family_dir, stem)
            structure_files = find_structure_files(family_dir, stem)
            if not structure_files:
                problems.append(Problem(family, stem, "PNG-backed case has no matching structure file"))
                continue

            structure_file = choose_classification_file(structure_files)
            try:
                classified = classify_structure(structure_file, symprec, angle_tolerance)
            except Exception as exc:
                problems.append(Problem(family, stem, f"classification failed: {exc}"))
                continue

            sg_folder = f"{classified['sg_number']}-{safe_sg_symbol(classified['sg_symbol'])}"
            target_folder = output_root / family / classified["hpkot_key"] / sg_folder
            files_to_copy = tuple([*structure_files, *pngs])
            cases.append(
                ClassifiedCase(
                    family=family,
                    case_stem=stem,
                    hpkot_key=classified["hpkot_key"],
                    sg_number=classified["sg_number"],
                    sg_symbol=classified["sg_symbol"],
                    laue_group=laue_group_from_point_group(classified["point_group"]),
                    target_folder=target_folder,
                    files_to_copy=files_to_copy,
                )
            )

    return cases, problems


def manifest_rows(cases: Iterable[ClassifiedCase]):
    for case in cases:
        yield {
            "family": case.family,
            "case_stem": case.case_stem,
            "hpkot_key": case.hpkot_key,
            "sg_number": case.sg_number,
            "sg_symbol": case.sg_symbol,
            "laue_group": case.laue_group,
            "manual_spin_laue": "",
            "target_folder": str(case.target_folder),
            "files_to_copy": ";".join(str(path) for path in case.files_to_copy),
        }


def write_manifest(path: Path, cases: list[ClassifiedCase]) -> None:
    fieldnames = [
        "family",
        "case_stem",
        "hpkot_key",
        "sg_number",
        "sg_symbol",
        "laue_group",
        "manual_spin_laue",
        "target_folder",
        "files_to_copy",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(manifest_rows(cases))


def write_problem_report(path: Path, problems: list[Problem]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=["family", "case_stem", "reason"])
        writer.writeheader()
        for problem in problems:
            writer.writerow(
                {
                    "family": problem.family,
                    "case_stem": problem.case_stem,
                    "reason": problem.reason,
                }
            )


def copy_cases(cases: list[ClassifiedCase], overwrite: bool) -> tuple[int, list[Problem]]:
    copied = 0
    problems: list[Problem] = []
    for case in cases:
        case.target_folder.mkdir(parents=True, exist_ok=True)
        for source in case.files_to_copy:
            destination = case.target_folder / source.name
            if destination.exists() and not overwrite:
                problems.append(
                    Problem(case.family, case.case_stem, f"destination exists, skipped: {destination}")
                )
                continue
            shutil.copy2(source, destination)
            copied += 1
    return copied, problems


def print_summary(
    cases: list[ClassifiedCase],
    problems: list[Problem],
    manifest: Path,
    problem_report: Path,
) -> None:
    print(f"Manifest: {manifest}")
    print(f"Problem report: {problem_report}")
    print(f"PNG-backed classified cases: {len(cases)}")
    print(f"Problems: {len(problems)}")

    by_class: dict[tuple[str, str], int] = {}
    for case in cases:
        by_class[(case.family, case.hpkot_key)] = by_class.get((case.family, case.hpkot_key), 0) + 1
    for (family, hpkot_key), count in sorted(by_class.items()):
        print(f"  {family}/{hpkot_key}: {count}")

    if problems:
        print("\nProblems:")
        for problem in problems:
            label = f"{problem.family}/{problem.case_stem}" if problem.case_stem else problem.family
            print(f"  {label}: {problem.reason}")


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, default=default_structure_root())
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--problem-report", type=Path)
    parser.add_argument("--copy", action="store_true", help="copy files after writing the manifest")
    parser.add_argument("--overwrite", action="store_true", help="allow replacing files in the output tree")
    parser.add_argument("--strict", action="store_true", help="return a nonzero exit code when unmatched cases are found")
    parser.add_argument("--symprec", type=float, default=1e-5)
    parser.add_argument("--angle-tolerance", type=float, default=-1.0)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    source_root = args.source_root.resolve()
    output_root = (args.output_root or default_output_root(source_root)).resolve()
    manifest = (args.manifest or output_root / "tested_cases_manifest.csv").resolve()
    problem_report = (
        args.problem_report or output_root / "tested_cases_unmatched.csv"
    ).resolve()

    cases, problems = discover_cases(
        source_root=source_root,
        output_root=output_root,
        symprec=args.symprec,
        angle_tolerance=args.angle_tolerance,
    )
    write_manifest(manifest, cases)
    write_problem_report(problem_report, problems)
    print_summary(cases, problems, manifest, problem_report)

    if args.copy:
        copied, copy_problems = copy_cases(cases, overwrite=args.overwrite)
        print(f"\nCopied files: {copied}")
        if copy_problems:
            print(f"Copy skips/problems: {len(copy_problems)}")
            for problem in copy_problems:
                print(f"  {problem.family}/{problem.case_stem}: {problem.reason}")
    else:
        print("\nDry-run only. Re-run with --copy to copy files.")

    return 1 if args.strict and problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
