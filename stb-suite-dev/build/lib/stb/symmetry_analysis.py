#!/usr/bin/env python3
import argparse
import numpy as np
import os
from ase.io import read
from ase.geometry import cell_to_cellpar
import spglib

def read_structure(filename, filetype):
    filetype = filetype.lower()
    if filetype == "cif":
        atoms = read(filename, format="cif")
    elif filetype == "poscar":
        atoms = read(filename, format="vasp")
    elif filetype == "siesta":
        atoms = read(filename, format="struct_out")
    else:
        raise ValueError(f"Unsupported file type: {filetype}")
    
    cell = atoms.get_cell()
    positions = atoms.get_scaled_positions()
    numbers = atoms.get_atomic_numbers()
    symbols = atoms.get_chemical_symbols()
    return atoms, cell, positions, numbers, symbols

def get_crystal_system(spacegroup_number):
    if spacegroup_number <= 2:
        return "Triclinic"
    elif spacegroup_number <= 15:
        return "Monoclinic"
    elif spacegroup_number <= 74:
        return "Orthorhombic"
    elif spacegroup_number <= 142:
        return "Tetragonal"
    elif spacegroup_number <= 167:
        return "Trigonal"
    elif spacegroup_number <= 194:
        return "Hexagonal"
    else:
        return "Cubic"

def analyze_symmetry(cell, positions, numbers, symprec=1e-3):
    structure = (cell, positions, numbers)
    dataset = spglib.get_symmetry_dataset(structure, symprec=symprec)
    symmetry = spglib.get_symmetry(structure, symprec=symprec)

    cellpar = cell_to_cellpar(cell)
    volume = np.linalg.det(cell)
    n_distinct = len(set(dataset.equivalent_atoms))

    return {
        "space_group_name": dataset.international,
        "space_group_number": dataset.number,
        "hall_symbol": dataset.hall,
        "point_group": dataset.pointgroup,
        "crystal_system": get_crystal_system(dataset.number),
        "cell_parameters": {
            "a": cellpar[0],
            "b": cellpar[1],
            "c": cellpar[2],
            "alpha": cellpar[3],
            "beta": cellpar[4],
            "gamma": cellpar[5]
        },
        "volume": volume,
        "wyckoff_positions": dataset.wyckoffs,
        "equivalent_atoms": dataset.equivalent_atoms,
        "n_distinct_sites": n_distinct,
        "n_symmetry_operations": len(symmetry["rotations"]),
        "symmetry_operations": [
            {
                "rotation": rot.tolist(),
                "translation": trans.tolist()
            }
            for rot, trans in zip(symmetry["rotations"], symmetry["translations"])
        ]
    }

def generate_txt_report(numbers, symbols, sym_data):
    lines = []
    lines.append("=== Crystal Symmetry Report ===\n")
    lines.append(f"Space group      : {sym_data['space_group_name']} (No. {sym_data['space_group_number']})")
    lines.append(f"Hall symbol      : {sym_data['hall_symbol']}")
    lines.append(f"Point group      : {sym_data['point_group']}")
    lines.append(f"Crystal system   : {sym_data['crystal_system']}\n")

    lines.append("== Cell Parameters ==")
    cp = sym_data["cell_parameters"]
    lines.append(f"a = {cp['a']:.6f} Å, b = {cp['b']:.6f} Å, c = {cp['c']:.6f} Å")
    lines.append(f"α = {cp['alpha']:.2f}°, β = {cp['beta']:.2f}°, γ = {cp['gamma']:.2f}°")
    lines.append(f"Cell volume      : {sym_data['volume']:.4f} Å³\n")

    lines.append("== Atomic Sites ==")
    for i, (Z, sym, wyck, eq) in enumerate(zip(numbers, symbols, sym_data["wyckoff_positions"], sym_data["equivalent_atoms"])):
        lines.append(f"Atom {i+1:>3}: Z={Z:<3} Symbol={sym:<2} Wyckoff={wyck}  Equivalent to atom {eq+1}")
    lines.append(f"\nNumber of symmetrically distinct sites: {sym_data['n_distinct_sites']}\n")

    lines.append("== Symmetry Operations ==")
    lines.append(f"Total operations: {sym_data['n_symmetry_operations']}\n")
    for i, op in enumerate(sym_data["symmetry_operations"]):
        lines.append(f"Operation {i+1}:")
        lines.append("Rotation matrix:")
        for row in op["rotation"]:
            lines.append("  " + "  ".join(f"{x:2d}" for x in row))
        lines.append("Translation vector: " + "  ".join(f"{x:.5f}" for x in op["translation"]))
        lines.append("")
    return "\n".join(lines)

def main():
    parser = argparse.ArgumentParser(
        description="Analyze crystal symmetry from a structure file using spglib.",
        epilog="Example:\n  python symmetry_analysis.py -i STRUCT_OUT -t siesta",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("-i", "--input", required=True, help="Input structure file (POSCAR, CIF, or SIESTA STRUCT_OUT)")
    parser.add_argument("-t", "--filetype", required=True, choices=["poscar", "cif", "siesta"],
                        help="File type: poscar, cif, siesta")
    parser.add_argument("--symprec", type=float, default=1e-3, help="Symmetry precision (default: 1e-3)")

    args = parser.parse_args()

    try:
        atoms, cell, positions, numbers, symbols = read_structure(args.input, args.filetype)
        sym_data = analyze_symmetry(cell, positions, numbers, symprec=args.symprec)

        txt_content = generate_txt_report(numbers, symbols, sym_data)
        with open("symmetry.dat", "w") as f:
            f.write(txt_content)
        print(txt_content)
        print("Output saved to symmetry.dat")

    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()

