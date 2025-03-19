#!/usr/bin/env python3

#################################################
# Siesta Tool Box - Structural                  #
# version 1.0.0                                 #
# University of Brasilia (UnB) - Brazil         #
# version 1.0.0 2025/02/05                      #
# Dr. Carlos M. O. Bastos                       #
#################################################


import argparse
import numpy as np
from pymatgen.core import Structure
from pymatgen.analysis.local_env import (
     JmolNN, MinimumDistanceNN, CrystalNN,
    BrunnerNNRelative, EconNN)
import warnings
import logging

# Configure the logger
logging.basicConfig(filename="warnings.log", level=logging.WARNING, format="%(message)s")

def warn_handler(message, category, filename, lineno, file=None, line=None):
    log_message = f"{category.__name__}: {message} (File: {filename}, Line: {lineno})"
    logging.warning(log_message)  # Save to log file
    print("\n⚠️ Warning detected! Check warnings.log for details.")


# Override the default warning handler
warnings.showwarning = warn_handler

__version__ = "1.0.0"

def compute_ecn(structure, mode, atoms_position=None):
    with open("structural_information.dat", "w") as f:

	# Lattice
        lattice = structure.lattice

        print("\nLattice parameters:")
        print(f"   a = {lattice.a:.3f} Å")
        print(f"   b = {lattice.b:.3f} Å")
        print(f"   c = {lattice.c:.3f} Å")
        print(f"   Alpha = {lattice.alpha:.2f}°")
        print(f"   Beta = {lattice.beta:.2f}°")
        print(f"   Gamma = {lattice.gamma:.2f}°")

        f.write("\nLattice parameters:\n")
        f.write(f"   a = {lattice.a:.3f} Å\n")
        f.write(f"   b = {lattice.b:.3f} Å\n")
        f.write(f"   c = {lattice.c:.3f} Å\n")
        f.write(f"   Alpha = {lattice.alpha:.2f}°\n")
        f.write(f"   Beta = {lattice.beta:.2f}°\n")
        f.write(f"   Gamma = {lattice.gamma:.2f}°\n")


        lattice_vectors = structure.lattice.matrix
        print("\nLattice vectors:")
        f.write("\nLattice vectors:\n")
        for i, vector in enumerate(lattice_vectors):
            print(f"   {'a_' + str(i+1)}: {vector[0]}   {vector[1]}   {vector[2]}")
            f.write(f"   {'a_' + str(i+1)}: {vector[0]}   {vector[1]}   {vector[2]}\n")




        # ECN Methods
        methods = {
            "JmolNN": JmolNN(),
            "MinDistNN": MinimumDistanceNN(),
            "CrystalNN": CrystalNN(),
            "BrunnerNN": BrunnerNNRelative(),
            "EconNN": EconNN()
        }
        ecn_results = {method: [] for method in methods}

        # Atomic positions
        pos_atomics = [[i+1, str(site.specie.symbol), site.coords] for i, site in enumerate(structure)]

        if mode == "mean":
            for i in range(len(structure)):
                for method_name, method in methods.items():
                    try:
                        ecn_results[method_name].append(method.get_cn(structure, i))
                    except:
                        ecn_results[method_name].append(None)

            ecn_avg = {method: np.nanmean([v for v in values if v is not None]) for method, values in ecn_results.items()}
            print("\nECN Average:")
            f.write("\nECN Average:\n")
            for method, value in ecn_avg.items():
                print(f"{method:15}: {value:.2f}")
                f.write(f"{method:15}: {value:.2f}\n")

        elif mode == "list" and atoms_position:
            for i in atoms_position:
                for method_name, method in methods.items():
                    try:
                        ecn_results[method_name].append(method.get_cn(structure, i-1))
                    except:
                        ecn_results[method_name].append(None)

            print("\nECN:")
            f.write("\nECN:\n")
            for i, atom_index in enumerate(atoms_position):
                print(f" Atom {pos_atomics[atom_index-1][0]}:")
                print(f"   Element: {pos_atomics[atom_index-1][1]}     Cartesian Position: {pos_atomics[atom_index-1][2]}")
                f.write(f" Atom {pos_atomics[atom_index-1][0]}:\n")
                f.write(f"   Element: {pos_atomics[atom_index-1][1]}     Cartesian Position: {pos_atomics[atom_index-1][2]}\n")
                for method, values in ecn_results.items():
                    print(f"      {method:15}: {values[i]}")
                    f.write(f"      {method:15}: {values[i]}\n")



        # Atomic positions
        pos_atomics=[]
        print("\nPositions of atoms:")
        f.write("\nPositions of atoms:\n")
        for i,site in enumerate(structure):
            print(f"index: {i+1} Element: {site.specie.symbol} cartesian position: {site.coords}")
            f.write(f"index: {i+1} Element: {site.specie.symbol} cartesian position: {site.coords}\n")
            pos_atomics.append([i+1,str(site.specie.symbol),site.coords])


def main():
    parser = argparse.ArgumentParser(description=f"Compute ECN from structure file. Version: {__version__}. WARNING: In this version the input file is cif or poscar format.")
    parser.add_argument("--file", required=True, help="Path to structure file. Only POSCAR and CIF are supported, but stb-translate can be used to convert file in poscar.")
    parser.add_argument("--mode", choices=["list", "mean"], required=True, help="Calculation mode: list or mean")
    parser.add_argument("--list", type=str, help="List of atom indices (comma-separated).Ex: [1,4,5,7] - Require for 'list' mode")

    args = parser.parse_args()

    if args.mode == "list" and not args.list:
        parser.error("--list is required when --mode is 'list'")

    atoms_position = list(map(int, args.list.strip('[]').split(','))) if args.list else None

    structure = Structure.from_file(args.file)
    compute_ecn(structure, args.mode, atoms_position)


    print("\n\n ###### ECN calculated! My atoms are more socially active than me. ######\n")


if __name__ == "__main__":
    main()
