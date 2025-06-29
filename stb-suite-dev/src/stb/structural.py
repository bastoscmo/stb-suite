#!/usr/bin/env python

#################################################
# Siesta Tool Box - Structural                  #
# version 1.0.0                                 #
# University of Brasilia (UnB) - Brazil         #
# version 1.0.0 2025/02/05                      #
# Dr. Carlos M. O. Bastos                       #
#################################################

VERSION = "1.5.1"

import os
import sys
import warnings
import subprocess
from time import sleep
import argparse
import textwrap
from typing import List, Dict
import argparse
import numpy as np
from pymatgen.core import Structure
from pymatgen.analysis.local_env import (
     JmolNN, MinimumDistanceNN, CrystalNN,
    BrunnerNNRelative, EconNN)
import warnings
import logging
from ase.io import read as ase_read
from pymatgen.io.ase import AseAtomsAdaptor

# ANSI colors for terminal output
COLORS = {
    'reset': '\033[0m',
    'cyan': '\033[96m',
    'blue': '\033[94m',
    'green': '\033[92m',
    'yellow': '\033[93m',
    'red': '\033[91m',
    'bold': '\033[1m',
    'underline': '\033[4m'
}

def color_text(text: str, color: str) -> str:
    """Return ANSI colored text."""
    return f"{COLORS[color]}{text}{COLORS['reset']}"

def show_intro() -> None:
    """Display toolbox introduction with logo."""
    os.system('cls' if os.name == 'nt' else 'clear')

    logo = color_text(r"""
.----------------.  .----------------.  .----------------.
| .--------------. || .--------------. || .--------------. |
| |    _______   | || |  _________   | || |   ______     | |
| |   /  ___  |  | || | |  _   _  |  | || |  |_   _ \    | |
| |  |  (__ \_|  | || | |_/ | | \_|  | || |    | |_) |   | |
| |   '.___`-.   | || |     | |      | || |    |  __'.   | |
| |  |`\____) |  | || |    _| |_     | || |   _| |__) |  | |
| |  |_______.'  | || |   |_____|    | || |  |_______/   | |
| |              | || |              | || |              | |
| '--------------' || '--------------' || '--------------' |
 '----------------'  '----------------'  '----------------'
 """, 'cyan')

    description = [
        "Siesta ToolBox Suite",
        "A comprehensive toolkit for SIESTA DFT simulations",
        f"Version {VERSION} | University of Brasilia - 2025",
        "Developed by Dr. Carlos M. O. Bastos"
    ]

    print(logo)
    print("\n" + "="*60)
    for line in description:
        print(line.center(60))
        sleep(0.2)
    print("="*60 + "\n")
    return

# Configure logger for warnings
logging.basicConfig(filename="warnings.log", level=logging.WARNING, format="%(message)s")

def warn_handler(message, category, filename, lineno, file=None, line=None):
    log_message = f"{category.__name__}: {message} (File: {filename}, Line: {lineno})"
    logging.warning(log_message)
    print("[WARNING] Warning detected! Check warnings.log for details.")

warnings.showwarning = warn_handler

__version__ = "1.0.0"

def compute_ecn(structure, mode, atoms_position=None):
    with open("structural_information.dat", "w") as f:

        # Lattice parameters
        lattice = structure.lattice
        print("\n[INFO] Lattice parameters:")
        print(f"   a = {lattice.a:.3f} Å")
        print(f"   b = {lattice.b:.3f} Å")
        print(f"   c = {lattice.c:.3f} Å")
        print(f"   Alpha = {lattice.alpha:.2f}°")
        print(f"   Beta = {lattice.beta:.2f}°")
        print(f"   Gamma = {lattice.gamma:.2f}°")

        f.write("\n[INFO] Lattice parameters:\n")
        f.write(f"   a = {lattice.a:.3f} Å\n")
        f.write(f"   b = {lattice.b:.3f} Å\n")
        f.write(f"   c = {lattice.c:.3f} Å\n")
        f.write(f"   Alpha = {lattice.alpha:.2f}°\n")
        f.write(f"   Beta = {lattice.beta:.2f}°\n")
        f.write(f"   Gamma = {lattice.gamma:.2f}°\n")

        # Lattice vectors
        lattice_vectors = structure.lattice.matrix
        print("\n[INFO] Writing Lattice vectors.")
        f.write("\nLattice vectors:\n")
        for i, vector in enumerate(lattice_vectors):
            f.write(f"   a_{i+1}: {vector[0]}   {vector[1]}   {vector[2]}\n")

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
            print("\n[INFO] Calculating the average ECN.")
            f.write("\nAverage ECN:\n")
            for method, value in ecn_avg.items():
                f.write(f"{method:15}: {value:.2f}\n")

        elif mode == "list" and atoms_position:
            for i in atoms_position:
                for method_name, method in methods.items():
                    try:
                        ecn_results[method_name].append(method.get_cn(structure, i-1))
                    except:
                        ecn_results[method_name].append(None)

            print("\n[INFO] Calculating the ECN for specified atoms.")
            f.write("\nECN for specified atoms:\n")
            for i, atom_index in enumerate(atoms_position):
                f.write(f" Atom {pos_atomics[atom_index-1][0]}:\n")
                f.write(f"   Element: {pos_atomics[atom_index-1][1]}     Cartesian Position: {pos_atomics[atom_index-1][2]}\n")
                for method, values in ecn_results.items():
                    print(f"      {method:15}: {values[i]}")
                    f.write(f"      {method:15}: {values[i]}\n")

        # Average bond distance calculation using CrystalNN
        print("\n[INFO] Calculating average bond distance...")
        f.write("\n[INFO] Average bond distance:\n")

        cnn = CrystalNN()
        distances = []

        indices = range(len(structure)) if mode == "mean" else [i - 1 for i in atoms_position]

        for i in indices:
            try:
                neighbors = cnn.get_nn_info(structure, i)
                for neighbor in neighbors:
                    dist = neighbor['site'].distance(structure[i])
                    distances.append(dist)
            except Exception as e:
                print(f"[WARNING] Failed to compute distances for atom {i+1}: {e}")

        if distances:
            avg_distance = np.mean(distances)
            print(f"   Average bond distance: {avg_distance:.4f} Å")
            f.write(f"   Average bond distance: {avg_distance:.4f} Å\n")
        else:
            print("[WARNING] No distances could be computed.")
            f.write("   No distances could be computed.\n")

        # Atomic positions
        pos_atomics = []
        f.write("\nAtomic positions:\n")
        for i, site in enumerate(structure):
            f.write(f"{i+1}  {site.specie.symbol} cartesian position: {site.coords}\n")
            pos_atomics.append([i+1, str(site.specie.symbol), site.coords])

def main():
    parser = argparse.ArgumentParser(description=f"Compute ECN from structure file. Version: {__version__}.")
    parser.add_argument("--file", required=True, help="Path to structure file.")
    parser.add_argument("--format", required=True, choices=["poscar", "cif", "siesta"], help="Input file format: poscar, cif, or siesta")
    parser.add_argument("--mode", choices=["list", "mean"], required=True, help="Calculation mode: list or mean")
    parser.add_argument("--list", type=str, help="List of atom indices (comma-separated). Example: [1,4,5,7] - Required for 'list' mode")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    args = parser.parse_args()

    if args.intro:
        show_intro()

    print("\n" + color_text("STRUCTURAL PROPERTIES:", 'bold'))
    print("-"*60)

    if args.mode == "list" and not args.list:
        parser.error("--list is required when --mode is 'list'")

    print("\n[INFO] Reading structure file...")
    atoms_position = list(map(int, args.list.strip('[]').split(','))) if args.list else None

    if args.format in ["poscar", "cif"]:
        structure = Structure.from_file(args.file)
    elif args.format == "siesta":
        atoms = ase_read(args.file, format="struct_out")
        structure = AseAtomsAdaptor.get_structure(atoms)
    else:
        raise ValueError("Unsupported file format. Use --format poscar, cif, or siesta.")

    compute_ecn(structure, args.mode, atoms_position)

    print("\n[INFO] Job complete!")
    print("\n"+"-"*60)

if __name__ == "__main__":
    main()

