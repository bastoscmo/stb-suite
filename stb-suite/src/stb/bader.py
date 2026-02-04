#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.9.0"

import os
import sys
import argparse
import warnings
import multiprocessing
from time import sleep

# --- Warnings Configuration ---
# Filter out specific deprecation warnings from external libraries (PyBader/PkgResources)
# to keep the terminal output clean for the user.
warnings.filterwarnings("ignore", category=UserWarning, module="pybader")
warnings.filterwarnings("ignore", category=DeprecationWarning)

# --- Library Imports with Error Handling ---
# We wrap imports in try/except blocks to provide clear instructions
# if the user is missing the required scientific libraries.


try:
    import sisl
except ImportError:
    print("\033[91m[CRITICAL] Library 'sisl' not found.\033[0m")
    sys.exit(1)

try:
    from pybader.interface import Bader as PyBaderCalc
except ImportError:
    try:
        import pybader
        PyBaderCalc = pybader.interface.Bader
    except Exception:
        print("\033[91m[CRITICAL] Could not initialize PyBader interface. Install setuptools.\033[0m")
        sys.exit(1)


# ================= ANSI COLORS =================
# Standard color palette for the STB Suite CLI
# Cores ANSI para terminal
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
    """Retorna texto formatado com cor ANSI"""
    return f"{COLORS[color]}{text}{COLORS['reset']}"

def show_intro() -> None:
    """Exibe a introdução estilizada da STB-SUITE"""
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

# ================= SCIENTIFIC DATA =================

# Default Valence Electron Dictionary (Z_val).
# IMPORTANT: This list represents the number of valence electrons typically
# used in standard pseudopotentials.
# NOTE: Transition metals (like Mo, W, Ti) may vary depending on whether
# semicore states are included in the pseudopotential generation.
VALENCE_DICT = {
    "H": 1.0,  "Li": 1.0, "Be": 2.0, "B": 3.0,  "C": 4.0,  "N": 5.0,  "O": 6.0,  "F": 7.0,
    "Na": 1.0, "Mg": 2.0, "Al": 3.0, "Si": 4.0, "P": 5.0,  "S": 6.0,  "Cl": 7.0,
    "K": 1.0,  "Ca": 2.0, "Ti": 4.0, "V": 5.0,  "Cr": 6.0, "Mn": 7.0, "Fe": 8.0, "Co": 9.0,
    "Ni": 10.0,"Cu": 11.0,"Zn": 12.0,"Ga": 3.0, "Ge": 4.0, "As": 5.0, "Se": 6.0, "Br": 7.0,
    "Mo": 14.0,"W": 14.0, "Au": 11.0 
}

def print_valence_warning():
    """Prints a warning to the user about checking their Pseudopotentials."""
    msg = (
        "[WARNING] Using DEFAULT valence electron counts.\n"
        "          Please verify if these match your Pseudopotentials (.psf/.psml)."
    )
    print(color_text(msg, 'yellow'))
    print("-" * 60)

# ================= MAIN LOGIC =================

def solve_bader(label, output_file=None, speed_mode='normal'):
    """
    Main function to execute the Bader Analysis workflow.
    
    Workflow:
    1. Read Siesta binary grid (.RHO) and Geometry (.XV) using sisl.
    2. Convert data to Gaussian Cube format (handles non-orthogonal cells).
    3. Run PyBader analysis on the Cube file.
    4. Normalize results (fix Angstrom vs Bohr unit mismatch).
    5. Print and save the report.
    """
    
    # Define filenames based on the SystemLabel
    file_rho = f"{label}.RHO"
    file_xv = f"{label}.XV"
    file_fdf = f"{label}.fdf"
    file_cube = f"{label}.cube"

    # Set default output filename if not provided
    if not output_file:
        output_file = f"{label}_BADER.txt"

    # Validation: Check if the grid file exists
    if not os.path.exists(file_rho):
        print(color_text(f"[ERROR] Grid file '{file_rho}' not found.", 'red'))
        return

    # Print initialization info
    print(f"[INFO] System: {color_text(label, 'bold')} | Mode: {color_text(speed_mode.upper(), 'cyan')}")
    print_valence_warning()

    # --- STEP 1: SISL (Reading and Conversion) ---
    try:
        print(f"1. [SISL] Reading geometry and charge density...")
        
        # Geometry Reading Priority:
        # 1. .XV file: Contains the final relaxed positions (Best choice).
        # 2. .fdf file: Contains input positions (Fallback).
        if os.path.exists(file_xv):
            geometry = sisl.get_sile(file_xv).read_geometry()
        elif os.path.exists(file_fdf):
            geometry = sisl.get_sile(file_fdf).read_geometry()
        else:
            print(color_text("[ERROR] No geometry file (.XV or .fdf) found.", 'red'))
            return

        # Read the Charge Density Grid (.RHO)
        rho_grid = sisl.get_sile(file_rho).read_grid()
        
        # Crucial Step: Attach the geometry object to the grid.
        # This ensures that when we write the Cube file, the atom positions
        # and unit cell vectors are correctly written in the header.
        rho_grid.set_geometry(geometry)
        
        # Export to Cube format (PyBader reads Cube natively)
        rho_grid.write(file_cube)
        
    except Exception as e:
        print(color_text(f"[ERROR] SISL processing failed: {e}", 'red'))
        return

    # --- STEP 2: PyBader (Calculation) ---
    # Auto-detect available CPU cores for parallel processing
    n_threads = multiprocessing.cpu_count()
    print(f"2. [PyBader] Starting calculation on {n_threads} threads...")
    
    try:
        # Load the temporary Cube file into PyBader
        bader_job = PyBaderCalc.from_file(file_cube, threads=n_threads)
        
        # Optimization: 'fast' mode
        # Uses 'minimum_distance' assignment instead of 'mass_weighted'.
        # Faster for large grids, but edges might be slightly less precise.
        if speed_mode == 'fast':
            print("   [FAST] Using 'minimum distance' refinement.")
            if hasattr(bader_job, 'refinement_method'):
                bader_job.refinement_method = 'minimum_distance' 
        
        # Execute the Bader Partitioning
        bader_job()
        
        # Retrieve raw electron populations per atom
        raw_populations = bader_job.atoms_charge
    except Exception as e:
        print(color_text(f"[ERROR] PyBader calculation failed: {e}", 'red'))
        return

    # --- STEP 3: Analysis and Normalization ---
    print("3. [Analysis] checking units...")
    
    total_theory = 0.0
    atoms_data = []

    # Map raw results to atom names and calculate theoretical valence
    for i, atom in enumerate(geometry.atoms):
        if i >= len(raw_populations): break
        sym = atom.symbol
        z_val = VALENCE_DICT.get(sym, 0.0)
        total_theory += z_val
        atoms_data.append({'id': i+1, 'sym': sym, 'z_val': z_val, 'pop_raw': raw_populations[i]})

    # Check for Unit Mismatch (Angstrom vs Bohr)
    # Siesta writes .cube in Angstroms, but standard Cube format expects Bohr.
    # This causes a cubic scaling error (~6.7x difference in total charge).
    total_raw = sum(raw_populations)
    
    # Avoid division by zero
    ratio = total_raw / total_theory if total_theory > 0 else 1.0
    
    # If the difference is > 10%, assume unit mismatch and calculate correction factor.
    correction_factor = 1.0 / ratio if abs(ratio - 1.0) > 0.1 else 1.0
    is_corrected = correction_factor != 1.0

    if is_corrected:
        print(color_text(f"   [INFO] Unit mismatch detected. Correction factor: {correction_factor:.4f}", 'cyan'))

    # --- STEP 4: Output Formatting ---
    out_lines = []
    out_lines.append(f"BADER CHARGE ANALYSIS REPORT - STB Suite v{VERSION}")
    out_lines.append(f"System: {label} | Speed Mode: {speed_mode}")
    out_lines.append("=" * 75)
    
    # Table Header
    header = f"{'Idx':<4} {'Elem':<5} {'Pop(e-)':<12} {'Z_val':<8} {'Net Charge':<12} {'State':<15}"
    out_lines.append(header)
    out_lines.append("-" * 75)
    
    # Print header to console
    print("\n" + header)
    print("-" * 75)

    total_final = 0.0

    # Process each atom for the final report
    for data in atoms_data:
        # Apply correction factor
        pop = data['pop_raw'] * correction_factor
        
        # Calculate Net Charge (Valence - Population)
        net = data['z_val'] - pop
        
        # Determine Oxidation State (Donor vs Acceptor)
        if net > 0.05:
            state, color = "Donor (+)", 'red'
        elif net < -0.05:
            state, color = "Acceptor (-)", 'blue'
        else:
            state, color = "Neutral", 'reset'

        # Formatted line for text file (no colors)
        line = f"{data['id']:<4} {data['sym']:<5} {pop:<12.4f} {data['z_val']:<8.2f} {net:<+12.4f} {state:<15}"
        out_lines.append(line)
        
        # Formatted line for terminal (with colors)
        print(f"{data['id']:<4} {data['sym']:<5} {pop:<12.4f} {data['z_val']:<8.2f} "
              f"{color_text(f'{net:<+12.4f}', 'bold')} {color_text(state, color)}")
        
        total_final += pop

    # Footer Statistics
    footer = [
        "-" * 75,
        f"Total Integrated: {total_final:.4f} (Target: {total_theory:.2f})",
        "=" * 75
    ]
    out_lines.extend(footer)
    for l in footer: print(l)

    # Save to file
    try:
        with open(output_file, "w") as f:
            f.write("\n".join(out_lines))
        print(f"\n{color_text('[OK]', 'green')} Results saved to: {color_text(output_file, 'bold')}")
    except IOError:
        print(color_text(f"[ERROR] Could not save file {output_file}", 'red'))

    # Optional: Clean up temporary cube file to save disk space
    # if os.path.exists(file_cube): os.remove(file_cube)

# ================= EXECUTION =================

def main():
    parser = argparse.ArgumentParser(description="STB Bader Tool")
    
    parser.add_argument("-l", "--label", required=True, 
                        help="The SystemLabel used in Siesta (e.g., 'siesta' or 'graphene')")
    parser.add_argument("-o", "--output", required=False, 
                        help="Name of the output text file.")
    parser.add_argument("--speed", choices=['normal', 'fast'], default='normal',
                        help="Optimization level. 'normal' (precise) or 'fast' (less refined edges).")
    parser.add_argument("--no-intro", dest="intro", action="store_false", 
                        help="Skip intro banner")

    args = parser.parse_args()

    if args.intro:
        show_intro()

    solve_bader(args.label, args.output, args.speed)
    
    print("\n" + "-" * 60)
    print(color_text("Electron counting is like accounting, but the currency is negative.\n", 'bold'))

if __name__ == "__main__":
    main()
