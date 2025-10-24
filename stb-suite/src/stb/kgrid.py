#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################
    
VERSION = "1.8.0"    

from time import sleep
import argparse
import numpy as np
import math
import sys
import os
# Try to import ASE, required only for .cif files
try:
    import ase.io
except ImportError:
    pass 

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


def parse_poscar(filename):
    """Reads a POSCAR file and returns the lattice vectors as a 3x3 numpy array."""
    with open(filename, 'r') as f:
        lines = [l.strip() for l in f if l.strip()]
    scale = float(lines[1])
    vecs = []
    for i in range(2, 5):
        parts = lines[i].split()
        vec = [float(p) for p in parts]
        vecs.append(vec)
    lattice = np.array(vecs) * scale
    return lattice

def parse_cif(filename):
    """Reads a CIF file using ASE and returns the lattice vectors."""
    try:
        atoms = ase.io.read(filename)
        lattice = atoms.get_cell()
        return lattice
    except NameError:
        print("Error: The 'ase' library is required to read .cif files.")
        print("Please install it using: pip install ase")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading CIF file '{filename}': {e}")
        sys.exit(1)

def parse_fhi(filename):
    """Reads a geometry.in (FHI-aims) file and returns the lattice vectors as a 3x3 numpy array."""
    vecs = []
    with open(filename, 'r') as f:
        for line in f:
            # Remove comments and whitespace
            cleaned_line = line.split('#', 1)[0].strip()
            
            if cleaned_line.startswith('lattice_vector'):
                parts = cleaned_line.split()
                if len(parts) >= 4:
                    try:
                        # Extract x, y, z values (indices 1, 2, 3)
                        vec = [float(parts[1]), float(parts[2]), float(parts[3])]
                        vecs.append(vec)
                    except ValueError:
                        print(f"Error: 'lattice_vector' line malformed in {filename}: {line}")
                        sys.exit(1)

    # Check if we found exactly 3 vectors
    if len(vecs) != 3:
        print(f"Error: Could not find 3 'lattice_vector' lines in file {filename}.")
        print(f"Found: {len(vecs)}")
        sys.exit(1)
        
    lattice = np.array(vecs)
    return lattice

def parse_fdf(filename):
    """Reads an .fdf (Siesta) file and returns the lattice vectors as a 3x3 numpy array."""
    in_block = False
    all_values = []
    lattice_constant = 1.0 # Default value

    with open(filename, 'r') as f:
        for line in f:
            # Remove comments (anything after '#') and strip whitespace
            cleaned_line = line.split('#', 1)[0].strip()
            
            if not cleaned_line: # Skip blank lines
                continue
                
            parts = cleaned_line.split()
            lower_line = cleaned_line.lower()

            # Look for LatticeConstant (can be anywhere)
            if lower_line.startswith('latticeconstant'):
                try:
                    lattice_constant = float(parts[1])
                except (IndexError, ValueError):
                    print(f"Error: 'LatticeConstant' line malformed in {filename}: {line}")
                    sys.exit(1)
                continue 

            # Look for the start of the block
            if lower_line == '%block latticevectors':
                in_block = True
                continue 

            # Look for the end of the block
            if lower_line == '%endblock latticevectors':
                in_block = False
                continue 

            # If we are inside the block, collect all numerical values
            if in_block:
                for part in parts:
                    try:
                        all_values.append(float(part))
                    except ValueError:
                        # Ignore non-numeric parts if any
                        pass
    
    # After reading the file, check if we have 9 values
    if len(all_values) != 9:
        print(f"Error: Expected 9 values in LatticeVectors block, but found {len(all_values)}.")
        sys.exit(1)
        
    # Convert list to 3x3 array and apply the constant
    lattice = np.array(all_values).reshape(3, 3)
    lattice = lattice * lattice_constant
    return lattice

def compute_monkhorts(cella, cellb, cellc, k_density):
    """Computes the reciprocal vectors and the number of Monkhorst-Pack divisions."""
    volume = np.dot(cella, np.cross(cellb, cellc))
    
    if abs(volume) < 1e-9:
        print("Error: Cell volume is zero. Check lattice vectors.")
        sys.exit(1)
        
    b1 = 2 * np.pi * np.cross(cellb, cellc) / volume
    b2 = 2 * np.pi * np.cross(cellc, cella) / volume
    b3 = 2 * np.pi * np.cross(cella, cellb) / volume

    lengths = [np.linalg.norm(b) for b in (b1, b2, b3)]
    divisions = [max(1, math.ceil(length / k_density)) for length in lengths]
    return divisions

def print_density_recommendation():
    """Prints a friendly k-point density recommendation table."""
    print("\n" + "="*65)
    print("              📐 K-Point Density Recommendation Guide              ")
    print("="*65)
    print("  Density (1/Å⁻¹)        Accuracy Level")
    print("  ---------------      --------------------------")
    print("  0.05 – 0.1           High precision")
    print("  0.10 – 0.30          Medium precision")
    print("  0.30 – 0.50          Low precision")
    print()
    print("  ⚠️  Tip: For most systems, a density between 0.2 and 0.3 is")
    print("     generally accurate enough while keeping cost reasonable.")
    print("="*65 + "\n")

def analyze_dimensionality(divisions):
    """Analyzes the computed grid to suggest the system's dimensionality."""
    ones_count = divisions.count(1)
    
    print("--- Dimensionality Analysis ---")
    if ones_count == 3:
        # Grid is [1, 1, 1]
        print("System appears to be 0D (e.g., a molecule).")
        print("A 1x1x1 grid (Gamma point) is typically sufficient.")
    elif ones_count == 2:
        # Grid is [N, 1, 1] or [1, N, 1] or [1, 1, N]
        print("System appears to be 1D (e.g., a nanotube or polymer).")
        print("The '1's in the grid correspond to the vacuum-padded directions.")
    elif ones_count == 1:
        # Grid is [N, M, 1] or [N, 1, M] or [1, N, M]
        print("System appears to be 2D (e.g., a slab or surface).")
        print("The '1' in the grid corresponds to the vacuum-padded direction.")
    else:
        # Grid is [N, M, P]
        print("System appears to be 3D (bulk material).")
    print("---------------------------------\n")

def main():
    parser = argparse.ArgumentParser(
        description="Compute the Monkhorst-Pack grid based on desired k-point density and a structure file."
    )
    parser.add_argument(
        "--density", "-d", type=float, required=True,
        help="Target k-point density (in 1/Å). Example: 0.03"
    )
    parser.add_argument(
        "--file", "-f", type=str, required=True,
        help="Path to the structure file."
    )
    
    parser.add_argument(
        "--type", "-t", type=str, required=True,
        # Changed 'geometry' to 'fhi'
        choices=['poscar', 'cif', 'fhi', 'fdf'], 
        help="Type of the structure file. Currently supports: 'poscar', 'cif', 'fhi', 'fdf'."
    )


    parser.add_argument("-v", "--version", action="version",
                        version=f"stb-kgrid {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()


    if args.intro == True:
        show_intro()

    print("\n" + color_text("Suggested Monkhorst-Pack k-grid from structure :", 'bold'))
    print("-"*60)

    args = parser.parse_args()

    filename = args.file
    file_type = args.type.lower()
    
    try:
        # --- Updated Decision Logic ---
        if file_type == 'cif':
            print(f"ℹ️  Reading file '{filename}' as type '{file_type}' (using ASE)...")
            lattice = parse_cif(filename)
            
        elif file_type == 'poscar':
            print(f"ℹ️  Reading file '{filename}' as type '{file_type}' (native method)...")
            lattice = parse_poscar(filename)

        elif file_type == 'fhi': # Changed from 'geometry'
            print(f"ℹ️  Reading file '{filename}' as type '{file_type}' (native method)...")
            lattice = parse_fhi(filename) # Renamed function call

        elif file_type == 'fdf':
            print(f"ℹ️  Reading file '{filename}' as type '{file_type}' (native method)...")
            lattice = parse_fdf(filename)
        # --- End of Update ---
            
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        return
    except Exception as e:
        print(f"An unexpected error occurred while reading the file: {e}")
        return

    # Calculate divisions
    divisions = compute_monkhorts(lattice[0], lattice[1], lattice[2], args.density)

    # Print recommendation table
    print_density_recommendation()
    
    # Print the suggested grid
    print(f"✅ Suggested Monkhorst-Pack grid: {divisions[0]} {divisions[1]} {divisions[2]}\n")
    
    # --- New feature ---
    # Analyze and print dimensionality
    analyze_dimensionality(divisions)
    # --- End of new feature ---

if __name__ == "__main__":
    main()
