#!/usr/bin/env python3

#################################################
# Siesta Tool Box - strain generator            #
# version 1.5.1                                 #
# UnB - 2025/04/22                              #
# Dr. Carlos M. O. Bastos                       #
# Modified for automatic strain type detection   #
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
import os
import argparse
import numpy as np



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


def determine_strain_type(direction):
    """Determine strain type based on direction input."""
    if len(direction) == 1:  # x, y, z
        return 'uniaxial'
    elif len(direction) == 2:  # xy, xz, yx, yz, zx, zy
        return 'biaxial'
    else:
        raise ValueError("[FAIL] Invalid direction. Use x, y, z for uniaxial or combinations like xy, yz for biaxial.")

def normalize_direction(direction):
    """Normalize direction input (e.g., yx -> xy)."""
    if len(direction) == 2:
        return ''.join(sorted(direction.lower()))
    return direction.lower()

def apply_cartesian_strain(lattice_vectors, strain, direction):
    """
    Apply uniaxial or biaxial strain in Cartesian coordinates.
    
    Args:
        lattice_vectors: 3x3 numpy array of lattice vectors
        strain: strain value (positive for tension, negative for compression)
        direction: strain direction (x, y, z, xy, xz, yz, etc.)
        
    Returns:
        Strained lattice vectors
    """
    # Normalize and validate direction
    direction = normalize_direction(direction)
    valid_directions = {'x', 'y', 'z', 'xy', 'xz', 'yz'}
    if direction not in valid_directions:
        raise ValueError(f"[FAIL] Invalid direction '{direction}'. Use x, y, z, xy, xz, or yz.")

    # Create strain tensor
    strain_tensor = np.zeros((3, 3))
    
    if len(direction) == 1:  # Uniaxial
        if direction == 'x':
            strain_tensor[0, 0] = strain
        elif direction == 'y':
            strain_tensor[1, 1] = strain
        elif direction == 'z':
            strain_tensor[2, 2] = strain
            
    else:  # Biaxial
        if 'x' in direction:
            strain_tensor[0, 0] = strain
        if 'y' in direction:
            strain_tensor[1, 1] = strain
        if 'z' in direction:
            strain_tensor[2, 2] = strain
    
    # Apply strain transformation: new_vec = (I + ε) · vec
    identity = np.eye(3)
    transformation = identity + strain_tensor
    
    # Transform each lattice vector
    strained_vectors = np.dot(transformation, lattice_vectors.T).T
    
    return strained_vectors

def main():
    parser = argparse.ArgumentParser(
        description="Applies strain in Cartesian coordinates to a SIESTA FDF file. "
                   "Type (uniaxial/biaxial) is inferred from direction. "
                   "IMPORTANT: atomic coordinates must be in fractional."
    )
    parser.add_argument("--file", required=True, help="Input FDF file.")
    parser.add_argument("--stdir", required=True, 
                       help="Direction of strain: x, y, z for uniaxial; xy, xz, yz, etc. for biaxial.")
    parser.add_argument("--stmin", type=float, default=0,
                       help="Minimum strain percentage (default: 0). Can be negative for compression.")
    parser.add_argument("--stmax", type=float, default=25,
                       help="Maximum strain percentage (default: 25).")
    parser.add_argument("--step", type=float, default=1,
                       help="Strain step percentage (default: 1).")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    args = parser.parse_args()

    # Validate strain range
    if args.stmin > args.stmax:
        raise ValueError("[FAIL] Minimum strain cannot be greater than maximum strain.")
    
    # Determine strain type from direction
    strain_type = determine_strain_type(args.stdir)
    norm_dir = normalize_direction(args.stdir)
    
    if args.intro == True:
        show_intro()

    print("\n" + color_text("STRAIN:", 'bold'))
    print("-"*60)
    
    print(f"[INFO] Detected strain type: {strain_type} (from direction '{args.stdir}')")
    print(f"[INFO] Strain direction: {norm_dir}")
    print(f"[INFO] Strain range: {args.stmin}% to {args.stmax}% with step {args.step}%")
    print("[INFO] Read file")
    print("[INFO] Generating FDF files with strain...")

    with open(args.file, "r") as f:
        lines = f.readlines()

    # Locate LatticeVectors section
    lv_start = None
    for i, line in enumerate(lines):
        if "LatticeVectors" in line:
            lv_start = i + 1
            break

    if lv_start is None:
        raise ValueError("[FAIL] LatticeVectors section not found in the FDF file.")

    lattice_vectors = []
    for i in range(3):
        vec = list(map(float, lines[lv_start + i].split()))
        lattice_vectors.append(np.array(vec))
    lattice_vectors = np.array(lattice_vectors)

    # Generate strained structures
#    strain_values = [i / 100 for i in range(args.stmin, args.stmax + 1, args.step)]
    strain_values = list(np.arange(args.stmin, args.stmax + args.step, args.step) / 100)


    for strain in strain_values:
        # Handle negative strain (compression) in folder name
        strain_prefix = "m" if strain < 0 else ""
        folder = f"strain_{norm_dir}_{strain_prefix}{abs((strain * 100)):.2f}"
        os.makedirs(folder, exist_ok=True)

        new_vectors = apply_cartesian_strain(lattice_vectors, strain, norm_dir)
        new_lv_lines = [f"  {'  '.join(f'{x:.6f}' for x in vec)}\n" for vec in new_vectors]
        new_lines = lines[:lv_start] + new_lv_lines + lines[lv_start + 3:]

        output_fdf = os.path.join(folder, args.file)
        with open(output_fdf, "w") as f:
            f.writelines(new_lines)
        print(f"[OK] Generated: {output_fdf} (strain: {strain*100:.1f}%)")
        
    print("[INFO] Complete job!") 
    print("\n"+"-"*60)
    print(color_text("This lattice has more tension than my last Zoom meeting.\n\n", 'bold'))        

if __name__ == "__main__":
    main()
