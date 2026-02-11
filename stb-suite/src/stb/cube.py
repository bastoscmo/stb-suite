#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.8.0"

import os
import sys
import argparse
import numpy as np
from time import sleep

# Try to import sisl
try:
    import sisl
except ImportError:
    print("\n\033[91m[CRITICAL ERROR] sisl library not found.\033[0m")
    print("Please install it using: pip install sisl")
    sys.exit(1)

# ANSI Colors for terminal
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
    """Returns text formatted with ANSI color."""
    return f"{COLORS[color]}{text}{COLORS['reset']}"

def show_intro() -> None:
    """Displays the stylized STB-SUITE introduction."""
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
        "Siesta ToolBox Suite - Grid Converter",
        "Professional Grid to Cube Converter",
        f"Version {VERSION} | University of Brasilia - 2025",
        "Developed by Dr. Carlos M. O. Bastos"
    ]

    print(logo)
    print("\n" + "="*60)
    for line in description:
        print(line.center(60))
        sleep(0.1)
    print("="*60 + "\n")
    return

def get_cell_matrix(lattice_obj):
    """
    Robustly retrieves the unit cell matrix from sisl object.
    Supports: .cell, .matrix, .vectors, .get_cell()
    """
    attributes = ['cell', 'matrix', 'vectors']
    
    for attr in attributes:
        if hasattr(lattice_obj, attr):
            return getattr(lattice_obj, attr)
            
    if hasattr(lattice_obj, 'get_cell') and callable(lattice_obj.get_cell):
        return lattice_obj.get_cell()
        
    try:
        return np.array(lattice_obj)
    except:
        return None

def main():
    parser = argparse.ArgumentParser(
        description="Convert SIESTA Grid files (VT, VH, RHO, BADER) to Gaussian Cube format.",
        epilog="Example usage:\n"
               "  stb_convert --label graphene --type RHO\n"
               "  stb_convert -l system -t VT -o potential.cube",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("-l", "--label", help="System Label (required).", required=True)
    parser.add_argument("-t", "--type", help="File type to convert.", choices=['VT', 'VH', 'RHO', 'BADER'], required=True)
    parser.add_argument("-o", "--output", help="Custom output filename. Default: label_type.cube", default=None)
    
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction.")
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")

    args = parser.parse_args()

    # Filenames
    input_file = f"{args.label}.{args.type}"
    xv_file = f"{args.label}.XV"
    output_file = args.output if args.output else f"{args.label}_{args.type}.cube"

    if args.intro:
        show_intro()

    print("\n" + color_text("GRID CONVERTER:", 'bold'))
    print("-" * 60)

    # --- 1. Check Input ---
    if not os.path.exists(input_file):
        print(f"{COLORS['red']}[ERROR] Input file '{input_file}' not found.{COLORS['reset']}")
        sys.exit(1)

    # --- 2. Load Grid ---
    print(f"\n[INFO] Reading SIESTA file: {color_text(input_file, 'cyan')}...")
    try:
        sile = sisl.get_sile(input_file)
        grid = sile.read_grid()
    except Exception as e:
        print(f"{COLORS['red']}[FATAL ERROR] Could not read grid: {e}{COLORS['reset']}")
        sys.exit(1)

    # --- 3. Load Geometry (Optional) ---
    if os.path.exists(xv_file):
        print(f"[INFO] Found geometry file: {color_text(xv_file, 'cyan')}. Mapping atoms...")
        try:
            geom = sisl.get_sile(xv_file).read_geometry()
            grid.geometry = geom
        except Exception as e:
            print(f"{COLORS['yellow']}[WARNING] Failed to read geometry from .XV: {e}{COLORS['reset']}")
    else:
        print(f"{COLORS['yellow']}[WARNING] No .XV file found. Atoms will be missing in the .cube output.{COLORS['reset']}")

    # --- 4. Metadata Display ---
    print(f"[INFO] Grid Metadata:")
    print(f"       Shape:   {grid.shape}")
    
    cell = get_cell_matrix(grid.lattice)
    if cell is not None:
        print(f"       Lattice: [{cell[0,0]:.2f}, {cell[1,1]:.2f}, {cell[2,2]:.2f}] Ang")
    
    if args.type.upper() == "RHO":
        # Approximate integral check
        total_val = grid.grid.sum() * grid.volume / (grid.shape[0] * grid.shape[1] * grid.shape[2])
        print(f"       Integral: {total_val:.4f} e (Check if close to total electrons)")

    # --- 5. Write Output ---
    print(f"[INFO] Converting to Cube format...")
    try:
        grid.write(output_file)
        print(f"[SUCCESS] File saved as: {color_text(output_file, 'green')}")
    except Exception as e:
        print(f"{COLORS['red']}[ERROR] Failed to write output file: {e}{COLORS['reset']}")
        sys.exit(1)

    print("\n[INFO] Complete job!") 
    print("\n"+"-"*60)
    # Funny closing phrase
    print(color_text("Conversion complete. It's hip to be square (cube).\n\n", 'bold'))

if __name__ == "__main__":
    main()
