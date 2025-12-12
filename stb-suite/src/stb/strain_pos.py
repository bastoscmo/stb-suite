#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.8.0"

import os
import sys
import re
import argparse
import numpy as np

# ANSI Colors for terminal (Same style as strain.py)
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
    """Returns text formatted with ANSI color codes."""
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
    
    print(logo)
    print(color_text(f"   Siesta Tool Box - Stress-Strain Post-Processing v{VERSION}", 'yellow'))
    print(color_text("   Developed by Dr. Carlos M. O. Bastos", 'blue'))
    print(color_text("   bastoscmo.github.io\n", 'blue'))

def parse_folder_name(folder_name):
    """
    Parses the folder name to extract direction and strain value.
    Expected format: strain_{dir}_{prefix}{value}
    Example: strain_x_m1.00 ('m' indicates minus/negative)
    """
    # Regex to capture direction, sign (optional), and value
    # Pattern: strain_ + (letters) + _ + (optional 'm') + (numbers.numbers)
    match = re.search(r"strain_([a-zA-Z0-9]+)_(m?)(\d+\.\d+)", folder_name)
    
    if match:
        direction = match.group(1)
        is_negative = match.group(2) == 'm'
        value_str = match.group(3)
        
        value = float(value_str)
        if is_negative:
            value = -value
            
        # The value in the folder name is percentage (e.g., 1.00 for 1%)
        return direction, value
    return None, None

def get_stress_from_file(filepath):
    """
    Reads the Siesta output file and extracts the LAST Voigt stress tensor.
    Returns a list [xx, yy, zz, yz, xz, xy] or None if not found.
    """
    last_stress = None
    
    try:
        with open(filepath, 'r') as f:
            for line in f:
                # Looks for the specific Voigt tensor line in kbar
                # Example: siesta: Stress tensor Voigt[x,y,z,yz,xz,xy] (kbar): ...
                if "Stress tensor Voigt" in line and "(kbar)" in line:
                    # Split by colon and take the part with numbers
                    parts = line.split(":")[-1].split()
                    if len(parts) >= 6:
                        # Convert to float and store (overwrites previous ones to get the final relaxed step)
                        last_stress = [float(x) for x in parts[:6]]
    except Exception as e:
        print(color_text(f"[ERROR] Failed to read {filepath}: {e}", 'red'))
        return None

    return last_stress

def main():
    parser = argparse.ArgumentParser(description="Extracts Stress-Strain curve data from Siesta outputs.")
    parser.add_argument("-f", "--file", type=str, required=True, help="Name of the Siesta output file inside the folders (e.g., calc.out).")
    parser.add_argument("-o", "--output", type=str, default="stress_strain_curve.dat", help="Name of the final data file (default: stress_strain_curve.dat).")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    
    args = parser.parse_args()

    if args.intro:
    	show_intro()
    print(color_text("-> Starting folder analysis...", 'green'))
    
    # Find all folders starting with 'strain_'
    all_items = os.listdir('.')
    strain_folders = [d for d in all_items if os.path.isdir(d) and d.startswith('strain_')]
    
    if not strain_folders:
        print(color_text("[WARNING] No 'strain_*' folders found in current directory.", 'red'))
        sys.exit(1)
        
    data_points = []
    
    print(f"   Found {len(strain_folders)} strain folders.")
    print("   Reading files and extracting tensors...\n")
    
    detected_direction = "Unknown"
    
    for folder in strain_folders:
        direction, strain_pct = parse_folder_name(folder)
        
        if direction is None:
            continue
            
        detected_direction = direction # Assumes all folders belong to the same experiment type
        
        filepath = os.path.join(folder, args.file)
        
        if not os.path.isfile(filepath):
            print(color_text(f"   [SKIP] File not found: {filepath}", 'yellow'))
            continue
            
        stress_voigt = get_stress_from_file(filepath)
        
        if stress_voigt:
            # Store: [Strain%, Sxx, Syy, Szz, Syz, Sxz, Sxy]
            # Siesta Voigt Stress: [x, y, z, yz, xz, xy]
            data_points.append([strain_pct] + stress_voigt)
            print(f"   Processed: {folder} -> Strain: {strain_pct:>6.2f}% | Sxx: {stress_voigt[0]:>8.2f} kBar")
        else:
            print(color_text(f"   [FAIL] Stress tensor not found in: {filepath}", 'red'))

    if not data_points:
        print(color_text("\n[ERROR] No valid data extracted.", 'red'))
        sys.exit(1)

    # Sort data by strain value (column 0)
    data_points.sort(key=lambda x: x[0])
    data_points = np.array(data_points)
    
    # Save file
    header = (f"Post-processing for Stress-Strain Curve\n"
              f"Direction detected: {detected_direction}\n"
              f"Columns:\n"
              f"1: Strain (%)   2: Strain (Frac)   "
              f"3: Sigma_xx (kBar)   4: Sigma_yy (kBar)   5: Sigma_zz (kBar)   "
              f"6: Sigma_yz (kBar)   7: Sigma_xz (kBar)   8: Sigma_xy (kBar)")
    
    # Prepare data for saving: Add fractional strain column
    # New order: Strain%, StrainFrac, xx, yy, zz, yz, xz, xy
    final_data = np.column_stack((data_points[:, 0], data_points[:, 0]/100.0, data_points[:, 1:]))
    
    np.savetxt(args.output, final_data, header=header, fmt='%12.6f')
    
    print(color_text(f"\n[SUCCESS] Data saved to: {args.output}", 'cyan'))
    print("-------------------------------------------------------")
    print(f"Total points: {len(data_points)}")
    print(f"Deformation type detected: {detected_direction.upper()}")
    print("-------------------------------------------------------")

if __name__ == "__main__":
    main()
