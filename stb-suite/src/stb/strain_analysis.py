#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
#     Mechanical Properties Analyzer            #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.9.0"

import os
import sys
import re
import argparse
import numpy as np
from scipy.stats import linregress
# Use trapezoid for compatibility with new Scipy
from scipy.integrate import trapezoid 

# ==========================================
#           UI / VISUALS
# ==========================================

COLORS = {
    'reset': '\033[0m',
    'cyan': '\033[96m',
    'blue': '\033[94m',
    'green': '\033[92m',
    'yellow': '\033[93m',
    'red': '\033[91m',
    'bold': '\033[1m',
    'underline': '\033[4m',
    'magenta': '\033[95m'
}

def color_text(text: str, color: str) -> str:
    """Returns text formatted with ANSI color codes."""
    return f"{COLORS.get(color, COLORS['reset'])}{text}{COLORS['reset']}"

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
    print(color_text(f"   Siesta Tool Box - Mechanical Analysis v{VERSION}", 'yellow'))
    print(color_text("   Developed by Dr. Carlos M. O. Bastos", 'blue'))
    print(color_text("   bastoscmo.github.io\n", 'blue'))

# ==========================================
#           CORE LOGIC
# ==========================================

def parse_folder_name(folder_name):
    """
    Parses the folder name to extract direction and strain value.
    Format: strain_{dir}_{prefix}{value} -> e.g., strain_xx_m1.00
    """
    match = re.search(r"strain_([a-zA-Z0-9]+)_(m?)(\d+\.\d+)", folder_name)
    
    if match:
        direction = match.group(1)
        is_negative = match.group(2) == 'm'
        value = float(match.group(3))
        if is_negative:
            value = -value
        return direction, value
    return None, None

def get_stress_from_file(filepath):
    """
    Reads the Siesta output file and extracts the LAST Voigt stress tensor.
    Returns [xx, yy, zz, yz, xz, xy] in kBar.
    """
    last_stress = None
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if "Stress tensor Voigt" in line and "(kbar)" in line:
                    parts = line.split(":")[-1].split()
                    if len(parts) >= 6:
                        last_stress = [float(x) for x in parts[:6]]
    except Exception:
        return None
    return last_stress

def get_lattice_z_vector(filepath):
    """Extracts the Z-vector magnitude from output for 2D height normalization."""
    z_len = 1.0
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
            for i in range(len(lines)-1, 0, -1):
                if "outcell: Unit cell vectors" in lines[i]:
                    vec_z = lines[i+3].split()
                    z_len = np.linalg.norm([float(x) for x in vec_z[:3]])
                    break
    except:
        pass
    return z_len

def calculate_yield_stress(strain, stress, modulus, offset=0.002):
    """
    Calculates 0.2% Offset Yield Strength.
    Returns (Yield Stress, Yield Strain).
    """
    # Offset line: Stress = E * (Strain - 0.002)
    offset_line = modulus * (strain - offset)
    
    idx = np.where(strain > 0)[0]
    if len(idx) < 2: return 0.0, 0.0
    
    diff = stress[idx] - offset_line[idx]
    for i in range(len(diff)-1):
        if diff[i] * diff[i+1] < 0:
            return stress[idx][i], strain[idx][i]
            
    return 0.0, 0.0

def analyze_mechanics(data, direction, is_2d=False, z_height=20.0):
    """
    Calculates Young's Modulus, UTS, Toughness.
    data: [Strain%, StrainFrac, Sxx, Syy, Szz, Syz, Sxz, Sxy] (kBar)
    """
    # Standard 3D: kBar -> GPa (1 kBar = 0.1 GPa)
    factor = 0.1 
    unit = "GPa"
    
    if is_2d:
        # Convert kBar to N/m using Z-height
        # 1 kBar = 1e8 Pa (N/m^2)
        # Z_height is in Angstrom (1e-10 m)
        # Stress_2D = Stress_kBar * 1e8 * Z * 1e-10
        # Stress_2D = Stress_kBar * Z * 0.01
        factor = 0.01 * z_height
        unit = "N/m"

    strain_pct = data[:, 0]
    strain_frac = data[:, 1]
    stress_raw = data[:, 2:8] # xx, yy, zz, yz, xz, xy
    
    # Identify Axial Component
    d_map = {'xx': 0, 'yy': 1, 'zz': 2, 'x': 0, 'y': 1, 'z': 2}
    ax_idx = d_map.get(direction.lower(), 0) # Default to xx if unknown
    
    sigma_axial = stress_raw[:, ax_idx] * factor
    
    # --- Young's Modulus (Stiffness) ---
    # Linear fit only in small strain region (-2% to 2%)
    limit_idx = np.where(abs(strain_pct) <= 2.0) 
    if len(limit_idx[0]) < 3: limit_idx = range(len(strain_pct))
    
    slope, intercept, r_val, _, _ = linregress(strain_frac[limit_idx], sigma_axial[limit_idx])
    youngs_modulus = slope
    
    # --- UTS (Ultimate Tensile Strength) ---
    uts = np.max(sigma_axial)
    
    # --- Toughness (Energy Density) ---
    toughness = trapezoid(sigma_axial, strain_frac)
    
    # --- Yield Strength (0.2% Offset) ---
    yield_stress, yield_strain = calculate_yield_stress(strain_frac, sigma_axial, youngs_modulus)
    
    return {
        'modulus': youngs_modulus,
        'uts': uts,
        'toughness': toughness,
        'yield': yield_stress,
        'r_squared': r_val**2,
        'unit': unit,
        'axial_stress': sigma_axial
    }

# ==========================================
#           MAIN
# ==========================================

def main():
    parser = argparse.ArgumentParser(description="Mechanical Properties from Stress-Strain Data.")
    parser.add_argument("-f", "--file", required=True, help="Siesta output file (e.g., calc.out).")
    parser.add_argument("-o", "--output", default="mechanical_curve.dat", help="Output raw data file.")
    parser.add_argument("--2d", dest="is2d", action="store_true", help="Enable 2D units (N/m).")
    parser.add_argument("--thickness", type=float, default=20.0, help="Vacuum height (Z) for 2D conversion (Angstrom).")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Hide intro.")
    
    args = parser.parse_args()

    if args.intro:
        show_intro()
        
    print(color_text("-> Scanning directories...", 'green'))
    
    folders = [d for d in os.listdir('.') if os.path.isdir(d) and d.startswith('strain_')]
    if not folders:
        sys.exit(color_text("[ERROR] No 'strain_*' folders found.", 'red'))
        
    data = []
    detected_dir = "xx"
    
    # Detect Z-height automatically if 2D
    if args.is2d and os.path.exists(os.path.join(folders[0], args.file)):
        z_auto = get_lattice_z_vector(os.path.join(folders[0], args.file))
        if z_auto > 1.0:
            print(f"{color_text('[INFO]', 'cyan')} Detected Cell Z-Height: {z_auto:.2f} Ang (used for N/m conversion)")
            args.thickness = z_auto

    print(f"   Found {len(folders)} strain steps.")
    
    for f in folders:
        d, val = parse_folder_name(f)
        if d: detected_dir = d
        
        fpath = os.path.join(f, args.file)
        if os.path.exists(fpath):
            stress = get_stress_from_file(fpath)
            if stress:
                row = [val, val/100.0] + stress
                data.append(row)
                print(f"   Reading {f}: Strain {val:>5.2f}% OK")
    
    if not data:
        sys.exit(color_text("[ERROR] No stress data found.", 'red'))
        
    # Sort and Convert
    data = np.array(sorted(data, key=lambda x: x[0]))
    
    # Perform Analysis
    results = analyze_mechanics(data, detected_dir, args.is2d, args.thickness)
    u = results['unit']
    
    # --- OUTPUT TO SCREEN ---
    print("\n" + "="*50)
    print(color_text("   MECHANICAL PROPERTIES REPORT", 'bold').center(50))
    print("="*50)
    print(f"Direction       : {detected_dir.upper()}")
    print(f"Dimensionality  : {'2D (Sheet)' if args.is2d else '3D (Bulk)'}")
    print("-" * 50)
    
    print(f"Young's Modulus : {color_text(f'{results['modulus']:.2f} {u}', 'green')} (R²={results['r_squared']:.4f})")
    print(f"UTS (Max Stress): {color_text(f'{results['uts']:.2f} {u}', 'red')}")
    print(f"Toughness       : {results['toughness']:.4f} {'J/m^2' if args.is2d else 'GJ/m^3'}")
    
    if results['yield'] > 0:
        print(f"Yield Str (0.2%): {results['yield']:.2f} {u}")
    else:
        print(f"Yield Str (0.2%): {color_text('Not detected (Linear)', 'yellow')}")
        
    print("="*50)
    
    # --- SAVE FILES ---
    
    # 1. Curve Data
    header = (f"Strain-Stress Curve | Dir: {detected_dir} | Unit: {u}\n"
              f"1:Strain(%) 2:Strain(Frac) 3:Stress_Axial({u})")
    
    out_data = np.column_stack((data[:, 0], data[:, 1], results['axial_stress']))
    np.savetxt(args.output, out_data, fmt="%12.6f", header=header)
    print(f"\n{color_text('[Saved]', 'cyan')} Curve data -> {args.output}")
    
    # 2. Report Text
    report_file = "mechanical_report.txt"
    with open(report_file, "w") as f:
        f.write("========================================\n")
        f.write("      MECHANICAL PROPERTIES REPORT      \n")
        f.write("========================================\n")
        f.write(f"Direction: {detected_dir}\n")
        f.write(f"Young's Modulus: {results['modulus']:.4f} {u}\n")
        f.write(f"UTS:             {results['uts']:.4f} {u}\n")
        f.write(f"Yield (0.2%):    {results['yield']:.4f} {u}\n")
        f.write(f"Toughness:       {results['toughness']:.4f}\n")
        
    print(f"{color_text('[Saved]', 'cyan')} Summary    -> {report_file}\n")

if __name__ == "__main__":
    main()
