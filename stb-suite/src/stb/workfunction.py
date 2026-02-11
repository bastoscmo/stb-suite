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
import numpy as np
import matplotlib.pyplot as plt
from time import sleep

# Try to import sisl
try:
    import sisl
except ImportError:
    print("\n\033[91m[CRITICAL ERROR] sisl library not found.\033[0m")
    print("Please install it using: pip install sisl")
    sys.exit(1)

# Constants
Ry2eV = 13.605693123

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
        "Siesta ToolBox Suite - Work Function",
        "A comprehensive toolkit for SIESTA DFT simulations",
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

def get_fermi_robust(filepath):
    """
    Extracts Fermi Energy from SIESTA output file (.out).
    Strategies:
    1. Final Summary ("siesta: Fermi = ...")
    2. Classic Format ("Fermi energy: ...")
    3. SCF Table (column 'Ef' or 'Ef(eV)')
    """
    E_f = None
    scf_ef_idx = None
    
    try:
        with open(filepath, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if not parts: continue
                
                # Strategy 1: Explicit Summary (SIESTA 5.x)
                if "Fermi" in parts and "=" in parts:
                    try:
                        eq_idx = parts.index("=")
                        if eq_idx + 1 < len(parts):
                            E_f = float(parts[eq_idx + 1])
                    except ValueError:
                        pass
                
                # Strategy 2: Classic Format
                elif "Fermi" in line and "energy:" in line:
                    try:
                        colon_idx = -1
                        for i, p in enumerate(parts):
                            if p.endswith(':'):
                                colon_idx = i
                        if colon_idx != -1 and colon_idx + 1 < len(parts):
                            val = parts[colon_idx + 1].replace('eV','')
                            E_f = float(val)
                    except ValueError:
                        pass

                # Strategy 3: SCF Table
                if "iscf" in line and ("Ef" in line or "Ef(eV)" in line):
                    for i, p in enumerate(parts):
                        if "Ef" in p:
                            scf_ef_idx = i
                            break
                
                elif scf_ef_idx is not None and parts[0] == 'scf:':
                    try:
                        data_idx = scf_ef_idx + 1
                        if data_idx < len(parts):
                            E_f = float(parts[data_idx])
                    except ValueError:
                        pass
        return E_f
    except Exception as e:
        print(f"{COLORS['red']}[ERROR] Parsing file {filepath}: {e}{COLORS['reset']}")
        return None

def read_grid_data(grid_file):
    """Reads the grid using sisl and converts to eV."""
    try:
        grid = sisl.get_sile(grid_file).read_grid()
        potential_ev = grid.grid * Ry2eV
        return grid, potential_ev
    except Exception as e:
        print(f"{COLORS['red']}[ERROR] Reading grid: {e}{COLORS['reset']}")
        sys.exit(1)

def calculate_planar_avg(grid, potential_ev, axis):
    """Calculates the planar average along the specified axis."""
    if axis == 2:   # Z (Normal)
        avg_axes = (0, 1)
        L = grid.lattice.length[2]
        N = grid.shape[2]
    elif axis == 1: # Y
        avg_axes = (0, 2)
        L = grid.lattice.length[1]
        N = grid.shape[1]
    elif axis == 0: # X
        avg_axes = (1, 2)
        L = grid.lattice.length[0]
        N = grid.shape[0]

    v_planar = np.mean(potential_ev, axis=avg_axes)
    z_vals = np.linspace(0, L, N)
    return z_vals, v_planar

def find_vacuum_level(v_planar):
    """Heuristic to find vacuum level by locating the flattest region."""
    grad = np.abs(np.gradient(v_planar))
    n_points = max(1, int(len(grad)*0.15)) # 15% of points
    flat_indices = np.argsort(grad)[:n_points]
    
    E_vac = np.mean(v_planar[flat_indices])
    vac_std = np.std(v_planar[flat_indices])
    
    return E_vac, vac_std

def write_gnuplot_wf(z_vals, v_planar, E_f, E_vac, label):
    """Generates a .gplot file and a data file for plotting."""
    data_filename = "workfunction_data.dat"
    
    # Write data file
    with open(data_filename, 'w') as f:
        f.write("# Distance(Ang)  Potential(eV)\n")
        for z, v in zip(z_vals, v_planar):
            f.write(f"{z:.6f}     {v:.6f}\n")

    # Write Gnuplot script
    fileout = []
    fileout.append('# Set terminal and output\n')
    fileout.append('set terminal pdfcairo enhanced font "Arial,20" size 8,6\n')
    fileout.append(f'set output "{label}_WF.pdf"\n')
    fileout.append('\n')
    fileout.append('# Labels and Styles\n')
    fileout.append('set ylabel "Potential (eV)" font "Arial,22"\n')
    fileout.append('set xlabel "Position (Angstrom)" font "Arial,22"\n')
    fileout.append('set grid xtics ytics lt 0 lw 1 lc rgb "#bbbbbb"\n')
    fileout.append(f'set title "Work Function: {label}"\n')
    fileout.append('\n')
    
    # Arrows and Labels for WF
    mid_z = z_vals[-1] / 2
    wf_val = E_vac - E_f
    
    fileout.append(f'# Fermi and Vacuum Levels\n')
    fileout.append(f'set label "Ef = {E_f:.2f} eV" at graph 0.05, graph 0.1 textcolor rgb "forest-green"\n')
    fileout.append(f'set label "Evac = {E_vac:.2f} eV" at graph 0.05, graph 0.15 textcolor rgb "red"\n')
    fileout.append(f'set label "WF = {wf_val:.3f} eV" at {mid_z}, {(E_f+E_vac)/2} center font ",14" \n')
    
    fileout.append(f'set arrow from graph 0, first {E_f} to graph 1, first {E_f} nohead dt 2 lc rgb "forest-green" lw 2\n')
    fileout.append(f'set arrow from graph 0, first {E_vac} to graph 1, first {E_vac} nohead dt 2 lc rgb "red" lw 2\n')
    fileout.append(f'set arrow from {mid_z}, {E_f} to {mid_z}, {E_vac} heads lc rgb "black" lw 2\n')
    
    fileout.append('\n')
    fileout.append(f'plot "{data_filename}" using 1:2 with lines lw 3 lc rgb "navy" title "Planar Avg"\n')

    with open('workfunction.gplot', 'w') as file:
        file.writelines(fileout)
    
    print(f"[INFO] Gnuplot script saved to 'workfunction.gplot'")
    print(f"[INFO] Data saved to '{data_filename}'")

def plot_matplotlib(z_vals, v_planar, E_f, E_vac, label, wf_val, axis):
    """Generates the Matplotlib preview."""
    try:
        plt.figure(figsize=(8, 6))
        plt.plot(z_vals, v_planar, label='Planar Avg Potential', color='navy', linewidth=2)
        plt.axhline(y=E_vac, color='firebrick', linestyle='--', label=f'$E_{{vac}}$ = {E_vac:.2f}')
        plt.axhline(y=E_f, color='forestgreen', linestyle='--', label=f'$E_F$ = {E_f:.2f}')
        
        mid_z = z_vals[-1] / 2
        
        # Arrow and Text
        plt.annotate('', xy=(mid_z, E_f), xytext=(mid_z, E_vac),
                     arrowprops=dict(arrowstyle='<->', color='black', lw=1.5))
        
        plt.text(mid_z * 1.05, (E_f + E_vac)/2, rf"$\Phi$ = {wf_val:.2f} eV", 
                 fontsize=12, verticalalignment='center', backgroundcolor='white')

        plt.xlabel(rf"Position along Axis {axis} ($\AA$)", fontsize=14)
        plt.ylabel("Potential (eV)", fontsize=14)
        plt.title(f"Work Function: {label}", fontsize=16)
        plt.legend(loc='best')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()
    except Exception as e:
        print(f"{COLORS['yellow']}[WARNING] Could not create interactive plot: {e}{COLORS['reset']}")

def main():
    parser = argparse.ArgumentParser(
        description="Calculate Work Function from SIESTA output.",
        epilog="Example usage:\n"
               "  stb_wf --label graphene --axis 2\n"
               "  stb_wf --label slab --fermi -4.5 --no-plot",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("-l", "--label", help="System Label (required).", required=True)
    parser.add_argument("-f", "--file", help="Output file (.out). Default: label.out", default=None)
    parser.add_argument("-g", "--grid", help="Potential grid file (.VT recommended). Default: label.VT", default=None)
    parser.add_argument("-z", "--axis", type=int, default=2, help="Axis normal to surface (0=x, 1=y, 2=z). Default: 2 (z)")
    
    parser.add_argument("--fermi", type=float, help="Manually force Fermi Energy (eV).", default=None)
    parser.add_argument("--no-plot", action="store_true", help="Disable automatic plotting.")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction.")
    
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")

    args = parser.parse_args()

    # Filenames
    out_file = args.file if args.file else f"{args.label}.out"
    grid_file = args.grid if args.grid else f"{args.label}.VT" 

    if args.intro:
        show_intro()

    print("\n" + color_text("WORK FUNCTION CALCULATOR:", 'bold'))
    print("-" * 60)

    # --- 1. Fermi Energy ---
    print(f"\n[INFO] Detecting Fermi Energy...")
    E_f = args.fermi 
    
    if E_f is None:
        if os.path.exists(out_file):
            E_f = get_fermi_robust(out_file)
        else:
            print(f"{COLORS['yellow']}[WARNING] File {out_file} not found. Cannot auto-detect Fermi.{COLORS['reset']}")

    if E_f is None:
        print("\n" + color_text("FATAL ERROR: Could not determine Fermi Energy.", 'red'))
        print("Please run with manual value found in your logs:")
        print(f"  python {sys.argv[0]} -l {args.label} --fermi -3.635")
        sys.exit(1)
        
    print(f"[INFO] Fermi Energy (Ef): {E_f:.6f} eV")

    # --- 2. Grid Reading ---
    if not os.path.exists(grid_file):
        print(f"{COLORS['red']}[ERROR] Grid file '{grid_file}' not found.{COLORS['reset']}")
        sys.exit(1)

    print(f"[INFO] Reading potential from {grid_file}...")
    grid, potential_ev = read_grid_data(grid_file)

    # --- 3. Planar Average ---
    print(f"[INFO] Calculating planar average along axis {args.axis}...")
    z_vals, v_planar = calculate_planar_avg(grid, potential_ev, args.axis)

    # --- 4. Vacuum Level & Work Function ---
    E_vac, vac_std = find_vacuum_level(v_planar)
    WF = E_vac - E_f

    # --- 5. Reporting ---
    print("-" * 40)
    print(color_text(f"RESULTS for {args.label}:", 'cyan'))
    print(f"  Fermi Level    = {E_f:8.4f} eV")
    print(f"  Vacuum Level   = {E_vac:8.4f} eV")
    print(color_text(f"  Work Function  = {WF:8.4f} eV", 'green'))
    print("-" * 40)

    if vac_std > 0.05:
        print(f"{COLORS['yellow']}[WARNING] Vacuum plateau is noisy (std={vac_std:.3f} eV). Results might be inaccurate.{COLORS['reset']}")

    # --- 6. Output Files & Plotting ---
    print(f"[INFO] Writing output files...")
    write_gnuplot_wf(z_vals, v_planar, E_f, E_vac, args.label)

    if not args.no_plot:
        plot_matplotlib(z_vals, v_planar, E_f, E_vac, args.label, WF, args.axis)

    print("\n[INFO] Complete job!") 
    print("\n"+"-"*60)
    print(color_text("Work Function calculated. Now I need a vacation.\n\n", 'bold'))

if __name__ == "__main__":
    main()
