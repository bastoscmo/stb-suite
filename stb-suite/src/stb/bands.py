#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.8.0"

import os
import sys
import warnings
import subprocess
from time import sleep
import argparse
import textwrap
from typing import List, Dict
import numpy as np
import re
import argparse
import matplotlib.pyplot as plt


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

def plot_gnuplot(high_sym):
######################### PDF Plot
    #verify gamma symbol and change for symbol
    for i in range(len(high_sym)):
        clean_text = re.sub(r'[^a-zA-Z\s]', '', high_sym[i][1]).lower()
        words = clean_text.split()
        if "gamma" in words:
            high_sym[i][1]="'{/Symbol G}'"
    # Plot
    fileout=[]
    fileout.append('# Set terminal and output\n')
    fileout.append('set terminal pdfcairo enhanced font "Arial,25" size 8,8\n')
    fileout.append('set output "energy_bands.pdf"\n')
    fileout.append('\n')
    fileout.append('# Title and axis labels\n')
    fileout.append('#set title "Advanced Data Visualization Example" font ",16"\n')
    fileout.append('set ylabel "Energy (eV)" font "Arial,28" #offset -0.1,0\n')
    fileout.append('set xlabel "" font "Arial,28"\n')
    fileout.append('\n')
    fileout.append('# Ranges and scales\n')
    fileout.append('#set xrange [0:20]\n')
    fileout.append('set yrange [-20:20]\n')
    fileout.append('#set logscale y 10  # Set log scale for y-axis\n')
    fileout.append('\n')
    fileout.append('# Grid settings\n')
    fileout.append('#set grid xtics ytics mxtics mytics\n')
    fileout.append('#set mxtics 5  # Minor grid on x-axis\n')
    fileout.append('#set mytics 5  # Minor grid on y-axis\n')
    fileout.append('\n')
    fileout.append('# Tics and formatting\n')
    xticsl="set xtics ("
    for i in range(len(high_sym)):
        xticsl=xticsl+f' "{high_sym[i][1]}" {float(high_sym[i][0])} , '
    xticsl=xticsl[:-2]
    xticsl=xticsl+')  font "Arial,28"\n'
    fileout.append(xticsl)
    fileout.append('set ytics format "%.1f" font "Arial,28" \n')
    fileout.append('set grid xtics front')
    fileout.append('#set no key \n')
    fileout.append('\n')
    fileout.append('# Annotations\n')
    fileout.append('#set label "Critical point" at 10, 100 front font ",12" textcolor rgb "red"\n')
    fileout.append(f'set arrow from {high_sym[0][1]},0.0 to {high_sym[-1][1]} ,0.0 nohead dt 2 lc rgb "dark-gray" lw 4 back\n')
    for i in range(len(high_sym)-2):
        fileout.append(f'set arrow from {float(high_sym[i+1][0])} ,graph 0 to {float(high_sym[i+1][0])},graph 1 nohead dt 2 lc rgb "dark-gray" lw 4 back\n')
    fileout.append('\n')
    fileout.append('# Color palette and colorbar (useful for pm3d plots)\n')
    fileout.append('#set palette defined (0 "blue", 1 "green", 2 "yellow", 3 "red")\n')
    fileout.append('#set colorbox vertical user origin 0.9, 0.2 size 0.03, 0.6\n')
    fileout.append('#set cblabel "Colorbar (units)" font ",12"\n')
    fileout.append('\n')
    fileout.append('# Multiplot (e.g., additional smaller plots in the same figure)\n')
    fileout.append('#set multiplot layout 1,1 title "Multiplot with Color Palette" font ",14"\n')
    fileout.append('#unset multiplot\n')
    fileout.append('\n')
    fileout.append('# Line styles\n')
    fileout.append('set style line 1 lc rgb "#1f77b4" lt 1 lw 3 pt 7 ps 1.5   # Solid blue line with points\n')
    fileout.append('#set style line 2 lc rgb "#ff7f0e" lt 2 lw 2 pt 5 ps 1.5   # Dashed orange line\n')
    fileout.append('#set style line 3 lc rgb "#2ca02c" lt 3 lw 2               # Solid green line\n')
    fileout.append('\n')
    fileout.append('# Plot commands\n')
    fileout.append(f'plot "bands_gnuplot.dat" using 1:2 with lines ls 1 title "" \n')
    with open('bands.gplot', 'w') as file:
        file.writelines(fileout)
    return

def read_data(file_path = "siesta.bands"):
    # Init vectors and dictionary
    high_sym=[]
    dic_bands={}
    # Read file
    with open(file_path, "r") as f:
        fermi_energy = float(f.readlines()[0].split()[0])

    with open(file_path, "r") as f:
        lines = f.readlines()[2:]
        for i, line in enumerate(lines):
            if i >1:
                if line.startswith("     "):  #Identify blank space
                    dic_bands[key].append(line.split())
                elif len(line.split())<=2:
                    high_sym.append(line.split())
                else:
                    key=float(line.split()[0])
                    dic_bands.setdefault(key,[])
                    dic_bands[key].append(line.split()[1:])
    # convert the dictionary to float
    for key in dic_bands:
        dic_bands[key]=np.array([elem for sublista in dic_bands[key] for elem in sublista],dtype=float)
    return fermi_energy,high_sym, dic_bands

def write_gnuplot_bands(dic_bands):
    # define initial key
    key_init = next(iter(dic_bands))
    # write the file
    file_name="bands_gnuplot.dat"
    with open(file_name, 'w') as file:
        for i in range(len(dic_bands[key_init])):
            for key in dic_bands.keys():
                file.write(f"{key}     {dic_bands[key][i]}\n")
            file.write("\n")
    return

def cbm_vbm(fermi_energy,high_sym,dic_bands):
    # Start value as infinity
    vbm = -np.inf
    cbm = np.inf
    below_fermi=[]
    above_fermi=[]
    for band in dic_bands.values():
        below_fermi = np.append(below_fermi,band[band <= fermi_energy])
        above_fermi = np.append(above_fermi,band[band > fermi_energy])
    if len(below_fermi) > 0:
        vbm = max(vbm, np.nanmax(below_fermi))
    if len(above_fermi) > 0:
        cbm = min(cbm, np.nanmin(above_fermi))
    band_gap = cbm - vbm if cbm > vbm else 0.0  # Avoid negative values
    print(f"[INFO] Fermi: {fermi_energy} \n[INFO] VBM: {vbm:.6f} \n[INFO] CBM: {cbm:.6f}\n[INFO] Band Gap in lines: {band_gap:.6f}")
    return vbm,cbm

def shift_bands(dic, val):
    return {k: [v - val for v in list] for k, list in dic.items()}

def plot(dic,custom_ticks):
    # Organize the data
    x_values = sorted(dic.keys())
    num_lines = len(next(iter(dic.values())))
    y_series = [[] for _ in range(num_lines)]
    for x in x_values:
        for i in range(num_lines):
            y_series[i].append(dic[x][i])
    # Organize the high symmetries points
    tick_positions = [float(t[0]) for t in custom_ticks]
    tick_labels = [t[1] for t in custom_ticks]
    # plot
    plt.figure(figsize=(8, 6))
    for i, y_vals in enumerate(y_series):
        plt.xticks(tick_positions, tick_labels)
        plt.plot(x_values, y_vals,color='blue')
    # vertical lines in High Symmetries
    for pos in tick_positions:
        plt.axvline(x=pos, color='gray', linestyle='--', linewidth=1)
    # plot limits
    plt.ylim(-20, 20)
    plt.ylabel("Energy")
    plt.grid(True)
    plt.show()


def main():
    parser = argparse.ArgumentParser(
        description="Process the band structure data.",
        epilog="Example usage:\n"
               "  stb_bands --file siesta.bands --shift fermi\n"
               "  stb_bands --file siesta.bands --shift manual --manual-value 0.5",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("--file",  dest="input_file", type=str, required=True,
                        help="Path to the input file containing band structure data (e.g., siesta.bands).")

    parser.add_argument("--shift", type=str, choices=["vbm", "cbm", "fermi", "manual"], required=True,
                        help="Reference energy shift:\n"
                             "  - 'vbm'    : Valence Band Maximum\n"
                             "  - 'cbm'    : Conduction Band Minimum\n"
                             "  - 'fermi'  : Fermi level\n"
                             "  - 'manual' : Custom shift value (requires --manual-value).")

    parser.add_argument("--manual-value", type=float,
                        help="Custom energy shift value (required if --shift manual is used).")

    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()


    # verify the manual value
    if args.shift == "manual" and args.manual_value is None:
        parser.error("--manual-value is required when --shift is set to 'manual'.")
    
    if args.intro == True:
        show_intro()

    print("\n" + color_text("BANDS:", 'bold'))
    print("-"*60)

    # Condition to shift the band structure
    print("\n[INFO] Read file ...")
    fermi_energy,high_sym,dic_bands=read_data(args.input_file)
    print("[INFO] Calculate VBM and CBM ...")
    vbm,cbm=cbm_vbm(fermi_energy,high_sym,dic_bands)

    if args.shift == "vbm":
        rshift = vbm
    elif args.shift == "cbm":
        rshift = cbm
    elif args.shift == "fermi":
        rshift = fermi_energy
    elif args.shift == "manual":
        rshift = args.manual_value
    print("[INFO] Write files...")
    print("[WARNING] \n")
    
    write_gnuplot_bands(shift_bands(dic_bands,rshift))
    plot(shift_bands(dic_bands,rshift),high_sym)
    plot_gnuplot(high_sym)
    
    print("\n[INFO] Complete job!") 
    print("\n"+"-"*60)
    print(color_text("Bands found! But still no sign of Metallica.\n\n", 'bold'))

if __name__ == "__main__":
    main()
