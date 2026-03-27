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
warnings.filterwarnings("ignore", category=UserWarning, module="pybader")
warnings.filterwarnings("ignore", category=DeprecationWarning)

# --- Library Imports ---
#try:
import sisl
#except ImportError:
#    print("\033[91m[CRITICAL] Library 'sisl' not found.\033[0m")
#    sys.exit(1)

#try:
from pybader.interface import Bader as PyBaderCalc
#except ImportError:
#    try:
import pybader
PyBaderCalc = pybader.interface.Bader
#    except Exception:
#        print("\033[91m[CRITICAL] Could not initialize PyBader interface. Install setuptools.\033[0m")
#        sys.exit(1)


# ================= ANSI COLORS =================
COLORS = {
    'reset': '\033[0m', 'cyan': '\033[96m', 'blue': '\033[94m',
    'green': '\033[92m', 'yellow': '\033[93m', 'red': '\033[91m',
    'bold': '\033[1m', 'underline': '\033[4m',
    'bg_red': '\033[41m', 'white': '\033[97m' # High visibility
}

def color_text(text: str, color: str) -> str:
    return f"{COLORS[color]}{text}{COLORS['reset']}"

def show_intro() -> None:
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
    print("\n" + "="*60)
    print("Siesta ToolBox Suite - Bader Analysis".center(60))
    print(f"Version {VERSION} | University of Brasilia".center(60))
    print("="*60 + "\n")

# ================= SCIENTIFIC DATA =================

# FALLBACK DICTIONARY
# Contains standard chemical valence (s+p for main group, s+d for transition).
# NOTE: DFT pseudopotentials (especially with semi-core states) may differ.
#       This is used ONLY if the .out file cannot be parsed.
FALLBACK_VALENCE = {
    # Period 1
    "H": 1.0,  "He": 2.0,
    # Period 2
    "Li": 1.0, "Be": 2.0, "B": 3.0,  "C": 4.0,  "N": 5.0,  "O": 6.0,  "F": 7.0,  "Ne": 8.0,
    # Period 3
    "Na": 1.0, "Mg": 2.0, "Al": 3.0, "Si": 4.0, "P": 5.0,  "S": 6.0,  "Cl": 7.0, "Ar": 8.0,
    # Period 4
    "K": 1.0,  "Ca": 2.0, "Sc": 3.0, "Ti": 4.0, "V": 5.0,  "Cr": 6.0, "Mn": 7.0, "Fe": 8.0,
    "Co": 9.0, "Ni": 10.0,"Cu": 11.0,"Zn": 12.0,"Ga": 3.0, "Ge": 4.0, "As": 5.0, "Se": 6.0,
    "Br": 7.0, "Kr": 8.0,
    # Period 5
    "Rb": 1.0, "Sr": 2.0, "Y": 3.0,  "Zr": 4.0, "Nb": 5.0, "Mo": 6.0, "Tc": 7.0, "Ru": 8.0,
    "Rh": 9.0, "Pd": 10.0,"Ag": 11.0,"Cd": 12.0,"In": 3.0, "Sn": 4.0, "Sb": 5.0, "Te": 6.0,
    "I": 7.0,  "Xe": 8.0,
    # Period 6 (Lanthanides usually 3, but can vary in DFT)
    "Cs": 1.0, "Ba": 2.0,
    "La": 3.0, "Ce": 4.0, "Pr": 3.0, "Nd": 3.0, "Pm": 3.0, "Sm": 3.0, "Eu": 2.0, "Gd": 3.0,
    "Tb": 3.0, "Dy": 3.0, "Ho": 3.0, "Er": 3.0, "Tm": 3.0, "Yb": 2.0, "Lu": 3.0,
    "Hf": 4.0, "Ta": 5.0, "W": 6.0,  "Re": 7.0, "Os": 8.0, "Ir": 9.0, "Pt": 10.0,"Au": 11.0,
    "Hg": 12.0,"Tl": 3.0, "Pb": 4.0, "Bi": 5.0, "Po": 6.0, "At": 7.0, "Rn": 8.0,
    # Period 7 (Actinides)
    "Fr": 1.0, "Ra": 2.0,
    "Ac": 3.0, "Th": 4.0, "Pa": 5.0, "U": 6.0,  "Np": 5.0, "Pu": 4.0, "Am": 3.0, "Cm": 3.0,
    "Bk": 3.0, "Cf": 3.0, "Es": 3.0, "Fm": 3.0, "Md": 3.0, "No": 2.0, "Lr": 3.0,
    "Rf": 4.0, "Db": 5.0, "Sg": 6.0, "Bh": 7.0, "Hs": 8.0, "Mt": 9.0, "Ds": 10.0,"Rg": 11.0,
    "Cn": 12.0,"Nh": 3.0, "Fl": 4.0, "Mc": 5.0, "Lv": 6.0, "Ts": 7.0, "Og": 8.0
}

def get_zval_from_output(label, override_path=None):
    """
    Reads Z_val from the .out file. Supports Siesta 4.x and 5.x formats.
    """
    
    if override_path:
        target_file = override_path
        if not os.path.exists(target_file):
            print(color_text(f"   [ERROR] Reference file '{target_file}' (via --ref) not found.", 'red'))
            return None
        print(f"   [INFO] Using external reference file: {color_text(target_file, 'bold')}")
    else:
        target_file = f"{label}.out"
        if not os.path.exists(target_file):
            return None

    dynamic_valence = {}
    current_label = None
    
    try:
        with open(target_file, 'r') as f:
            for line in f:
                # ---------------- LABEL DETECTION ----------------
                
                # PRIORITY 1: Siesta 5 Standard ("atom: Called for C")
                # This is the safest method as it appears right before the processing block
                if "atom: Called for" in line:
                    parts = line.split()
                    try:
                        # Find "for" and get the next element
                        idx = parts.index("for")
                        current_label = parts[idx+1]
                        # Clean cleanup (e.g. "C(Z=6)" -> "C")
                        current_label = current_label.split('(')[0]
                    except ValueError:
                        pass
                
                # PRIORITY 2: Old Standard ("Species number: ... Label: C")
                elif "Species number:" in line and "Label:" in line:
                    parts = line.split()
                    try:
                        lbl_candidate = parts[-1]
                        current_label = lbl_candidate
                    except IndexError:
                        continue

                # ---------------- ZVAL DETECTION ----------------
                
                # Check if we have an active Label
                if current_label:
                    zval = None
                    
                    # CASE 1: Standard Vna ("Vna: chval, zval: ...")
                    if "Vna: chval, zval:" in line:
                        parts = line.split()
                        try:
                            zval = float(parts[-1])
                        except ValueError: pass

                    # CASE 2: PSML/Pseudopotential Generation
                    # "Valence charge in pseudo generation:    4.00000"
                    elif "Valence charge in pseudo generation:" in line:
                        parts = line.split()
                        try:
                            zval = float(parts[-1])
                        except ValueError: pass
                    
                    if zval is not None:
                        dynamic_valence[current_label] = zval
                        
    except Exception as e:
        print(color_text(f"[WARN] Error reading {target_file}: {e}", 'yellow'))
        return None
        
    return dynamic_valence if dynamic_valence else None

# ================= MAIN LOGIC =================

def solve_bader(label, output_file=None, speed_mode='normal', ref_file=None):
    
    file_rho = f"{label}.RHO"
    file_xv = f"{label}.XV"
    file_fdf = f"{label}.fdf"
    file_cube = f"{label}.cube"
    
    if not output_file:
        output_file = f"{label}_BADER.txt"

    if not os.path.exists(file_rho):
        print(color_text(f"[ERROR] Grid file '{file_rho}' not found.", 'red'))
        return

    print(f"[INFO] System: {color_text(label, 'bold')} | Mode: {color_text(speed_mode.upper(), 'cyan')}")

    # --- Z_val Setup ---
    print(f"0. [Setup] Configuring valence charges...")
    
    detected_valence = get_zval_from_output(label, override_path=ref_file)
    
    if detected_valence:
        print(f"   {color_text('[SUCCESS]', 'green')} Z_vals detected: {detected_valence}")
        valence_source = detected_valence
        source_name = f"Siesta Output ({ref_file if ref_file else label+'.out'})"
    else:
        # --- HIGH VISIBILITY WARNING BLOCK ---
        warn_msg = [
            "\n" + "!" * 65,
            "[WARNING] .out FILE NOT FOUND OR UNREADABLE! USING DEFAULT VALUES.",
            "Please ensure your pseudopotential Z_val matches standard values.",
            "Use --ref to point to a valid .out file.",
            "!" * 65 + "\n"
        ]
        # Print with RED Background and WHITE text
        for line in warn_msg:
            print(f"{COLORS['bg_red']}{COLORS['white']}{COLORS['bold']}{line.center(65)}{COLORS['reset']}")

        valence_source = FALLBACK_VALENCE
        source_name = "Hardcoded Defaults (FALLBACK)"

    # --- STEP 1: SISL ---
    try:
        print(f"1. [SISL] Reading geometry and charge density...")
        
        if os.path.exists(file_xv):
            geometry = sisl.get_sile(file_xv).read_geometry()
        elif os.path.exists(file_fdf):
            geometry = sisl.get_sile(file_fdf).read_geometry()
        else:
            print(color_text("[ERROR] No geometry file (.XV or .fdf) found.", 'red'))
            return

        rho_grid = sisl.get_sile(file_rho).read_grid()
        rho_grid.set_geometry(geometry)
        rho_grid.write(file_cube)
        
    except Exception as e:
        print(color_text(f"[ERROR] SISL processing failed: {e}", 'red'))
        return

    # --- STEP 2: PyBader ---
    n_threads = multiprocessing.cpu_count()
    print(f"2. [PyBader] Starting calculation on {n_threads} threads...")
    
    try:
        bader_job = PyBaderCalc.from_file(file_cube, threads=n_threads)
        if speed_mode == 'fast':
            if hasattr(bader_job, 'refinement_method'):
                bader_job.refinement_method = 'minimum_distance' 
        bader_job()
        raw_populations = bader_job.atoms_charge
    except Exception as e:
        print(color_text(f"[ERROR] PyBader calculation failed: {e}", 'red'))
        return

    # --- STEP 3: Analysis ---
    print("3. [Analysis] compiling results...")
    
    total_theory = 0.0
    atoms_data = []

    for i, atom in enumerate(geometry.atoms):
        if i >= len(raw_populations): break
        sym = atom.symbol
        
        z_val = valence_source.get(sym, 0.0)
        
        if z_val == 0.0 and sym not in valence_source:
             print(color_text(f"   [WARN] Element {sym} not found in Z_val dictionary!", 'red'))

        total_theory += z_val
        atoms_data.append({'id': i+1, 'sym': sym, 'z_val': z_val, 'pop_raw': raw_populations[i]})

    # Unit Correction
    total_raw = sum(raw_populations)
    ratio = total_raw / total_theory if total_theory > 0 else 1.0
    correction_factor = 1.0 / ratio if abs(ratio - 1.0) > 0.1 else 1.0
    
    if correction_factor != 1.0:
        print(color_text(f"   [INFO] Unit mismatch corrected. Factor: {correction_factor:.4f}", 'cyan'))

    # --- STEP 4: Output ---
    out_lines = []
    out_lines.append(f"BADER CHARGE ANALYSIS REPORT - STB Suite v{VERSION}")
    out_lines.append(f"System: {label}")
    out_lines.append(f"Z_val Source: {source_name}")
    out_lines.append("=" * 75)
    
    header = f"{'Idx':<4} {'Elem':<5} {'Pop(e-)':<12} {'Z_val':<8} {'Net Charge':<12} {'State':<15}"
    out_lines.append(header)
    out_lines.append("-" * 75)
    
    print("\n" + header)
    print("-" * 75)

    total_final = 0.0

    for data in atoms_data:
        pop = data['pop_raw'] * correction_factor
        net = data['z_val'] - pop
        
        if net > 0.05: state, color = "Donor (+)", 'red'
        elif net < -0.05: state, color = "Acceptor (-)", 'blue'
        else: state, color = "Neutral", 'reset'

        line = f"{data['id']:<4} {data['sym']:<5} {pop:<12.4f} {data['z_val']:<8.2f} {net:<+12.4f} {state:<15}"
        out_lines.append(line)
        print(f"{data['id']:<4} {data['sym']:<5} {pop:<12.4f} {data['z_val']:<8.2f} "
              f"{color_text(f'{net:<+12.4f}', 'bold')} {color_text(state, color)}")
        
        total_final += pop

    footer = [
        "-" * 75,
        f"Total Integrated: {total_final:.4f} (Target: {total_theory:.2f})",
        "=" * 75
    ]
    out_lines.extend(footer)
    for l in footer: print(l)

    try:
        with open(output_file, "w") as f:
            f.write("\n".join(out_lines))
        print(f"\n{color_text('[OK]', 'green')} Results saved to: {color_text(output_file, 'bold')}")
    except IOError:
        print(color_text(f"[ERROR] Could not save file {output_file}", 'red'))

# ================= EXECUTION =================

def main():
    parser = argparse.ArgumentParser(description="STB Bader Tool")
    parser.add_argument("-l", "--label", required=True, help="SystemLabel used in Siesta")
    parser.add_argument("-o", "--output", required=False, help="Output filename")
    parser.add_argument("--ref", required=False, default=None, help="Path to a specific .out file (Overrides default)")
    parser.add_argument("--speed", choices=['normal', 'fast'], default='normal', help="Optimization level")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Skip intro")

    args = parser.parse_args()

    if args.intro:
        show_intro()

    solve_bader(args.label, args.output, args.speed, args.ref)
    
    print("\n" + "-" * 60)
    print(color_text("Electron counting is like accounting, but the currency is negative.\n", 'bold'))

if __name__ == "__main__":
    main()
    

