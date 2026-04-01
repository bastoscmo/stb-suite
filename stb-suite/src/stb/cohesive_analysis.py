#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################
    
VERSION = "1.3.1"    

import os
import sys
import argparse
from time import sleep

# ANSI colors for terminal
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
    """Returns ANSI formatted text"""
    return f"{COLORS[color]}{text}{COLORS['reset']}"
                          
def show_intro() -> None: 
    """Displays the stylized STB-SUITE intro"""
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
        "Cohesive Energy Post-Processing",
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

def get_atom_counts(fdf_path):
    """Parses structure.fdf to count the number of atoms for each chemical species"""
    species_map = {}
    counts = {}
    
    in_species_block = False
    in_coords_block = False
    
    try:
        with open(fdf_path, 'r') as f:
            for line in f:
                cleaned_line = line.split('#', 1)[0].strip()
                if not cleaned_line: continue
                
                lower_line = cleaned_line.lower()
                
                # Block toggles
                if lower_line == '%block chemicalspecieslabel':
                    in_species_block = True; continue
                if lower_line == '%endblock chemicalspecieslabel':
                    in_species_block = False; continue
                if lower_line == '%block atomiccoordinatesandatomicspecies':
                    in_coords_block = True; continue
                if lower_line == '%endblock atomiccoordinatesandatomicspecies':
                    in_coords_block = False; continue

                # Parse species mapping (ID -> Symbol)
                if in_species_block:
                    parts = cleaned_line.split()
                    if len(parts) >= 3:
                        species_map[parts[0]] = parts[2]
                        counts[parts[2]] = 0 # Initialize count

                # Parse coordinates to count atoms
                if in_coords_block:
                    parts = cleaned_line.split()
                    if len(parts) >= 4:
                        specie_id = parts[3]
                        if specie_id in species_map:
                            symbol = species_map[specie_id]
                            counts[symbol] += 1
                            
        return counts
    except FileNotFoundError:
        print(color_text(f"[ERROR] Structure file '{fdf_path}' not found.", 'red'))
        sys.exit(1)

def extract_final_energy(folder_path, out_filename):
    """Finds the specific output file in the folder and extracts the final FreeEng"""
    target_file = os.path.join(folder_path, out_filename)
    
    if not os.path.exists(target_file):
        return None
    
    energy = None
    
    with open(target_file, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            # Search ignoring exact spaces
            if "siesta: FreeEng" in line and "=" in line:
                try:
                    # Get everything after '=' and extract the first number
                    val_str = line.split('=')[1].split()[0]
                    energy = float(val_str)
                except (IndexError, ValueError):
                    pass
        
    return energy

def main():
    parser = argparse.ArgumentParser(
        description="Process cohesive energy results from SIESTA calculations.",
        epilog="Example usage:\n"
               "  stb_cohesive_analysis -o calc.out\n"
               "  stb_cohesive_analysis -o calc.out -d /path/to/results",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("-o", "--out", dest="out_file", type=str, required=True, 
                        help="Name of the SIESTA output file (e.g., calc.out)")
    
    parser.add_argument("-d", "--dir", dest="dir_path", type=str, default=".", 
                        help="Path to the folder containing 'structure' and 'atoms' directories (default: current directory)")
    
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")

    parser.add_argument("--no-intro", dest="intro", action="store_false", 
                        help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro == True:
        show_intro()

    print("\n" + color_text("COHESIVE ANALYSIS:", 'bold'))
    print("-" * 60)
    
    root_dir = args.dir_path
    struct_dir = os.path.join(root_dir, "structure")
    atoms_dir = os.path.join(root_dir, "atoms")
    
    print("\n[INFO] Checking directories ...")
    if not os.path.exists(struct_dir) or not os.path.exists(atoms_dir):
        print(color_text(f"[ERROR] Required directories 'structure' and/or 'atoms' not found in '{root_dir}'.", 'red'))
        sys.exit(1)

    # 1. Parse Structure for Atom Counts
    struct_fdf = os.path.join(struct_dir, "structure.fdf")
    print("[INFO] Read structure file ...")
    atom_counts = get_atom_counts(struct_fdf)
    
    if not atom_counts:
        print(color_text("[ERROR] Could not extract atom counts. Check your structure.fdf", 'red'))
        sys.exit(1)

    total_atoms = sum(atom_counts.values())
    print(f"[INFO] Total atoms in cell: {total_atoms}")
    for sym, count in atom_counts.items():
        print(f"       -> {sym}: {count} atoms")

    # 2. Extract Energies
    print("\n[INFO] Extracting energies ...")
    errors = False
    
    # Bulk Energy
    e_bulk = extract_final_energy(struct_dir, args.out_file)
    if e_bulk is None:
        print(color_text(f"[WARNING] Could not find '{args.out_file}' or finished calculation for the full structure.", 'yellow'))
        errors = True
    else:
        print(f"[INFO] Energy (Full Structure): {e_bulk:12.6f} eV")

    # Isolated Atoms Energy
    isolated_energies = {}
    sum_isolated_energy = 0.0
    
    for sym, count in atom_counts.items():
        sym_dir = os.path.join(atoms_dir, sym)
        e_atom = extract_final_energy(sym_dir, args.out_file)
        
        if e_atom is None:
            print(color_text(f"[WARNING] Could not find '{args.out_file}' or results for isolated atom: {sym}", 'yellow'))
            errors = True
        else:
            isolated_energies[sym] = e_atom
            sum_isolated_energy += (e_atom * count)
            print(f"[INFO] Energy (Isolated {sym}):   {e_atom:12.6f} eV")

    # 3. Calculate Cohesive Energy
    if errors:
        print(color_text("\n[ERROR] Cannot calculate cohesive energy because some calculations are missing or incomplete.", 'red'))
        sys.exit(1)

	# Cohesive Energy Formula: E_coh = (E_bulk - Sum(E_iso)) / N_atoms
	# Negative value means a bound, stable structure.
    print("\n[INFO] Calculating Cohesive Energy ...")
    e_coh_total =  e_bulk - sum_isolated_energy
    e_coh_per_atom = e_coh_total / total_atoms

    print("\n[INFO] FINAL RESULTS:")
    print(f"       Sum of Isolated Atoms: {sum_isolated_energy:12.6f} eV")
    print(f"       Bulk Structure Energy: {e_bulk:12.6f} eV")
    print(f"       Total Cohesive Energy:     {e_coh_total:10.4f} eV")
    print(f"       Cohesive Energy per Atom:  {e_coh_per_atom:10.4f} eV/atom")

    # 4. Save results to a .dat file
    print("\n[INFO] Write files ...")
    out_file_path = os.path.join(root_dir, "cohesive_results.dat")
    try:
        with open(out_file_path, 'w') as f_out:
            f_out.write("==================================================\n")
            f_out.write("             COHESIVE ENERGY RESULTS              \n")
            f_out.write("==================================================\n\n")
            
            f_out.write(f"Total atoms in cell: {total_atoms}\n")
            for sym, count in atom_counts.items():
                f_out.write(f"  -> {sym}: {count} atoms\n")
            f_out.write("-" * 50 + "\n")
            
            f_out.write(f"Energy (Full Structure): {e_bulk:16.6f} eV\n")
            for sym, count in atom_counts.items():
                f_out.write(f"Energy (Isolated {sym}):   {isolated_energies[sym]:16.6f} eV\n")
            f_out.write("-" * 50 + "\n")
            
            f_out.write(f"Sum of Isolated Atoms:   {sum_isolated_energy:16.6f} eV\n\n")
            f_out.write(f"Total Cohesive Energy:       {e_coh_total:12.4f} eV\n")
            f_out.write(f"Cohesive Energy per Atom:    {e_coh_per_atom:12.4f} eV/atom\n")
            f_out.write("\n==================================================\n")
            
        print(f"[INFO] Results saved to: {out_file_path}")
    except Exception as e:
        print(color_text(f"[ERROR] Failed to save results to file: {e}", 'red'))

    print("\n[INFO] Complete job!") 
    print("\n"+"-"*60)
    print(color_text("Cohesive energy analysis complete!\n\n", 'bold'))

if __name__ == "__main__":
    main()
