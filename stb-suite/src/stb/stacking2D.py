#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

"""
STB Monolayer Stacker (ZSL Heteroepitaxy Mode)
Reads two Siesta .fdf files, finds a commensurate supercell using the ZSL algorithm,
and stacks them into a van der Waals heterostructure.
"""

VERSION = "1.8.0"

import os
import sys
import warnings
import argparse
import textwrap
import numpy as np
from time import sleep
from pymatgen.core import Structure, Lattice
from pymatgen.core.operations import SymmOp
from pymatgen.analysis.interfaces.zsl import ZSLGenerator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

# Suppress Pymatgen warnings for cleaner CLI output
warnings.filterwarnings("ignore")

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
    """Returns text formatted with ANSI colors"""
    return f"{COLORS[color]}{text}{COLORS['reset']}"

def show_intro() -> None:
    """Displays the stylized STB-SUITE introduction"""
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
        "Siesta ToolBox Suite - Monolayer Stacker",
        "ZSL Heteroepitaxy Algorithm & Geometric Control",
        f"Version {VERSION} | University of Brasilia - 2026",
        "Developed by Dr. Carlos M. O. Bastos"
    ]

    print(logo)
    print("\n" + "="*60)
    for line in description:
        print(line.center(60))
        sleep(0.2)
    print("="*60 + "\n")
    return

def parse_fdf_to_pymatgen(filepath):
    """Parses a Siesta .fdf file and returns a Pymatgen Structure object."""
    lattice_vectors = []
    species_dict = {}
    coords = []
    atomic_species = []
    
    in_species = False
    in_lattice = False
    in_coords = False
    
    try:
        with open(filepath, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if not parts or line.startswith('#'): continue
                    
                if '%block ChemicalSpeciesLabel' in line: in_species = True; continue
                elif '%endblock ChemicalSpeciesLabel' in line: in_species = False; continue
                elif '%block LatticeVectors' in line: in_lattice = True; continue
                elif '%endblock LatticeVectors' in line: in_lattice = False; continue
                elif '%block AtomicCoordinatesAndAtomicSpecies' in line: in_coords = True; continue
                elif '%endblock AtomicCoordinatesAndAtomicSpecies' in line: in_coords = False; continue
                
                if in_species: species_dict[parts[0]] = parts[2]
                elif in_lattice: lattice_vectors.append([float(x) for x in parts[:3]])
                elif in_coords:
                    coords.append([float(x) for x in parts[:3]])
                    atomic_species.append(species_dict[parts[3]])
                    
        if not lattice_vectors or not coords:
            raise ValueError(f"Missing lattice vectors or coordinates in {filepath}.")

        lattice = Lattice(lattice_vectors)
        return Structure(lattice, atomic_species, coords, coords_are_cartesian=False)

    except Exception as e:
        print(color_text(f"[ERROR] Error parsing {filepath}: {e}", 'red'))
        sys.exit(1)

def stack_heterostructure(layer1, layer2, gaps, target_vacuum, max_area, max_strain, shift_x=0.0, shift_y=0.0, twist_angle=0.0, strain_mode='top', match_id=0, interactive=False, batch_sym=False):
    """Matches lattices, handles selection, and generates heterostructures."""
    
    # 1. APPLY MATHEMATICALLY EXACT ROTATION
    if twist_angle != 0.0:
        print(color_text(f"[INFO] Applying initial twist angle of {twist_angle}° to Layer 2...", 'cyan'))
        op = SymmOp.from_axis_angle_and_translation([0, 0, 1], twist_angle)
        layer2.apply_operation(op, fractional=False)

    print(color_text("[INFO] Searching for commensurate supercells (ZSL Algorithm)...", 'cyan'))
    zsl = ZSLGenerator(max_area=max_area, max_length_tol=max_strain, max_angle_tol=0.1) 
    
    raw_matches = list(zsl(layer1.lattice.matrix[:2], layer2.lattice.matrix[:2]))
    if not raw_matches:
        print(color_text(f"\n[ERROR] No lattice match found! Try increasing --max_area or --max_strain.", 'red'))
        sys.exit(1)

    evaluated_matches = []
    for m in raw_matches:
        t_mat1 = np.eye(3); t_mat1[:2, :2] = m.film_transformation
        t_mat2 = np.eye(3); t_mat2[:2, :2] = m.substrate_transformation
        
        t_l1 = layer1.copy(); t_l1.make_supercell(t_mat1)
        t_l2 = layer2.copy(); t_l2.make_supercell(t_mat2)
        
        # Calculate Linear Strain
        strain_a = abs(t_l1.lattice.a / t_l2.lattice.a - 1.0) * 100
        strain_b = abs(t_l1.lattice.b / t_l2.lattice.b - 1.0) * 100
        max_len_strain = max(strain_a, strain_b)
        
        # Calculate Angular Strain (Un-rotation penalty)
        vec1 = t_l1.lattice.matrix[0]
        vec2 = t_l2.lattice.matrix[0]
        cos_angle = np.clip(np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2)), -1.0, 1.0)
        angle_diff = np.degrees(np.arccos(cos_angle))
        
        evaluated_matches.append({
            'match': m,
            'strain': max_len_strain,
            'angle_strain': angle_diff,
            'area': m.match_area
        })

    evaluated_matches.sort(key=lambda x: (x['strain'], x['angle_strain']))

    if interactive:
        display_limit = 15
        print("\n" + color_text(f"--- ZSL Commensurate Supercells Found ({len(evaluated_matches)} matches) ---", 'bold'))
        header = f"{'ID':<4} | {'Area (Å²)':<10} | {'Strain (%)':<10} | {'Ang. Strain (°)':<15} | {'Matrix L1'}"
        print(color_text(header, 'blue'))
        print("-" * 75)
        
        for i, data in enumerate(evaluated_matches[:display_limit]):
            f_mat = np.round(data['match'].film_transformation).astype(int)
            mat_str = f"[{f_mat[0][0]:2d}, {f_mat[0][1]:2d}], [{f_mat[1][0]:2d}, {f_mat[1][1]:2d}]"
            row_text = f"{i:<4} | {data['area']:<10.2f} | {data['strain']:<10.2f} | {data['angle_strain']:<15.2f} | {mat_str}"
            # Highlight best matches
            if data['strain'] < 1.0 and data['angle_strain'] < 1.0:
                print(color_text(row_text, 'green'))
            else:
                print(row_text)
            
        while True:
            try:
                user_input = input(color_text(f"\nSelect a Match ID (0 to {len(evaluated_matches)-1}) or 'q' to quit: ", 'yellow'))
                if user_input.lower() == 'q': sys.exit(0)
                selected_id = int(user_input)
                if 0 <= selected_id < len(evaluated_matches):
                    match_id = selected_id
                    break
            except ValueError:
                pass
    else:
        if match_id >= len(evaluated_matches) or match_id < 0:
            print(color_text(f"[ERROR] Selected match_id {match_id} is out of range.", 'red'))
            sys.exit(1)
        
    best_match_data = evaluated_matches[match_id]
    best_match = best_match_data['match']
    
    t_mat1 = np.eye(3); t_mat1[:2, :2] = best_match.film_transformation
    t_mat2 = np.eye(3); t_mat2[:2, :2] = best_match.substrate_transformation
    
    print(color_text(f"\n[INFO] Selected Match ID {match_id} | Area: {best_match_data['area']:.2f} Å² | Angular Strain: {best_match_data['angle_strain']:.2f}°", 'green'))
    if twist_angle != 0.0 and best_match_data['angle_strain'] > 1.0:
        print(color_text(f"[WARNING] High Angular Strain! This match will 'un-rotate' your structure by ~{best_match_data['angle_strain']:.1f}° to force PBC fit.", 'yellow'))
        print(color_text(f"          Increase --max_area to find a true Moiré supercell with 0° Angular Strain.", 'yellow'))

    # Define Shifts
    shifts_to_run = {}
    if batch_sym:
        print(color_text("[INFO] Batch Symmetry Mode active. Generating high-symmetry configurations...", 'cyan'))
        if shift_x != 0.0 or shift_y != 0.0:
            print(color_text(f"  -> [WARNING] Custom shifts (-tx {shift_x}, -ty {shift_y}) are IGNORED in --batch_sym mode.", 'yellow'))
        if twist_angle != 0.0:
            print(color_text(f"  -> [WARNING] Using --batch_sym with a twist angle ({twist_angle}°). Shifting may not yield true AA/AB symmetries due to Moiré patterns.", 'yellow'))

        gamma = layer1.lattice.angles[2]
        is_hex = (abs(gamma - 120.0) < 2.0 or abs(gamma - 60.0) < 2.0) and abs(layer1.lattice.a - layer1.lattice.b) < 0.1

        shifts_to_run = {"AA": [0.0, 0.0], "Bridge_X": [0.5, 0.0], "Bridge_Y": [0.0, 0.5], "Center": [0.5, 0.5]}
        if is_hex:
            print(color_text("  -> [INFO] Hexagonal lattice detected. Adding Hex_AB and Hex_BA stackings.", 'green'))
            shifts_to_run.update({"Hex_AB": [1/3, 2/3], "Hex_BA": [2/3, 1/3]})
    else:
        shifts_to_run = {"Custom": [shift_x, shift_y]}

    results = {}
    print(color_text("[INFO] Building heterostructures...", 'cyan'))
    for name, shift in shifts_to_run.items():
        for gap in gaps:
            layer1_supercell = layer1.copy(); layer1_supercell.make_supercell(t_mat1)
            layer2_supercell = layer2.copy()
            layer2_supercell.translate_sites(range(len(layer2_supercell)), [shift[0], shift[1], 0.0])
            layer2_supercell.make_supercell(t_mat2)

            # Strain logic
            if strain_mode == 'top':
                base_lattice = layer1_supercell.lattice
            elif strain_mode == 'bottom':
                base_lattice = layer2_supercell.lattice
            else:
                base_lattice = Lattice((layer1_supercell.lattice.matrix + layer2_supercell.lattice.matrix) / 2.0)

            max_strain_val = max(abs(base_lattice.a / layer1_supercell.lattice.a - 1.0), abs(base_lattice.a / layer2_supercell.lattice.a - 1.0))

            # 2. FLAWLESS MAPPING OF ROTATED ATOMS TO THE NEW LATTICE BOX
            l1_frac = [s.frac_coords for s in layer1_supercell]
            l1_cart_strained = [base_lattice.get_cartesian_coords(f) for f in l1_frac]
            
            l2_frac = [s.frac_coords for s in layer2_supercell]
            l2_cart_strained = [base_lattice.get_cartesian_coords(f) for f in l2_frac]

            z_max_l1 = max([c[2] for c in l1_cart_strained])
            z_min_l1 = min([c[2] for c in l1_cart_strained])
            z_min_l2 = min([c[2] for c in l2_cart_strained])
            
            z_shift_cart = z_max_l1 - z_min_l2 + gap
            l2_cart_shifted = [c + np.array([0, 0, z_shift_cart]) for c in l2_cart_strained]
            
            calc_vacuum = target_vacuum if target_vacuum is not None else layer1_supercell.lattice.c - (z_max_l1 - z_min_l1)
            new_z_max = max([c[2] for c in l2_cart_shifted])
            new_c_length = (new_z_max - z_min_l1) + calc_vacuum

            new_matrix = base_lattice.matrix.copy()
            new_matrix[2] = [0, 0, new_c_length]
            new_lattice = Lattice(new_matrix)

            all_species = [s.specie for s in layer1_supercell] + [s.specie for s in layer2_supercell]
            all_coords_cart = l1_cart_strained + l2_cart_shifted

            hetero = Structure(new_lattice, all_species, all_coords_cart, coords_are_cartesian=True)

            res_name = f"{name}_gap_{gap:.2f}" if len(gaps) > 1 else name
            results[res_name] = (hetero, max_strain_val)

    return results

def analyze_and_report_symmetry(layer1, layer2, hetero, filename="symmetry_report.txt", symprec=0.01, print_stdout=True):
    """Analyzes symmetry for both layers and the heterostructure, printing a formatted table."""
    def get_details(struct):
        try:
            sga = SpacegroupAnalyzer(struct, symprec=symprec)
            return {
                "Space Group": f"{sga.get_space_group_symbol()} ({sga.get_space_group_number()})",
                "Point Group": sga.get_point_group_symbol(),
                "Crystal System": str(sga.get_crystal_system()).title(),
                "Hall Symbol": sga.get_hall(),
            }
        except Exception as e:
            return {"Error": f"Analysis failed ({e})"}

    l1_info = get_details(layer1)
    l2_info = get_details(layer2)
    het_info = get_details(hetero)

    # Define table properties
    properties = ["Crystal System", "Space Group", "Point Group", "Hall Symbol"]
    
    # Build the formatted table as a single string
    table_lines = []
    separator = "=" * 72
    thin_separator = "-" * 72
    
    table_lines.append(separator)
    table_lines.append(f"  DETAILED SYMMETRY ANALYSIS (Tolerance: {symprec} Å)")
    table_lines.append(separator)
    
    # Table Header
    header = f"{'Property':<17} | {'Layer 1':<16} | {'Layer 2':<16} | {'Heterostructure':<16}"
    table_lines.append(header)
    table_lines.append(thin_separator)
    
    # Table Rows
    if "Error" in l1_info or "Error" in l2_info or "Error" in het_info:
        table_lines.append("  [!] Error encountered during symmetry analysis.")
        table_lines.append(f"  L1: {l1_info.get('Error', 'OK')} | L2: {l2_info.get('Error', 'OK')} | Het: {het_info.get('Error', 'OK')}")
    else:
        for prop in properties:
            l1_val = l1_info.get(prop, "N/A")
            l2_val = l2_info.get(prop, "N/A")
            het_val = het_info.get(prop, "N/A")
            row = f"{prop:<17} | {str(l1_val):<16} | {str(l2_val):<16} | {str(het_val):<16}"
            table_lines.append(row)
            
    table_lines.append(separator)
    report_text = "\n".join(table_lines)

    # Print to Terminal with STB Colors
    if print_stdout: 
        print(f"\n{color_text(separator, 'blue')}")
        print(color_text(f"  DETAILED SYMMETRY ANALYSIS (Tolerance: {symprec} Å)", 'bold'))
        print(color_text(separator, 'blue'))
        print(color_text(header, 'cyan'))
        print(thin_separator)
        
        if "Error" not in table_lines[5]: # Check if error occurred
            for prop in properties:
                l1_val = l1_info.get(prop, "N/A")
                l2_val = l2_info.get(prop, "N/A")
                het_val = het_info.get(prop, "N/A")
                
                # Highlight the Heterostructure column in yellow to make it stand out
                het_colored = color_text(f"{str(het_val):<16}", 'yellow')
                row_colored = f"{prop:<17} | {str(l1_val):<16} | {str(l2_val):<16} | {het_colored}"
                print(row_colored)
        else:
            print(table_lines[5])
            print(table_lines[6])
            
        print(color_text(separator, 'blue'))
        print(color_text(f"[INFO] Symmetry report saved to: {filename}", 'green'))
    
    # Save to text file
    with open(filename, 'w') as f: 
        f.write(report_text + "\n")

def export_to_fdf(structure, filename="stacked_structure.fdf"):
    with open(filename, 'w') as f:
        f.write("# Generated by STB Monolayer Stacker\n\n")
        
        unique_species = list(set([site.specie.symbol for site in structure]))
        num_atoms = len(structure)
        
        f.write(f"NumberOfSpecies    {len(unique_species)}\n")
        f.write(f"NumberOfAtoms      {num_atoms}\n\n")
        
        f.write("%block ChemicalSpeciesLabel\n")
        from pymatgen.core.periodic_table import Element
        for i, symbol in enumerate(unique_species, start=1):
            z = Element(symbol).Z
            f.write(f" {i}   {z}   {symbol}\n")
        f.write("%endblock ChemicalSpeciesLabel\n\n")
        
        f.write("LatticeConstant 1.00 Ang\n\n")
        f.write("AtomicCoordinatesFormat  Fractional\n\n")
        
        f.write("%block LatticeVectors\n")
        for vec in structure.lattice.matrix:
            f.write(f" {vec[0]:.8f}   {vec[1]:.8f}   {vec[2]:.8f}\n")
        f.write("%endblock LatticeVectors\n\n")
        
        f.write("%block AtomicCoordinatesAndAtomicSpecies\n")
        species_id_map = {symbol: i for i, symbol in enumerate(unique_species, start=1)}
        for site in structure:
            specie_id = species_id_map[site.specie.symbol]
            f.write(f"  {site.frac_coords[0]:.8f}   {site.frac_coords[1]:.8f}   {site.frac_coords[2]:.8f}   {specie_id}\n")
        f.write("%endblock AtomicCoordinatesAndAtomicSpecies\n")

    print(color_text(f"[INFO] Structure exported to: {filename}", 'green'))

def main():
    parser = argparse.ArgumentParser(
        description="STB Monolayer Stacking Tool with Geometric Control & Symmetry Analysis.",
        epilog="Example usage:\n"
               "  stb_stacking -l1 bottom.fdf -l2 top.fdf -i\n"
               "  stb_stacking -l1 bottom.fdf -l2 top.fdf --batch_sym -g 3.2",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument("-l1", "--layer1", required=True, help="Path to the bottom monolayer .fdf file")
    parser.add_argument("-l2", "--layer2", required=True, help="Path to the top monolayer .fdf file")
    
    parser.add_argument("-a", "--max_area", type=float, default=150.0, help="Maximum supercell area (default: 150.0)")
    parser.add_argument("-s", "--max_strain", type=float, default=0.05, help="Maximum allowed strain fraction (default: 0.05)")
    
    match_group = parser.add_mutually_exclusive_group()
    match_group.add_argument("-i", "--interactive", action="store_true", help="Interactive mode: shows table to select match ID")
    match_group.add_argument("-id", "--match_id", type=int, default=0, help="Directly select a specific match ID (default: 0)")
    
    parser.add_argument("--batch_sym", action="store_true", help="Generate all high-symmetry stackings automatically")
    
    gap_group = parser.add_mutually_exclusive_group()
    gap_group.add_argument("-g", "--gap", type=float, default=3.2, help="Van der Waals gap in Angstroms (default: 3.2)")
    gap_group.add_argument("--gap_range", type=float, nargs=3, metavar=('START', 'END', 'STEPS'), help="Range for energy curve. Ex: 2.5 4.5 11")
    
    parser.add_argument("-v", "--vacuum", type=float, default=None, help="Target vacuum space in Angstroms. Inherits L1 by default.")
    parser.add_argument("-t", "--twist", type=float, default=0.0, help="Twist angle of layer 2 in degrees (default: 0.0)")
    parser.add_argument("-tx", "--shift_x", type=float, default=0.0, help="Fractional shift of layer 2 in X (default: 0.0)")
    parser.add_argument("-ty", "--shift_y", type=float, default=0.0, help="Fractional shift of layer 2 in Y (default: 0.0)")
    parser.add_argument("-sm", "--strain_mode", choices=['top', 'bottom', 'sym'], default='top', help="Strain distribution mode (default: top)")
    
    parser.add_argument("-o", "--output", default="stacked_structure.fdf", help="Output .fdf base file name")
    parser.add_argument("--sym_out", default="symmetry_report.txt", help="Output text file for symmetry analysis")
    parser.add_argument("-sp", "--symprec", type=float, default=0.01, help="Symmetry tolerance in Angstroms (default: 0.01)")
    
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro()

    print("\n" + color_text("STACKING PROCESS:", 'bold'))
    print("-" * 60)

    print(color_text("[INFO] Parsing FDF files...", 'cyan'))
    layer1 = parse_fdf_to_pymatgen(args.layer1)
    layer2 = parse_fdf_to_pymatgen(args.layer2)

    if args.gap_range:
        gaps = np.linspace(args.gap_range[0], args.gap_range[1], int(args.gap_range[2]))
        print(color_text(f"[INFO] Energy Curve Mode detected. Generating {len(gaps)} distance points.", 'cyan'))
    else:
        gaps = [args.gap]

    results = stack_heterostructure(
        layer1, layer2, gaps=gaps, target_vacuum=args.vacuum,
        max_area=args.max_area, max_strain=args.max_strain,
        shift_x=args.shift_x, shift_y=args.shift_y, twist_angle=args.twist,
        strain_mode=args.strain_mode, match_id=args.match_id,
        interactive=args.interactive, batch_sym=args.batch_sym
    )

    base_out_name = args.output[:-4] if args.output.endswith('.fdf') else args.output
    base_sym_name = args.sym_out[:-4] if args.sym_out.endswith('.txt') else args.sym_out

    is_multi_output = args.batch_sym or len(gaps) > 1

    for name, (hetero, applied_strain) in results.items():
        print("\n" + color_text("=" * 60, 'blue'))
        print(color_text(f"  Final Heterostructure Summary: [{name}]", 'bold'))
        print(color_text("=" * 60, 'blue'))
        
        # Extrai a simetria rapidamente para o sumário
        try:
            sga_het = SpacegroupAnalyzer(hetero, symprec=args.symprec)
            het_sg = f"{sga_het.get_space_group_symbol()} ({sga_het.get_space_group_number()})"
            het_pg = sga_het.get_point_group_symbol()
            sym_str = f"{het_sg} | Point Group: {het_pg}"
        except Exception as e:
            sym_str = f"Analysis failed ({e})"

        # Formatação com padding para alinhamento vertical perfeito
        print(f"{'Formula':<26}: {color_text(hetero.composition.reduced_formula, 'yellow')}")
        print(f"{'Total Atoms':<26}: {color_text(str(len(hetero)), 'yellow')}")
        print(f"{'Lattice C Parameter':<26}: {color_text(f'{hetero.lattice.c:.2f} Å', 'yellow')} (Vacuum Corrected)")
        print(f"{'Max Applied Linear Strain':<26}: {color_text(f'{applied_strain:.2%}', 'yellow')}")
        print(f"{'Hetero Symmetry':<26}: {color_text(sym_str, 'cyan')}")
        print(color_text("-" * 60, 'blue'))
        
        if is_multi_output:
            out_filename = f"{base_out_name}_{name}.fdf"
            sym_filename = f"{base_sym_name}_{name}.txt"
        else:
            out_filename = args.output
            sym_filename = args.sym_out
        
        export_to_fdf(hetero, out_filename)
        
        # Evita poluir o terminal com várias tabelas grandes no modo batch
        print_to_stdout = not is_multi_output
        analyze_and_report_symmetry(layer1, layer2, hetero, filename=sym_filename, symprec=args.symprec, print_stdout=print_to_stdout)

    print("\n" + color_text("[INFO] Complete job!", 'green')) 
    print("-" * 60)
    print(color_text("Stacking successful! Let's hope your supercell doesn't break the cluster.\n", 'bold'))

if __name__ == "__main__":
    main()
