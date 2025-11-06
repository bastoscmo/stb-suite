#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################
    
VERSION = "1.9.0"  

import os
import sys
import subprocess
from time import sleep
import argparse
import textwrap
from typing import List, Dict, Callable

try:
    import readline
    readline.parse_and_bind("tab: complete")
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
    """Returns text formatted with ANSI color"""
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
        "Siesta ToolBox Suite",
        "A comprehensive toolkit for SIESTA DFT simulations",
        f"Version {VERSION} | University of Brasilia - 2025",
        "Developed by Dr. Carlos M. O. Bastos"
    ]

    print(logo)
    print("\n" + "="*60)
    for line in description:
        print(line.center(60))
        sleep(0.2) # Mantive seu sleep original
    print("="*60 + "\n")

def show_main_menu() -> None:
    """Displays the main category menu"""
    print("\n" + color_text("STB-SUITE Main Menu:", 'bold'))
    print("-"*60)
    print(f"{color_text('1.', 'yellow')} {color_text('Calculation (Preparation)', 'blue')}\n    Tools to set up new calculations (inputs, k-grids, etc.)\n")
    print(f"{color_text('2.', 'yellow')} {color_text('Analysis (Post-processing)', 'blue')}\n    Tools to analyze simulation results (bands, DOS, structures)\n")
    print(f"{color_text('3.', 'yellow')} {color_text('Utilities & Interfaces', 'blue')}\n    Helper tools for file management and conversion\n")
    print(f"{color_text('0.', 'yellow')} {color_text('Exit', 'red')}")
    print("-"*60)

def show_sub_menu(title: str, tools_dict: Dict) -> None:
    """Displays a sub-menu for a specific tool category"""
    print("\n" + "="*60)
    print(color_text(f"--- {title} ---", 'cyan').center(68))
    print("="*60 + "\n")
    
    for key, info in tools_dict.items():
        menu_title = color_text(info['title'], 'blue')
        desc = textwrap.fill(info['description'], width=55, subsequent_indent='    ')
        print(f"{color_text(str(key)+'.', 'yellow')} {menu_title}\n    {desc}\n")
    
    print(f"{color_text('0.', 'yellow')} {color_text('Back to Main Menu', 'red')}")
    print("-"*60)

def run_tool(tool_name: str, args: List[str]) -> None:
    """Executes a suite tool as a subprocess"""
    try:
        cmd = [f"{tool_name}"] + args
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(color_text(f"\nError running {tool_name}: {e}", 'red'))
    except FileNotFoundError:
        print(color_text(f"\nTool {tool_name} not found!", 'red'))
        print(color_text(f"Make sure {tool_name} is in your system's PATH.", 'yellow'))
    input("\nPress Enter to continue...")

def get_input(prompt: str, color: str = 'green') -> str:
    """
    Gets user input with a colored prompt.
    Graças ao 'import readline', esta função agora suporta Tab-completion!
    """
    return input(color_text(prompt, color))

def get_float_input(prompt: str, default: float = None) -> float:
    """Gets a float number from the user"""
    while True:
        try:
            # get_input() já tem o tab-completion
            value_str = get_input(prompt)
            if value_str == "" and default is not None:
                return default
            return float(value_str)
        except ValueError:
            print(color_text("Please enter a valid number", 'red'))

def get_int_input(prompt: str, default: int = None) -> int:
    """Gets an integer number from the user"""
    while True:
        try:
            # get_input() já tem o tab-completion
            value_str = get_input(prompt)
            if value_str == "" and default is not None:
                return default
            return int(value_str)
        except ValueError:
            print(color_text("Please enter a valid integer", 'red'))

# ==========================================================
# TOOL FUNCTIONS
# ==========================================================

def run_input_generator() -> None:
    """Interface for the Input File Generator (stb-inputfile)"""
    print("\n" + "="*60)
    print(color_text("INPUT FILE GENERATOR (stb-inputfile)", 'bold').center(60))
    print("="*60 + "\n")
    
    # --- Validação do ficheiro de entrada ---
    # Esta linha agora terá Tab-completion!
    input_file = get_input("Input structure file (e.g., struct.fdf): ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input structure file: ")
    
    # --- PONTO 1: Menu Numérico para o tipo de cálculo ---
    mode_list = [
        'total_energy', 'total_energy+d3',
        'relax', 'relax+d3',
        'aimd', 'aimd+d3',
        'bands', 'bands+d3'
    ]
    
    print(f"\n{color_text('Available calculation modes:', 'yellow')}")
    for i, mode in enumerate(mode_list, 1):
        print(f"  {color_text(str(i)+'.', 'yellow')} {mode}")

    choice = 0
    max_choice = len(mode_list)
    
    while not (1 <= choice <= max_choice):
        choice = get_int_input(f"\nSelect calculation mode (1-{max_choice}): ")
        if not (1 <= choice <= max_choice):
            print(color_text(f"Invalid choice! Please select between 1 and {max_choice}.", 'red'))
            
    calc_type = mode_list[choice - 1]
    print(f"Selected mode: {color_text(calc_type, 'cyan')}") 

    # --- PONTO 2: Validação do caminho do Pseudopotencial ---
    
    args = [
        input_file, 
        "--type", calc_type,
        "--no-intro"
    ]
    
    while True:
        # Esta linha agora terá Tab-completion!
        pp_path_input = get_input("\nPseudopotentials path (optional, press Enter to skip): ")
        
        if not pp_path_input.strip():
            print(color_text("Skipping pseudopotential path.", 'yellow'))
            break 
        
        pp_path = os.path.expanduser(pp_path_input)
        
        if os.path.isdir(pp_path):
            args.extend(["--pp-path", pp_path])
            print(color_text(f"Using PP path: {pp_path}", 'green'))
            break
        else:
            print(color_text(f"Path not found: '{pp_path}'", 'red'))
            print(color_text("Please enter a valid path or press Enter to skip.", 'yellow'))
    
    run_tool("stb-inputfile", args)

def run_kgrid_generator() -> None:
    """Interface for the K-Grid Generator (stb-kgrid)"""
    print("\n" + "="*60)
    print(color_text("K-GRID GENERATOR (stb-kgrid)", 'bold').center(60))
    print("="*60 + "\n")
    
    # Esta linha agora terá Tab-completion!
    input_file = get_input("Input structure file (fdf/poscar/cif/fhi): ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input structure file: ")
    
    type_list = ['fdf', 'poscar', 'cif', 'fhi']
    print(f"\n{color_text('Available types:', 'yellow')} {', '.join(type_list)}")
    file_type = get_input("Structure file type: ").lower()
    while file_type not in type_list:
        print(color_text("Invalid type!", 'red'))
        file_type = get_input("Structure file type: ").lower()
    
    density = get_float_input("\nK-point density (e.g., 0.2): ")
    while density <= 0:
        print(color_text("Density must be a positive number!", 'red'))
        density = get_float_input("K-point density (e.g., 0.2): ")
    
    args = [
        "--file", input_file,
        "--type", file_type,
        "--density", str(density),
        "--no-intro"
    ]
    
    run_tool("stb-kgrid", args)

def run_kpath_generator() -> None:
    """Interface for the K-Path Generator (stb-kpath)"""
    print("\n" + "="*60)
    print(color_text("K-PATH GENERATOR (stb-kpath)", 'bold').center(60))
    print("="*60 + "\n")
    
    # Esta linha agora terá Tab-completion!
    input_file = get_input("Input structure file (fdf/poscar): ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input structure file: ")

    type_list = ['fdf', 'poscar']
    print(f"\n{color_text('Available types:', 'yellow')} {', '.join(type_list)}")
    file_type = get_input("Structure file type (fdf/poscar): ").lower()
    while file_type not in type_list:
        print(color_text(f"Invalid type! Must be one of {type_list}", 'red'))
        file_type = get_input("Structure file type (fdf/poscar): ").lower()
    
    precision = get_float_input("\nSymmetry precision (default: 0.01): ", 0.01)

    args = [
        "--file", input_file,
        "--type", file_type,
        "--prec", str(precision),
        "--no-intro"
    ]
    
    run_tool("stb-kpath", args)

def run_dos_parser() -> None:
    """Interface for the PDOS XML Parser (stb-dos)"""
    print("\n" + "="*60)
    print(color_text("PDOS XML PARSER (stb-dos)", 'bold').center(60))
    print("="*60 + "\n")
    
    # Esta linha agora terá Tab-completion!
    input_file = get_input("Input PDOS.xml file: ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input PDOS.xml file: ")

    type_list = ['total', 'atom', 'species']
    print(f"\n{color_text('Available types:', 'yellow')} {', '.join(type_list)}")
    dos_types_str = get_input(f"DOS types (space-separated, default: total atom species): ")
    if not dos_types_str.strip():
        dos_types = ['total', 'atom', 'species']
    else:
        dos_types = dos_types_str.split()
        if not all(t in type_list for t in dos_types):
            print(color_text("Input contains invalid types. Using default.", 'yellow'))
            dos_types = ['total', 'atom', 'species']

    shift = get_input("\nEnergy shift ('fermi', '0.0', or a number, default: fermi): ")
    if not shift.strip():
        shift = 'fermi'
    
    print(f"\n{color_text('Select projection mode:', 'yellow')}")
    print(f"  {color_text('1.', 'yellow')} l (s, p, d, f) [Default]")
    print(f"  {color_text('2.', 'yellow')} ml (s, px, py, pz, dxy...)")

    choice = 0
    while not (1 <= choice <= 2):
        # Usamos get_int_input com default = 1
        choice = get_int_input(f"\nSelect mode (1-2) [default: 1]: ", 1) 
        if not (1 <= choice <= 2):
            print(color_text(f"Invalid choice! Please select 1 or 2.", 'red'))
            
    projection_mode = 'l' if choice == 1 else 'ml'
    print(f"Selected mode: {color_text(projection_mode, 'cyan')}")
    
    
    args = [
        input_file, # Positional argument
        "--shift", shift,
        "--type"
    ]
   
    args.extend(dos_types)
    args.extend(["--projection", projection_mode])
    args.append("--no-intro")
    
    run_tool("stb-dos", args)

def run_strain_generator() -> None:
    """Interface for the Strain Generator (stb-strain)"""
    print("\n" + "="*60)
    print(color_text("STRAIN GENERATOR (stb-strain)", 'bold').center(60))
    print("="*60 + "\n")
    
    # Esta linha agora terá Tab-completion!
    input_file = get_input("Input FDF file: ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input FDF file: ")
    
    direction = get_input("Strain direction (x,y,z,xy,xz,yz): ").lower()
    while not all(c in 'xyz' for c in direction) or len(direction) not in (1, 2):
        print(color_text("Invalid direction! Use x,y,z for uniaxial or xy,xz,yz for biaxial", 'red'))
        direction = get_input("Strain direction (x,y,z,xy,xz,yz): ").lower()
    
    stmin = get_float_input("Minimum strain % (default 0): ", 0.0)
    stmax = get_float_input("Maximum strain % (default 25): ", 25.0)
    while stmax <= stmin:
        print(color_text("Maximum strain must be greater than minimum strain!", 'red'))
        stmax = get_float_input("Maximum strain % (default 25): ", 25.0)
    
    step = get_float_input("Step % (default 1): ", 1.0)
    while step <= 0:
        print(color_text("Step must be positive!", 'red'))
        step = get_float_input("Step % (default 1): ", 1.0)
    
    args = [
        "--file", input_file,
        "--stdir", direction,
        "--stmin", str(stmin),
        "--stmax", str(stmax),
        "--step", str(step),
        "--no-intro"
    ]
    
    run_tool("stb-strain", args)

def run_bands_analyzer() -> None:
    """Interface for the Bands Analyzer (stb-bands)"""
    print("\n" + "="*60)
    print(color_text("BANDS ANALYZER (stb-bands)", 'bold').center(60))
    print("="*60 + "\n")
    
    # Esta linha agora terá Tab-completion!
    input_file = get_input("Input bands file (e.g., siesta.bands): ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input bands file: ")
    
    shift_options = {
        '1': ('vbm', "Valence Band Maximum"),
        '2': ('cbm', "Conduction Band Minimum"),
        '3': ('fermi', "Fermi level"),
        '4': ('manual', "Custom value")
    }
    
    print("\nEnergy reference options:")
    for key, (_, desc) in shift_options.items():
        print(f" {color_text(key, 'yellow')}. {desc}")
    
    choice = get_input("\nSelect reference (1-4): ")
    while choice not in shift_options:
        print(color_text("Invalid choice!", 'red'))
        choice = get_input("Select reference (1-4): ")
    
    shift_type, _ = shift_options[choice]
    args = ["--file", input_file, "--shift", shift_type, "--no-intro"]
    
    if shift_type == "manual":
        manual_value = get_float_input("Enter custom shift value: ")
        args.extend(["--manual-value", str(manual_value)])
    
    run_tool("stb-bands", args)

def run_dos_convolution() -> None:
    """Interface for the DOS Processor (Convolution) (stb-convdos)"""
    print("\n" + "="*60)
    print(color_text("DOS PROCESSOR (CONVOLUTION) (stb-convdos)", 'bold').center(60))
    print("="*60 + "\n")
    
    # Esta linha agora terá Tab-completion!
    input_file = get_input("Input DOS file (e.g., dos_total.dat): ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input DOS file: ")
    
    # Esta linha agora terá Tab-completion!
    out_file = get_input("Output file (default: dos_filtered.dat): ", 'green') or "dos_filtered.dat"
    size = get_int_input("Gaussian mask size (default: 11): ", 11)
    sigma = get_float_input("Standard deviation (default: 1.0): ", 1.0)
    
    args = [
        "--file", input_file,
        "--size", str(size),
        "--sigma", str(sigma),
        "--out", out_file,
        "--no-intro"
    ]
    
    run_tool("stb-convdos", args)

def run_structure_analyzer() -> None:
    """Interface for the Structure Analyzer (stb-structural)"""
    print("\n" + "="*60)
    print(color_text("STRUCTURE ANALYZER (stb-structural)", 'bold').center(60))
    print("="*60 + "\n")
    
    # Esta linha agora terá Tab-completion!
    input_file = get_input("Input structure file (CIF/POSCAR/SIESTA): ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input structure file: ")
    format_type = get_input("Format file (cif/poscar/siesta): ")
    while format_type not in ['cif','poscar','siesta'] :
        print(color_text("File type is not available!", 'red'))
        format_type = get_input("Format file (cif/poscar/siesta): ")
    
    mode = get_input("Analysis mode (list/mean): ").lower()
    while mode not in ['list', 'mean']:
        print(color_text("Invalid mode! Choose 'list' or 'mean'", 'red'))
        mode = get_input("Analysis mode (list/mean): ").lower()
    
    args = ["--file", input_file, "--mode", mode,"--format",format_type,"--no-intro"]
    
    if mode == "list":
        atom_list = get_input("Enter atom indices (comma-separated, e.g. 1,4,5): ")
        args.extend(["--list", f"[{atom_list}]"])
    
    run_tool("stb-structural", args)

def run_file_translator() -> None:
    """Interface for the File Translator (stb-translate)"""
    print("\n" + "="*60)
    print(color_text("FILE TRANSLATOR (stb-translate)", 'bold').center(60))
    print("="*60 + "\n")
    
    # Formatos suportados (o 'cif' foi adicionado à lista de saída)
    input_formats = ['fdf','poscar', 'cif', 'siesta', 'xyz', 'fhi', 'dftb', 'xsf']
    output_formats = ['cif', 'xyz', 'poscar', 'fdf', 'dftb', 'xsf', 'fhi'] # Adicionei 'cif' aqui
    
    input_file = get_input("Input file path: ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input file path: ")

    print(f"\n{color_text('Supported input formats:', 'yellow')}")
    for i, fmt in enumerate(input_formats, 1):
        print(f"  {color_text(str(i)+'.', 'yellow')} {fmt}")

    choice_in = 0
    max_in = len(input_formats)
    while not (1 <= choice_in <= max_in):
        choice_in = get_int_input(f"\nSelect input format (1-{max_in}): ")
        if not (1 <= choice_in <= max_in):
            print(color_text(f"Invalid choice! Please select between 1 and {max_in}.", 'red'))
    
    in_format = input_formats[choice_in - 1]
    print(f"Selected input format: {color_text(in_format, 'cyan')}")

    out_file = get_input("\nOutput file path: ")
    
    print(f"\n{color_text('Supported output formats:', 'yellow')}")
    for i, fmt in enumerate(output_formats, 1):
        print(f"  {color_text(str(i)+'.', 'yellow')} {fmt}")

    choice_out = 0
    max_out = len(output_formats)
    while not (1 <= choice_out <= max_out):
        choice_out = get_int_input(f"\nSelect output format (1-{max_out}): ")
        if not (1 <= choice_out <= max_out):
            print(color_text(f"Invalid choice! Please select between 1 and {max_out}.", 'red'))
            
    out_format = output_formats[choice_out - 1]
    print(f"Selected output format: {color_text(out_format, 'cyan')}")

    # ##### NOVO BLOCO: Seleção do Formato de Coordenadas #####
    print(f"\n{color_text('Select output coordinate format:', 'yellow')}")
    print(f"  {color_text('1.', 'yellow')} Cartesian (Angstroms)")
    print(f"  {color_text('2.', 'yellow')} Direct (Fractional)")
    print(f"  {color_text('3.', 'yellow')} Default (Use input format or output's default)")

    coord_choice = 0
    # Usamos default=3 para que pressionar Enter selecione a opção "Default"
    while not (1 <= coord_choice <= 3):
        coord_choice = get_int_input(f"\nSelect format (1-3) [default: 3]: ", 3) 
        if not (1 <= coord_choice <= 3):
            print(color_text(f"Invalid choice! Please select between 1 and 3.", 'red'))

    coord_format_value = None # Valor a ser passado para o argumento
    
    if coord_choice == 1:
        coord_format_value = "cartesian"
        print(f"Selected coordinate format: {color_text('Cartesian', 'cyan')}")
    elif coord_choice == 2:
        coord_format_value = "direct"
        print(f"Selected coordinate format: {color_text('Direct', 'cyan')}")
    else:
        # coord_format_value permanece None
        print(f"Selected coordinate format: {color_text('Default', 'cyan')}")
    # ##### FIM DO NOVO BLOCO #####

    # Construção dos argumentos base
    args = [
        "--in-format", in_format,
        "--in-file", input_file,
        "--out-format", out_format,
        "--out-file", out_file,
        "--no-intro"
    ]
    
    if coord_format_value:
        args.extend(["--coord-format", coord_format_value])


    if in_format == "xyz":
        print(color_text("\nXYZ format requires a separate lattice file.", 'yellow'))
        # Esta linha agora terá Tab-completion!
        lattice_file = get_input("Lattice vectors file (required for XYZ): ")
        while not os.path.isfile(lattice_file):
            print(color_text("File not found!", 'red'))
            lattice_file = get_input("Lattice vectors file: ")
        args.extend(["--lattice", lattice_file])
    
    run_tool("stb-translate", args)




def run_clean_tool() -> None:
    """Interactive interface for the Clean Files tool (stb-clean)"""
    print("\n" + "="*60)
    print(color_text("CLEAN FILES TOOL (stb-clean)", 'bold').center(60))
    print("="*60 + "\n")

    # Esta linha agora terá Tab-completion!
    path = get_input("Directory to clean (default: current): ").strip()
    if path == "":
        path = "."
    while not os.path.isdir(path):
        print(color_text("Directory not found!", 'red'))
        path = get_input("Enter a valid directory: ").strip()

    default_exts = ['.psml', '.psf', '.fdf', '.sh']
    print(f"\nExtensions to keep (space-separated, default: {' '.join(default_exts)}):")
    ext_input = get_input("Extensions: ").strip()
    if ext_input:
        extensions = ext_input.split()
    else:
        extensions = default_exts

    confirm_choice = get_input("Skip confirmation and delete directly? [y/N]: ").strip().lower()
    no_confirm = confirm_choice == 'y'

    dry_run_choice = get_input("Perform a dry run (show what would be deleted)? [y/N]: ").strip().lower()
    dry_run = dry_run_choice == 'y'

    args = ["--path", path, "--keep"] + extensions
    if no_confirm:
        args.append("--no-confirm")
    if dry_run:
        args.append("--dry-run")
    
    args.append("--no-intro")

    print()
    run_tool("stb-clean", args)

    if not dry_run:
        print("\n" + color_text("Cleanup complete. Your folder is now cleaner than my browser history.", "green"))

def run_symmetry_analyzer() -> None:
    """Interface for the Symmetry Analyzer (stb-symmetry)"""
    print("\n" + "="*60)
    print(color_text("SYMMETRY ANALYZER (stb-symmetry)", 'bold').center(60))
    print("="*60 + "\n")
    
    # Esta linha agora terá Tab-completion!
    input_file = get_input("Input structure file (CIF/POSCAR/SIESTA): ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input structure file: ")
    
    file_type = get_input("File type (poscar/cif/siesta): ").lower()
    while file_type not in ['poscar', 'cif', 'siesta']:
        print(color_text("Invalid file type! Use poscar/cif/siesta", 'red'))
        file_type = get_input("File type: ").lower()
    
    args = ["--input", input_file, "--filetype", file_type, "--no-intro"]
    run_tool("stb-symmetry", args)

def run_wantibexos_interface() -> None:
    """Interface for the Wantibexos (stb-siesta2wtb)"""
    print("\n" + "="*60)
    print(color_text("WANTIBEXOS INTERFACE (stb-siesta2wtb)", 'bold').center(60))
    print("="*60 + "\n")
    
    # Esta linha agora terá Tab-completion!
    input_file = get_input("Input FDF file: ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input FDF file: ")
    
    # Esta linha agora terá Tab-completion!
    output_file = get_input("SIESTA output file (optional): ")
    fermi = get_input("Manual Fermi level (optional): ")
    
    args = ["--input", input_file]
    if output_file:
        # Validação extra para o ficheiro opcional
        if os.path.isfile(output_file):
            args.extend(["--output", output_file])
        else:
            print(color_text(f"Warning: Output file '{output_file}' not found, skipping.", 'yellow'))
    if fermi:
        args.extend(["--fermi-level", fermi])
    
    args.append("--no-intro")
    
    run_tool("stb-siesta2wtb", args)

# ==========================================================
# SUB-MENU LOGIC
# ==========================================================

# Define the tool dictionaries
PREPARATION_TOOLS = {
    1: {'title': "Input File Generator (stb-inputfile)",
        'description': "Create a 'calc.fdf' input file from a structure file.",
        'func': run_input_generator},
    2: {'title': "K-Grid Generator (stb-kgrid)",
        'description': "Suggest a Monkhorst-Pack grid (k-points) based on desired density.",
        'func': run_kgrid_generator},
    3: {'title': "K-Path Generator (stb-kpath)",
        'description': "Generate a high-symmetry k-path for band structure calculations.",
        'func': run_kpath_generator},
    4: {'title': "Strain Generator (stb-strain)",
        'description': "Generate strained structures for calculations.",
        'func': run_strain_generator},
}

ANALYSIS_TOOLS = {
    1: {'title': "Bands Analyzer (stb-bands)",
        'description': "Analyze .bands files and calculate band gaps.",
        'func': run_bands_analyzer},
    2: {'title': "PDOS XML Parser (stb-dos)",
        'description': "Extract data from PDOS.xml by total, atom, and species.",
        'func': run_dos_parser},
    3: {'title': "DOS Processor (Convolution) (stb-convdos)", 
        'description': "Apply Gaussian convolution to Density of States (DOS) files.",
        'func': run_dos_convolution},
    4: {'title': "Structure Analyzer (stb-structural)", 
        'description': "Calculate ECN and analyze structural properties.",
        'func': run_structure_analyzer},
    5: {'title': "Symmetry Analyzer (stb-symmetry)",
        'description': "Analyze the symmetry of crystal structures.",
        'func': run_symmetry_analyzer},
}

UTILITY_TOOLS = {
    1: {'title': "File Translator (stb-translate)",
        'description': "Convert between file formats (CIF, POSCAR, fdf, xyz...).",
        'func': run_file_translator},
    2: {'title': "Clean File Tools (stb-clean)",
        'description': "Clean the directory of calculation files (except essential ones).",
        'func': run_clean_tool},
    3: {'title': "Wantibexos Interface (stb-siesta2wtb)",
        'description': "Convert SIESTA Hamiltonian to Wantibexos format.",
        'func': run_wantibexos_interface},
}

def run_sub_menu(title: str, tools_dict: Dict) -> None:
    """Handles the logic for showing and running a sub-menu"""
    while True:
        show_sub_menu(title, tools_dict)
        try:
            choice_str = get_input(f"\nSelect an option (0-{len(tools_dict)}): ")
            
            if choice_str == '0':
                break # Go back to the main menu
            
            try:
                choice = int(choice_str)
            except ValueError:
                choice = float(choice_str)

            if choice in tools_dict:
                tools_dict[choice]['func']() # Run the selected tool
            else:
                print(color_text(f"\nInvalid choice! Please select between 0 and {len(tools_dict)}.", 'red'))
                sleep(1)
                
        except ValueError:
            print(color_text("\nPlease enter a valid number!", 'red'))
            sleep(1)
        except KeyboardInterrupt:
            break # Go back to the main menu

# ==========================================================
# MAIN FUNCTION
# ==========================================================

def main():
    """Main function to run the STB-SUITE interface"""
    show_intro()
    
    while True:
        show_main_menu()
        
        try:
            choice = get_input("\nSelect an option (0-3): ")
            
            if choice == '1':
                run_sub_menu("Calculation (Preparation)", PREPARATION_TOOLS)
            elif choice == '2':
                run_sub_menu("Analysis (Post-processing)", ANALYSIS_TOOLS)
            elif choice == '3':
                run_sub_menu("Utilities & Interfaces", UTILITY_TOOLS)
            elif choice == '0':
                print(color_text("\nThank you for using STB-SUITE!", 'cyan'))
                break
            else:
                print(color_text("\nInvalid choice! Please select between 0 and 3.", 'red'))
                sleep(1)
                
        except ValueError:
            print(color_text("\nPlease enter a valid number!", 'red'))
            sleep(1)
        except KeyboardInterrupt:
            print(color_text("\n\nOperation cancelled by user.", 'yellow'))
            break

if __name__ == "__main__":
    main()
