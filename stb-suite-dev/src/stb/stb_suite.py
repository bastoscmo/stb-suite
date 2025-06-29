#!/usr/bin/env python3

#################################################
# STB-SUITE - Siesta ToolBox Suite              #
# version 1.5.1                                 #
# UnB - 2025                                    #
# Dr. Carlos M. O. Bastos                       #
# Unified Beautiful Interface                   #
#################################################

import os
import sys
import subprocess
from time import sleep
import argparse
import textwrap
from typing import List, Dict

VERSION = "1.5.7"

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

def show_menu(options: List[Dict[str, str]]) -> None:
    """Exibe um menu bonito com opções numeradas"""
    print("\n" + color_text("STB-SUITE Main Menu:", 'bold'))
    print("-"*60)
    
    for i, option in enumerate(options, 1):
        title = color_text(option['title'], 'blue')
        desc = textwrap.fill(option['description'], width=55, subsequent_indent='    ')
        print(f"{color_text(str(i)+'.', 'yellow')} {title}\n    {desc}\n")
    
    print(f"{color_text('0.', 'yellow')} {color_text('Exit', 'red')}")
    print("-"*60)

def run_tool(tool_name: str, args: List[str]) -> None:
    """Executa uma ferramenta da suíte como subprocesso"""
    try:
        cmd = [f"{tool_name}"] + args
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(color_text(f"\nError running {tool_name}: {e}", 'red'))
    except FileNotFoundError:
        print(color_text(f"\nTool {tool_name}.py not found!", 'red'))
    input("\nPress Enter to continue...")

def get_input(prompt: str, color: str = 'green') -> str:
    """Obtém entrada do usuário com prompt colorido"""
    return input(color_text(prompt, color))

def get_float_input(prompt: str, default: float = None) -> float:
    """Obtém um número float do usuário"""
    while True:
        try:
            value = get_input(prompt)
            return float(value) if value else default
        except ValueError:
            print(color_text("Please enter a valid number", 'red'))

def get_int_input(prompt: str, default: int = None) -> int:
    """Obtém um número inteiro do usuário"""
    while True:
        try:
            value = get_input(prompt)
            return int(value) if value else default
        except ValueError:
            print(color_text("Please enter a valid integer", 'red'))

def run_strain_generator() -> None:
    """Interface para o Strain Generator"""
    print("\n" + "="*60)
    print(color_text("STRAIN GENERATOR", 'bold').center(60))
    print("="*60 + "\n")
    
    # Get parameters
    input_file = get_input("Input FDF file: ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input FDF file: ")
    
    direction = get_input("Strain direction (x,y,z,xy,xz,yz): ").lower()
    while not all(c in 'xyz' for c in direction) or len(direction) not in (1, 2):
        print(color_text("Invalid direction! Use x,y,z for uniaxial or xy,xz,yz for biaxial", 'red'))
        direction = get_input("Strain direction (x,y,z,xy,xz,yz): ").lower()
    
    stmin = get_float_input("Minimum strain % (default 0): ", 0)
    stmax = get_float_input("Maximum strain % (default 25): ", 25)
    while stmax <= stmin:
        print(color_text("Maximum strain must be greater than minimum strain!", 'red'))
        stmax = get_float_input("Maximum strain % (default 25): ", 25)
    
    step = get_float_input("Step % (default 1): ", 1)
    while step <= 0:
        print(color_text("Step must be positive!", 'red'))
        step = get_float_input("Step % (default 1): ", 1)
    
    # Build command
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
    """Interface para o Bands Analyzer"""
    print("\n" + "="*60)
    print(color_text("BANDS ANALYZER", 'bold').center(60))
    print("="*60 + "\n")
    
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

def run_dos_processor() -> None:
    """Interface para o DOS Processor"""
    print("\n" + "="*60)
    print(color_text("DOS PROCESSOR", 'bold').center(60))
    print("="*60 + "\n")
    
    input_file = get_input("Input DOS file: ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input DOS file: ")
    
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
    """Interface para o Structure Analyzer"""
    print("\n" + "="*60)
    print(color_text("STRUCTURE ANALYZER", 'bold').center(60))
    print("="*60 + "\n")
    
    input_file = get_input("Input structure file (CIF/POSCAR/SIESTA): ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input structure file: ")
    format_type = get_input("Format file (cif/poscar/siesta): ")
    while format_type not in ['cif','poscar','siesta'] :
        print(color_text("File type are no avaliable!", 'red'))
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
    """Interface para o File Translator"""
    print("\n" + "="*60)
    print(color_text("FILE TRANSLATOR", 'bold').center(60))
    print("="*60 + "\n")
    
    input_formats = ['poscar', 'cif', 'siesta', 'xyz', 'fhi', 'dftb', 'xsf']
    output_formats = ['xyz', 'poscar', 'fdf', 'dftb', 'xsf', 'fhi']
    
    print("Supported input formats: " + ", ".join(input_formats))
    print("Supported output formats: " + ", ".join(output_formats))
    
    input_file = get_input("\nInput file path: ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input file path: ")
    
    in_format = get_input("Input format: ").lower()
    while in_format not in input_formats:
        print(color_text("Invalid input format!", 'red'))
        in_format = get_input("Input format: ").lower()
    
    out_format = get_input("Output format: ").lower()
    while out_format not in output_formats:
        print(color_text("Invalid output format!", 'red'))
        out_format = get_input("Output format: ").lower()
    
    out_file = get_input("Output file path: ")
    
    args = [
        "--in-format", in_format,
        "--in-file", input_file,
        "--out-format", out_format,
        "--out-file", out_file,
        "--no-intro"
    ]
    
    if in_format == "xyz":
        lattice_file = get_input("Lattice vectors file (required for XYZ): ")
        args.extend(["--lattice", lattice_file])
    
    run_tool("stb-translate", args)

def run_clean_tool() -> None:
    """Interactive interface for the Clean Files tool"""
    print("\n" + "="*60)
    print(color_text("CLEAN FILES TOOL", 'bold').center(60))
    print("="*60 + "\n")

    path = get_input("Directory to clean (default: current): ").strip()
    if path == "":
        path = "."
    while not os.path.isdir(path):
        print(color_text("Directory not found!", 'red'))
        path = get_input("Enter a valid directory: ").strip()

    # Ask for extensions to keep
    default_exts = ['.psml', '.psf', '.fdf', '.sh']
    print("\nExtensions to keep (space-separated, default: .psml .psf .fdf .sh):")
    ext_input = get_input("Extensions: ").strip()
    if ext_input:
        extensions = ext_input.split()
    else:
        extensions = default_exts

    # Confirmation preference
    confirm_choice = get_input("Skip confirmation and delete directly? [y/N]: ").strip().lower()
    no_confirm = confirm_choice == 'y'

    # Dry run option
    dry_run_choice = get_input("Perform a dry run (show what would be deleted)? [y/N]: ").strip().lower()
    dry_run = dry_run_choice == 'y'

    # Prepare arguments
    args = ["--path", path, "--keep"] + extensions
    if no_confirm:
        args.append("--no-confirm")
    if dry_run:
        args.append("--dry-run")

    print()
    run_tool("stb-clean", args)

    # Funny outro if not dry run
    if not dry_run:
        print("\n" + color_text("Cleanup complete. Your folder is now cleaner than my browser history.", "green"))

def run_symmetry_analyzer() -> None:
    """Interface para o Symmetry Analyzer"""
    print("\n" + "="*60)
    print(color_text("SYMMETRY ANALYZER", 'bold').center(60))
    print("="*60 + "\n")
    
    input_file = get_input("Input structure file (CIF/POSCAR/SIESTA): ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input structure file: ")
    
    file_type = get_input("File type (poscar/cif/siesta): ").lower()
    while file_type not in ['poscar', 'cif', 'siesta']:
        print(color_text("Invalid file type! Use poscar/cif/siesta", 'red'))
        file_type = get_input("File type: ").lower()
    
    args = ["--input", input_file, "--filetype", file_type]
    run_tool("stb-symmetry", args)

def run_wantibexos_interface() -> None:
    """Interface para o Wantibexos"""
    print("\n" + "="*60)
    print(color_text("WANTIBEXOS INTERFACE", 'bold').center(60))
    print("="*60 + "\n")
    
    input_file = get_input("Input FDF file: ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input FDF file: ")
    
    output_file = get_input("SIESTA output file (optional): ")
    fermi = get_input("Manual Fermi level (optional): ")
    
    args = ["--input", input_file]
    if output_file:
        args.extend(["--output", output_file])
    if fermi:
        args.extend(["--fermi-level", fermi])
    
    run_tool("stb-siesta2wtb", args)

def main():
    show_intro()
    
    menu_options = [
        {
            'title': "Strain Generator",
            'description': "Generate strained structures for SIESTA calculations"
        },
        {
            'title': "Bands Analyzer",
            'description': "Analyze band structure results and calculate band gaps"
        },
        {
            'title': "DOS Processor",
            'description': "Apply Gaussian convolution to Density of States"
        },
        {
            'title': "Structure Analyzer",
            'description': "Calculate ECN and analyze crystal structures"
        },
        {
            'title': "File Translator",
            'description': "Convert between different file formats for DFT calculations"
        },
        {
            'title': "Clean File Tools",
            'description': "Remove all files except those with specified extensions."
        },
        {
            'title': "Symmetry Analyzer",
            'description': "Analyze crystal symmetry"
        },
        {
            'title': "Wantibexos Interface",
            'description': "Convert SIESTA Hamiltonian to Wantibexos format"
        }
    ]
    
    tool_functions = {
        1: run_strain_generator,
        2: run_bands_analyzer,
        3: run_dos_processor,
        4: run_structure_analyzer,
        5: run_file_translator,
        6: run_clean_tool,
        7: run_symmetry_analyzer,
        8: run_wantibexos_interface
    }
    
    while True:
        show_menu(menu_options)
        
        try:
            choice = get_input("\nSelect an option (0-8): ")
            if choice == '0':
                print(color_text("\nThank you for using STB-SUITE!", 'cyan'))
                break
            
            choice = int(choice)
            if choice in tool_functions:
                tool_functions[choice]()
        except ValueError:
            print(color_text("\nPlease enter a valid number!", 'red'))
            sleep(1)
        except KeyboardInterrupt:
            print(color_text("\n\nOperation cancelled by user.", 'yellow'))
            break

if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="STB-SUITE - Unified Interface")
    parser.add_argument("--tool", type=int, help="Directly launch a specific tool (1-8)")
    args = parser.parse_args()
    
    if args.tool:
        # Modo de execução direta
        tool_functions = {
            1: run_strain_generator,
            2: run_bands_analyzer,
            3: run_dos_processor,
            4: run_structure_analyzer,
            5: run_file_translator,
            6: run_clean_tool,
            7: run_symmetry_analyzer,
            8: run_wantibexos_interface
        }
        
        if args.tool in tool_functions:
            show_intro()
            tool_functions[args.tool]()
        else:
            print(color_text("Invalid tool number! Use 1-8.", 'red'))
    else:
        # Modo interativo
        main()
