#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.8.0"

import os
import sys
import shutil
import argparse
from time import sleep
import glob
from phonopy import Phonopy
from phonopy.interface.siesta import read_siesta, write_siesta

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

def get_required_pseudos(symbols: list, pseudo_dir: str):
    """
    Verifica se os pseudopotenciais existem no diretório fornecido.
    Retorna a lista de caminhos encontrados e a lista de elementos ausentes.
    """
    unique_elements = set(symbols)
    found_pseudos = []
    missing_elements = []

    for element in unique_elements:
        psf_path = os.path.join(pseudo_dir, f"{element}.psf")
        psml_path = os.path.join(pseudo_dir, f"{element}.psml")

        if os.path.exists(psf_path):
            found_pseudos.append(psf_path)
        elif os.path.exists(psml_path):
            found_pseudos.append(psml_path)
        else:
            missing_elements.append(element)

    return found_pseudos, missing_elements

def main():
    parser = argparse.ArgumentParser(
        description="Automate SIESTA phonon displacement folders with Phonopy.",
        epilog="Example usage:\n"
               "  stb_phonon --structure structure.fdf --calc calc.fdf --dim 2 2 2\n"
               "  stb_phonon --pseudo-dir ~/pseudos --distance 0.015",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("-dim", type=int, nargs=3, default=[2, 2, 2], 
                        help="Supercell dimensions (default: 2 2 2)")
    
    parser.add_argument("-d", "--distance", type=float, default=0.01, 
                        help="Displacement distance in Angstroms (default: 0.01)")
    
    parser.add_argument("-s", "--structure", type=str, default="structure.fdf", 
                        help="Input structure file containing the unit cell (default: structure.fdf)")
    
    parser.add_argument("-c", "--calc", type=str, default="calc.fdf", 
                        help="Input calculation parameters file (default: calc.fdf)")
    
    parser.add_argument("-p", "--pseudo-dir", type=str, default=".", 
                        help="Directory containing the pseudopotentials (.psf or .psml) (default: current directory)")
    
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    
    parser.add_argument("--no-intro", dest="intro", action="store_false", 
                        help="Do not show the introduction")

    args = parser.parse_args()

    # Exibe a introdução se a flag não for acionada
    if args.intro:
        show_intro()

    print("\n" + color_text("PHONON CALCULATION (VIA PHONOPY):", 'bold'))
    print("-"*60)

    # 1. Validação de arquivos de entrada
    print("\n[INFO] Validating input files ...")
    if not os.path.exists(args.structure):
        print(color_text(f"[ERROR] Structure file '{args.structure}' not found in the current directory.", 'red'))
        sys.exit(1)
    if not os.path.exists(args.calc):
        print(color_text(f"[ERROR] Calculation file '{args.calc}' not found in the current directory.", 'red'))
        sys.exit(1)

    # 2. Leitura da estrutura
    print(f"[INFO] Reading unit cell from '{args.structure}' ...")
    try:
        unitcell = read_siesta(args.structure)
    except Exception as e:
        print(color_text(f"[ERROR] Failed to read {args.structure}. Make sure it's properly formatted.\nDetails: {e}", 'red'))
        sys.exit(1)

    # 3. Extração e validação de Pseudopotenciais
    symbols = unitcell.symbols
    unique_elements = list(set(symbols))
    print(f"[INFO] Elements found in unit cell: {', '.join(unique_elements)}")
    
    print(f"[INFO] Searching for pseudopotentials in '{args.pseudo_dir}' ...")
    pseudos_to_copy, missing = get_required_pseudos(unique_elements, args.pseudo_dir)

    if missing:
        print(color_text(f"\n[CRITICAL ERROR] Missing pseudopotentials for the following elements: {', '.join(missing)}", 'red'))
        print(color_text(f"Action required: Please add the necessary '{missing[0]}.psf' or '{missing[0]}.psml' files into the '{args.pseudo_dir}' directory and rerun the script.", 'yellow'))
        sys.exit(1)
        
    print(f"[INFO] Found all required pseudopotentials: {', '.join([os.path.basename(p) for p in pseudos_to_copy])}")

    # 4. Inicialização do Phonopy
    print(f"[INFO] Generating supercell {args.dim} with {args.distance} Å displacements ...")
    supercell_matrix = [
        [args.dim[0], 0, 0], 
        [0, args.dim[1], 0], 
        [0, 0, args.dim[2]]
    ]
    
    phonon = Phonopy(unitcell, supercell_matrix=supercell_matrix)
    phonon.generate_displacements(distance=args.distance)
    supercells = phonon.supercells_with_displacements

    # 5. Criação dos diretórios e cópia
    output_root = "phonon_runs"
    os.makedirs(output_root, exist_ok=True)
    
    print(f"[INFO] Building {len(supercells)} displacement folders in '{output_root}' ...")

    for i, scell in enumerate(supercells):
        if scell is None:
            continue
            
        folder_name = os.path.join(output_root, f"disp-{i+1:03d}")
        os.makedirs(folder_name, exist_ok=True)
        
        # A. Escreve a supercélula deslocada
        disp_struct_path = os.path.join(folder_name, args.structure)
        write_siesta(disp_struct_path, scell)
        
        # B. Copia o calc.fdf
        shutil.copy(args.calc, os.path.join(folder_name, args.calc))
        
        # C. Copia apenas os pseudopotenciais exigidos
        for pseudo_path in pseudos_to_copy:
            pseudo_filename = os.path.basename(pseudo_path)
            shutil.copy(pseudo_path, os.path.join(folder_name, pseudo_filename))

    # 6. Salvar metadados do Phonopy
    yaml_path = os.path.join(output_root, "phonopy_disp.yaml")
    phonon.save(yaml_path)
    print(f"[INFO] Saved Phonopy metadata to '{yaml_path}'")
    
    print("\n[INFO] Complete job!") 
    print("\n"+"-"*60)
    print(color_text("Phonon folders ready! Let the atoms shake, rattle and roll.\n\n", 'bold'))

if __name__ == "__main__":
    main()
