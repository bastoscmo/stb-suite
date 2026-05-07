#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.8.0"

import os
import sys
import glob
import subprocess
import argparse
from time import sleep
try:
    import phonopy
except ImportError:
    print("\n[ERROR] Phonopy is not installed. Please install it using: pip install phonopy")
    sys.exit(1)

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

def main():
    parser = argparse.ArgumentParser(
        description="Post-processing for SIESTA phonon calculations using Phonopy.",
        epilog="Example usage:\n"
               "  phonons_post -dir phonon_runs -l siesta -m 20 20 20\n"
               "  phonons_post --tmin 0 --tmax 800 --tstep 5",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("-dir", "--directory", type=str, default="phonon_runs", 
                        help="Directory containing the displacement folders (default: phonon_runs)")
    parser.add_argument("-l", "--label", type=str, default="siesta", 
                        help="SystemLabel used in calc.fdf (default: siesta)")
    parser.add_argument("-m", "--mesh", type=int, nargs=3, default=[20, 20, 20], 
                        help="Q-point mesh for thermal properties (default: 20 20 20)")
    parser.add_argument("--tmin", type=float, default=0.0, help="Minimum temperature in K (default: 0)")
    parser.add_argument("--tmax", type=float, default=1000.0, help="Maximum temperature in K (default: 1000)")
    parser.add_argument("--tstep", type=float, default=10.0, help="Temperature step in K (default: 10)")
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro()

    print("\n" + color_text("PHONON POST-PROCESSING:", 'bold'))
    print("-" * 60)

    phonon_dir = args.directory
    system_label = args.label

    # 1. Validação do Diretório
    print(f"\n[INFO] Validating phonon directory '{phonon_dir}' ...")
    if not os.path.exists(phonon_dir):
        print(color_text(f"[ERROR] Directory '{phonon_dir}' not found.", 'red'))
        sys.exit(1)

    yaml_file = os.path.join(phonon_dir, "phonopy_disp.yaml")
    if not os.path.exists(yaml_file):
        print(color_text(f"[ERROR] '{yaml_file}' not found. Did you run the displacement generator?", 'red'))
        sys.exit(1)

    # 2. Extração de Forças (Criando FORCE_SETS)
    print(f"[INFO] Extracting forces from .FA files with label '{system_label}' ...")
    fa_files_abs = sorted(glob.glob(os.path.join(phonon_dir, "disp-*", f"{system_label}.FA")))
    
    if not fa_files_abs:
        print(color_text(f"[ERROR] No {system_label}.FA files found in {phonon_dir}/disp-*", 'red'))
        print(color_text("Make sure SIESTA calculations finished successfully.", 'yellow'))
        sys.exit(1)

    # Pegamos os caminhos relativos ao phonon_dir para rodar o comando lá dentro
    fa_files_rel = [os.path.relpath(f, phonon_dir) for f in fa_files_abs]
    print(f"[INFO] Found {len(fa_files_rel)} force files. Generating FORCE_SETS...")
    
    cmd = ["phonopy", "--siesta", "-f"] + fa_files_rel
    
    try:
        # Roda o comando de extração de forças DENTRO da pasta phonon_runs
        subprocess.check_call(cmd, cwd=phonon_dir, stdout=subprocess.DEVNULL)
        print(color_text("[SUCCESS] FORCE_SETS generated successfully.", 'green'))
    except subprocess.CalledProcessError as e:
        print(color_text(f"[ERROR] Failed to generate FORCE_SETS. Phonopy error: {e}", 'red'))
        sys.exit(1)

    # 3. Propriedades Térmicas usando a API Python do Phonopy
    print(f"\n[INFO] Initializing Phonopy API and loading FORCE_SETS ...")
    try:
        # Como o yaml e o FORCE_SETS estão na mesma pasta, precisamos avisar o phonopy ou mudar o dir de execução
        os.chdir(phonon_dir)
        phonon = phonopy.load("phonopy_disp.yaml")
    except Exception as e:
        print(color_text(f"[ERROR] Could not load phonopy data: {e}", 'red'))
        sys.exit(1)

    print(f"[INFO] Running thermal properties calculation ...")
    print(f"       -> Q-Mesh: {args.mesh}")
    print(f"       -> Temperature Range: {args.tmin} K to {args.tmax} K (step: {args.tstep} K)")

    phonon.run_mesh(args.mesh)
    phonon.run_thermal_properties(t_min=args.tmin, t_max=args.tmax, t_step=args.tstep)

    # 4. Salvando Gráficos e Dados
    print("[INFO] Exporting results ...")
    
    tp_plot = phonon.plot_thermal_properties()
    plot_filename = "thermal_properties.png"
    tp_plot.savefig(plot_filename, dpi=300)
    print(color_text(f" -> Saved plot as '{os.path.join(phonon_dir, plot_filename)}'", 'cyan'))
    
    tp_dict = phonon.get_thermal_properties_dict()
    temperatures = tp_dict['temperatures']
    free_energy = tp_dict['free_energy']
    entropy = tp_dict['entropy']
    heat_capacity = tp_dict['heat_capacity']
    
    dat_filename = "thermal_properties.dat"
    with open(dat_filename, "w") as f:
        f.write("# T(K)       FreeEnergy(kJ/mol)  Entropy(J/K/mol)  HeatCapacity(J/K/mol)\n")
        for i in range(len(temperatures)):
            f.write(f"{temperatures[i]:10.2f} {free_energy[i]:18.6f} {entropy[i]:18.6f} {heat_capacity[i]:18.6f}\n")
            
    print(color_text(f" -> Saved raw data as '{os.path.join(phonon_dir, dat_filename)}'", 'cyan'))
    
    # Retorna ao diretório original
    os.chdir("..")

    print("\n" + "-" * 60)
    print(color_text("Post-processing complete! Results are in your phonon directory.\n", 'bold'))

if __name__ == "__main__":
    main()
