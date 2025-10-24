#!/usr/bin/env python
"""
SIESTA ToolBox - HAM_WANTIBEXOS Interface
Convert SIESTA Hamiltonian to Wantibexos code
"""

import os
import sys
import re
import argparse
import textwrap
from time import sleep
import numpy as np
import sisl
from typing import Tuple

VERSION = "1.8.0"
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
    """Format text with ANSI color codes"""
    return f"{COLORS[color]}{text}{COLORS['reset']}"

def show_intro() -> None:
    """Display STB-SUITE introduction banner"""
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
        "SIESTA-to-Wantibexos Hamiltonian Interface",
        f"Version {VERSION} | University of Brasilia - 2025",
        ""
    ]
    
    print(logo)
    print("\n" + "="*60)
    for line in description:
        print(line.center(60))
        sleep(0.1)
    print("="*60 + "\n")

def parse_arguments() -> argparse.Namespace:
    """Handle command-line arguments with argparse"""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''
        Convert SIESTA Hamiltonian to TB format for wantibexos
        Requires SIESTA calculation with SaveHS = .true. and SaveRho = .true.
        '''))
    
    parser.add_argument('-i', '--input', required=True,
                        help='Input FDF file from SIESTA calculation')
    parser.add_argument('-o', '--output', default=None,
                        help='SIESTA output file for Fermi energy extraction')
    parser.add_argument('-f', '--fermi-level', type=float, default=None,
                        help='Manual Fermi energy setting (optional)')
    parser.add_argument('-q', '--quiet', action='store_false',
                        help='Suppress terminal animations and colors')
    parser.add_argument("-v", "--version", action="version",
                        version=f"stb-translate {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    
    return parser.parse_args()

def get_system_label(fdf_path: str) -> str:
    """Extract SystemLabel from SIESTA .fdf file with fallback to 'siesta'"""
    # NOTE: This function isn't strictly necessary if sisl can
    # find the .HSX file from the .fdf, but we'll keep it.
    default_label = "siesta"
    try:
        with open(fdf_path, 'r') as f:
            for line in f:
                # Strip comments and whitespace
                line = re.sub(r'#.*', '', line).strip()
                if not line:
                    continue
                
                # Split into key-value pairs
                key_value = re.split(r'\s*=\s*|\s+', line, maxsplit=1)
                if key_value[0].lower() == 'systemlabel':
                    return key_value[1].split()[0] if len(key_value) > 1 else default_label
    except Exception as e:
        pass
    return default_label

def get_fermi_energy(output_file: str) -> float:
    """Extract Fermi energy from SIESTA output"""
    try:
        with open(output_file, 'r') as f:
            for line in reversed(f.readlines()):
                if 'Fermi =' in line:
                    return float(line.split()[-1])
    except (FileNotFoundError, IndexError) as e:
        raise RuntimeError(f"Fermi energy extraction failed: {str(e)}")

from typing import Tuple

def spin_mapping(spin_obj) -> Tuple[str, int]:
    """Map SIESTA spin types to Wantibexos conventions (code, multiplicity)."""
    # This function was correct. It looks for the keyword
    # in the full string (e.g., 'unpolarized' in '<Spin.unpolarized>')
    spin_str = str(spin_obj).lower()
    mapping = {
        'unpolarized': ('NP', 1),
        'polarized': ('SP', 2),
        'non-colinear': ('SOC', 2),
        'spin-orbit': ('SOC', 2)
    }
    for key in mapping:
        if key in spin_str:
            return mapping[key]
    raise ValueError(f"Unsupported spin type in: '{spin_str}'")



def format_float(value: float) -> str:
    """Format floating point numbers for TB files"""
    return f"{value:12.6f}"


def write_basis(hamiltonian, spin_suffix: str, spin_factor: int) -> None:
    """Write basis set information to file"""
    filename = f"basis_set-{spin_suffix}.basis"
    with open(filename, 'w') as f:
        f.write("# bindex aspecie ax ay az l m spin\n")
        
        for spin in range(spin_factor):
            io = 0
            for ia, atom in enumerate(hamiltonian.geometry.atoms):
                xyz = hamiltonian.geometry.xyz[ia, :]
                for orbital in atom:
                    spin_val = 1 if spin == 0 else -1
                    if spin_suffix == 'SOC':
                        spin_val = 1 - 2*spin  # 1 for up, -1 for down
                    f.write(f"{io+1 + spin*hamiltonian.no:6d} {atom.tag:6} "
                            f"{format_float(xyz[0])} {format_float(xyz[1])} "
                            f"{format_float(xyz[2])} {orbital.l:2} "
                            f"{orbital.m:2} {spin_val:3}\n")
                    io += 1

def write_hamiltonian(hamiltonian, spin_suffix: str, fermi: float) -> None:
    """Write Hamiltonian data in TB format"""
    filename = f"tb-{spin_suffix}.ham"
    nbasis = hamiltonian.no
    ncell = np.prod(hamiltonian.nsc)

    with open(filename, 'w') as f:
        # Header section
        f.write(f"{spin_suffix}\n")
        f.write(f"{format_float(0.0)}\n")  # Scissors operator
        f.write(f"{format_float(fermi)}\n")
        
        # Cell vectors
        for vec in hamiltonian.cell:
            f.write(" ".join([format_float(v) for v in vec]) + "\n")
        
        # Dimensions
        f.write(f"{hamiltonian.no * (2 if spin_suffix != 'NP' else 1)}\n")
        f.write(f"{ncell}\n")
        f.write("# rcell_x rcell_y rcell_z i j ReH ImH S\n")
        
        # Hamiltonian data handling based on spin type
        if spin_suffix == 'NP':
            # Unpolarized case
            for icell in range(ncell):
                for i in range(nbasis):
                    for j in range(nbasis):
                        sc_index = hamiltonian.geometry.o2sc(j + icell*nbasis)
                        H = hamiltonian[i, j + icell*nbasis][0]
                        S = hamiltonian.S[i, j + icell*nbasis]
                        
                        line = (f"{' '.join(format_float(v) for v in sc_index)} "
                                f"{i+1} {j+1} {format_float(H)} "
                                f"{format_float(0.0)} {format_float(S)}\n")
                        f.write(line)
        
        elif spin_suffix == 'SP':
            # Spin polarized case
            S = np.zeros((ncell+1, 2*nbasis+2, 2*nbasis+2))
            reH = np.zeros((ncell+1, 2*nbasis+2, 2*nbasis+2))
            imH = np.zeros((ncell+1, 2*nbasis+2, 2*nbasis+2)) # Although not used here, it's good for consistency
            a = np.zeros((ncell+1, 2*nbasis+2))
            b = np.zeros((ncell+1, 2*nbasis+2))
            c = np.zeros((ncell+1, 2*nbasis+2))

            # Populate matrices
            for icell in range(ncell):
                for i in range(nbasis):
                    for j in range(nbasis):
                        # Up-up block
                        reH[icell+1, i+1, j+1] = hamiltonian[i, j + icell*nbasis][0]
                        S[icell+1, i+1, j+1] = hamiltonian.S[i, j + icell*nbasis]
                        
                        # Down-down block
                        reH[icell+1, i+1+nbasis, j+1+nbasis] = hamiltonian[i, j + icell*nbasis][1]
                        S[icell+1, i+1+nbasis, j+1+nbasis] = hamiltonian.S[i, j + icell*nbasis]
                        
                        # Cell coordinates
                        sc = hamiltonian.geometry.o2sc(j + icell*nbasis)
                        a[icell+1, j+1] = sc[0]
                        b[icell+1, j+1] = sc[1]
                        c[icell+1, j+1] = sc[2]
                        a[icell+1, j+1+nbasis] = sc[0]
                        b[icell+1, j+1+nbasis] = sc[1]
                        c[icell+1, j+1+nbasis] = sc[2]

            # Write to file
            for icell in range(1, ncell+1):
                for i in range(1, 2*nbasis+1):
                    for j in range(1, 2*nbasis+1):
                        line = (f"{format_float(a[icell,j])} {format_float(b[icell,j])} "
                                f"{format_float(c[icell,j])} {i} {j} "
                                f"{format_float(reH[icell,i,j])} {format_float(imH[icell,i,j])} " # ImH will be 0.0
                                f"{format_float(S[icell,i,j])}\n")
                        f.write(line)
        
        elif spin_suffix == 'SOC':
            # Non-colinear and spin-orbit cases
            S = np.zeros((ncell+1, 2*nbasis+2, 2*nbasis+2))
            reH = np.zeros((ncell+1, 2*nbasis+2, 2*nbasis+2))
            imH = np.zeros((ncell+1, 2*nbasis+2, 2*nbasis+2))
            a = np.zeros((ncell+1, 2*nbasis+2))
            b = np.zeros((ncell+1, 2*nbasis+2))
            c = np.zeros((ncell+1, 2*nbasis+2))

            for icell in range(ncell):
                for i in range(nbasis):
                    for j in range(nbasis):
                        sc = hamiltonian.geometry.o2sc(j + icell*nbasis)
                        idx = icell+1
                        
                        # Common cell coordinates
                        a[idx, j+1] = sc[0]
                        b[idx, j+1] = sc[1]
                        c[idx, j+1] = sc[2]
                        a[idx, j+1+nbasis] = sc[0]
                        b[idx, j+1+nbasis] = sc[1]
                        c[idx, j+1+nbasis] = sc[2]

                        # Matrix elements
                        H = hamiltonian[i, j + icell*nbasis]
                        S_val = hamiltonian.S[i, j + icell*nbasis]
                        
                        # Up-up block
                        reH[idx, i+1, j+1] = H[0]
                        S[idx, i+1, j+1] = S_val
                        
                        # Down-down block
                        reH[idx, i+1+nbasis, j+1+nbasis] = H[1]
                        S[idx, i+1+nbasis, j+1+nbasis] = S_val
                        
                        # CORRECTION APPLIED HERE
                        # siesta2wtb.py shows that H[4] and H[5] are the imaginary parts
                        # of the diagonal blocks (up-up, dn-dn) for SOC.
                        # Non-colinear (NC) doesn't have H[4] or H[5].
                        
                        # sisl returns 4 elements for NC (H[0..3])
                        # and 8 elements for SOC (H[0..7])
                        if len(H) > 4:  # Spin-orbit (SOC) case
                            
                            # Up-up block (imaginary)
                            imH[idx, i+1, j+1] = H[4] # <- THIS LINE WAS ADDED/CORRECTED

                            # Up-down block (real and imaginary)
                            reH[idx, i+1, j+1+nbasis] = H[2]
                            imH[idx, i+1, j+1+nbasis] = H[3]
                            
                            # Down-up block (real and imaginary)
                            reH[idx, i+1+nbasis, j+1] = H[6]
                            imH[idx, i+1+nbasis, j+1] = H[7]
                            
                            # Down-down block (imaginary)
                            imH[idx, i+1+nbasis, j+1+nbasis] = H[5]

                        elif len(H) > 2:  # Non-colinear (NC) case
                            # (Has 4 elements, H[0], H[1], H[2], H[3])
                            
                            # Up-down block (real)
                            reH[idx, i+1, j+1+nbasis] = H[2]
                            
                            # Down-up block (real)
                            reH[idx, i+1+nbasis, j+1] = H[3]

            # Write to file
            for icell in range(1, ncell+1):
                for i in range(1, 2*nbasis+1):
                    for j in range(1, 2*nbasis+1):
                        line = (f"{format_float(a[icell,j])} {format_float(b[icell,j])} "
                                f"{format_float(c[icell,j])} {i} {j} "
                                f"{format_float(reH[icell,i,j])} "
                                f"{format_float(imH[icell,i,j])} "
                                f"{format_float(S[icell,i,j])}\n")
                        f.write(line)



def main():
    """Main workflow controller"""

    args = parse_arguments()
    
    # Hides the intro if --quiet is used
    # Note the quiet logic is inverted (action='store_false')
    # So args.quiet will be 'False' if -q is passed
    show_banner = args.quiet 
    
    if show_banner:
        show_intro()
        print("\n" + color_text("SIESTA-WANTIBEXOS Interface:", 'bold'))
        print("-"*60)
        print(color_text("[INFO] Initializing Hamiltonian processing...", 'yellow'))
    
    try:
        # Get system label (although not used directly by sisl.get_sile)
        system_label = get_system_label(args.input)

        # Load SIESTA Hamiltonian
        
        # CORRECTION 1: Use args.input to read geometry, not "calc.fdf"
        print(f"[INFO] Reading geometry from: {args.input}")
        geom = sisl.get_sile(args.input).read_geometry()
        print("[OK] Read geometry")
        
        print(f"[INFO] Reading Hamiltonian from: {args.input}")
        ham = sisl.get_sile(args.input).read_hamiltonian(geometry=geom)
        print("[OK] Read Hamiltonian from SIESTA")

        # CORRECTION 2: Pass the entire str(ham.spin) string to the function
        # The spin_mapping function is already prepared for this.
        spin_type_str = str(ham.spin)
        spin_suffix, spin_factor = spin_mapping(spin_type_str)
        print(f"[OK] Determine the spin type calculation: {spin_suffix} (Spin factor: {spin_factor})")
        
        # Get Fermi energy
        fermi = 0.0
        if args.fermi_level is not None:
            fermi = args.fermi_level
            print(f"[OK] Using manual Fermi Energy ( {fermi} eV )")
        elif args.output is not None:
            print(f"[INFO] Reading Fermi energy from: {args.output}")
            fermi = get_fermi_energy(args.output)
            print(f"[OK] Read the Fermi Energy ( {fermi} eV )")
        else:
            print(f"[WARN] No output file (-o) or manual Fermi level (-f) provided. Using Fermi = 0.0 eV")

        # Generate output files
        print("[INFO] Writing the basis... wait!")
        write_basis(ham, spin_suffix, spin_factor)
        print(f"[OK] Wrote the basis: basis_set-{spin_suffix}.basis")
        
        print("[INFO] Writing the Hamiltonian... wait!")
        # CORRECTION 3 is inside the write_hamiltonian function
        write_hamiltonian(ham, spin_suffix, fermi)
        print(f"[OK] Wrote the Hamiltonian: tb-{spin_suffix}.ham")
        
    except Exception as e:
        print(color_text(f"\nCritical error: {str(e)}", 'red'), file=sys.stderr)
        sys.exit(1) # Exit with error code

    if show_banner:
        print(color_text("\n[INFO] Successfully generated TB files.", 'green'))
    
    print("\n[INFO] Complete job!")
    print("\n"+"-"*60)
    if show_banner:
        print(color_text("Saving Hamiltonians before they collapse their own wavefunctions.\n\n", 'bold'))

if __name__ == "__main__":
    main()
