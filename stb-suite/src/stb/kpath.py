#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.8.0"

from time import sleep
from pymatgen.core import Structure, Lattice
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import sys
import os
import argparse 

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


# --- Color Class for UI ---
class Cores:
    """Class to store ANSI color codes for the terminal."""
    RESET = '\033[0m'
    BOLD = '\033[1m'
    RED = '\033[91m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    CYAN = '\033[96m'
# --- End of Color Class ---

def parse_fdf_to_structure(fdf_file="struct.fdf"):
    """
    Reads a SIESTA .fdf file and converts it into a 
    pymatgen Structure object.
    
    This parser is now more robust against comments.
    """
    
    try:
        with open(fdf_file, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"{Cores.RED}Error: File '{fdf_file}' not found.{Cores.RESET}")
        return None

    # Variables to store the data
    lattice_vectors = []
    species_map = {}
    atom_coords = []
    atom_species_names = []
    
    lattice_constant = 1.0 # Default
    coords_format = "Fractional" # Default
    
    # Flags for reading blocks
    in_species_block = False
    in_lattice_block = False
    in_coords_block = False
    
    for line in lines:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        
        # Convert line to lowercase for robust detection
        lower_line = line.lower()
        
        # Block detection
        if lower_line.startswith('%block chemicalspecieslabel'):
            in_species_block = True
            continue
        elif lower_line.startswith('%endblock chemicalspecieslabel'):
            in_species_block = False
            continue
        elif lower_line.startswith('%block latticevectors'):
            in_lattice_block = True
            continue
        elif lower_line.startswith('%endblock latticevectors'):
            in_lattice_block = False
            continue
        elif lower_line.startswith('%block atomiccoordinatesandatomicspecies'):
            in_coords_block = True
            continue
        elif lower_line.startswith('%endblock atomiccoordinatesandatomicspecies'):
            in_coords_block = False
            continue
            
        # Reading simple parameters
        if not (in_species_block or in_lattice_block or in_coords_block):
            if lower_line.startswith('latticeconstant'):
                lattice_constant = float(line.split()[1])
            if lower_line.startswith('atomiccoordinatesformat'):
                coords_format = line.split()[1]
            continue
            
        # --- PARSER (HANDLING COMMENTS) ---
        try:
            if in_species_block:
                parts = line.split()
                # Format: 1 50 Sn [optional comment]
                species_map[parts[0]] = parts[2] # Maps index (string '1') to name ('Sn')
                
            elif in_lattice_block:
                parts = line.split()
                # Format: 4.88... 0.0... 0.0... [optional comment]
                # Takes only the first 3 parts
                lattice_vectors.append([float(p) for p in parts[0:3]])
                
            elif in_coords_block:
                parts = line.split()
                if len(parts) < 4:
                    continue # Skip invalid lines
                
                # Format: x y z species_index [optional comment]
                coords = [float(parts[0]), float(parts[1]), float(parts[2])]
                atom_coords.append(coords)
                
                species_index = parts[3] # ex: '1'
                species_name = species_map.get(species_index) # ex: 'Sn'
                if species_name:
                    atom_species_names.append(species_name)
                else:
                    print(f"{Cores.RED}Error: Species index '{species_index}' not found in the map.{Cores.RESET}")
                    return None
                    
        except (ValueError, IndexError) as e:
            print(f"{Cores.YELLOW}Warning: Skipping malformed line in FDF: {line} ({e}){Cores.RESET}")
            continue
        # --- END OF PARSER ---

    # Check if the coordinate format is as expected
    is_cartesian = (coords_format.lower() != "fractional")

    # Create the Lattice object (assuming vectors are in Angstrom)
    if not lattice_vectors:
        print(f"{Cores.RED}Error: No lattice vectors found in file.{Cores.RESET}")
        return None
        
    lattice = Lattice(lattice_vectors)
    
    # Create the Structure object
    if not atom_species_names or not atom_coords:
        print(f"{Cores.RED}Error: No atoms found in file.{Cores.RESET}")
        return None
        
    structure = Structure(lattice, atom_species_names, atom_coords, coords_are_cartesian=is_cartesian)
    
    return structure


def write_siesta_kpath_file(kpoints_dict, path_segments, num_points=50, output_filename="kpath_bs.fdf"):
    """
    Writes the k-path to a SIESTA-formatted FDF file for band structure.
    
    This function now writes ALL suggested path segments, even if disjointed.
    """
    
    # 1. Determine the full sequence of k-points in order
    if not path_segments:
        print(f"  {Cores.YELLOW}Warning: No k-path segments found. File will not be written.{Cores.RESET}")
        return
        
    # --- START OF CORRECTION ---
    # Rebuilds the 'path_sequence' based on the correct data structure
    # (which is a list of paths, e.g: [['A', 'B', 'C'], ['D', 'E']])
    
    path_sequence = []
    
    for i, segment_list in enumerate(path_segments): # ex: segment_list = ['\Gamma', 'Y', 'H']
        if not segment_list:
            continue # Skip if segment is empty
            
        if i == 0:
            # For the first path, simply add all points
            path_sequence.extend(segment_list)
        else:
            # For subsequent paths
            last_point_in_sequence = path_sequence[-1]
            first_point_of_new_segment = segment_list[0]
            
            if last_point_in_sequence == first_point_of_new_segment:
                # The path is continuous (e.g., ends in Z, starts in Z)
                # Add all points, *except* the first one (which is duplicated)
                path_sequence.extend(segment_list[1:])
            else:
                # The path is disjoint (it's a "jump", e.g., ends in H_1, starts in M)
                print(f"  {Cores.CYAN}Note: Disjointed path detected (jump from {last_point_in_sequence} to {first_point_of_new_segment}). Including full segment.{Cores.RESET}")
                # Add ALL points from the new segment
                path_sequence.extend(segment_list)
    
    # --- END OF CORRECTION ---
    
    if not path_sequence:
        print(f"  {Cores.RED}Error: Could not determine path sequence.{Cores.RESET}")
        return
        
    # This line will now print the full path, including jumps
    path_str_display = '-'.join([r'\Gamma' if p == 'GAMMA' else p for p in path_sequence])
    print(f"  {Cores.CYAN}SIESTA Path to be written:{Cores.RESET} {path_str_display}")

    try:
        with open(output_filename, 'w') as f:
            f.write("### BANDS\n")
            f.write(" BandLinesScale  ReciprocalLatticeVectors\n\n")
            f.write("%block BandLines\n")
            
            # 2. Write the first point (with '1' point)
            first_label = path_sequence[0]
            first_coords = kpoints_dict[first_label]
            coord_str = " ".join([f"{c:14.10f}" for c in first_coords])
            # Fix 'GAMMA' to '\Gamma' label in file
            f_label = r'\Gamma' if first_label == 'GAMMA' else first_label
            f.write(f"1   {coord_str}   {f_label}\n")
            
            # 3. Write the rest of the points
            for label in path_sequence[1:]:
                coords = kpoints_dict[label]
                coord_str = " ".join([f"{c:14.10f}" for c in coords])
                # Fix 'GAMMA' to '\Gamma' label in file
                f_label = r'\Gamma' if label == 'GAMMA' else label
                # Use the number of points (ex: 50) for all following segments
                f.write(f"{num_points}   {coord_str}   {f_label}\n")
                
            f.write("%endblock BandLines\n")
        
        print(f"  {Cores.GREEN}Success:{Cores.RESET} SIESTA file '{Cores.BOLD}{output_filename}{Cores.RESET}' has been created.")
    
    except KeyError as e:
        # This error happens if a label in the path (ex: '\Gamma') has a 
        # different name in the dictionary (ex: 'GAMMA')
        
        # Tries to substitute '\Gamma' with 'GAMMA' if 'GAMMA' exists
        if str(e) == r"'\Gamma'" and 'GAMMA' in kpoints_dict:
            print(f"  {Cores.YELLOW}Warning: Found '\\Gamma' label, attempting to use 'GAMMA' internally.{Cores.RESET}")
            # Recreate the sequence, substituting \Gamma with GAMMA
            path_sequence_fixed = ['GAMMA' if label == r'\Gamma' else label for label in path_sequence]
            # Try to write the file AGAIN
            write_siesta_kpath_file_fixed(kpoints_dict, path_sequence_fixed, num_points, output_filename)
        else:
             print(f"  {Cores.RED}Error writing file: K-point label {e} found in path but not in k-points list.{Cores.RESET}")
             print(f"  {Cores.RED}Available labels: {list(kpoints_dict.keys())}{Cores.RESET}")

    except Exception as e:
        print(f"  {Cores.RED}Error writing {output_filename}: {e}{Cores.RESET}")

# --- HELPER FUNCTION TO FIX LABEL NAMES ---
def write_siesta_kpath_file_fixed(kpoints_dict, path_sequence, num_points, output_filename):
    """
    Emergency "helper" function to write the file if the KeyError
    for '\Gamma' vs 'GAMMA' occurs.
    """
    print(f"  {Cores.CYAN}Retrying file write with 'GAMMA' label...{Cores.RESET}")
    try:
        with open(output_filename, 'w') as f:
            f.write("### BANDS\n")
            f.write(" BandLinesScale  ReciprocalLatticeVectors\n\n")
            f.write("%block BandLines\n")
            
            first_label = path_sequence[0]
            first_coords = kpoints_dict[first_label]
            coord_str = " ".join([f"{c:14.10f}" for c in first_coords])
            f_label = r'\Gamma' if first_label == 'GAMMA' else first_label
            f.write(f"1   {coord_str}   {f_label}\n")
            
            for label in path_sequence[1:]:
                coords = kpoints_dict[label]
                coord_str = " ".join([f"{c:14.10f}" for c in coords])
                f_label = r'\Gamma' if label == 'GAMMA' else label
                f.write(f"{num_points}   {coord_str}   {f_label}\n")
                
            f.write("%endblock BandLines\n")
        
        print(f"  {Cores.GREEN}Success:{Cores.RESET} SIESTA file '{Cores.BOLD}{output_filename}{Cores.RESET}' has been created (with label fix).")
    except Exception as e:
        print(f"  {Cores.RED}Error during retry: {e}{Cores.RESET}")
# --- END OF SIESTA WRITE FUNCTIONS ---


# --- NEW FUNCTION FOR VASP KPOINTS ---
def write_vasp_kpoints_file(kpoints_dict, path_segments, num_points=50, output_filename="KPOINTS"):
    """
    Writes the k-path to a VASP-formatted KPOINTS file for band structure.
    
    This function handles disjointed paths by adding a blank line.
    """
    
    # --- Fix \Gamma vs GAMMA issue ---
    # Create a copy of the dictionary and paths to avoid modifying the original
    kpoints_dict_fixed = kpoints_dict.copy()
    
    # Pymatgen sometimes uses '\Gamma' in paths but 'GAMMA' in the kpoints dict.
    # We replace '\Gamma' with 'GAMMA' everywhere to be consistent.
    if r'\Gamma' in kpoints_dict and 'GAMMA' not in kpoints_dict:
        kpoints_dict_fixed['GAMMA'] = kpoints_dict_fixed.pop(r'\Gamma')
        print(f"  {Cores.CYAN}Note: Standardizing '\\Gamma' to 'GAMMA' internally.{Cores.RESET}")
    elif r'\Gamma' not in kpoints_dict and 'GAMMA' in kpoints_dict:
         pass # Already in 'GAMMA' format
    elif r'\Gamma' in kpoints_dict and 'GAMMA' in kpoints_dict:
         # This case is ambiguous, but we'll prefer GAMMA
         kpoints_dict_fixed.pop(r'\Gamma', None)
         print(f"  {Cores.YELLOW}Warning: Both '\\Gamma' and 'GAMMA' found, using 'GAMMA'.{Cores.RESET}")

    # Fix the path segments to use 'GAMMA'
    path_segments_fixed = []
    for segment in path_segments:
        path_segments_fixed.append(['GAMMA' if label == r'\Gamma' else label for label in segment])
    # --- End of fix ---

    try:
        with open(output_filename, 'w') as f:
            f.write(f"K-Path for Band Structure calculation\n")
            f.write(f" {num_points}\n")
            f.write("Line-mode\n")
            f.write("Reciprocal\n")
            
            # path_segments is like [['A', 'B', 'C'], ['D', 'E']]
            for i, segment_list in enumerate(path_segments_fixed):
                
                # Add a blank line *between* disjointed paths
                # (e.g., between C and D)
                if i > 0:
                    f.write("\n") 
                    # Fix 'GAMMA' to '\Gamma' for display
                    label_from = r'\Gamma' if path_segments_fixed[i-1][-1] == 'GAMMA' else path_segments_fixed[i-1][-1]
                    label_to = r'\Gamma' if segment_list[0] == 'GAMMA' else segment_list[0]
                    print(f"  {Cores.CYAN}Note: Disjointed path detected (jump from {label_from} to {label_to}).{Cores.RESET}")

                # Write all segments in the current continuous path
                # (e.g., A-B, then B-C)
                for j in range(len(segment_list) - 1):
                    start_label = segment_list[j]
                    end_label = segment_list[j+1]
                    
                    start_coords = kpoints_dict_fixed[start_label]
                    end_coords = kpoints_dict_fixed[end_label]
                    
                    # Fix 'GAMMA' to '\Gamma' label in file
                    f_start_label = r'\Gamma' if start_label == 'GAMMA' else start_label
                    f_end_label = r'\Gamma' if end_label == 'GAMMA' else end_label
                    
                    # Write start point
                    f.write(f"  {start_coords[0]:.10f}  {start_coords[1]:.10f}  {start_coords[2]:.10f}   ! {f_start_label}\n")
                    # Write end point
                    f.write(f"  {end_coords[0]:.10f}  {end_coords[1]:.10f}  {end_coords[2]:.10f}   ! {f_end_label}\n")
                    
                    # Add a blank line to separate this segment (e.g., A-B)
                    # from the next (e.g., B-C)
                    f.write("\n")

        print(f"  {Cores.GREEN}Success:{Cores.RESET} VASP file '{Cores.BOLD}{output_filename}{Cores.RESET}' has been created.")

    except KeyError as e:
         print(f"  {Cores.RED}Error writing file: K-point label {e} found in path but not in k-points list.{Cores.RESET}")
         print(f"  {Cores.RED}Available labels: {list(kpoints_dict_fixed.keys())}{Cores.RESET}")
    except Exception as e:
        print(f"  {Cores.RED}Error writing {output_filename}: {e}{Cores.RESET}")
# --- END OF VASP WRITE FUNCTION ---


def get_kpath_from_structure(structure, symprec=0.01, write_fdf_file=False, write_kpoints_file=False): 
    """
    Takes a Pymatgen Structure object and prints
    its space group and high-symmetry k-path.
    
    If write_fdf_file is True, it also writes the k-path to 'kpath_bs.fdf'.
    If write_kpoints_file is True, it also writes the k-path to 'KPOINTS'.
    """
    try:
        
        # We pass the structure to HighSymmKpath, which finds the
        # primitive cell itself.
        kpath = HighSymmKpath(structure, symprec=symprec) 
        
        # Use SpacegroupAnalyzer on the primitive cell (kpath.prim)
        analyzer = SpacegroupAnalyzer(kpath.prim)
        
        space_group_symbol = analyzer.get_space_group_symbol()
        space_group_number = analyzer.get_space_group_number()
        bravais_lattice = analyzer.get_lattice_type()

        # --- IMPROVED OUTPUT ---
        print(f"{Cores.BOLD}--- 1. Structure Analysis ---{Cores.RESET}")
        print(f"  {Cores.CYAN}Formula:{Cores.RESET} {structure.composition.reduced_formula}")
        print(f"  {Cores.CYAN}Space Group:{Cores.RESET} {Cores.BOLD}{space_group_symbol} (No. {space_group_number}){Cores.RESET}")
        print(f"  {Cores.CYAN}Bravais Lattice:{Cores.RESET} {bravais_lattice}")
        print(f"  {Cores.YELLOW}Using precision (symprec):{Cores.RESET} {symprec}")


        print(f"\n{Cores.BOLD}--- 2. High-Symmetry K-Points (Fractional Coordinates) ---{Cores.RESET}")
        kpoints = kpath.kpath['kpoints']
        for label, coords in kpoints.items():
            # Format coordinates
            coord_str = ", ".join([f"{c:8.5f}" for c in coords])
            # Fix 'GAMMA' to '\Gamma' for display
            display_label = r'\Gamma' if label == 'GAMMA' else label
            print(f"  {Cores.GREEN}{display_label:<5}:{Cores.RESET} ({coord_str})")

        print(f"\n{Cores.BOLD}--- 3. Suggested K-Path ---{Cores.RESET}")
        path = kpath.kpath['path']
        path_segments_str_list = []
        
        # Transform the list of lists into a readable string
        # ex: [['A', 'B'], ['C', 'D']] -> "A-B | C-D"
        for segment_list in path:
             # Fix 'GAMMA' to '\Gamma' for display
            display_segment = [r'\Gamma' if label == 'GAMMA' else label for label in segment_list]
            path_segments_str_list.append("-".join(display_segment))
        
        print(f"  Path: {Cores.GREEN}{Cores.BOLD}{' | '.join(path_segments_str_list)}{Cores.RESET}")
        
        # --- ADDED REFERENCE ---
        print(f"\n{Cores.CYAN}Methodology:{Cores.RESET}")
        print(f"  The path follows the Setyawan & Curtarolo (2010) convention,")
        print(f"  as implemented by the Pymatgen library.")
        print(f"  {Cores.YELLOW}Reference:{Cores.RESET}")
        print(f"    Setyawan, W., & Curtarolo, S. (2010).")
        print(f"    High-throughput electronic band structure calculations: Challenges and tools.")
        print(f"    Computational Materials Science, 49(2), 299-312.")
        print(f"    DOI: 10.1016/j.commatsci.2010.05.010")
        # --- END OF REFERENCE ---
        
        num_points_per_segment = 50 
        
        # --- FILE GENERATION ---
        print(f"\n{Cores.BOLD}--- 4. Output File Generation ---{Cores.RESET}")
        
        # Write the SIESTA k-path file if requested
        if write_fdf_file:
            write_siesta_kpath_file(
                kpath.kpath['kpoints'], 
                kpath.kpath['path'], # Pass the original list of lists
                num_points_per_segment, 
                "kpath_bs.fdf" # Output filename
            )
        
        # Write the VASP k-path file if requested
        elif write_kpoints_file:
            write_vasp_kpoints_file(
                kpath.kpath['kpoints'], 
                kpath.kpath['path'], # Pass the original list of lists
                num_points_per_segment, 
                "KPOINTS" # Output filename
            )
        # --- END OF GENERATION ---

    except Exception as e:
        print(f"\n{Cores.RED}An error occurred while analyzing the structure: {e}{Cores.RESET}")
        print(f"{Cores.YELLOW}It might be necessary to adjust the 'symprec' parameter (current: {symprec}).{Cores.RESET}")
        print(f"{Cores.YELLOW}This often happens if the input structure is distorted or the FDF/POSCAR parse failed.{Cores.RESET}")


def main():
    
    # 1. Configure ArgParse
    parser = argparse.ArgumentParser(
        description=f"""{Cores.BOLD}Finds the high-symmetry k-path for FDF (SIESTA) or POSCAR (VASP) files.{Cores.RESET}
Uses the Setyawan & Curtarolo (2010) methodology via Pymatgen.""",
        formatter_class=argparse.RawTextHelpFormatter 
    )
    
    # Argument 1: Filename (now a required flag)
    parser.add_argument(
        "-f", "--file",
        dest="filename", # Stores the value in 'args.filename'
        type=str,
        required=True,  # Makes this flag mandatory
        help="The structure file name (e.g., struct.fdf, POSCAR, CONTCAR)."
    )
    
    # Argument 2: File type (optional)
    parser.add_argument(
        "-t", "--type",
        type=str,
        choices=['fdf', 'poscar'], 
        default=None, 
        help="Specifies the file type.\n"
             "  'fdf':    For SIESTA FDF files.\n"
             "  'poscar': For VASP POSCAR/CONTCAR files.\n"
             "If omitted, the script will try to guess from the filename."
    )
    
    # Argument 3: Precision (optional)
    parser.add_argument(
        "-p", "--prec",
        type=float,
        default=0.01,
        help="Sets the precision (symprec) for symmetry analysis.\n"
             "Default: 0.01"
    )
    
    parser.add_argument("-v", "--version", action="version",
                        version=f"stb-kpath {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()


    if args.intro == True:
        show_intro()

    print("\n" + color_text("Suggested k-path from structure (only bulk):", 'bold'))
    print("-"*60)


    # Read the arguments provided on the command line
    args = parser.parse_args()
    
    # 2. Process the arguments
    
    filename = args.filename
    precision = args.prec
    file_type = args.type
    
    # Check if the file exists
    if not os.path.exists(filename):
        print(f"{Cores.RED}Error: File '{filename}' not found.{Cores.RESET}")
        sys.exit(1) 

    # Try to guess the type if not provided
    if file_type is None:
        print(f"{Cores.YELLOW}File type not specified, trying to guess...{Cores.RESET}")
        if filename.lower().endswith('.fdf'):
            file_type = 'fdf'
        # Common VASP filenames
        elif 'poscar' in filename.lower() or 'contcar' in filename.lower() or filename.lower().endswith('.vasp'):
            file_type = 'poscar'
        else:
            print(f"{Cores.RED}Error: Could not guess file type for '{filename}'.{Cores.RESET}")
            print(f"Please use the -t 'fdf' or -t 'poscar' option.")
            sys.exit(1)
            
    print(f"{Cores.CYAN}Reading file: {filename} (Type: {file_type}){Cores.RESET}\n")
    
    # 3. Load the Structure object
    
    structure_obj = None
    
    try:
        if file_type == 'fdf':
            # Use your parsing function for FDF
            structure_obj = parse_fdf_to_structure(filename)
            
        elif file_type == 'poscar':
            # Use Pymatgen's native reader for POSCAR
            structure_obj = Structure.from_file(filename)
            
    except Exception as e:
        print(f"{Cores.RED}An error occurred while READING the file '{filename}':{Cores.RESET}")
        print(f"{e}")
        sys.exit(1)
        
    
    # 4. Run the analysis (if loading was successful)
    if structure_obj:
        
        # Determine which output file to write
        create_fdf_output = (file_type == 'fdf')
        create_kpoints_output = (file_type == 'poscar')
        
        # Get the k-path from the Structure object, passing the flags
        get_kpath_from_structure(
            structure_obj, 
            symprec=precision, 
            write_fdf_file=create_fdf_output,
            write_kpoints_file=create_kpoints_output
        )
        
    else:
        print(f"{Cores.RED}Error: Could not create the structure object from the file.{Cores.RESET}")
        
if __name__ == "__main__":
    main()
