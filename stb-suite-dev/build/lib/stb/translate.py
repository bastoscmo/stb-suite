#!/usr/bin/env python

#################################################
# Siesta Tool Box - Translate                   #
# version 1.0.1                                 #
# UnB - 2025/02/05                              #
# Dr. Carlos M. O. Bastos                       #
#################################################

VERSION = "1.5.10"

import os
import sys
import warnings
import subprocess
from time import sleep
import argparse
import textwrap
from typing import List, Dict
import argparse
from pymatgen.core.periodic_table import Element
from pymatgen.core import SiteCollection, Structure
import numpy as np


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

# Import libraries:

# Supported formats
INPUT_FORMATS = {"poscar", "cif", "siesta", "xyz", "fhi", "dftb", "xsf"}
OUTPUT_FORMATS = {"cif","xyz", "poscar", "fdf", "dftb", "xsf", "fhi"}


# Dictionary with all periodic elements table
def periodic_table():
    element_atomicnumber = {
        "H": 1,     # Hidrogênio
        "He": 2,    # Hélio
        "Li": 3,    # Lítio
        "Be": 4,    # Berílio
        "B": 5,     # Boro
        "C": 6,     # Carbono
        "N": 7,     # Nitrogênio
        "O": 8,     # Oxigênio
        "F": 9,     # Flúor
        "Ne": 10,   # Neônio
        "Na": 11,   # Sódio
        "Mg": 12,   # Magnésio
        "Al": 13,   # Alumínio
        "Si": 14,   # Silício
        "P": 15,    # Fósforo
        "S": 16,    # Enxofre
        "Cl": 17,   # Cloro
        "Ar": 18,   # Argônio
        "K": 19,    # Potássio
        "Ca": 20,   # Cálcio
        "Sc": 21,   # Escândio
        "Ti": 22,   # Titânio
        "V": 23,    # Vanádio
        "Cr": 24,   # Cromo
        "Mn": 25,   # Manganês
        "Fe": 26,   # Ferro
        "Co": 27,   # Cobalto
        "Ni": 28,   # Níquel
        "Cu": 29,   # Cobre
        "Zn": 30,   # Zinco
        "Ga": 31,   # Gálio
        "Ge": 32,   # Germanio
        "As": 33,   # Arsênio
        "Se": 34,   # Selênio
        "Br": 35,   # Bromo
        "Kr": 36,   # Kriptônio
        "Rb": 37,   # Rubídio
        "Sr": 38,   # Estrôncio
        "Y": 39,    # Ítrio
        "Zr": 40,   # Zircônio
        "Nb": 41,   # Nióbio
        "Mo": 42,   # Molibdênio
        "Tc": 43,   # Tecnécio
        "Ru": 44,   # Rutênio
        "Rh": 45,   # Ródio
        "Pd": 46,   # Paládio
        "Ag": 47,   # Prata
        "Cd": 48,   # Cádmio
        "In": 49,   # Índio
        "Sn": 50,   # Estanho
        "Sb": 51,   # Antimônio
        "Te": 52,   # Telúrio
        "I": 53,    # Iodo
        "Xe": 54,   # Xenônio
        "Cs": 55,   # Césio
        "Ba": 56,   # Bário
        "La": 57,   # Lantânio
        "Ce": 58,   # Cério
        "Pr": 59,   # Praseodímio
        "Nd": 60,   # Neodímio
        "Pm": 61,   # Promécio
        "Sm": 62,   # Samário
        "Eu": 63,   # Európio
        "Gd": 64,   # Gadolínio
        "Tb": 65,   # Térbio
        "Dy": 66,   # Disprósio
        "Ho": 67,   # Holmium
        "Er": 68,   # Érbio
        "Tm": 69,   # Túlio
        "Yb": 70,   # Itérbio
        "Lu": 71,   # Lutécio
        "Hf": 72,   # Háfnio
        "Ta": 73,   # Tântalo
        "W": 74,    # Wolframio
        "Re": 75,   # Rênio
        "Os": 76,   # Ósmio
        "Ir": 77,   # Irídio
        "Pt": 78,   # Platina
        "Au": 79,   # Ouro
        "Hg": 80,   # Mercúrio
        "Tl": 81,   # Tálio
        "Pb": 82,   # Chumbo
        "Bi": 83,   # Bismuto
        "Po": 84,   # Polônio
        "At": 85,   # Astato
        "Rn": 86,   # Radônio
        "Fr": 87,   # Frâncio
        "Ra": 88,   # Radônio
        "Ac": 89,   # Actínio
        "Th": 90,   # Tório
        "Pa": 91,   # Protactínio
        "U": 92,    # Urânio
        "Np": 93,   # Netúnio
        "Pu": 94,   # Plutônio
        "Am": 95,   # Amerício
        "Cm": 96,   # Curió
        "Bk": 97,   # Berquélio
        "Cf": 98,   # Califórnio
        "Es": 99,   # Einstênio
        "Fm": 100,  # Férmio
        "Md": 101,  # Mendelevio
        "No": 102,  # Nobelio
        "Lr": 103,  # Laurêncio
        "Rf": 104,  # Rutherfórdio
        "Db": 105,  # Dúbnio
        "Sg": 106,  # Seabórgio
        "Bh": 107,  # Bóhrio
        "Hs": 108,  # Hassio
        "Mt": 109,  # Meitnério
        "Ds": 110,  # Darmstádio
        "Rg": 111,  # Roentgênio
        "Cn": 112,  # Copernício
        "Nh": 113,  # Nihônio
        "Fl": 114,  # Fleróvio
        "Mc": 115,  # Moscóvio
        "Lv": 116,  # Livermório
        "Ts": 117,  # Tenesso
        "Og": 118   # Oganessônio
    }

    # Dictionary with all atomic numbers table (invert Dictionary)
    atomicnumber_element = {str(v): k for k, v in element_atomicnumber.items()}
    return element_atomicnumber, atomicnumber_element


def readfile(filedata):
    with open(filedata, 'r') as fil:
        data = [line.split() for line in fil if line.strip()
                ]
    return data


def dic_atoms_position(atomsposition):
    dic_atomspos = {}
    for elt in atomsposition:
        if elt[0] not in dic_atomspos:
            dic_atomspos[elt[0]] = []
        dic_atomspos[elt[0]].append([elt[1], elt[2], elt[3]])
    return dic_atomspos


############# Extract Data Functions ######################

# This Function define the number of atoms types for xyz file
def getatomsandvectors_xyz(dataxyz, latticedata):
    element, atomicnumber = periodic_table()
    xyz = readfile(dataxyz)[2:]
    lattice = readfile(latticedata)
    # Lattice
    latticeparameter = lattice[0][0]
    typevectors = 'Cartesian'
    vectors = [lattice[1], lattice[2], lattice[3]]
    # Atoms position
    atoms = []
    j = 1
    getatoms = []
    for i in range(len(xyz)):
        if xyz[i][0] not in atoms:
            atoms.append(xyz[i][0])
            getatoms.append([j, element[xyz[i][0]], xyz[i][0]])
            j = j+1
    atomsposition = xyz
    for elem in getatoms:
        cont = 0
        for elxyz in xyz:
            if elem[2] == elxyz[0]:
                cont = cont+1
        elem.append(str(cont))
    atomic_position = dic_atoms_position(atomsposition)
    return typevectors, latticeparameter, vectors, getatoms, atomic_position


# This Function define the number of atoms types for vasp file
def getatomsandvectors_vasp(poscar):
    element, atomicnumber = periodic_table()
    datavasp = readfile(poscar)
    latticeparameter = datavasp[1][0]
    vectors = datavasp[2:5]
    typevectors = datavasp[7][0]
    getatoms = []
    for i in range(len(datavasp[5])):
        getatoms.append([i+1, element[datavasp[5][i]],datavasp[5][i], datavasp[6][i]])
    atomsposition = []
    cont = 8 
    for el in getatoms:
        for i in range(int(el[3])):
            atomsposition.append([el[2], datavasp[cont][0], datavasp[cont][1], datavasp[cont][2]])
            cont = cont+1
    atomic_position = dic_atoms_position(atomsposition)
    return typevectors, latticeparameter, vectors, getatoms, atomic_position


# This Function define the data for cif file
def getatomsandvectors_cif(input_cif):
    element, atomicnumber = periodic_table()
    structure = Structure.from_file(input_cif)
    atom_data = [(i+1, site.species_string, site.coords)
                 for i, site in enumerate(structure.sites)]
    getatoms = []
    dicatoms = {}
    icont = 1
    for atom in atom_data:
        if atom[1] not in dicatoms:
            dicatoms[atom[1]] = []
        dicatoms[atom[1]].append([f"{atom[2][0]:.8f}",
                                  f"{atom[2][1]:.8f}",
                                  f"{atom[2][2]:.8f}"
                                  ])
    for line in dicatoms:
        getatoms.append([f"{icont}",
                         f"{element[line]}",
                         f"{line}",
                         f"{len(dicatoms[line])}"])
    atomic_position = dicatoms
    typevectors = 'Cartesian'
    latticeparameter = '1.00'
    vectors = []
    for line in structure.lattice.matrix:
        vectors.append([f"{line[0]:.8f}", f"{line[1]:.8f}", f"{line[2]:.8f}"])
    coord = structure.cart_coords
    return typevectors, latticeparameter, vectors, getatoms, atomic_position


# This Function define the data for cif file
def getatomsandvectors_fhi(input_fhi):
    element, atomicnumber = periodic_table()
    datafhi = readfile(input_fhi)
    vectors = []
    atomdata = {}
    latticeparameter = "1.00"
    for lines in datafhi:
        if lines[0] == "atom_frac":
            atomdata.setdefault(lines[4], [])
            atomdata[lines[4]].append([
                f"{float(lines[1]):.8f}",
                f"{float(lines[2]):.8f}",
                f"{float(lines[3]):.8f}"])
            typevectors = 'Direct'
        elif lines[0] == "atom":
            atomdata.setdefault(lines[4], [])
            atomdata[lines[4]].append([
                f"{float(lines[1]):.8f}",
                f"{float(lines[2]):.8f}",
                f"{float(lines[3]):.8f}"])
            typevectors = 'Cartesian'
        elif lines[0] == 'lattice_vector':
            vectors.append([lines[1], lines[2], lines[3]])
    getatoms = []
    icont = 1
    for el in atomdata:
        getatoms.append([str(icont), str(element[el]),
                        el, str(len(atomdata[el]))])
        icont = icont+1
    atomic_position = atomdata

    return typevectors, latticeparameter, vectors, getatoms, atomic_position


def getatomsandvectors_siesta(input_siesta):
    element, atomicnumber = periodic_table()
    datasiesta = readfile(input_siesta)
    latticeparameter = "1.00"
    typevectors = "Direct"
    vectors = datasiesta[:3]
    dicatoms = {}
    for pos in datasiesta[4:]:
        if pos[0] not in dicatoms:
            dicatoms[pos[0]] = []
        dicatoms[pos[0]].append([pos[0], pos[1], atomicnumber[pos[1]]])
    getatoms = []
    for atoms in dicatoms:
        getatoms.append([dicatoms[atoms[0]][0][0],
                        dicatoms[atoms[0]][0][1],
                        dicatoms[atoms[0]][0][2],
                        str(len(dicatoms[atoms[0]]))])
    atomic_position = {}
    for position in datasiesta[4:]:
        if atomicnumber[position[1]] not in atomic_position:
            atomic_position[atomicnumber[position[1]]] = []
        atomic_position[atomicnumber[position[1]]].append(
            [position[2], position[3], position[4]])
    return typevectors, latticeparameter, vectors, getatoms, atomic_position


def getatomsandvectors_dftb(input_dftb):
    element, atomicnumber = periodic_table()
    datadftb = readfile(input_dftb)
    vectors = []
    getatoms = []
    atomic_position = {}
    dic_data = {}
    latticeparameter = "1.00"
    if datadftb[0][1] == 'S':
        typevectors = 'Cartesian'
    elif datadftb[0][1] == 'F':
        typevectors = 'Direct'
    else:
        print("[FAIL] Type of coordinate not define: Only S or F are accepted")
        exit()
    for i in range(1, 4):
        vectors.append([f"{float(datadftb[-i][0]):.8f}",
                        f"{float(datadftb[-i][1]):.8f}",
                        f"{float(datadftb[-i][2]):.8f}"])
    icont = 0
    dic_data = {f"{i + 1}": elem for i, elem in enumerate(datadftb[1])}
    for line in datadftb[2:-4]:
        atomic_position.setdefault(dic_data[line[1]], [])
        atomic_position[dic_data[line[1]]].append(
            [f"{float(line[2]):.8f}", f"{float(line[3]):.8f}", f"{float(line[4]):.8f}"])
    icont = 1
    for el in atomic_position:
        getatoms.append([f"{icont}",
                         f"{element[el]}",
                         f"{el}",
                         f"{len(atomic_position[el])}"])
        icont = icont+1

    return typevectors, latticeparameter, vectors, getatoms, atomic_position


def getatomsandvectors_xsf(input_xsf):
    print("[WARNING] Only for PRIMVEC format")
    element, atomicnumber = periodic_table()
    dataxsf = readfile(input_xsf)
    vectors = []
    getatoms = []
    atomic_position = {}
    dic_data = {}
    latticeparameter = "1.00"
    typevectors = 'Cartesian'
    dataxsf = [v for v in dataxsf if not str(v[0]).startswith("#")]
    for j in range(len(dataxsf)):
        if dataxsf[j][0] == 'PRIMVEC':
            for i in range(1, 4):
                vectors.append([f"{float(dataxsf[j+i][0]):.8f}",
                                f"{float(dataxsf[j+i][1]):.8f}",
                                f"{float(dataxsf[j+i][2]):.8f}"])
    for j in range(len(dataxsf)):
        if dataxsf[j][0] == 'PRIMCOORD':
            na = dataxsf[j+1][0]
            for i in range(int(na)):
                atomic_position.setdefault(atomicnumber[dataxsf[j+i+2][0]], [])
                atomic_position[atomicnumber[dataxsf[j+i+2][0]]].append(
                    [f"{float(dataxsf[j+i+2][1]):.8f}",
                     f"{float(dataxsf[j+i+2][2]):.8f}",
                     f"{float(dataxsf[j+i+2][3]):.8f}"])
    icont = 1
    for el in atomic_position:
        getatoms.append([f"{icont}",
                         f"{element[el]}",
                         f"{el}",
                         f"{len(atomic_position[el])}"])
        icont = icont+1

    return typevectors, latticeparameter, vectors, getatoms, atomic_position

###################### Write functions ################


def writefilefdf(typevectors, latticeparameter, vectors, getatoms, atomsposition, outfilename):
    numberofatoms = 0
    for lin in getatoms:
        numberofatoms = numberofatoms+int(lin[3])
    outfile = []
    outfile.append(
        '# automatic create translation file using sstranslate\n\n')
    outfile.append(f"NumberOfSpecies    {len(getatoms)}")
    outfile.append(f"NumberofAtoms      {numberofatoms}\n\n")
    outfile.append("%block ChemicalSpeciesLabel")
    for atoms in getatoms:
        outfile.append(f" {atoms[0]}   {atoms[1]}   {atoms[2]}")
    outfile.append("%endblock ChemicalSpeciesLabel \n")
    outfile.append(f"LatticeConstant {latticeparameter} Ang \n")
    if typevectors == 'Direct':
        outfile.append("AtomicCoordinatesFormat  Fractional \n\n")
    if typevectors == 'Cartesian':
        outfile.append("AtomicCoordinatesFormat  Ang\n\n")
    outfile.append("%block LatticeVectors")
    for lin in vectors:
        outfile.append(f" {lin[0]}   {lin[1]}   {lin[2]} ")
    outfile.append("%endblock LatticeVectors\n\n")
    outfile.append("%block AtomicCoordinatesAndAtomicSpecies")
    for elem in getatoms:
        for position in atomsposition[elem[2]]:
            outfile.append(
                f"  {position[0]}   {position[1]}   {position[2]}   {elem[0]}  ")
    outfile.append("%endblock AtomicCoordinatesAndAtomicSpecies")
    np.savetxt(outfilename, outfile, fmt='%s')
    return


def writefileposcar(typevectors, latticeparameter, vectors, getatoms, atomsposition, outfilename):
    outfile = []
    outfile.append(
        '# automatic create using sstranslate')
    outfile.append(f"{latticeparameter}")
    for lin in vectors:
        outfile.append(f"{lin[0]}   {lin[1]}   {lin[2]} ")
    lineatoms = ""
    linenatoms = ""
    for atoms in getatoms:
        lineatoms = lineatoms + f"{atoms[2]}   "
        linenatoms = linenatoms+f"{atoms[3]}   "
    outfile.append(f"{lineatoms}")
    outfile.append(f"{linenatoms}")
    if typevectors == 'Direct':
        outfile.append("Direct")
    if typevectors == 'Cartesian':
        outfile.append("Cartesian")
    for elem in getatoms:
        for position in atomsposition[elem[2]]:
            outfile.append(
                f"{position[0]}   {position[1]}   {position[2]}")
    np.savetxt(outfilename, outfile, fmt='%s')
    return

    # Calcula ângulos entre vetores
def angle(u, v):
        cos_theta = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
        return np.degrees(np.arccos(np.clip(cos_theta, -1.0, 1.0)))
        

def writefilecif(typevectors, latticeparameter, vectors, getatoms, atomsposition, outfilename):
    """
    Escreve um arquivo CIF no formato padrão.
    """
    vectors = np.array(vectors, dtype=float) * float(latticeparameter)

    outfile = []

    outfile.append("data_generated")
    outfile.append("_symmetry_space_group_name_H-M   'P 1'")
    outfile.append("_symmetry_Int_Tables_number      1")
    outfile.append("_cell_length_a    {:.8f}".format(np.linalg.norm(vectors[0])))
    outfile.append("_cell_length_b    {:.8f}".format(np.linalg.norm(vectors[1])))
    outfile.append("_cell_length_c    {:.8f}".format(np.linalg.norm(vectors[2])))

    alpha = angle(vectors[1], vectors[2])
    beta = angle(vectors[0], vectors[2])
    gamma = angle(vectors[0], vectors[1])

    outfile.append("_cell_angle_alpha  {:.8f}".format(alpha))
    outfile.append("_cell_angle_beta   {:.8f}".format(beta))
    outfile.append("_cell_angle_gamma  {:.8f}".format(gamma))
    outfile.append(" ")

    outfile.append("loop_")
    outfile.append("_symmetry_equiv_pos_as_xyz")
    outfile.append("  'x, y, z'")
    outfile.append(" ")

    outfile.append("loop_")
    outfile.append("_atom_site_label")
    outfile.append("_atom_site_type_symbol")
    outfile.append("_atom_site_fract_x")
    outfile.append("_atom_site_fract_y")
    outfile.append("_atom_site_fract_z")

    # Conversão de coordenadas cartesianas para fracionárias, se necessário
    inv_lattice = np.linalg.inv(vectors)

    for elem in getatoms:
        for pos in atomsposition[elem[2]]:
            if typevectors.lower() == 'direct':
                xf, yf, zf = map(float, pos)
            elif typevectors.lower() == 'cartesian':
                cart = np.array([float(pos[0]), float(pos[1]), float(pos[2])])
                fract = np.dot(inv_lattice, cart)
                xf, yf, zf = fract
            outfile.append(
                f"{elem[2]}   {elem[2]}   {xf:.8f}   {yf:.8f}   {zf:.8f}"
            )

    np.savetxt(outfilename, outfile, fmt='%s')
    return




def writefilexyz(typevectors, latticeparameter, vectors, getatoms, atomsposition, outfilename):
    outfile = []
    sum = 0
    comment = ""
    for el in getatoms:
        comment = comment+f"{el[2]}{el[3]} "
        sum = sum+int(el[3])
    outfile.append(f"{sum}")
    outfile.append(comment)
    vectors = np.array(vectors, dtype="float")*float(latticeparameter)
    if typevectors == 'Cartesian':
        for elem in getatoms:
            for position in atomsposition[elem[2]]:
                outfile.append(
                    f"{elem[2]}   {position[0]}   {position[1]}   {position[2]}")
    elif typevectors == 'Direct':
        for elem in getatoms:
            for position in atomsposition[elem[2]]:
                vcart = np.dot(np.array(vectors, dtype="float"),
                               np.array(position, dtype='float'))
                outfile.append(
                    f"{elem[2]}   {float(vcart[0]):.8f}   {float(vcart[1]):.8f}   {float(vcart[2]):.8f}")
    np.savetxt(outfilename, outfile, fmt='%s')
    return


def writefiledftb(typevectors, latticeparameter, vectors, getatoms, atomsposition, outfilename):
    outfile = []
    sum = 0
    comment = ""
    for el in getatoms:
        comment = comment+f"{el[2]}{el[3]} "
        sum = sum+int(el[3])
    if typevectors == 'Cartesian':
        outfile.append(f"{sum}   S")
    elif typevectors == 'Direct':
        outfile.append(f"{sum}   F")
    atoms = ""
    for i in range(len(getatoms)):
        atoms = atoms + f"{getatoms[i][2]}   "
    outfile.append(atoms)

    icont = 1
    for elem in getatoms:
        for position in atomsposition[elem[2]]:
            outfile.append(
                f"    {icont}  {elem[0]}  {position[0]}   {position[1]}   {position[2]}")
            icont = icont+1
    outfile.append(f"    0.00000000  0.00000000 0.00000000")
    vectors = np.array(vectors, dtype="float")*float(latticeparameter)
    for i in range(3):
        outfile.append(f"    {vectors[i][0]:.8f}   {vectors[i][1]:.8f}   {vectors[i][2]:.8f}")
    np.savetxt(outfilename, outfile, fmt='%s')
    return


def writefilexsf(typevectors, latticeparameter, vectors, getatoms, atomsposition, outfilename):
    outfile = []
    sum = 0
    comment = ""
    for el in getatoms:
        comment = comment+f"{el[2]}{el[3]} "
        sum = sum+int(el[3])
    vectors = np.array(vectors, dtype="float")*float(latticeparameter)
    outfile.append('# create by stb-translate ')
    outfile.append(f'# {comment}\n')
    outfile.append('CRYSTAL')
    outfile.append(f"PRIMVEC")
    for i in range(3):
        outfile.append(f"    {vectors[i][0]:.8f}   {vectors[i][1]:.8f}   {vectors[i][2]:.8f}")
    if typevectors == 'Cartesian':
        outfile.append(f"PRIMCOORD")
        outfile.append(f"{sum}  1")
        for elem in getatoms:
            for position in atomsposition[elem[2]]:
                outfile.append(
                    f"    {elem[1]}   {position[0]}   {position[1]}   {position[2]}")
    elif typevectors == 'Direct':
        outfile.append(f"PRIMCOORD")
        outfile.append(f"{sum}  1")
        for elem in getatoms:
            for position in atomsposition[elem[2]]:
                vcart = np.dot(np.array(vectors, dtype="float"),np.array(position, dtype='float'))
                outfile.append(f"{elem[1]}   {float(vcart[0]):.8f}   {float(vcart[1]):.8f}   {float(vcart[2]):.8f}")
    np.savetxt(outfilename, outfile, fmt='%s')
    return


def writefilefhi(typevectors, latticeparameter, vectors, getatoms, atomsposition, outfilename):
    outfile = []
    sum = 0
    comment = ""
    for el in getatoms:
        comment = comment+f"{el[2]}{el[3]} "
        sum = sum+int(el[3])
    vectors = np.array(vectors, dtype="float")*float(latticeparameter)
    for vector in vectors:
        outfile.append(f"lattice_vector   {vector[0]:.8f}   {vector[1]:.8f}   {vector[2]:.8f}")
    outfile.append(" ")
    if typevectors == 'Cartesian':
        for elem in getatoms:
            for position in atomsposition[elem[2]]:
                outfile.append(
                    f"atom   {position[0]}   {position[1]}   {position[2]}   {elem[2]}")
    elif typevectors == 'Direct':
        for elem in getatoms:
            for position in atomsposition[elem[2]]:
                outfile.append(
                    f"atom_frac   {position[0]}   {position[1]}   {position[2]}   {elem[2]}")
    np.savetxt(outfilename, outfile, fmt='%s')
    return


def main():

    parser = argparse.ArgumentParser(
        description="File format converter using stb-translate."
    )

    parser.add_argument("-if", "--in-format", required=True, choices=INPUT_FORMATS,
                        help="Input file format (options: poscar, cif, siesta, xyz, fhi, dftb, xsf)")
    parser.add_argument("-i", "--in-file", required=True,
                        help="Path to the input file")
    parser.add_argument("-of", "--out-format", required=True, choices=OUTPUT_FORMATS,
                        help="Output file format (options: cif , xyz, poscar, fdf, dftb, xsf, fhi)")
    parser.add_argument("-o", "--out-file", required=True,
                        help="Path to the output file")
    parser.add_argument(
        "--lattice", help="Lattice vectors file, required only for XYZ output")
    parser.add_argument("-v", "--version", action="version",
                        version=f"stb-translate {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()


    if args.intro == True:
        show_intro()

    print("\n" + color_text("TRANSLATE:", 'bold'))
    print("-"*60)




    # Validate lattice parameter requirement
    if args.in_format == "xyz" and not args.lattice:
        parser.error(
            "The --lattice argument is required when input format is XYZ.")

    print(f"\n[INFO] Converting {args.in_file} ({args.in_format}) to {args.out_file} ({args.out_format})...")

    if args.out_format == "xyz":
        print(f"[INFO] Lattice vector file: {args.lattice}")

    match (args.in_format):
        case "poscar":
            typevectors, latticeparameter, vectors, getatoms, atomsposition = getatomsandvectors_vasp(
                args.in_file)
        case "cif":
            typevectors, latticeparameter, vectors, getatoms, atomsposition = getatomsandvectors_cif(
                args.in_file)
        case "siesta":
            typevectors, latticeparameter, vectors, getatoms, atomsposition = getatomsandvectors_siesta(
                args.in_file)
        case "xyz":
            typevectors, latticeparameter, vectors, getatoms, atomsposition = getatomsandvectors_xyz(
                args.in_file, args.lattice)
        case "fhi":
            typevectors, latticeparameter, vectors, getatoms, atomsposition = getatomsandvectors_fhi(
                args.in_file)
        case "dftb":
            typevectors, latticeparameter, vectors, getatoms, atomsposition = getatomsandvectors_dftb(
                args.in_file)
        case "xsf":
            typevectors, latticeparameter, vectors, getatoms, atomsposition = getatomsandvectors_xsf(
                args.in_file)
        

    print(f"[OK] Read the file {args.in_file} ({args.in_format})")

    match (args.out_format):
        case "xyz":
            writefilexyz(typevectors, latticeparameter, vectors,
                         getatoms, atomsposition, args.out_file)
        case "poscar":
            writefileposcar(typevectors, latticeparameter, vectors,
                            getatoms, atomsposition, args.out_file)
        case "fdf":
            writefilefdf(typevectors, latticeparameter, vectors,
                         getatoms, atomsposition, args.out_file)
        case "dftb":
            writefiledftb(typevectors, latticeparameter, vectors,
                          getatoms, atomsposition, args.out_file)
        case "xsf":
            writefilexsf(typevectors, latticeparameter, vectors,
                         getatoms, atomsposition, args.out_file)
        case "fhi":
            writefilefhi(typevectors, latticeparameter, vectors,
                         getatoms, atomsposition, args.out_file)

        case "cif":
            writefilecif(typevectors, latticeparameter, vectors,
                         getatoms, atomsposition, args.out_file)

    print(f"[OK] Writing the file {args.out_file} ({args.out_format})")
    
    print("[INFO] Complete job!") 
    print("\n"+"-"*60)
    print(color_text("Converting input files is 10% coding, 90% crying.\n\n", 'bold'))

if __name__ == "__main__":
    main()
