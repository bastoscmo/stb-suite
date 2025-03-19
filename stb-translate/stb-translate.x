#!/usr/bin/env python3

#################################################
# Siesta Tool Box - Translate                   #
# version 1.0.1                                 #
# UnB - 2025/02/05                              #
# Dr. Carlos M. O. Bastos                       #
#################################################

VERSION = "1.0.1"

# Import libraries
from lib import sstranslate
import sys
import argparse

# Supported formats
INPUT_FORMATS = {"poscar", "cif", "siesta", "xyz", "fhi", "dftb", "xsf"}
OUTPUT_FORMATS = {"xyz", "poscar", "fdf", "dftb", "xsf", "fhi"}

def main():
    parser = argparse.ArgumentParser(
        description="File format converter using stb-translate."
    )

    parser.add_argument("-if", "--in-format", required=True, choices=INPUT_FORMATS,
                        help="Input file format (options: poscar, cif, siesta, xyz, fhi, dftb, xsf)")
    parser.add_argument("-i", "--in-file", required=True, help="Path to the input file")
    parser.add_argument("-of", "--out-format", required=True, choices=OUTPUT_FORMATS,
                        help="Output file format (options: xyz, poscar, fdf, dftb, xsf, fhi)")
    parser.add_argument("-o", "--out-file", required=True, help="Path to the output file")
    parser.add_argument("--lattice", help="Lattice vectors file, required only for XYZ output")
    parser.add_argument("-v", "--version", action="version", version=f"stb-translate {VERSION}")

    args = parser.parse_args()

    #Validate lattice parameter requirement
    if args.in_format == "xyz" and not args.lattice:
        parser.error("The --lattice argument is required when input format is XYZ.")

    print(f"Converting {args.in_file} ({args.in_format}) to {args.out_file} ({args.out_format})...")

    if args.out_format == "xyz":
        print(f"Lattice vector file: {args.lattice}")


    match (args.in_format):
        case "poscar":
            typevectors, latticeparameter, vectors, getatoms, atomsposition =sstranslate.getatomsandvectors_vasp(args.in_file)
        case "cif":
            typevectors, latticeparameter, vectors, getatoms, atomsposition =sstranslate.getatomsandvectors_cif(args.in_file)
        case "siesta":
            typevectors, latticeparameter, vectors, getatoms, atomsposition =sstranslate.getatomsandvectors_siesta(args.in_file)
        case "xyz":
            typevectors, latticeparameter, vectors, getatoms, atomsposition =sstranslate.getatomsandvectors_xyz(args.in_file,args.lattice)
        case "fhi":
            typevectors, latticeparameter, vectors, getatoms, atomsposition =sstranslate.getatomsandvectors_fhi(args.in_file)
        case "dftb":
            typevectors, latticeparameter, vectors, getatoms, atomsposition =sstranslate.getatomsandvectors_dftb(args.in_file)
        case "xsf":
            typevectors, latticeparameter, vectors, getatoms, atomsposition =sstranslate.getatomsandvectors_xsf(args.in_file)

    match (args.out_format):
        case "xyz":
            sstranslate.writefilexyz(typevectors, latticeparameter, vectors, getatoms, atomsposition,args.out_file)
        case "poscar":
            sstranslate.writefileposcar(typevectors, latticeparameter, vectors, getatoms, atomsposition,args.out_file)
        case "fdf":
            sstranslate.writefilefdf(typevectors, latticeparameter, vectors, getatoms, atomsposition,args.out_file)
        case "dftb":
            sstranslate.writefiledftb(typevectors, latticeparameter, vectors, getatoms, atomsposition,args.out_file)
        case "xsf":
            sstranslate.writefilexsf(typevectors, latticeparameter, vectors, getatoms, atomsposition,args.out_file)
        case "fhi":
            sstranslate.writefilefhi(typevectors, latticeparameter, vectors, getatoms, atomsposition,args.out_file)

if __name__ == "__main__":
    main()
