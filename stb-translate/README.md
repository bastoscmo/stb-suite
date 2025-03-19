# stb-translate

stb-translate is a command-line tool for converting atomic structure files between different formats, including POSCAR (VASP), CIF, SIESTA (OUT_STRUCTURE), XYZ, FHI, DFTB (GEN), and XSF.

## Features
- Convert atomic structure files between multiple formats
- Supports both input and output format validation
- Requires lattice parameters when converting to XYZ format
- Easy-to-use command-line interface

## Requirements

This project depends on the following Python libraries:

- `numpy`: For numerical array manipulation and scientific calculations.
- `pymatgen`: For analysis and manipulation of crystalline structures.

## Installing Dependencies

To install the dependencies, use the following command:

```bash
pip install -r requirements.txt
```

## Installation
Ensure you have Python 3 installed. Clone this repository and navigate to the script's directory:

```sh
 git clone https://github.com/bastoscmo/stb-translate.git          
 cd stb-translate
```

## Usage
### Basic Command Format
```sh
python stb-translate.x --in-format FORMAT --in-file INPUT_FILE --out-format FORMAT --out-file OUTPUT_FILE [--lattice LATTICE-VECTORS_FILE]
```

### Arguments
| Argument | Description |
|----------|-------------|
| `-if`, `--in-format` | Input file format (poscar, cif, siesta, xyz, fhi, dftb, xsf) |
| `-i`, `--in-file` | Path to the input file |
| `-of`, `--out-format` | Output file format (xyz, poscar, fdf, dftb, xsf, fhi) |
| `-o`, `--out-file` | Path to the output file |
| `--lattice` | file with lattice vectors **Required only when output format is XYZ** |
| `-v`, `--version` | Display the current version |
| `-h`, `--help` | Show help message |

### Examples
#### Convert CIF to POSCAR
```sh
python stb-translate.x -if cif -i structure.cif -of poscar -o structure.POSCAR
```
#### Convert POSCAR to XYZ (with lattice parameters)
```sh
python stb-translate.x -if poscar -i structure.POSCAR -of xyz -o structure.xyz --lattice lattice.dat
```

#### Convert SIESTA to DFTB (GEN)
```sh
python stb-translate.x -if siesta -i structure.OUT_STRUCTURE -of dftb -o structure.gen
```

## Error Handling
- If an unsupported format is provided, the script will return an error.
- If `--lattice` is missing when outputting XYZ format, the script will prompt for it.

## License
This project is licensed under the MIT License.

## Contributing
Feel free to submit issues or pull requests for improvements!

## Contact
For any questions or suggestions, open an issue or contact the repository owner.

