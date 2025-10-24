# 📦 STB-SUITE – Siesta Toolbox Suite

**Version:** 1.8.0
**Author:** Dr. Carlos M. O. Bastos – University of Brasília (UnB) – 2025
**License:** MIT
**Compatible with:** Python ≥3.9, Linux/macOS/Windows

---

## 🔍 Description

**STB-SUITE** (Siesta Toolbox Suite) is a comprehensive collection of command-line tools developed to assist users of the **SIESTA** DFT code with common tasks in computational materials physics. The suite now includes robust tools for **calculation preparation** (input file generation, k-grids, and k-paths), **results analysis** (bands, DOS, structure, and symmetry), and **utilities** (format conversion and workspace cleanup).

It provides a unified and intuitive interface that streamlines the workflow of researchers in materials science.

---

## 🚀 Features

Version 1.8.0 significantly expands the suite's capabilities, organizing the functionalities into three main categories:

### 🔸 Calculation (Preparation)
- ✅ **Input File Generator (`stb-inputfile`)**
  Automatically generate FDF input files (`calc.fdf`) from structure files, including suggestions for calculation type (total energy, relaxation, AIMD, bands) and pseudopotential paths.
- ✅ **K-Grid Generator (`stb-kgrid`)**
  Suggest and automatically generate Monkhorst-Pack K-point grids based on the desired K-point density to optimize convergence.
- ✅ **K-Path Generator (`stb-kpath`)**
  Generate high-symmetry paths for band structure calculations (bands) from structure files.
- ✅ **Strain Generator (`stb-strain`)**
  Automatically generate supercells under uniaxial or biaxial strain in Cartesian coordinates.

### 🔸 Analysis (Post-processing)
- ✅ **Bands Analyzer (`stb-bands`)**
  Process and visualize SIESTA band structures, calculate band gaps, and customize the energy reference (VBM, CBM, Fermi, or manual value).
- ✅ **PDOS XML Parser (`stb-dos`)**
  Extract and analyze Projected Density of States (PDOS) data from the `PDOS.xml` file, with options for total DOS, by atom, and by species.
- ✅ **DOS Processor (Convolution) (`stb-convdos`)**
  Apply Gaussian smoothing (convolution) to DOS data to improve visualization.
- ✅ **Structural Analyzer (`stb-structural`)**
  Calculate lattice parameters, nearest neighbor analysis (ECN), and coordination numbers using multiple algorithms.
- ✅ **Symmetry Analyzer (`stb-symmetry`)**
  Extract space group, crystal system, point group, Wyckoff positions, and symmetry operations from structures.

### 🔸 Utilities & Interfaces
- ✅ **File Translator (`stb-translate`)**
  Seamlessly convert between structure file formats: CIF, POSCAR, XYZ, FDF (Siesta), DFTB, FHI-aims, and XSF.
- ✅ **Wantibexos Interface (`stb-siesta2wtb`)**
  Convert SIESTA Hamiltonians to the tight-binding format compatible with **Wantibexos**.
- ✅ **Clean File Tools (`stb-clean`)**
  Utility to automatically clean directories by removing temporary and unnecessary calculation files.
- ✅ **Unified GUI-like Terminal Interface**
  A user-friendly terminal interface for selecting tools and parameters interactively.

---

## 🧠 Requirements

- Python ≥3.9
- Conda (recommended)

### Python Dependencies (automatically installed via Conda or `pip`):
- `numpy`
- `matplotlib`
- `ase`
- `pymatgen`
- `spglib`
- `sisl`
- `argparse` (built-in)

---

## 📦 Installation

### 🔸 Using Conda (recommended)
```bash
conda install bastoscmo::stb_suite
```

### 🔸 Manual Installation (from GitHub)
```bash
git clone https://github.com/bastoscmo/stb-suite.git
cd stb-suite
pip install .
```

---

## 🔧 Usage

### ✨ Start the main suite:
```bash
stb-suite
```
An interactive terminal menu will guide you through the available tools.

---

## 🛠️ Individual Commands

| Command                  | Category                    | Description                                                              |
|--------------------------|-----------------------------|------------------------------------------------------------------------|
| `stb-inputfile`          | Preparation                 | Generates the `calc.fdf` input file from a structure.                  |
| `stb-kgrid`              | Preparation                 | Generates the K-point grid (Monkhorst-Pack).                           |
| `stb-kpath`              | Preparation                 | Generates the high-symmetry path for band structure calculation.         |
| `stb-strain`             | Preparation                 | Generates structures under strain.                                     |
| `stb-bands`              | Analysis                    | Analyzes `.bands` files and calculates band gaps.                      |
| `stb-dos`                | Analysis                    | Analyzes the `PDOS.xml` file (DOS and PDOS).                           |
| `stb-convdos`            | Analysis                    | Applies Gaussian convolution to DOS data.                              |
| `stb-structural`         | Analysis                    | Analyzes structural properties (ECN, coordination).                    |
| `stb-symmetry`           | Analysis                    | Analyzes the symmetry of crystal structures.                           |
| `stb-translate`          | Utility                     | Structure file format converter.                                       |
| `stb-siesta2wtb`         | Utility                     | Converts SIESTA Hamiltonian to Wantibexos.                             |
| `stb-clean`              | Utility                     | Cleans temporary and calculation files in a directory.                 |

---

## 📑 Examples

### ▶️ Input File Generation
```bash
stb-inputfile structure.fdf --type relax+d3 --pp-path /path/to/pseudopotentials
```

### ▶️ K-Grid Generation
```bash
stb-kgrid --file structure.cif --type cif --density 0.15
```

### ▶️ K-Path Generation
```bash
stb-kpath --file structure.poscar --type poscar --prec 0.001
```

### ▶️ Bands Analysis
```bash
stb-bands --file siesta.bands --shift vbm --plot
```

### ▶️ PDOS Analysis
```bash
stb-dos PDOS.xml --shift fermi --type total atom
```

### ▶️ File Conversion
```bash
stb-translate --in-format cif --in-file structure.cif --out-format poscar --out-file POSCAR
```

### ▶️ Directory Cleanup
```bash
stb-clean --path . --keep .fdf .psml .sh --dry-run
```

---

## 📚 Documentation

- The complete documentation for each tool is available through the `--help` argument:
```bash
stb-toolname --help
```

---

## ✍️ Citation

If you use **STB-SUITE** in your work, please cite:

> Comming son

Optionally, link to the GitHub repository.

---

## 🤝 Contributions

Contributions, issues, and feature requests are welcome. Feel free to check [issues](https://github.com/username/stb-suite/issues) or submit a pull request.

---

## ⚖️ License

This project is licensed under the MIT License – see the [LICENSE](LICENSE) file for details.

---

## ❤️ Acknowledgments

- Developed at the **LCCMat - Institute of Physics - University of Brasília (UnB)**.
- Thanks to the SIESTA development team and the open-source community.
