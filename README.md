
# 📦 STB-SUITE – Siesta Toolbox Suite

**Version:** 1.5.10  
**Author:** Dr. Carlos M. O. Bastos – University of Brasília (UnB) – 2025  
**License:** MIT   
**Compatible with:** Python ≥3.9, Linux/macOS/Windows

---

## 🔍 Description

**STB-SUITE** (Siesta Toolbox Suite) is a comprehensive collection of command-line tools developed to assist users of the **SIESTA** DFT code with common tasks such as structural manipulation, symmetry analysis, electronic structure post-processing (bands and DOS), strain generation, and file format conversion.

It provides an intuitive and unified interface that streamlines the workflow of computational materials science researchers.

---

## 🚀 Features

- ✅ **Band Structure Analysis**  
  Process and visualize SIESTA band structures with customizable energy references.

- ✅ **Density of States Convolution**  
  Apply Gaussian smoothing to DOS data to improve visualization.

- ✅ **Strain Generator**  
  Automatically generate supercells under uniaxial or biaxial strain in Cartesian coordinates.

- ✅ **Structural Analysis**  
  Compute lattice parameters, nearest neighbor analysis (ECN), and coordination numbers using multiple algorithms.

- ✅ **Symmetry Analysis**  
  Extract space group, crystal system, point group, Wyckoff positions, and symmetry operations.

- ✅ **File Format Conversion**  
  Seamlessly convert between CIF, POSCAR, XYZ, FDF (Siesta), DFTB, FHI-aims, and XSF formats.

- ✅ **Interface to Wantibexos**  
  Convert SIESTA Hamiltonians to tight-binding format compatible with **Wantibexos**.

- ✅ **Clean Workspace Utility**  
  Automatically clean directories by removing temporary and unnecessary files.

- ✅ **Unified GUI-like Terminal Interface**  
  A user-friendly terminal interface for selecting tools and parameters interactively.

---

## 🧠 Requirements

- Python ≥3.9  
- Conda (recommended)  

### Python Dependencies (installed automatically via Conda or `pip`):
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
git clone https://github.com/username/stb-suite.git
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

| Command                  | Description                          |
|--------------------------|--------------------------------------|
| `stb-bands`              | Band structure analysis             |
| `stb-convdos`            | DOS convolution (Gaussian)          |
| `stb-strain`             | Strain generator                    |
| `stb-structural`         | Structural analysis                 |
| `stb-symmetry`           | Symmetry analysis                   |
| `stb-translate`          | File format converter               |
| `stb-wantibexos`         | SIESTA → Wantibexos converter       |
| `stb-clean`              | Clean files in directories          |

---

## 📑 Examples

### ▶️ Band structure processing
```bash
stb-bands --file siesta.bands --shift fermi
```

### ▶️ Apply Gaussian convolution to DOS
```bash
stb-convdos --file siesta.DOS --size 11 --sigma 0.2 --out dos_filtered.dat
```

### ▶️ Generate strained structures
```bash
stb-strain --file structure.fdf --stdir xy --stmin -5 --stmax 5 --step 1
```

### ▶️ Structural analysis
```bash
stb-structural --file structure.cif --format cif --mode mean
```

### ▶️ Symmetry analysis
```bash
stb-symmetry --file structure.cif --format cif
```

### ▶️ File conversion
```bash
stb-translate --file structure.cif --in cif --out poscar
```

### ▶️ Convert to Wantibexos
```bash
stb-wantibexos --input siesta.fdf --output siesta.out
```

### ▶️ Clean directory
```bash
stb-clean --keep .fdf .psml .sh
```

---

## 📚 Documentation

- Full documentation: **(Coming Soon)**  
- Each tool provides help via:
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

