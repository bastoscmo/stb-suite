# stb-structural.x

**stb-structural.x** is a command-line tool designed to compute the Effective Coordination Number (ECN) of crystalline structures using `pymatgen`. It supports two operation modes: **list** (for specific atoms) and **mean** (for the entire structure).

## 📌 Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/your-username/stb-structural.git
   cd stb-structural
   ```
2. Install dependencies (a `conda` or `venv` environment is recommended):
   ```bash
   pip install pymatgen numpy
   ```
3. Make the script executable:
   ```bash
   chmod +x stb-structural.x
   ```

## 🚀 Usage

Run the program with the following syntax:
```bash
stb-structural.x --file FILE --mode (list or mean) --list [1,2,3,4,5]
```

### **Options**:
| Option          | Description |
|---------------|------------------------------------------------|
| `--file FILE` | Input structure file (POSCAR, CIF, etc.). |
| `--mode list` | Computes ECN for specific atoms. Requires `--list`. |
| `--mode mean` | Computes the average ECN for the entire structure. |
| `--list [1,2,3]` | List of atom indices for ECN calculation (used with `list` mode). |
| `--version`   | Displays the program version. |

## 📊 Examples

### **1️⃣ Compute the average ECN of the entire structure**
```bash
stb-structural.x --file structure.cif --mode mean
```

### **2️⃣ Compute ECN for specific atoms (list mode)**
```bash
stb-structural.x --file structure.cif --mode list --list [1,5,10]
```

### **3️⃣ Display the program version**
```bash
stb-structural.x --version
```

## 📜 Program Output
Results are printed to the terminal and saved in the file `structural_information.dat`. The output includes:
- **Lattice parameters** (lattice vectors, angles, etc.).
- **Atomic positions** in Cartesian coordinates.
- **ECN values** for different coordination number methods.

## 🔗 Dependencies
- Python 3.x
- `pymatgen`
- `numpy`

## 🛠️ Author
Developed by **Dr. Carlos M. O. Bastos** – University of Brasília (UnB).
