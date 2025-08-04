

# 🧬 DockMate

**DockMate** is an all-in-one Python package for automated protein-ligand docking using AutoDock Vina. It wraps protein and ligand preparation, pocket detection, and docking in a simple, scriptable interface—perfect for batch screening and rapid prototyping.

## 🚀 Features

* Protein preparation using PDBFixer or MGLTools
* Ligand preparation and format conversion
* Pocket detection via P2Rank
* Fully automated docking with AutoDock Vina
* Visualization and fetching utilities
* Support for multi-ligand and multi-pocket docking

## 📦 Installation

```bash
git clone https://github.com/yourusername/dockmate.git
cd dockmate
pip install .
```

> ⚠️ Make sure you’re on **Windows** with dependencies (Open Babel, MGLTools, P2Rank, Vina) in the `bin/` folder as bundled.

## 📁 Directory Structure

```
dockmate/
├── bin/                # Includes binaries (vina, obabel, mgltools, p2rank)
├── utils/              # Helper utilities: fetch, clean, visualize
├── protein.py          # Protein preparation
├── ligand.py           # Ligand preparation
├── find_pocket.py      # Pocket detection
├── docking_engine.py   # Vina docking engine
├── __init__.py       
```

---

## 🧪 Usage

### ✅ Basic API

```python
from dockmate import Protein, Ligand, PocketFinder, VinaDocking
from dockmate.utils import fetch_pdb, fetch_sdf, clean_temp_folder, view_molecule
```

---

## 🧬 Example 1: Single Ligand Docking

```python
from dockmate import Protein, Ligand, PocketFinder, VinaDocking
from dockmate.utils import view_molecule, clean_temp_folder, fetch_pdb, fetch_sdf

clean_temp_folder()
fetch_pdb("1HVR")

protein = Protein("1HVR.pdb")
protein.prepare()
protein.save_pdbqt()

pf = PocketFinder(protein.pdb_path)
pockets = pf.run()
pf.save_report()

fetch_sdf(cid="338")
ligand = Ligand("338.sdf")
ligand.prepare(minimize=None)
ligand.save_pdbqt()

vina = VinaDocking(
    receptor=protein.pdbqt_path,
    ligand=ligand.pdbqt_path,
    center=pockets[0]['center'],
    size=(40, 40, 40),
    cpu=4
)
vina.run()
vina.save_results("vina_results/1HVR_338_pocket1_docking")

```

---

## 🔁 Example 2: Multi-Ligand Docking

```python
from pathlib import Path
from dockmate import Protein, Ligand, PocketFinder, VinaDocking
from dockmate.utils import clean_temp_folder, fetch_sdf

clean_temp_folder()

protein = Protein("1HVR.pdb")
protein.prepare()
protein.save_pdbqt()

pf = PocketFinder(protein.pdb_path)
pockets = pf.run()

# Download ligands
fetch_sdf("338", save_dir="ligands")
fetch_sdf("5288826", save_dir="ligands")

ligand_files = sorted(Path("ligands").glob("*.sdf"))

for ligand_file in ligand_files:
    ligand = Ligand(str(ligand_file))
    ligand.prepare(minimize=None)
    ligand.save_pdbqt()

    for pocket in pockets:
        vina = VinaDocking(
            receptor=protein.pdbqt_path,
            ligand=ligand.pdbqt_path,
            center=pocket["center"],
            size=(40, 40, 40),
            exhaustiveness=8,
            num_modes=9
        )
        vina.run()
        vina.save_results(f"vina_results/{protein.pdbqt_path.stem}_{ligand.pdbqt_path.stem}_pocket_{pocket['rank']}_docking")

```

---

## 🧰 Utility Functions

```python
fetch_pdb("1HVR")               # Downloads PDB file
fetch_sdf("338")                # Downloads SDF file
view_molecule("path/to/file")   # Visualize molecule in browser
clean_temp_folder()             # Clear temporary workspace
```

---

## 📄 License

MIT License. See `LICENSE` for details.
