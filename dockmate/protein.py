import os
import subprocess
from pathlib import Path
from typing import Union, Optional
import shutil

from pdbfixer import PDBFixer
from openmm.app import PDBFile

# Constants
MGLTOOLS_PATH = (Path(__file__).parent / "bin" / "mgltools").resolve()
MGL_PYTHON_EXE = (MGLTOOLS_PATH / "python.exe").resolve()
PREPARE_RECEPTOR_SCRIPT = (
    MGLTOOLS_PATH / "Lib" / "site-packages" /
    "AutoDockTools" / "Utilities24" / "prepare_receptor4.py"
).resolve()

OBABEL_EXE = (Path(__file__).parent / "bin" /
              "OpenBabel-3.1.1" / "obabel.exe").resolve()
TEMP_DIR = (Path(__file__).parent / "temp").resolve()
TEMP_DIR.mkdir(parents=True, exist_ok=True)


class Protein:
    """Handles protein preparation for docking using PDBFixer, Open Babel, and MGLTools.

    This class automates the preprocessing of protein structures by:
    - Converting various input formats to PDB using Open Babel
    - Fixing missing residues and atoms using PDBFixer
    - Optionally removing water molecules and adding hydrogens
    - Converting the fixed PDB file to PDBQT format using MGLTools
    """

    SUPPORTED_INPUTS = {".pdb", ".mol2", ".sdf", ".pdbqt", ".ent", ".xyz"}

    def __init__(self, file_path: Union[str, Path]):
        """
        Initialize the Protein object.

        Args:
            file_path (str | Path): Path to the input protein file.

        Raises:
            FileNotFoundError: If the file does not exist.
            ValueError: If the file format is unsupported.
        """
        self.file_path = Path(file_path).resolve()
        self.pdb_path: Optional[Path] = None
        self.pdbqt_path: Optional[Path] = None

        if not self.file_path.is_file():
            raise FileNotFoundError(
                f"❌ Protein file not found: {self.file_path}")

        self.ext = self.file_path.suffix.lower()
        if self.ext not in self.SUPPORTED_INPUTS:
            raise ValueError(
                f"❌ Unsupported file format '{self.ext}'. Supported formats: {self.SUPPORTED_INPUTS}")

    def prepare(
        self,
        pH: float = 7.4,
        add_hydrogens: bool = True,
        remove_water: bool = True,
        add_charges: bool = True  # currently unused, can be implemented if needed
    ) -> None:
        """
        Prepares the protein for docking by performing the following steps:

        1. Converts the input file to PDB format using Open Babel (if needed).
        2. Fixes the structure using PDBFixer (adds missing atoms/residues, removes heterogens).
        3. Optionally adds hydrogens at the given pH.
        4. Saves the fixed structure to a temporary .pdb file.
        5. Converts the fixed PDB to PDBQT format using MGLTools.

        Args:
            pH (float): The pH at which to add hydrogens. Defaults to 7.4.
            add_hydrogens (bool): Whether to add hydrogens. Defaults to True.
            remove_water (bool): Whether to remove water molecules. Defaults to True.
            add_charges (bool): Placeholder for adding charges (not implemented).
        """
        # Convert to .pdb if needed
        if self.ext != ".pdb":
            self.pdb_path = TEMP_DIR / f"{self.file_path.stem}.pdb"
            cmd = [str(OBABEL_EXE), str(self.file_path),
                   "-O", str(self.pdb_path), "--gen3d"]
            subprocess.run(cmd, check=True)
        else:
            self.pdb_path = self.file_path

        # Fix structure using PDBFixer
        fixer = PDBFixer(filename=str(self.pdb_path))
        fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.removeHeterogens(keepWater=not remove_water)

        if add_hydrogens:
            fixer.addMissingHydrogens(pH)

        # Save fixed PDB
        fixed_pdb_path = TEMP_DIR / f"{self.file_path.stem}_fixed.pdb"
        with open(fixed_pdb_path, "w") as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f)

        self.pdb_path = fixed_pdb_path

        # Convert to PDBQT using MGLTools
        output_pdbqt = TEMP_DIR / f"{self.file_path.stem}.pdbqt"
        cmd = [
            str(MGL_PYTHON_EXE),
            str(PREPARE_RECEPTOR_SCRIPT),
            "-r", str(self.pdb_path),
            "-o", str(output_pdbqt)
        ]
        if add_hydrogens:
            cmd += ["-A", "hydrogens"]
        if remove_water:
            cmd += ["-U", "waters"]

        subprocess.run(cmd, check=True)
        self.pdbqt_path = output_pdbqt

    def save_pdbqt(self, save_path: Union[str, Path] = ".") -> None:
        """
        Saves the prepared PDBQT file to the specified destination.

        Args:
            save_path (str | Path): File path or directory where the PDBQT file will be saved.
                                    If a directory is provided, the file is saved with its original name.

        Raises:
            RuntimeError: If the PDBQT file has not been generated via prepare().
        """
        if self.pdbqt_path is None or not self.pdbqt_path.exists():
            raise RuntimeError("❌ Protein not prepared. Run prepare() first.")

        save_path = Path(save_path).resolve()
        save_path.parent.mkdir(parents=True, exist_ok=True)

        shutil.copy2(self.pdbqt_path, save_path)
