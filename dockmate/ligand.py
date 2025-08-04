import shutil
import subprocess
from pathlib import Path
from typing import Optional, Union

# === CONFIGURATION ===
MGLTOOLS_PATH = (Path(__file__).parent / "bin" / "mgltools").resolve()
MGL_PYTHON_EXE = (MGLTOOLS_PATH / "python.exe").resolve()
PREPARE_LIGAND_SCRIPT = (MGLTOOLS_PATH / "Lib" / "site-packages" /
                         "AutoDockTools" / "Utilities24" / "prepare_ligand4.py").resolve()
OBABEL_EXE = (Path(__file__).parent / "bin" /
              "OpenBabel-3.1.1" / "obabel.exe").resolve()

TEMP_DIR = (Path(__file__).parent / "temp").resolve()
TEMP_DIR.mkdir(exist_ok=True)


class Ligand:
    """
    Ligand preparation pipeline using Open Babel and MGLTools.

    This class supports automatic conversion, 3D generation, optional energy minimization, 
    and conversion to the PDBQT format required by AutoDock Vina.

    Supported input formats: mol2, sdf, pdb, mol, smi  
    Supported forcefields: mmff94, mmff94s, uff, gaff
    """

    SUPPORTED_INPUTS = {"mol2", "sdf", "pdb", "mol", "smi"}
    SUPPORTED_FORCEFIELDS = {"mmff94", "mmff94s", "uff", "gaff"}

    def __init__(self, file_path: Union[str, Path]):
        self.file_path = Path(file_path).resolve()
        self.mol2_path: Optional[Path] = None
        self.pdbqt_path: Optional[Path] = None

        if not self.file_path.is_file():
            raise FileNotFoundError(
                f"❌ Ligand file not found: {self.file_path}")

        ext = self.file_path.suffix.lower().lstrip(".")
        if ext not in self.SUPPORTED_INPUTS:
            raise ValueError(
                f"❌ Unsupported file format '.{ext}'. Supported formats: {self.SUPPORTED_INPUTS}")
        self.input_format = ext

    def prepare(
        self,
        minimize: Optional[str] = None,
        add_hydrogens: bool = True,
        remove_nphs: bool = True,
        remove_lps: bool = True,
        remove_waters: bool = True,
    ) -> Path:
        """
        Prepare the ligand by converting to MOL2, optionally minimizing energy, 
        and generating a final PDBQT file using MGLTools.

        Args:
            minimize (str, optional): Forcefield to use for energy minimization 
                ("mmff94", "mmff94s", "uff", or "gaff"). If None, no minimization is performed.
            add_hydrogens (bool): Whether to add hydrogens during MGLTools processing. Default is True.
            remove_nphs (bool): Remove non-polar hydrogens using MGLTools. Default is True.
            remove_lps (bool): Remove lone pairs using MGLTools. Default is True.
            remove_waters (bool): Remove water molecules using MGLTools. Default is True.

        Returns:
            Path: Path to the generated PDBQT file.

        Raises:
            ValueError: If an unsupported forcefield or input format is provided.
            RuntimeError: If MGLTools fails to generate the PDBQT file.
        """

        # === Step 1: Convert + Gen3D + Minimize to MOL2 ===
        self.mol2_path = TEMP_DIR / f"{self.file_path.stem}.mol2"
        cmd = [
            str(OBABEL_EXE), "-i", self.input_format, str(self.file_path),
            "-o", "mol2", "-O", str(self.mol2_path),
            "--gen3d", "-h"
        ]

        if minimize:
            forcefield = minimize.lower()
            if forcefield not in self.SUPPORTED_FORCEFIELDS:
                raise ValueError(
                    f"❌ Unsupported forcefield '{forcefield}'. Supported: {self.SUPPORTED_FORCEFIELDS}")
            cmd += ["--minimize", "--ff", forcefield]

        subprocess.run(cmd, check=True)

        # === Step 2: MGLTools to PDBQT ===
        pdbqt_filename = f"{self.mol2_path.stem}.pdbqt"
        mgl_cmd = [
            str(MGL_PYTHON_EXE), str(PREPARE_LIGAND_SCRIPT),
            "-l", self.mol2_path.name, "-o", pdbqt_filename
        ]

        if add_hydrogens:
            mgl_cmd += ["-A", "hydrogens"]

        remove_flags = []
        if remove_nphs:
            remove_flags.append("nphs")
        if remove_lps:
            remove_flags.append("lps")
        if remove_waters:
            remove_flags.append("waters")
        if remove_flags:
            mgl_cmd += ["-U", "_".join(remove_flags)]

        result = subprocess.run(
            mgl_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            cwd=TEMP_DIR
        )

        if result.returncode != 0:
            print("❌ MGLTools failed.")
            print("STDOUT:\n", result.stdout)
            print("STDERR:\n", result.stderr)
            raise RuntimeError("❌ MGLTools ligand preparation failed.")

        self.pdbqt_path = TEMP_DIR / pdbqt_filename

    def save_pdbqt(self, save_path: Union[str, Path] = ".") -> None:
        """
        Save the prepared PDBQT file to the specified location.

        Args:
            save_path (str | Path): Destination file or directory where the PDBQT file should be saved. 
                If a directory is given, the original filename is preserved.

        Raises:
            RuntimeError: If prepare() has not been called or the PDBQT file does not exist.
        """
        if self.pdbqt_path is None or not self.pdbqt_path.exists():
            raise RuntimeError("❌ Ligand not prepared. Run prepare() first.")

        save_path = Path(save_path).resolve()
        save_path.parent.mkdir(parents=True, exist_ok=True)

        shutil.copy2(self.pdbqt_path, save_path)
