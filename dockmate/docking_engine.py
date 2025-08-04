import subprocess
from pathlib import Path
from typing import Optional, Union
import pandas as pd
import shutil
import os

# Path to AutoDock Vina executable
VINA_PATH = (Path(__file__).parent / "bin" / "Vina" / "vina.exe").resolve()


class VinaDocking:
    """
    A wrapper class for running AutoDock Vina to perform molecular docking and extract results.

    Attributes:
        receptor (Path): Path to the receptor file in PDBQT format.
        ligand (Path): Path to the ligand file in PDBQT format.
        center (tuple): Grid center coordinates as (x, y, z).
        size (tuple): Grid box dimensions as (x, y, z).
        exhaustiveness (int): Vina exhaustiveness setting controlling search thoroughness.
        num_modes (int): Number of output binding modes to generate.
        cpu (int): Number of CPU cores to use for docking.
        seed (Optional[int]): Optional random seed for reproducibility.
    """

    def __init__(
        self,
        receptor: Union[str, Path],
        ligand: Union[str, Path],
        center: tuple[float, float, float],
        size: tuple[float, float, float],
        exhaustiveness: int = 8,
        num_modes: int = 9,
        cpu: int = os.cpu_count(),
        seed: Optional[int] = None,
    ):
        self.receptor = Path(receptor).resolve()
        self.ligand = Path(ligand).resolve()
        self._validate_inputs(center, size)

        if not self.receptor.is_file():
            raise FileNotFoundError(
                f"❌ Receptor file not found: {self.receptor}")
        if not self.ligand.is_file():
            raise FileNotFoundError(f"❌ Ligand file not found: {self.ligand}")

        self.center = center
        self.size = size
        self.exhaustiveness = exhaustiveness
        self.num_modes = num_modes
        self.cpu = cpu
        self.seed = seed

        self.base_name = f"{self.receptor.stem}_{self.ligand.stem}_docked"

        # Temporary results directory
        TEMP_RESULTS_DIR = (Path(__file__).parent /
                            "temp" / "vina_results").resolve()
        TEMP_RESULTS_DIR.mkdir(parents=True, exist_ok=True)

        # Output files
        self.output_pdbqt = TEMP_RESULTS_DIR / f"{self.base_name}.pdbqt"
        self.output_log = TEMP_RESULTS_DIR / f"{self.base_name}.log"
        self.output_csv = TEMP_RESULTS_DIR / f"{self.base_name}.csv"

        self._vina_output: Optional[str] = None
        self._docking_df: Optional[pd.DataFrame] = None

    def _validate_inputs(self, center, size):
        if not (isinstance(center, tuple) and len(center) == 3):
            raise ValueError("⚠️ 'center' must be a 3-tuple of floats.")
        if not (isinstance(size, tuple) and len(size) == 3):
            raise ValueError("⚠️ 'size' must be a 3-tuple of floats.")
        if any(not isinstance(v, (float, int)) for v in center + size):
            raise TypeError(
                "⚠️ Grid center and size values must be float or int.")

    def run(self) -> pd.DataFrame:
        """
        Executes AutoDock Vina with the specified parameters.

        Returns:
            pd.DataFrame: A DataFrame containing the docking results with columns:
                - mode: Binding mode index.
                - affinity_kcal_mol: Binding affinity (kcal/mol).
                - rmsd_lb: Lower bound of RMSD from best mode.
                - rmsd_ub: Upper bound of RMSD from best mode.

        Raises:
            RuntimeError: If Vina execution fails or produces no output.
        """
        cmd = [
            str(VINA_PATH),
            "--receptor", str(self.receptor),
            "--ligand", str(self.ligand),
            "--center_x", str(self.center[0]),
            "--center_y", str(self.center[1]),
            "--center_z", str(self.center[2]),
            "--size_x", str(self.size[0]),
            "--size_y", str(self.size[1]),
            "--size_z", str(self.size[2]),
            "--out", str(self.output_pdbqt),
            "--exhaustiveness", str(self.exhaustiveness),
            "--num_modes", str(self.num_modes),
            "--cpu", str(self.cpu),
        ]

        if self.seed is not None:
            cmd += ["--seed", str(self.seed)]

        try:
            result = subprocess.run(
                cmd, capture_output=True, text=True, check=True)
            self._vina_output = result.stdout

            if self._vina_output:
                with open(self.output_log, "w") as log_file:
                    log_file.write(self._vina_output)

        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"❌ AutoDock Vina failed:\n{e.stderr or e.stdout}")

        self._docking_df = self._parse_output()
        self._docking_df.to_csv(self.output_csv, index=False)
        return self._docking_df

    def _parse_output(self) -> pd.DataFrame:
        """
        Parses the textual output from AutoDock Vina and extracts docking mode data.

        Returns:
            pd.DataFrame: A DataFrame containing parsed results with the columns:
                - mode
                - affinity_kcal_mol
                - rmsd_lb
                - rmsd_ub

        Raises:
            RuntimeError: If Vina output has not been generated.
            ValueError: If no valid docking results are found.
        """
        if self._vina_output is None:
            raise RuntimeError(
                "❌ No Vina output to parse. Please run docking first.")

        lines = self._vina_output.splitlines()
        data = []
        header_found = False

        for line in lines:
            if line.strip().startswith("mode"):
                header_found = True
                continue
            if header_found and line.strip() and not line.startswith("-"):
                tokens = line.split()
                if len(tokens) == 4:
                    try:
                        mode, affinity, rmsd_lb, rmsd_ub = tokens
                        data.append((int(mode), float(affinity),
                                    float(rmsd_lb), float(rmsd_ub)))
                    except ValueError:
                        continue  # Skip malformed lines

        if not data:
            raise ValueError("❌ No docking results found in Vina output.")

        return pd.DataFrame(
            data, columns=["mode", "affinity_kcal_mol", "rmsd_lb", "rmsd_ub"]
        )

    def save_results(self, save_dir: Union[str, Path] = Path("./vina_results")) -> Path:
        """
        Copies the docking result files (.pdbqt, .log, .csv) to the specified directory.

        Args:
            save_dir (Union[str, Path], optional): Destination directory to save the output files.
                Defaults to './vina_results'.

        Returns:
            Path: The full resolved path of the target output folder.

        Raises:
            RuntimeError: If any result file is missing (docking not run or failed).
        """
        # Check if results exist before proceeding
        if not all(file.exists() for file in [self.output_pdbqt, self.output_log, self.output_csv]):
            raise RuntimeError(
                "❌ Docking results are missing. Please run docking before saving results.")

        save_dir = Path(save_dir).resolve()
        save_dir.mkdir(parents=True, exist_ok=True)

        for file in [self.output_pdbqt, self.output_log, self.output_csv]:
            shutil.copy(file, save_dir / file.name)
