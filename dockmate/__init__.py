"""
DockMate: All-in-one Protein-Ligand Docking Package
Integrates MGLTools, P2Rank, and AutoDock Vina
"""

__version__ = "1.0.0"
__author__ = "DockMate Team"

from .protein import Protein
from .ligand import Ligand
from .docking_engine import VinaDocking
from .pocket_finder import PocketFinder


__all__ = [
    "Protein",
    "Ligand",
    "VinaDocking",
    "PocketFinder",

]
