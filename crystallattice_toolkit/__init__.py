# crystallattice_toolkit/__init__.py

"""
CrystalLatticeToolkit
=====================

A small library for:
 - Real-space lattice construction (`Lattice`)
 - Space-group symmetry (`SpaceGroup`)
 - Reciprocal-space tools (`ReciprocalLattice`)
 - Fermi-surface reconstruction (`FermiSurface`)
"""

__version__ = "0.1.0"

# Core classes
from .lattice import Lattice
from .symmetry import SpaceGroup
from .reciprocal import ReciprocalLattice
from .fermi import FermiSurface

# Convenience imports (optional)
__all__ = [
    "Lattice",
    "SpaceGroup",
    "ReciprocalLattice",
    "FermiSurface",
    "__version__",
]
