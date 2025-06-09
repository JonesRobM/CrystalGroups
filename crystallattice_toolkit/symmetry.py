"""
Module for handling space-group symmetry via spglib.
"""
import numpy as np
import spglib
from .lattice import Lattice

class SpaceGroup:
    """
    Wrap spglib to extract symmetry operations for a given lattice.

    Attributes:
        spg_number (int): International space-group number (1–230).
        lattice (Lattice): Real-space lattice object.
    """
    def __init__(
        self,
        spg_number: int,
        lattice: Lattice,
        symprec: float = 1e-5
    ):
        """
        Initialise SpaceGroup.

        Args:
            spg_number: integer between 1 and 230.
            lattice: Lattice instance with unit cell definition.
            symprec: float symmetry tolerance.
        """
        self.spg_number = spg_number
        self.lattice = lattice
        self.symprec = symprec

    def get_symmetry_operations(self) -> tuple[np.ndarray, np.ndarray]:
        """
        Compute rotation matrices and translation vectors.

        Returns:
            rotations: M×3×3 array of integer rotation matrices.
            translations: M×3 array of fractional translations.
        """
        cell = (
            self.lattice.lattice_vectors,
            self.lattice.positions_frac,
            self.lattice.species
        )
        sym_data = spglib.get_symmetry(cell, symprec=self.symprec)
        rotations = np.array(sym_data['rotations'], dtype=int)
        translations = np.array(sym_data['translations'], dtype=float)
        return rotations, translations

    def get_spacegroup_symbol(self) -> str:
        """
        Return Hermann–Mauguin symbol for the space group.
        """
        cell = (
            self.lattice.lattice_vectors,
            self.lattice.positions_frac,
            self.lattice.species
        )
        dataset = spglib.get_spacegroup(cell, symprec=self.symprec)
        return dataset