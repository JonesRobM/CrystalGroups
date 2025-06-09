"""
Module for building real‐space crystal lattices.
"""
import numpy as np
from scipy.spatial import cKDTree

class Lattice:
    """
    Represents a crystal lattice in Cartesian coordinates.

    Attributes:
        lattice_vectors (np.ndarray): 3×3 array of basis vectors (rows are vectors).
        species (list[str]): List of atomic species labels.
        positions_frac (np.ndarray): N×3 array of fractional coordinates (0≤u<1).
    """
    def __init__(
        self,
        lattice_vectors: np.ndarray,
        species: list,
        positions_frac: np.ndarray
    ):
        """
        Initialise Lattice.

        Args:
            lattice_vectors: array_like of shape (3,3), real-space basis vectors as rows.
            species: list of length N of atomic symbols.
            positions_frac: array_like of shape (N,3), fractional coordinates in [0,1).
        """
        self.lattice_vectors = np.asarray(lattice_vectors, dtype=float)
        self.species = list(species)
        self.positions_frac = np.asarray(positions_frac, dtype=float)
        if self.positions_frac.ndim != 2 or self.positions_frac.shape[1] != 3:
            raise ValueError("positions_frac must be N×3 array")
        if len(self.species) != self.positions_frac.shape[0]:
            raise ValueError("Length of species list must match number of positions")

    def get_cartesian_positions(self) -> np.ndarray:
        """
        Convert fractional coordinates to Cartesian positions.

        Returns:
            cartesian (np.ndarray): N×3 array of positions in Å.
        """
        return self.positions_frac.dot(self.lattice_vectors)

    def get_neighbors(
        self,
        cutoff: float
    ) -> list[tuple[int,int,float]]:
        """
        Find atomic pairs within a given distance cutoff (periodic images included).

        Args:
            cutoff: float, distance threshold in same units as lattice_vectors.

        Returns:
            neighbors: list of tuples (i, j, distance) with i<j.
        """
        # Build supercell points within one extra shell
        # Generate 3^3 neighbour images
        shifts = np.array([[i, j, k]
                           for i in (-1, 0, 1)
                           for j in (-1, 0, 1)
                           for k in (-1, 0, 1)])
        # Expand positions
        frac_images = (self.positions_frac[None, :, :] + shifts[:, None, :])  # (27, N, 3)
        frac_images = frac_images.reshape(-1, 3)
        cart_images = frac_images.dot(self.lattice_vectors)

        tree = cKDTree(cart_images)
        pairs = tree.query_pairs(r=cutoff, output_type='set')

        N = len(self.positions_frac)
        neighbors = []
        for idx1, idx2 in pairs:
            # Map flat indices back to original atoms
            i = idx1 % N
            j = idx2 % N
            if i < j:
                d = np.linalg.norm(cart_images[idx1] - cart_images[idx2])
                neighbors.append((i, j, d))
        return neighbors

    def make_supercell(
        self,
        multipliers: tuple[int, int, int]
    ) -> "Lattice":
        """
        Generate a supercell by replicating the unit cell.

        Args:
            multipliers: (n_a, n_b, n_c) replication counts along each lattice vector.

        Returns:
            supercell (Lattice): new Lattice object.
        """
        na, nb, nc = multipliers
        shifts = np.array([[i, j, k]
                           for i in range(na)
                           for j in range(nb)
                           for k in range(nc)])
        new_frac = (self.positions_frac[None, :, :] + shifts[:, None, :])
        new_frac = new_frac.reshape(-1, 3) / np.array([na, nb, nc])
        new_species = self.species * len(shifts)
        new_vectors = np.vstack([
            self.lattice_vectors[0] * na,
            self.lattice_vectors[1] * nb,
            self.lattice_vectors[2] * nc
        ])
        return Lattice(new_vectors, new_species, new_frac)