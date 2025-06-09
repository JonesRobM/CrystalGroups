"""
Reciprocal lattice and first Brillouin zone computation.
"""
import numpy as np
from scipy.spatial import ConvexHull
from .lattice import Lattice

class ReciprocalLattice:
    """
    Represents the reciprocal-space lattice of a real-space Lattice.
    """
    def __init__(self, real_lattice: Lattice):
        """
        Compute reciprocal basis vectors.
        """
        a1, a2, a3 = real_lattice.lattice_vectors
        vol = np.dot(a1, np.cross(a2, a3))
        b1 = 2*np.pi * np.cross(a2, a3) / vol
        b2 = 2*np.pi * np.cross(a3, a1) / vol
        b3 = 2*np.pi * np.cross(a1, a2) / vol
        self.basis = np.vstack([b1, b2, b3])

    def get_reciprocal_points(self, shell: int=1) -> np.ndarray:
        """
        Generate reciprocal lattice points within given shell.
        """
        shifts = [np.array([i,j,k]) for i in range(-shell,shell+1)
                  for j in range(-shell,shell+1)
                  for k in range(-shell,shell+1)]
        pts = np.array([s.dot(self.basis) for s in shifts])
        return pts

    def first_brillouin_zone(self) -> np.ndarray:
        """
        Compute vertices of the first Brillouin zone via convex hull of Wigner–Seitz cell.

        Returns:
            vertices: P×3 array of BZ vertex coordinates.
        """
        pts = self.get_reciprocal_points(shell=1)
        # Voronoi cell at origin approximated by points within one shell
        hull = ConvexHull(pts)
        vertices = pts[hull.vertices]
        return vertices