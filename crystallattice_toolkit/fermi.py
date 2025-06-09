"""
Fermi-surface reconstruction from electronic band data.
"""
import numpy as np
from scipy.interpolate import griddata
from skimage.measure import marching_cubes
from .reciprocal import ReciprocalLattice

class FermiSurface:
    """
    Generate a Fermi isosurface in k-space.
    """
    def __init__(
        self,
        reciprocal_lattice: ReciprocalLattice,
        kpoints: np.ndarray,
        energies: np.ndarray
    ):
        """
        Args:
            reciprocal_lattice: ReciprocalLattice instance.
            kpoints: M×3 array of k-point coordinates.
            energies: M×B array (or M,) of band energies.
        """
        self.reciprocal = reciprocal_lattice
        self.kpoints = kpoints
        self.energies = energies

    def interpolate_on_grid(
        self,
        resolution: int = 50
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Interpolate energy on a regular 3D grid.

        Returns:
            grid_x, grid_y, grid_z: 3D coordinate arrays.
            grid_e: energy values on grid.
        """
        # Define grid in reciprocal cell
        verts = self.reciprocal.basis
        # Create meshgrid in fractional coordinates
        lin = np.linspace(0, 1, resolution)
        fx, fy, fz = np.meshgrid(lin, lin, lin, indexing='ij')
        pts_frac = np.vstack([fx.ravel(), fy.ravel(), fz.ravel()]).T
        pts_cart = pts_frac.dot(verts)
        # Interpolate
        grid_e = griddata(self.kpoints, self.energies, pts_cart, method='linear')
        grid_e = grid_e.reshape((resolution, resolution, resolution))
        return fx, fy, fz, grid_e

    def mesh_isosurface(
        self,
        fermi_level: float,
        resolution: int = 50
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Extract mesh (verts, faces) of constant-energy surface at Fermi level.
        """
        fx, fy, fz, grid_e = self.interpolate_on_grid(resolution)
        verts, faces, normals, values = marching_cubes(grid_e, level=fermi_level)
        # Convert fractional vertices to Cartesian
        verts_cart = np.dot(verts / (resolution-1), self.reciprocal.basis)
        return verts_cart, faces, normals