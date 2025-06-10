import numpy as np
import matplotlib.pyplot as plt

def plot_atoms(lattice, species=None, ax=None, show_cell=True, atom_size=300, colours=None):
    """
    Plots atomic positions in Cartesian coordinates.

    Parameters
    ----------
    lattice : Lattice
        An instance of the Lattice class.
    species : list[str], optional
        A list of chemical species for each atom. Must match the length of lattice.positions.
    ax : matplotlib.axes._subplots.Axes3DSubplot, optional
        An existing 3D axis. If None, a new figure and axis are created.
    show_cell : bool
        Whether to draw the unit cell box.
    atom_size : float
        Size of the atom markers.
    colours : dict[str, str], optional
        Mapping of chemical species to matplotlib colour codes.
    """
    cart_coords = lattice.get_cartesian_positions()
    lattice_vectors = lattice.lattice_vectors

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

    if species is None:
        species = ['X'] * len(cart_coords)

    if colours is None:
        colours = {}
    unique_species = sorted(set(species))
    default_palette = plt.cm.tab10.colors
    for i, sp in enumerate(unique_species):
        if sp not in colours:
            colours[sp] = default_palette[i % len(default_palette)]

    for coord, sp in zip(cart_coords, species):
        ax.scatter(*coord, s=atom_size, c=[colours[sp]], label=sp)

    if show_cell:
        draw_unit_cell(ax, lattice_vectors)

    # Aesthetic tweaks
    ax.set_box_aspect([1, 1, 1])
    ax.set_xlabel('x (Å)')
    ax.set_ylabel('y (Å)')
    ax.set_zlabel('z (Å)')
    ax.legend()
    plt.tight_layout()
    return ax

def draw_unit_cell(ax, lattice_vectors):
    """
    Draws the edges of the unit cell given lattice vectors.
    """
    origin = np.zeros(3)
    a1, a2, a3 = lattice_vectors
    corners = [
        origin,
        a1, a2, a3,
        a1 + a2, a1 + a3, a2 + a3,
        a1 + a2 + a3
    ]
    corners = np.array(corners)

    # Define edges by index pairs
    edges = [
        (0,1), (0,2), (0,3),
        (1,4), (1,5), (2,4), (2,6),
        (3,5), (3,6),
        (4,7), (5,7), (6,7)
    ]
    for start, end in edges:
        xs, ys, zs = zip(corners[start], corners[end])
        ax.plot(xs, ys, zs, color='k', linewidth=0.5, alpha=0.7)
