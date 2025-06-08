## Project: CrystalLatticeToolkit

### Overview
A small Python library for constructing crystal lattices from space-group specifications, plotting atomic positions in Cartesian space, evaluating the first Brillouin zone, and reconstructing Fermi surfaces where applicable.

### Dependencies
```text
numpy
scipy
spglib  # for space-group operations (Le Bail & Digne, 2002)
pymatgen  # structure manipulation and reciprocal-lattice computations
matplotlib  # plotting
mayavi or plotly  # optional 3D visualisation of Brillouin zone and Fermi surface
``` 

### Repository Structure
```
crystallattice_toolkit/
├── README.md
├── setup.py
├── requirements.txt
├── crystallattice_toolkit/
│   ├── __init__.py
│   ├── lattice.py        # define unit cell, atomic positions
│   ├── symmetry.py       # wrappers around spglib for space groups
│   ├── reciprocal.py     # compute reciprocal lattice, BZ
│   ├── fermi.py          # Fermi-surface reconstruction
│   └── plot/
│       ├── __init__.py
│       ├── plot_atoms.py
│       ├── plot_bz.py
│       └── plot_fermi.py
└── examples/
    ├── generate_fcc.py
    └── plot_bz_bcc.py
```
