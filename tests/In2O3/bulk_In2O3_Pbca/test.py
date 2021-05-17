import os.path
import lattice

# The primitive unit cell was used for this
properties = lattice.do_lattice(
    os.path.abspath(os.path.dirname(__file__)), "orthorhombic"
)
