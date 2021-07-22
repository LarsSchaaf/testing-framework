"""
This is a generic relax script that takes two arguments: .xyz file to relax + index

Using the mEspresso class + the relaxation alorithm for the vacancy


Usefull links for relaxations:

- Relaxing also the cell: Exp Cell Filter

https://wiki.fysik.dtu.dk/ase/ase/constraints.html#the-expcellfilter-class

- Fixing symetry: Here is an example of using these tools to demonstrate the difference
 between minimising a perturbed bcc Al cell with and without symmetry-preservation.
 Since bcc is unstable with respect to fcc with a Lennard Jones model, the unsymmetrised
 case relaxes to fcc, while the constraint keeps the original symmetry.

    (https://wiki.fysik.dtu.dk/ase/ase/constraints.html#the-fixsymmetry-class)
    https://gitlab.com/ase/ase/-/blob/master/doc/ase/fix_symmetry_example.py
"""

from gapalysis.espresso.calculator import calculator as calc_qe_dft


# ----- For Relaxation
import ase
import ase.io
import os
import numpy as np
from ase.constraints import StrainFilter, ExpCellFilter

from ase.optimize.precon import PreconLBFGS
from ase.spacegroup.symmetrize import FixSymmetry, check_symmetry
from gapalysis.scripts.ase.run_relax import print_status, log_prop, write_frame


import logging
from gapalysis.utils.logging_utils import initialize_logger
import sys


# ---- Setting Variables

MAX_STEPS = 30
FMAX = 0.005
FNAME = sys.argv[1]
INDEX = sys.argv[2]
OPT_CELL = False
FIX_SYMMETRY = True

END_FNAME = (
    os.path.splitext(os.path.split(FNAME)[1])[0] + "_relaxed_qe.xyz"
)  # Currently in run_dir, eg. /data/test/tu.pdf -> ./tu_relaxed_qe.xyz

OUT_FNAME = (
    os.path.splitext(os.path.split(FNAME)[1])[0] + "_relaxing_traj_qe.xyz"
)  # .xyz file where itteration of relaxations are saved

# --- Relaxation algorithm

if __name__ == "__main__":

    # Logging
    log = logging.getLogger("QE-relax")
    initialize_logger(log, log_name="QE-relax")
    log.info("Logging setup was successful.")

    log.info(f"END_FNAME: {END_FNAME}")
    log.info(f"OUT_FNAME: {OUT_FNAME}")

    # Config
    base_tags = ["MD", "slab", "relax_qe"]

    config = {
        "name": "Relax QE - Generic",
        "in_fname": FNAME,
        "tags": [*base_tags, "itt1", "<50"],
    }

    # Load Atoms
    ats = ase.io.read(FNAME, INDEX)
    atoms = ats.copy()
    atoms.wrap()

    log.info(f"Initial Symmetry: {check_symmetry(atoms, verbose=False)}")

    atoms.set_calculator(calc_qe_dft)

    if FIX_SYMMETRY:
        log.info("Fixing symmetry")
        atoms.set_constraint(FixSymmetry(atoms))

    if OPT_CELL:
        # log.info("Optimise lattice parameters using StrainFilter")
        # atoms = StrainFilter(atoms.copy())
        log.info("Optimise lattice parameters using ExpCellFilter")
        atoms_run = ExpCellFilter(atoms)
    else:
        # Else same as atoms object.
        # This is done because ExpCellFilter is not an atoms object so doesn't have
        # .write property for example

        atoms_run = atoms

    # https://wiki.fysik.dtu.dk/ase/tutorials/lattice_constant.html

    # Setup + Run Dynacics
    dyn = PreconLBFGS(atoms_run, restart=OUT_FNAME[:-4] + ".pckl")
    dyn.attach(write_frame, interval=1, dyn=dyn, fname=OUT_FNAME)
    dyn.run(fmax=FMAX, steps=MAX_STEPS)

    # Process Results
    tenergy = [atoms.get_total_energy()]
    relaxed_ats = atoms.copy()

    log.info(f"Final Symmetry: {check_symmetry(atoms, verbose=False)}")
    log.info(f"Saving end Energy / Writing relaxed structure to file {END_FNAME}")
    np.savetxt("total_energies.txt", tenergy)
    ase.io.write(END_FNAME, relaxed_ats)

else:
    print("Relax script imported rather than run! __name__ !=__main__ ")
