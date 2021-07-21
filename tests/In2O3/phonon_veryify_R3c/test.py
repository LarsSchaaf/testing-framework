from testingframework.share.evaluate import evaluate_all_xyz
import os
from gapalysis.phonons.quick_implementation import do_phonons  # , analyze_phonons
import numpy as np
import ase

xyz_file_dir = os.path.abspath(os.path.dirname(__file__))  # All in testing directory
properties0 = evaluate_all_xyz(xyz_file_dir, compare=True)

to_dir = os.path.abspath(os.getcwd())
bulk_fname = os.path.join(xyz_file_dir, "R3c.inp")

properties1, ph_init = do_phonons(
    bulk_fname=bulk_fname,
    bulk_struct_tests=["bulk_In2O3_R3c"],
    n_supercell=2,
    band_paths=[["G", "L", "B1", "B", "Z", "G", "X", "Q", "F", "P1", "Z", "L", "P"]],
    dx=0.02 * ase.units.Bohr,
    args={"SETUP": True, "CALCULATE": True, "PROCESS": True},
    to_dir="",
    load_force_fnames=os.path.join(to_dir, "R3c_phonons_supercell-00{}.xyz"),
    load_force_format=np.arange(5, dtype=int) + 1,
    un_displ_fname=os.path.join(to_dir, "R3c_phonons_supercell.xyz"),
)

properties = {"force_compare": properties0, "phonons": properties1}
