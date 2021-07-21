import sys

import numpy as np
import phonopy
import ase.units
import ase.io
import ase
import os
from glob import glob

from testingframework.share.utilities import model_test_root, evaluate


def evaluate_all_xyz(dir, compare=False):
    """
    Takes all .xyz files inside the test dir and evaluates them using the calculator.

    To Do: Check wether we want force, energies, stesses ect.

    Args:
        dir (str): directory where .xyz files are stored
        compare (bool): wether to compare calculated values with values in .xyz file.
            This is usefull, when the saved .xyz file is a DFT calculation.
    Returns:
        dict: properites. Energy and forces for each file stored in json file
    """

    # Fname pattern
    fname_p = os.path.join(dir, "*.xyz")
    fnames = glob(fname_p)
    properties = {}

    # to_dir = model_test_root()
    to_dir = os.path.abspath(os.getcwd())
    print(f"To Dir: {to_dir}")

    for fname in fnames:

        name = os.path.splitext(os.path.split(fname)[-1])[0]

        ats = ase.io.read(fname)
        ats_org = ats.copy()
        evaluate(ats)

        fs = ats.get_forces()
        pe = ats.get_potential_energy()
        xyz = ats.get_positions()

        properties[name] = {
            "xyz": xyz.tolist(),
            "PE": pe,
            "all_forces": fs.tolist(),
        }

        ase.io.write(os.path.join(to_dir, name + ".xyz"), ats)

        if compare:

            df = np.linalg.norm(fs - ats_org.arrays["forces"], axis=1)

            properties[name].update(
                {"df": df.tolist(), "dPE": pe - ats_org.info["energy"]}
            )

    return properties
