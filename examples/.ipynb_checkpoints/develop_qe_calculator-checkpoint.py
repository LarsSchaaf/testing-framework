"""
Python script for running QE with the github calculator

This was originally created by Sandip on Feb 2021. The code requires the espresso
.esspresso to be installed. This was modified to work on womble0 using the GridEngine

The script can be called by

`python qe_simple.py traj.xyz 2`

where 2 is the configuraiton in traj.xyz that should be calculated

 - still need to do out directory
"""

from espresso.espresso import Espresso  # Espresso Library
from ase import units  # Ase Units
from espresso.siteconfig import SiteConfig
from ase.io import read, write
from ase.io.espresso import kspacing_to_grid
import sys
import os
import numpy as np
import ase
import ase.io


# ---- Setting Variables

# Environements
os.environ["SCRATCH"] = "/home/lls34/trash/scratch/"

atoms = ase.io.read(
    "/home/lls34/GitHub/01_PhD/PhD_Code/submodules/testing-framework/testing-framework/examples/In_mp-1055994_primitive.cif"
)
atoms.center(vacuum=5)
os.environ[
    "SGE_O_WORKDIR"
] = "/home/lls34/GitHub/01_PhD/PhD_Code/submodules/testing-framework/testing-framework/examples"
os.environ["NSLOTS"] = "1"


PSP_PATH = (
    "/home/lls34/programs/QuantumExpresso/pseudopotentials/use_pot"  # Pseudo Potentails
)

# Sceduling Properties
site = SiteConfig("GE", scratchenv="SCRATCH")
site = SiteConfig(scheduler=None)


# Atom specific properties

atoms = atoms

# info = {"uid": atoms.info["uid"], "type": atoms.info["type"]}

prefix = "test"  # "QE_{}_0".format(info["uid"])

j = 0
# Make sure unique output file
while os.path.exists(os.path.join(prefix)) or os.path.exists(
    os.path.join(prefix, ".xyz")
):
    prefix = "QE_{}_{}".format(info["uid"], j)
    j += 1

    if j > 1000:
        raise ValueError("While loop ran too often")
print(f"New prefix {prefix}")

# Now Properites


kspacing = 0.25 / (2.0 * np.pi)

kpts = kspacing_to_grid(atoms, spacing=0.25 / (2.0 * np.pi))
xc = "PBE"
dipole = {"status": False}
calcstress = True
spinpol = True

convergence = {
    "energy": 1e-7,
    "mixing": 0.5,
    "nmix": 10,
    "maxsteps": 500,
    "diag": "david",
}

# if "Slab" in info["type"]:
#     calcstress = False
#     print("type {}: dipole correction on".format(info["type"]))
#     dipole = {"status": True}
#     kpts[2] = 1

print(f"using kpts: {kpts}")

magmom_dict = {"In": 10, "O": 4, "C": 3, "H": 1, "N": 3}
magmoms = []

for i in range(len(atoms)):
    magmoms.append(magmom_dict[atoms[i].symbol])

atoms.set_initial_magnetic_moments(magmoms)

# if "isolated_atom" in info["type"]:
#     calcstress = False
#     spinpol = False
#     kpts = "gamma"
#     atoms.arrays.pop("initial_magmoms")

calc = Espresso(
    # pseudopotentials=pseudopotentials,
    pw=65 * units.Ry,  # plane-wave cut-off in eV, defaults to 350.0
    dw=8 * 65 * units.Ry,  # default 10*pw
    xc=xc,
    kpts=kpts,
    calculation="scf",
    psppath=PSP_PATH,
    parflags="-npool 1",
    nbands=-20,
    occupations="smearing",
    sigma=0.1,
    smearing="gaussian",
    convergence=convergence,
    dipole=dipole,
    calcstress=calcstress,
    spinpol=spinpol,
    # vdw_corr='D3',
    output={
        "removesave": False,
        "avoidio": False,
        "removewf": False,
        "wf_collect": True,
    },
    outdir=prefix,
    site=site,
    txt="pw.out",
)

try:
    atoms.set_calculator(calc)
    print(f"PE: {atoms.get_potential_energy()}")
    # atoms.info = info
    write("{}.xyz".format(prefix), [atoms], format="extxyz")
    atoms = read("{}.xyz".format(prefix))
    # info = atoms.info
    if "stress" in atoms.info.keys():
        info["virial"] = -atoms.info["stress"] / atoms.get_volume()
        atoms.info = info
        print("writing new")
        write("{}.xyz".format(prefix), [atoms], format="extxyz")

except ImportError as e:
    open("{}.FAIL".format(prefix), "w").close()
    print(f"Error: {prefix}, {e}")
