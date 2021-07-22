from espresso.espresso import Espresso
from ase.optimize import QuasiNewton
from ase import units
from espresso.siteconfig import SiteConfig
from glob import glob
from ase.optimize import LBFGS
from ase.io import read, write
from ase.io.espresso import kspacing_to_grid
import sys, os
import numpy as np

# cwd = os.getcwd()
# psppath = os.path.join(cwd, '../psp')
psppath = "/gpfs/nobackup/projects/quantumchemistry/rom_cq/SOAPML/data/In2O3/espresso_calculations/psp_copy/"
# psppath='/gpfs/nobackup/projects/quantumchemistry/rom_cq/SOAPML/data/In2O3/espresso_calculations/psp/'
traj = read(sys.argv[-1], index=":")
# traj = read('/gpfs/nobackup/projects/quantumchemistry/rom_cq/SOAPML/data/In2O3/espresso_calculations/LARS_inF_1/inF_1_traj.xyz', index=':')

site = SiteConfig("PBS", scratchenv="SCRATCH")
# site=SiteConfig(scheduler=None)


for atoms in traj:

    info = {"uid": atoms.info["uid"], "type": atoms.info["type"]}

    prefix = "{}_{}".format(info["type"], info["uid"])
    print(prefix)
    if len(glob("./" + prefix + "*.pbs")) > 0:
        print("skipping", prefix)
        continue

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

    if "slab+ads" in info["type"]:
        calcstress = False
        print("type {}: dipole correction on".format(info["type"]))
        dipole = {"status": True}
        kpts[2] = 1

    print("using kpts:", kpts)

    magmom_dict = {"In": 10, "O": 4, "C": 3, "H": 1, "N": 3}
    magmoms = []
    for i in range(len(atoms)):

        magmoms.append(magmom_dict[atoms[i].symbol])

    atoms.set_initial_magnetic_moments(magmoms)

    if "isolated_atom" in info["type"]:
        calcstress = False
        spinpol = False
        kpts = "gamma"
        atoms.arrays.pop("initial_magmoms")

    calc = Espresso(
        pw=65 * units.Ry,  # plane-wave cut-off in eV, defaults to 350.0
        dw=8 * 65 * units.Ry,  # default 10*pw
        xc=xc,
        kpts=kpts,
        calculation="scf",
        psppath=psppath,
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
            "removesave": True,
            "avoidio": True,
            "removewf": True,
            "wf_collect": False,
        },
        outdir=prefix,
        site=site,
        txt="pw.out",
    )

    try:
        atoms.set_calculator(calc)
        print(atoms.get_potential_energy())
        atoms.info = info
        write("{}.xyz".format(prefix), [atoms], format="extxyz")
        atoms = read("{}.xyz".format(prefix))
        info = atoms.info
        if "stress" in atoms.info.keys():
            info["virial"] = -atoms.info["stress"] / atoms.get_volume()
            atoms.info = info
            print("writing new")
            write("{}.xyz".format(prefix), [atoms], format="extxyz")

    except:
        open("{}.FAIL".format(prefix), "w").close()  ######
        print("error", prefix)
