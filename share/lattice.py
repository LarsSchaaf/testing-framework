import numpy as np
from utilities import relax_config, run_root
import ase.io, sys, os.path
from ase.optimize.precon import PreconLBFGS
import matscipy.elasticity
from ase.units import GPa


import os
import logging
import logging.handlers
import pathlib

DATA_DIR = "/data/lls34/"

LOG_DIR = pathlib.Path(os.path.join(DATA_DIR, "logs/"))


def initialize_logger(logger, log_name=None):
    """
    Initialize the given `logger` to log to files in the given `log_name` directory.
    If log_name is None, then no logs are stored.
    """
    # Create Main Log directory
    if os.path.exists(DATA_DIR) and not os.path.exists(LOG_DIR):
        os.makedirs(LOG_DIR, exist_ok=True)

    # Set up general logger and formatting
    logger.setLevel(logging.DEBUG)  # listen to everything
    formatter = logging.Formatter(
        "%(asctime)s [%(name)-10s] [%(threadName)-12s] %(message)s [in %(pathname)s:%(lineno)d]"
    )

    # If file_name provided, add file handler for logs
    if log_name:
        if not os.path.exists(LOG_DIR / f"{log_name}"):
            os.makedirs(LOG_DIR / f"{log_name}", exist_ok=True)

        file_handler = logging.handlers.RotatingFileHandler(
            LOG_DIR / f"{log_name}/{log_name}.log", maxBytes=1000240, backupCount=5
        )
        file_handler.setFormatter(formatter)
        file_handler.setLevel(logging.DEBUG)
        logger.addHandler(file_handler)

    # Add console handler for logs:
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.DEBUG)
    console_formatter = logging.Formatter("[%(name)-10s] %(message)s")
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)


# logging.basicConfig(
#     filename=__file__[-3] + ".txt",
#     filemode="a",
#     format="%(asctime)s %(msecs)d- %(process)d-%(levelname)s - %(message)s",
#     datefmt="%d-%b-%y %H:%M:%S %p",
#     level=logging.DEBUG,
# )

log = logging.getLogger(__name__)
initialize_logger(log, log_name=__name__)

# Calculate B for Hexagonal, Tetragonal & Trigonal
# DOI: 10.1103/PhysRevB.77.104118 Eq. (29)
def HTT_B(c11, c33, c12, c13):
    numerator = c33 * (c11 + c12) - 2 * c13 ** 2
    denominator = c11 + c12 + 2 * c33 - 4 * c13

    return numerator / denominator


# Mater Trans v 53 p 1247 (2012), Eq. 4-10
def VRH_B(c11, c33, c12, c13, c44, c66):
    Bv = (1.0 / 9.0) * (2.0 * (c11 + c12) + c33 + 4.0 * c13)

    M = c11 + c12 + 2.0 * c33 - 4.0 * c13
    Csq = (c11 + c12) * c33 - 2.0 * c13 ** 2
    Br = Csq / M

    return (Bv + Br) / 2.0


def calc_E_vs_V(
    bulk, dV=0.025, n_steps=(-10, 10), tol=1e-2, method="lbfgs"
):  # hack tol to deal with Te C2/m
    """
    Run the E vs V analysis

    Changes the cell size isometrically to scan for different volumes and scales
    the possitions accordingly. Each structure is then relaxed. The trajectory,
    unrelaxed and relaxed .xyz files are saved in the run_folder

    Args:
        bulk (ase.atoms): the starting atomic positions
        dV (float, optional): Increment of volumes V0*dV. Defaults to 0.025.
        n_steps (tuple, optional): Number of dV steps in each direction.
            Defaults to (-10, 10).
        tol (float, optional): relaxation tolerance. Defaults to 1e-2.
        method (str, optional): relaxation method. Defaults to "lbfgs".

    Returns:
        list: results of E vs V investigation (Volume/num_at, Energy/num_at, Stress)
    """
    import model

    log.debug(f"Starting calc_E_vs_V for {model.name}")

    V0 = bulk.get_volume()
    dV *= V0
    E_vs_V = []

    scaled_bulk = bulk.copy()

    # Firstly Reduce Volume
    log.debug("Runs for decreasing Volume: ")

    for i in range(0, n_steps[0] - 1, -1):
        log.info(f"Decreased volume {i}/{n_steps[0]}")
        V_cur = scaled_bulk.get_volume()

        # Increase volume with scaled atoms
        scaled_bulk.set_cell(
            scaled_bulk.get_cell() * ((V0 + i * dV) / V_cur) ** (1.0 / 3.0),
            scale_atoms=True,
        )
        # ase.io.write(sys.stdout, scaled_bulk, format="extxyz")
        # print("trying to relax i", i)
        try:
            if hasattr(model, "fix_cell_dependence"):
                model.fix_cell_dependence(scaled_bulk)

            # Writing unlrelaxed
            to_fname = run_root + "-E_vs_V_%03d-unrelaxed.xyz" % i
            log.debug(f"    - writing unrelaxed {to_fname}")
            ase.io.write(
                to_fname, scaled_bulk, format="extxyz",
            )

            # Relax config
            log.debug(f"    - About to relax using relax_config() call")
            scaled_bulk = relax_config(
                scaled_bulk,
                relax_pos=True,
                relax_cell=True,
                tol=tol,
                max_steps=200,
                save_traj=True,
                constant_volume=True,
                method=method,
                refine_symmetry_tol=1.0e-1,
                keep_symmetry=True,
                config_label="E_vs_V_%03d" % i,
                from_base_model=True,
                save_config=True,
            )
        except Exception as e:
            msg = "WARNING: failed config in calc_E_vs_V", str(e)
            log.error(msg)
            print(msg)
            sys.exit(1)  #### NB
            break

        # ase.io.write(sys.stdout, scaled_bulk, format="extxyz")

        log.debug(f"    - Appending E_vs_V dictionary")
        E_vs_V.insert(
            0,
            (
                scaled_bulk.get_volume() / len(scaled_bulk),
                scaled_bulk.get_potential_energy(force_consistent=True) / len(bulk),
                list(scaled_bulk.get_stress()),
            ),
        )

    # Increase Volume

    log.debug("Runs for increasing Volume: ")
    scaled_bulk = bulk.copy()
    for i in range(1, n_steps[1] + 1):
        log.debug(f"    Increased volume {i}/{n_steps[1]}")
        V_cur = scaled_bulk.get_volume()
        scaled_bulk.set_cell(
            scaled_bulk.get_cell() * ((V0 + i * dV) / V_cur) ** (1.0 / 3.0),
            scale_atoms=True,
        )
        ase.io.write(sys.stdout, scaled_bulk, format="extxyz")
        print("trying to relax i", i)
        try:
            if hasattr(model, "fix_cell_dependence"):
                model.fix_cell_dependence(scaled_bulk)
            ase.io.write(
                run_root + "-E_vs_V_%02d-unrelaxed.xyz" % i,
                scaled_bulk,
                format="extxyz",
            )
            scaled_bulk = relax_config(
                scaled_bulk,
                relax_pos=True,
                relax_cell=True,
                tol=tol,
                max_steps=200,
                save_traj=True,
                constant_volume=True,
                method=method,
                refine_symmetry_tol=1.0e-1,
                keep_symmetry=True,
                config_label="E_vs_V_%02d" % i,
                from_base_model=True,
                save_config=True,
            )
        except Exception as e:
            print("failed", str(e))
            break
        ase.io.write(sys.stdout, scaled_bulk, format="extxyz")
        E_vs_V.append(
            (
                scaled_bulk.get_volume() / len(scaled_bulk),
                scaled_bulk.get_potential_energy(force_consistent=True) / len(bulk),
                list(scaled_bulk.get_stress()),
            )
        )

    if hasattr(model, "fix_cell_dependence"):
        model.fix_cell_dependence()

    return E_vs_V


def do_lattice(
    test_dir,
    lattice_type,
    dV=0.025,
    n_steps=(-10, 10),
    tol=1.0e-2,
    method="lbfgs",
    applied_P=0.0,
):

    import model

    bulk = ase.io.read(test_dir + "/bulk.xyz", format="extxyz")

    results_dict = {}

    print("relax bulk")
    # relax the initial unit cell and atomic positions
    (orig_cell, new_cell) = (None, None)
    while (
        new_cell is None
        or np.max(np.abs(np.dot(np.linalg.inv(new_cell), orig_cell) - np.eye(3))) > 0.05
    ):
        if hasattr(model, "fix_cell_dependence"):
            model.fix_cell_dependence(bulk)
        orig_cell = bulk.get_cell()
        bulk = relax_config(
            bulk,
            relax_pos=True,
            relax_cell=True,
            tol=tol,
            save_traj=True,
            method=method,
            refine_symmetry_tol=1.0e-2,
            keep_symmetry=True,
            config_label="bulk",
            from_base_model=True,
            save_config=True,
            applied_P=applied_P,
        )
        new_cell = bulk.get_cell()
        if hasattr(model, "fix_cell_dependence"):
            model.fix_cell_dependence()
        else:
            break

    print("final relaxed bulk")
    ase.io.write(sys.stdout, bulk, format="extxyz")
    ase.io.write(os.path.join("..", run_root + "-relaxed.xyz"), bulk, format="extxyz")

    print("calculating E vs. V")
    E_vs_V = calc_E_vs_V(bulk, dV=dV, n_steps=n_steps, tol=tol)
    results_dict.update({"E_vs_V": E_vs_V})

    print("calculating elastic constants")

    if hasattr(model, "fix_cell_dependence"):
        model.fix_cell_dependence(bulk)

    opt = lambda atoms, **kwargs: PreconLBFGS(atoms, **kwargs)
    if lattice_type == "cubic":
        elastic_consts = matscipy.elasticity.fit_elastic_constants(
            bulk, symmetry="cubic", optimizer=opt, logfile=sys.stdout
        )
        c11 = elastic_consts[0][0, 0] / GPa
        c12 = elastic_consts[0][0, 1] / GPa
        c44 = elastic_consts[0][3, 3] / GPa
        results_dict.update(
            {"c11": c11, "c12": c12, "c44": c44, "B": (c11 + 2.0 * c12) / 3.0}
        )
    elif lattice_type == "orthorhombic":
        elastic_consts = matscipy.elasticity.fit_elastic_constants(
            bulk, optimizer=opt, logfile=sys.stdout
        )
        c11 = elastic_consts[0][0, 0] / GPa
        c22 = elastic_consts[0][1, 1] / GPa
        c33 = elastic_consts[0][2, 2] / GPa
        c12 = elastic_consts[0][0, 1] / GPa
        c13 = elastic_consts[0][0, 2] / GPa
        c23 = elastic_consts[0][1, 2] / GPa
        c44 = elastic_consts[0][3, 3] / GPa
        c55 = elastic_consts[0][4, 4] / GPa
        c66 = elastic_consts[0][5, 5] / GPa
        results_dict.update(
            {
                "c11": c11,
                "c22": c22,
                "c33": c33,
                "c12": c12,
                "c13": c13,
                "c23": c23,
                "c44": c44,
                "c55": c55,
                "c66": c66,
            }
        )
    elif lattice_type == "tetragonal":
        elastic_consts = matscipy.elasticity.fit_elastic_constants(
            bulk, symmetry="tetragonal_high", optimizer=opt, logfile=sys.stdout
        )
        c11 = elastic_consts[0][0, 0] / GPa
        c33 = elastic_consts[0][2, 2] / GPa
        c12 = elastic_consts[0][0, 1] / GPa
        c13 = elastic_consts[0][0, 2] / GPa
        c44 = elastic_consts[0][3, 3] / GPa
        c66 = elastic_consts[0][5, 5] / GPa
        results_dict.update(
            {
                "c11": c11,
                "c33": c33,
                "c12": c12,
                "c13": c13,
                "c44": c44,
                "c66": c66,
                "B": VRH_B(c11, c33, c12, c13, c44, c66),
            }
        )
    elif lattice_type == "hexagonal":
        # Need to check if hexagonal structures are truly trigonal_high
        # symmetry=triginal_high not hexagonal until matscipy is debugged
        elastic_consts = matscipy.elasticity.fit_elastic_constants(
            bulk, symmetry="trigonal_high", optimizer=opt, logfile=sys.stdout
        )
        c11 = elastic_consts[0][0, 0] / GPa
        c33 = elastic_consts[0][2, 2] / GPa
        c12 = elastic_consts[0][0, 1] / GPa
        c13 = elastic_consts[0][0, 2] / GPa
        c44 = elastic_consts[0][3, 3] / GPa
        c14 = elastic_consts[0][0, 3] / GPa
        c15 = elastic_consts[0][0, 4] / GPa
        c25 = elastic_consts[0][1, 4] / GPa
        c66 = elastic_consts[0][5, 5] / GPa
        results_dict.update(
            {
                "c11": c11,
                "c33": c33,
                "c12": c12,
                "c13": c13,
                "c44": c44,
                "c14": c14,
                "c15": c15,
                "c25": c25,
                "c66": c66,
                "B": HTT_B(c11, c33, c12, c13),
            }
        )
    elif lattice_type == "trigonal":
        elastic_consts = matscipy.elasticity.fit_elastic_constants(
            bulk, symmetry="trigonal_high", optimizer=opt, logfile=sys.stdout
        )
        c11 = elastic_consts[0][0, 0] / GPa
        c33 = elastic_consts[0][2, 2] / GPa
        c12 = elastic_consts[0][0, 1] / GPa
        c13 = elastic_consts[0][0, 2] / GPa
        c44 = elastic_consts[0][3, 3] / GPa
        c14 = elastic_consts[0][0, 3] / GPa
        c15 = elastic_consts[0][0, 4] / GPa
        c25 = elastic_consts[0][1, 4] / GPa
        c66 = elastic_consts[0][5, 5] / GPa
        results_dict.update(
            {
                "c11": c11,
                "c33": c33,
                "c12": c12,
                "c13": c13,
                "c44": c44,
                "c14": c14,
                "c15": c15,
                "c25": c25,
                "c66": c66,
                "B": HTT_B(c11, c33, c12, c13),
            }
        )

    if hasattr(model, "fix_cell_dependence"):
        model.fix_cell_dependence()

    return results_dict
