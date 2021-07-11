"""
Taken from the testing-framework: phonons.py

"""
import sys

# from utilities import *
import numpy as np
import phonopy
import ase.units
import ase.io
import ase
import os

from utilities import run_root, evaluate

FILE_LABEL = "xyz"


def do_phonons(
    bulk_struct_tests,
    n_supercell,
    band_paths=None,
    dx=0.01,
    args={"SETUP": True, "CALCULATE": True, "PROCESS": True},
):
    """
    Run Calculations to determine phonon dispersion

    This funciton uses phonopy to:
        1. Setup : Create all the required displacements
        2. Perform Calculations
        3. Record results in properties dir

    Args:
        bulk_struct_tests (list): list of strings, containining the name of the bulk strucutre, such as 'bulk_In2O3_Ia3'
        n_supercell (int): size of supercell, eg. 2 corresponds to [2 2 2] supercell
        band_paths (list, optional): list of special points to use for band path. Defaults to None.
        dx (float, optional): Size of displacements for phonopy. Defaults to 0.01.
        args (dict, optional): which part of analysis to perform. Defaults to {"SETUP": True, "CALCULATE": True, "PROCESS": True}.

    Raises:
        RuntimeError:
        RuntimeError:

    Returns:
        properties: dictionary of resulting properties containining {??} keys
    """

    # Check if band_paths same length as bulk_struc_tests
    if band_paths is not None and len(band_paths) != len(bulk_struct_tests):
        raise RuntimeError(
            "got {} bulk structs but different {} band paths".format(
                len(bulk_struct_tests), len(band_paths)
            )
        )
    properties = {}
    for bulk_i, bulk_struct_test in enumerate(bulk_struct_tests):

        to_dir = run_root
        at0 = get_relaxed_bulk(bulk_struct_test)

        # magnetic moments could change the symmetry, ignored here for now
        phonopy_atoms = phonopy.structure.atoms.PhonopyAtoms(
            symbols=at0.get_chemical_symbols(),
            scaled_positions=at0.get_scaled_positions(),
            masses=at0.get_masses(),
            cell=at0.get_cell(),
        )

        phonons = phonopy.Phonopy(
            phonopy_atoms, np.diag([n_supercell] * 3), factor=phonopy.units.VaspToTHz
        )
        phonons.generate_displacements(distance=dx)

        # convert from chosen Phonopy units (THz) to cm^-1
        THz_per_invcm = ase.units._c * 1.0e-10

        ####################################################################################################
        if args["SETUP"]:
            sys.stderr.write("SETUP\n")

            displ_supercells = phonons.get_supercells_with_displacements()

            at0_sc = at0 * n_supercell
            # reorder at0_sc to match order from phonopy
            at = ase.Atoms(
                pbc=True,
                cell=displ_supercells[0].get_cell(),
                positions=displ_supercells[0].get_positions(),
                numbers=displ_supercells[0].get_atomic_numbers(),
            )
            matched_pos = np.zeros(at.positions.shape)
            mapping = [-1] * len(at)
            at0_sc_scaled_pos = at0_sc.get_scaled_positions()
            at_scaled_pos = at.get_scaled_positions()
            for at0_sc_i in range(len(at)):
                scaled_dists = at0_sc_scaled_pos[at0_sc_i] - at_scaled_pos
                scaled_dists -= np.round(scaled_dists)
                closest_i = np.argmin(np.linalg.norm(scaled_dists, axis=1))
                matched_pos[at0_sc_i] = at.positions[closest_i]
                mapping[closest_i] = at0_sc_i
            if -1 in mapping:
                raise RuntimeError("Failed to map orig and displaced atom positions")
            at0_sc = at0_sc[mapping]

            to_fname = os.path.join(to_dir, "UNDISPL.{}".format(FILE_LABEL))
            ase.io.write(to_fname, at0_sc)

            # create displaced cells
            sys.stderr.write(
                "Creating {} displacements\n".format(len(displ_supercells))
            )
            at_sets = []
            for (displ_i, displ_config) in enumerate(displ_supercells):
                at = ase.Atoms(
                    pbc=True,
                    cell=displ_config.get_cell(),
                    positions=displ_config.get_positions(),
                    numbers=displ_config.get_atomic_numbers(),
                )
                at_sets.append(at)
                to_fname = os.path.join(
                    to_dir, "DISPL_{}.{}".format(displ_i, FILE_LABEL)
                )
                ase.io.write(to_fname, at)

        ####################################################################################################
        all_forces = []
        if args["CALCULATE"]:
            sys.stderr.write("CALCULATE\n")

            sys.stderr.write("Calculating for {} displacements\n".format(len(at_sets)))

            evaluate(at0_sc)
            f0 = at0_sc.get_forces()

            to_fname = os.path.join(to_dir, "FORCES.UNDISPL.{}".format("txt"))
            np.savetxt(to_fname, f0)

            for (displ_i, at) in enumerate(at_sets):
                at0_sc.set_positions(at.positions)
                at0_sc.set_cell(at.get_cell())
                sys.stderr.write("start evaluate displacement {}\n".format(displ_i))
                evaluate(at0_sc)
                f = at0_sc.get_forces()
                to_fname = os.path.join(
                    to_dir, "CALC_DISPL_{}.{}".format(displ_i, FILE_LABEL)
                )
                print("To force calculated fname:")
                ase.io.write(to_fname, at0_sc)

                all_forces.append(f)

                to_fname = os.path.join(
                    to_dir, "FORCES.DISPL_{}.{}".format(displ_i, "txt")
                )
                np.savetxt(to_fname, all_forces[-1])

        ####################################################################################################
        if args["PROCESS"]:
            # This has been adjusted to read previously command line run dirs
            sys.stderr.write("PROCESS\n")

            Nat = all_forces[0].shape[0]
            for f in all_forces:
                f -= f0
                f -= np.outer(np.ones(Nat), np.sum(f, axis=0) / Nat)
            all_forces = np.asarray(all_forces)

            properties[bulk_struct_test] = {
                "dx": dx,
                "all_forces": all_forces.tolist(),
                "symb": at0.get_chemical_symbols(),
                "scaled_pos": at0.get_scaled_positions().tolist(),
                "m": at0.get_masses().tolist(),
                "c": at0.get_cell().tolist(),
                "n_cell": np.diag([n_supercell] * 3).tolist(),
                "unit_factor": phonopy.units.VaspToTHz,
                "band_path": band_paths[bulk_i],
            }

        ####################################################################################################

    return properties
