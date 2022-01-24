""" Horizontal Density of States plot

This script is based on the phonon analysis script writen by Cas + modified by lls34
Cas developped it for Si, lls34 used it for Indium Oxide.

Run using  `python ../../../scripts/analyze_phonons_vDOS.py -s In2O3`

"""

#!/usr/bin/env python
import numpy as np
from analyze_utils import *
import matplotlib

matplotlib.use("PDF")
from matplotlib import pyplot
import matplotlib.pyplot as plt
import phonopy
import ase.units
import sys, os

if len(sys.argv) < 4:
    TEST_NAME = "phonons_In2O3_R3c"
    TEST_NAME = "phonons_In2O3_Ia3"
else:
    TEST_NAME = sys.argv[3]
# TEST_NAME = "phonons_In2O3_Ia3_relaxed"


(args, models, tests, default_analysis_settings) = analyze_start([TEST_NAME])
data = read_properties(models, tests, args.test_set)
ref_model_name = "DFT_QE"  # default_analysis_settings["ref_model"]

print(f"Data keys: {data.keys()}")
print(f"RefModel: {data[ref_model_name]}")


def analyze_phonons(m_data):

    # Unit Conversion
    THz_per_invcm = ase.units._c * 1.0e-10
    at0 = ase.atoms.Atoms(
        symbols=m_data["symb"],
        scaled_positions=np.asarray(m_data["scaled_pos"]),
        masses=np.asarray(m_data["m"]),
        cell=np.asarray(m_data["c"]),
        pbc=[True] * 3,
    )

    # Create phonopy atoms based on .json file
    phonopy_atoms = phonopy.structure.atoms.PhonopyAtoms(
        symbols=m_data["symb"],
        scaled_positions=np.asarray(m_data["scaled_pos"]),
        masses=np.asarray(m_data["m"]),
        cell=np.asarray(m_data["c"]),
    )
    phonons = phonopy.Phonopy(
        phonopy_atoms, m_data["n_cell"], factor=m_data["unit_factor"]
    )
    phonons.generate_displacements(distance=m_data["dx"])
    phonons.produce_force_constants(
        forces=np.asarray(m_data["all_forces"]), calculate_full_force_constants=True
    )
    phonons.symmetrize_force_constants()

    # DOS
    n_dos_mesh = 6
    phonons.run_mesh([n_dos_mesh] * 3)
    phonons.run_total_dos()
    frequencies = phonons.get_total_dos_dict()["frequency_points"]
    PHdos = phonons.get_total_dos_dict()["total_dos"]
    # contains by columns the frequencies in cm^{-1} and the vDOS
    # the vDOS ins in units of "number of states/(unit cell x frequency[cm^{-1}])"
    # i.e. if you integrate the vDOS throughout frequency, it will be 3N, where N is the number of atoms in the unit cell
    m_data["DOS"] = {"freq": frequencies, "val": PHdos}

    # band path
    if "band_path" in m_data:
        band_path = m_data["band_path"]
        band_n = [20]
        lat = at0.get_cell().get_bravais_lattice()
        special_points = lat.get_special_points()
        print("special points:", special_points)
        path_pts = []
        path_labels = []
        pt = []

        # Cycle through all points in Band Path
        for s in " ".join(band_path).split():

            # Try and convert special points to floats
            # Can't use if == str, because they are all strings
            try:
                qi = float(s)
                pt.append(qi)

            # Except if can't convert str to floats
            except ValueError as e:
                # Check wether previous points where numbers if so through error
                if len(pt) > 0:
                    raise RuntimeError(f"got non-float after {len(pt)} floats")
                # Append path pts
                for p in s:
                    try:
                        path_pts.append(np.asarray(special_points[p]))
                        path_labels.append(p)
                    except KeyError:
                        raise RuntimeError(
                            "Failed to find special point {}, known {}".format(
                                p, special_points.keys()
                            )
                        )
            if len(pt) == 3:
                path_pts.append(np.asarray(pt))
                path_labels.append("{:.2f}_{:.2f}_{:.2f}".format(*pt))
                pt = []
        if len(band_n) == 1:
            band_n *= len(path_pts) - 1
        assert len(band_n) == len(path_pts) - 1
        q_pts = []
        q_pt_labels = []
        for (seg_i, (label, n, p0, p1)) in enumerate(
            zip(path_labels[0:-1], band_n, path_pts[0:-1], path_pts[1:])
        ):
            for p_i in range(n):
                x = float(p_i) / float(n)
                q_pts.append((1.0 - x) * p0 + x * p1)
                if p_i == 0:
                    q_pt_labels.append(label)
                else:
                    q_pt_labels.append(".")
        q_pts.append(path_pts[-1])
        q_pt_labels.append(path_labels[-1])
        phonons.run_band_structure([q_pts])
        bs = phonons.get_band_structure_dict()
        m_data["BAND_PATH"] = {
            "positions": bs["distances"][0],
            "frequencies": bs["frequencies"][0] / THz_per_invcm,
            "labels": q_pt_labels,
        }

    else:
        print("No band path found in data")


if ref_model_name in data and TEST_NAME in data[ref_model_name]:
    for (bulk_i, bulk_struct_test) in enumerate(data[ref_model_name][TEST_NAME]):
        print("analyze ref model", ref_model_name, bulk_struct_test)
        analyze_phonons(data[ref_model_name][TEST_NAME][bulk_struct_test])
bulk_struct_tests = []
for model_name in models:
    if model_name in data and TEST_NAME in data[model_name]:
        for bulk_struct_test in data[model_name][TEST_NAME]:
            if bulk_struct_test not in bulk_struct_tests:
                bulk_struct_tests.append(bulk_struct_test)

f, (ax_BP, ax_DOS) = plt.subplots(
    1, 2, gridspec_kw={"width_ratios": [5, 1]}, sharey=True, figsize=(6, 4)
)

# Selecting colors -> Change this to colorplot object
colors = ["red", "blue", "red", "green", "cyan"]

# Renaming for paper
D = {}
D["DFT_QE"] = "DFT"
D["GAP_it1_50s"] = "GAP"
D["pACE_B3_N4_13_rid_1.05_mlearn"] = "ACE(mlearn)"

# Add all others as just their name
D = {**dict(zip(models, models)), **D}

k = 0

# Renaming ticks dict
# D_tick_normal = dict(
#     zip(model_data["BAND_PATH"]["labels"], model_data["BAND_PATH"]["labels"])
# )
D_tick = {
    "G": "$\\Gamma$",
    ".": ".",
    "X": "X",
    "L": "L",
    "K": "K",
    "R": "R",
    "M": "M",
    "N": "N",
    "H": "H",
    "A": "A",
    "P": "P",
    "Z": "Z",
}

# models = ["DFT_QE"]  # ["GAP_inA_50s", "GAP_itB_1_a_2", "DFT_QE"]
for (j, model_name) in enumerate(models):
    if model_name == ref_model_name and len(models) > 1:
        continue
    print("analyze model", model_name)
    for (bulk_i, bulk_struct_test) in enumerate(bulk_struct_tests):
        if TEST_NAME not in data[model_name]:
            continue
        try:
            ref_model_data = data[ref_model_name][TEST_NAME][bulk_struct_test]
            THz_per_invcm = ase.units._c * 1.0e-10 * 100
            print(ref_model_data["DOS"]["val"])
            print(ref_model_data["DOS"]["freq"] / THz_per_invcm)  # THz_per_invcm)
            ax_DOS.plot(
                ref_model_data["DOS"]["val"],
                ref_model_data["DOS"]["freq"] / THz_per_invcm,
                color="black",
            )  # * THz_per_invcm)
        except:
            ref_model_data = None
        try:
            model_data = data[model_name][TEST_NAME][bulk_struct_test]
        except:
            continue
        print("analyze model-bulk", model_name, bulk_struct_test)
        analyze_phonons(model_data)
        THz_per_invcm = ase.units._c * 1.0e-10 * 100
        ax_DOS.plot(
            model_data["DOS"]["val"],
            ref_model_data["DOS"]["freq"] / THz_per_invcm,
            color=colors[j % len(colors)],
        )
        if ref_model_data is not None and "BAND_PATH" in ref_model_data:
            for i in range(ref_model_data["BAND_PATH"]["frequencies"].shape[1]):
                ax_BP.plot(
                    ref_model_data["BAND_PATH"]["positions"],
                    ref_model_data["BAND_PATH"]["frequencies"][:, i] * 0.01,
                    "-",
                    color="black",
                    label=D[ref_model_name] if k == 0 else None,
                )
                k += 1
        if "BAND_PATH" in model_data:
            for i in range(model_data["BAND_PATH"]["frequencies"].shape[1]):
                print(f"D keys: {D.keys()}")
                ax_BP.plot(
                    model_data["BAND_PATH"]["positions"],
                    model_data["BAND_PATH"]["frequencies"][:, i] * 0.01,
                    "-",
                    color=colors[j % len(colors)],
                    label=D[model_name] if i == 0 else None,
                )  ##############################
            ax_BP.set_xticks(
                [
                    p
                    for p, l in zip(
                        model_data["BAND_PATH"]["positions"],
                        model_data["BAND_PATH"]["labels"],
                    )
                    if l != "."
                ]
            )

            D_tick = {
                **dict(
                    zip(
                        model_data["BAND_PATH"]["labels"],
                        model_data["BAND_PATH"]["labels"],
                    )
                ),
                **D_tick,
            }
            ax_BP.set_xticklabels(
                [D_tick[l] for l in model_data["BAND_PATH"]["labels"] if l != "."]
            )
            print(model_data["BAND_PATH"]["labels"])
            ax_BP.set_xlim(
                model_data["BAND_PATH"]["positions"][0],
                model_data["BAND_PATH"]["positions"][-1],
            )
            # ax_BP.set_ylim(ylim)
            ylim = ax_DOS.get_ylim()
            if bulk_i == len(bulk_struct_tests) - 1:
                ylim = ax_BP.get_ylim()
                # print(ref_model_data['BAND_PATH']['labels'][1:-1])
                for p, l in zip(
                    ref_model_data["BAND_PATH"]["positions"][1:-1],
                    ref_model_data["BAND_PATH"]["labels"][1:-1],
                ):
                    if l != ".":
                        ax_BP.plot(
                            [p, p], ylim, "-", color="black", label=None, linewidth=0.5
                        )
# print(ref_model_data['DOS'].keys())
# THz_per_invcm = ase.units._c * 1.0e-10 * 100
# print(ref_model_data['DOS']['val'])
# print(ref_model_data['DOS']['freq'] / THz_per_invcm)# THz_per_invcm)
# ax_DOS.plot(ref_model_data['DOS']['val'], ref_model_data['DOS']['freq'] / THz_per_invcm, color="black")# * THz_per_invcm)
# ax_BP.set_xlabel("BZ point")
# ax_BP.set_ylim(ylim)
# ax_BP.set_title(bulk_struct_test)
ax_BP.legend(loc="upper right")
ax_BP.axhline(color="black", linewidth=0.5)
ax_BP.set_ylabel("Frequency (THz)")
ax_BP.set_xlabel("BZ point")
ax_BP.legend(loc="lower right")
ax_DOS.set_xlabel("DOS (THz⁻¹)")
# ax_DOS.set_xlabel("freq (cm$^{-1}$)")
# ax_DOS.set_ylabel("DOS (arb. units)")
# ylim = ax_DOS.get_ylim()
# ax_BP.set_ylim(-0.2, 4.0)
# ax_DOS.legend()
f.tight_layout()

for i in range(1000):
    fname = f"../plots/{TEST_NAME}_vDOS_{i}.pdf"
    if os.path.exists(fname):
        continue
    print(f"saved: {fname}")
    f.savefig(fname)
    break
# fig_BP.savefig("phonons_In2O3_Ia3_BAND_PATH-PAPER.pdf")
pyplot.clf()