from phonons import *

properties = do_phonons(
    ["bulk_In2O3_R3c"],
    n_supercell=2,
    band_paths=[["G", "L", "B1", "B", "Z", "G", "X", "Q", "F", "P1", "Z", "L", "P"]],
    dx=0.010583544211276823,
)  # GHNGPHPN from matrials project, before: GXKGL
