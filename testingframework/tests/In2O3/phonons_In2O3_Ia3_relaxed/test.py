from phonons import *

properties = do_phonons(
    ["bulk_In2O3_Ia3_relaxed"],
    n_supercell=1,
    band_paths=["GHNGPHPN"],
    dx=0.010583544211276823,
)  # GHNGPHPN from matrials project, before: GXKGL
