from phonons import *

properties = do_phonons(
    ["bulk_In2O3_Ia3"], n_supercell=2, band_paths=["GHNGPHPN"]
)  # GHNGPHPN from matrials project, before: GXKGL
