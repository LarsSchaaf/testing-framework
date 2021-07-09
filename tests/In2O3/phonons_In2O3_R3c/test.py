from phonons import *

properties = do_phonons(
    ["bulk_In2O3_R3c"], n_supercell=1, band_paths=["GLBZGXQFPZLP"]
)  # GHNGPHPN from matrials project, before: GXKGL
