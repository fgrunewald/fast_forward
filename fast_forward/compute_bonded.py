import argparse
import numpy as np
from numba import jit, njit, prange
import MDAnalysis as mda
from fast_forward.bonded_functions import NORMAL_FUNCS

def compute_value_for_interaction(universe, inter_type, valid_pairs):
    for idxs in valid_pairs:
        pair_pos = [ universe.trajectory.coordinate_array[:, pair, :] for pair in idxs]
    time_series = NORMAL_FUNCS[inter_type](*pair_pos)
    return time_series
