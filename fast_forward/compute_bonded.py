import argparse
import numpy as np
from numba import jit, njit, prange
import MDAnalysis as mda
from fast_forward.bonded_functions import NORMAL_FUNCS

def compute_value_for_interaction(universe, inter_type, valid_pairs):
    nframes = universe.trajectory.n_frames
    time_series = np.zeros((len(valid_pairs) * nframes))
    for idx, idxs in enumerate(valid_pairs):
        pair_pos = [ universe.trajectory.coordinate_array[:, pair, :] for pair in idxs]
        time_series[idx*nframes:(idx+1)*nframes] = NORMAL_FUNCS[inter_type](*pair_pos)
    return time_series
