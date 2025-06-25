import numpy as np
from fast_forward.bonded_functions import NORMAL_FUNCS

ARR_SHAPES = {'bonds': 1,
              'angles': 1,
              'dihedrals': 1,
              'virtual_sites3fd': 2,
              'virtual_sites3out': 3
              }


def compute_value_for_interaction(universe, inter_type, valid_pairs):
    nframes = universe.trajectory.n_frames
    # this is designed to make an empty array with the correct shape.
    # virtual_sitesn is missing, because we can't assume how many atoms are needed.
    # in this case, take the len of the first valid pairs entry
    time_series = np.zeros((len(valid_pairs) * nframes, ARR_SHAPES.get(inter_type, len(valid_pairs[0]))))
    for idx, idxs in enumerate(valid_pairs):
        pair_pos = [ universe.trajectory.coordinate_array[:, pair, :] for pair in idxs]
        time_series[idx*nframes:(idx+1)*nframes] = NORMAL_FUNCS[inter_type](*pair_pos)
    return time_series
