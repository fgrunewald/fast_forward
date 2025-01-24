
from fast_forward.compute_bonded import compute_value_for_interaction
import numpy as np

BINS_DICT = {"bonds": np.arange(0, 7, 0.01),
             "angles": np.arange(181),
             "dihedrals": np.arange(-180, 181)
             }

def interaction_distribution(u, inter_type, pair_idxs, group_name, prefix):
    time_series = compute_value_for_interaction(u, inter_type, pair_idxs)
    np.savetxt("{prefix}{name}_{inter_type}.dat".format(name=group_name,
                                                        inter_type=inter_type,
                                                        prefix=prefix),
               time_series)

    probs, edges = np.histogram(time_series, density=True, bins=BINS_DICT[inter_type])
    center_points = edges[:-1] + np.diff(edges)/2.
    distr = np.transpose((center_points, probs))
    np.savetxt("{prefix}{name}_{inter_type}_distr.dat".format(name=group_name,
                                                              inter_type=inter_type,
                                                              prefix=prefix),
               distr)

    return distr
