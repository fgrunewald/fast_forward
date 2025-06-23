
from fast_forward.compute_bonded import compute_value_for_interaction
import numpy as np

BINS_DICT = {"bonds": np.arange(0, 7, 0.01),
             "angles": np.arange(181),
             "dihedrals": np.arange(-180, 181)
             }

VIRTUALSITE_TYPES = ['virtual_sitesn', 'virtual_sites3']

def interaction_distribution(u, inter_type, pair_idxs, group_name, prefix, save):
    time_series = compute_value_for_interaction(u, inter_type, pair_idxs)
    if save:
        np.savetxt("{prefix}{name}_{inter_type}.dat".format(name=group_name,
                                                            inter_type=inter_type,
                                                            prefix=prefix),
                   time_series)
    if inter_type not in VIRTUALSITE_TYPES:
        probs, edges = np.histogram(time_series, density=True, bins=BINS_DICT[inter_type])
        center_points = edges[:-1] + np.diff(edges)/2.
        distr = np.transpose((center_points, probs))
        if save:
            np.savetxt("{prefix}{name}_{inter_type}_distr.dat".format(name=group_name,
                                                                      inter_type=inter_type,
                                                                      prefix=prefix),
                       distr)

        return distr
    else:
        return time_series