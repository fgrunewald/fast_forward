import numpy as np
from fast_forward.bonded_functions import NORMAL_FUNCS

BINS_DICT = {"bonds": np.arange(0, 7, 0.01),
             "angles": np.arange(181),
             "dihedrals": np.arange(-180, 181),
             "distances": np.arange(0,30,0.05) # might need adjustment to allow for larger molecules
             }

VIRTUALSITE_TYPES = {'virtual_sitesn': {'1' : 'virtual_sitesn'},
                     'virtual_sites3': {'2' : 'virtual_sites3fd',
                                        '4' : 'virtual_sites3out'}}
ARR_SHAPES = {'virtual_sitesn': 1,
              'virtual_sites3fd': 2,
              'virtual_sites3out': 3
              }

RECOGNISED_INTERACTIONS = ['bonds', 'angles', 'dihedrals', 'virtual_sitesn', 'virtual_sites3', 'distances']

def interaction_distribution(u, inter_type, pair_idxs, group_name="", prefix="", save=False):
    # i.e. if inter_type not one of virtual_sitesn, virtual_sites3fd, etc.
    if inter_type not in [v for subdict in VIRTUALSITE_TYPES.values() for v in subdict.values()]:
        nframes = u.trajectory.n_frames
        time_series = np.zeros((len(pair_idxs) * nframes))
        for idx, idxs in enumerate(pair_idxs):
            pair_pos = [ u.trajectory.coordinate_array[:, pair, :] for pair in idxs]
            time_series[idx*nframes:(idx+1)*nframes] = NORMAL_FUNCS[inter_type](*pair_pos)

        if save:
            np.savetxt("{prefix}{name}_{inter_type}.dat".format(name=group_name,
                                                                inter_type=inter_type,
                                                                prefix=prefix),
                       time_series)
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
        vs_fitted = np.zeros((len(pair_idxs), int(ARR_SHAPES[inter_type])))
        for idx, idxs in enumerate(pair_idxs):
            pair_pos = [ u.trajectory.coordinate_array[:, pair, :] for pair in idxs]
            vs_fitted[idx] = NORMAL_FUNCS[inter_type](*pair_pos)
        if save:
            np.savetxt("{prefix}{name}_{inter_type}.dat".format(name=group_name,
                                                                inter_type=inter_type,
                                                                prefix=prefix),
                       vs_fitted)
        return vs_fitted

