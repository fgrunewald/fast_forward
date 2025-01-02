# Copyright 2024 University of Groningen
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from lmfit import create_params, minimize
import numpy as np
from fast_forward.interaction_fit import interaction_fitter


def func(pars, x, data):
    vals = pars.valuesdict()
    a = vals['a']
    b = vals['b']

    i, j, k = data

    r_ij = j - i
    r_jk = k - j
    r_ijk = r_ij + a * r_jk

    r_ijk_norm = r_ijk / np.linalg.norm(r_ijk)
    r_predicted_vs = i + (b * r_ijk_norm)

    return x - r_predicted_vs


def _vs3fd(pair_list):
    arr1, arr2, arr3, arr4 = pair_list
    frames = arr1.shape[0]
    ab = np.zeros((frames, 2))
    for fdx in range(frames):
        fit_params = create_params(a=1, b=1)
        data = [arr2[fdx], arr3[fdx], arr4[fdx]]
        out = minimize(func, fit_params, args=(arr1[fdx],), kws={'data': data})
        ab[fdx] = [out.params["a"].value, out.params["b"].value]
    return ab


def vs_handler(interaction_type, universe, indices, atom_list, prefix, plot=False):
    if interaction_type == 'virtual_sites3':
        inds = [indices]
        nframes = universe.trajectory.n_frames
        time_series = np.zeros((len(inds) * nframes, 2))
        for idx, idxs in enumerate(inds):
            pair_pos = [universe.trajectory.coordinate_array[:, pair, :] for pair in idxs]
            time_series[idx * nframes:(idx + 1) * nframes] = _vs3fd(pair_pos)

        params_out = []
        for param, var in zip(time_series.T, ['a', 'b']):
            np.savetxt("{prefix}{name}_{inter_type}_{variable}.dat".format(name=atom_list,
                                                                           inter_type=interaction_type,
                                                                           prefix=prefix,
                                                                           variable=var),
                       time_series)

            probs, edges = np.histogram(param, density=True, bins=30)
            center_points = edges[:-1] + np.diff(edges) / 2.
            distr = np.transpose((center_points, probs))
            np.savetxt("{prefix}{name}_{inter_type}_{variable}_distr.dat".format(name=atom_list,
                                                                                 inter_type=interaction_type,
                                                                                 prefix=prefix,
                                                                                 variable=var),
                       distr)

            center, _ = interaction_fitter(distr, "virtual_sites3", atom_list+f'_{var}', plot=plot)
            params_out.append(center)

        par_inds = params_out, indices

    elif interaction_type == 'virtual_sitesn':
        par_inds = [], indices

    return par_inds