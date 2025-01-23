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

from fast_forward.compute_bonded import compute_value_for_interaction
import numpy as np

def interaction_distribution(u, inter_type, pair_idxs, group_name, prefix):
    time_series = compute_value_for_interaction(u, inter_type, pair_idxs)
    np.savetxt("{prefix}{name}_{inter_type}.dat".format(name=group_name,
                                                        inter_type=inter_type,
                                                        prefix=prefix),
               time_series)

    bins_dict = {"bonds": np.arange(time_series.min()-0.1, time_series.max()+0.1, 0.01),
                 "angles": np.arange(181),
                 "dihedrals": np.arange(-180, 181)
                 }

    probs, edges = np.histogram(time_series, density=True, bins=bins_dict[inter_type])
    center_points = edges[:-1] + np.diff(edges)/2.
    distr = np.transpose((center_points, probs))
    np.savetxt("{prefix}{name}_{inter_type}_distr.dat".format(name=group_name,
                                                              inter_type=inter_type,
                                                              prefix=prefix),
               distr)

    return distr
