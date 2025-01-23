# Copyright 2025 University of Groningen
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

import matplotlib.pyplot as plt
import numpy as np
import pickle

def make_distribution_plot(data, out, atom_list, interaction, plot_data=False):

    x = data.T[0]
    y = data.T[1]

    fig, ax = plt.subplots()

    ax.plot(x, y, c='#6970E0', label='Distribution')
    ax.plot(x, out.best_fit, c='#E06B69', label='Fit')

    curr_lims = ax.get_ylim()
    ax.set_ylim(0, curr_lims[1] + (curr_lims[1] * 0.1))

    ax.axvline(out.params['center'].value, c='#506155', label=f"Center = {out.params['center'].value: .2f}")
    ax.fill_between(np.linspace(out.params['center'].value - out.params['sigma'].value,
                                out.params['center'].value + out.params['sigma'].value,
                                100),
                    curr_lims[1] + (curr_lims[1] * 0.1),
                    color='#506155',
                    alpha=0.35,
                    label=f"Sigma = {out.params['sigma'].value: .2f}")

    ax.legend()
    ax.set_title(atom_list)

    if (interaction == "dihedrals") or (interaction == 'angles'):
        ax.set_xlabel('Angle')
        ax.set_xlim(-180, 180)
    else:
        ax.set_xlabel('Distance')

    if plot_data:
        data_out = {ax.get_xlabel(): x,
                    'simulated_distribution': y,
                    'fitted_distribution': out.best_fit,
                    'distribution_center': out.params['center'].value,
                    'distribution_sigma': out.params['sigma'].value,
                    'full_fit_data': out
                    }
        pickle.dump(data_out, open(atom_list + f'_{interaction}.p', 'wb'))

    fig.savefig(atom_list + f'_{interaction}.png')
    plt.close(fig)
