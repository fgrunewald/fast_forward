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
import numpy as np
from lmfit.models import GaussianModel
import matplotlib.pyplot as plt
from MDAnalysis.units import constants
from lmfit import create_params


def interaction_fitter(data, interaction, atom_list, T=310, plot=False):
    # this is not the boltzmann constant, but for some reason mda calls it this?
    R = constants['Boltzmann_constant']

    x = data.T[0]
    y = data.T[1]

    mod = GaussianModel()

    if interaction in ['angles', 'dihedrals']:
        pars = create_params(amplitude=dict(value=y.mean(), min=0),
                             center=dict(value=x.mean(), min=x.mean() - 20, max=x.mean() + 20),
                             sigma=dict(value=x.std(), min=x.std() / 2, max=x.std() * 1.5)
                             )
    else:
        pars = mod.guess(y, x=x)

    out = mod.fit(y, pars, x=x)

    center = np.round(out.params["center"].value, 2)

    # need this here because mdanalysis read gromacs coords in angstroms but need in nm.
    # can't convert earlier because otherwise the force constant goes stonks with small widths
    if interaction in ['bonds', 'constraints']:
        center = np.round(center / 10, 2)

    if (interaction == "dihedrals") or (interaction == 'angles'):
        sin_term = np.sin(np.deg2rad(np.round(out.params["center"].value, 2))) ** 2
        var = np.deg2rad(out.params["sigma"].value) ** 2
        sigma = np.round((R * T) / (sin_term * var), 2)
    else:
        sigma = np.round((R * T) / ((out.params["sigma"].value/10) ** 2), -1)

    if plot:
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

        fig.savefig(atom_list+f'_{interaction}.png')
        plt.close(fig)

    return center, sigma

