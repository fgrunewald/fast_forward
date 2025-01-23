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


def make_distribution_plot(x, y, out, atom_list, interaction):
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

    fig.savefig(atom_list + f'_{interaction}.png')
    plt.close(fig)

def _bonds_fitter(initial_center, initial_sigma, precision, R, T):
    # need this here because mdanalysis read gromacs coords in angstroms but need in nm.
    center = f'{np.round(initial_center / 10, precision):.{precision}f}'
    sigma = np.round((R * T) / ((initial_sigma / 10) ** 2), -1)
    return center, sigma

def _angles_fitter(initial_center, initial_sigma, precision, R, T):
    center = f'{np.round(initial_center, precision):.{precision}f}'
    sin_term = np.sin(np.deg2rad(float(center))) ** 2
    var = np.deg2rad(initial_sigma) ** 2
    sigma = np.round((R * T) / (sin_term * var), 2)
    return center, sigma

def _dihedrals_fitter(initial_center, initial_sigma, precision, R, T):

    initial_center -= 180
    center = ((initial_center + 180) % 360) - 180
    center = f'{np.round(center, precision):.{precision}f}'

    init = np.deg2rad(initial_sigma)
    sin = np.sin(init)
    cos = np.cos(init)
    circ = np.atan2(sin, cos)
    sigma = np.round((R * T) / (circ ** 2), 1)
    return center, sigma

def interaction_fitter(data, interaction, precision=3, T=310):

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

    fit_result = mod.fit(y, pars, x=x)

    center = fit_result.params["center"].value
    sigma = fit_result.params["sigma"].value

    func_dict = {'bonds': _bonds_fitter,
                 'constraints': _bonds_fitter,
                 'angles': _angles_fitter,
                 'dihedrals': _dihedrals_fitter,
                 }

    center, sigma = func_dict[interaction](center, sigma, precision, R, T)

    return center, sigma, fit_result

