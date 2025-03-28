

import numpy as np
from lmfit.models import GaussianModel
import matplotlib.pyplot as plt
from MDAnalysis.units import constants
from lmfit import create_params
from scipy import fft
from .interaction import Interaction


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

def _bonds_fitter(initial_center, initial_sigma, atoms, precision, R, T, convert_constraints=10000):
    # need this here because mdanalysis read gromacs coords in angstroms but need in nm.
    center = f'{np.round(initial_center / 10, precision):.{precision}f}'
    sigma = np.round((R * T) / ((initial_sigma / 10) ** 2), -1)
    # return center, sigma

    if sigma < convert_constraints:
        return Interaction(name='bonds', func_type=1,
                           location=center, force_constant=sigma, atoms=atoms[0])
    else:
        return Interaction(name='constraints', func_type=1,
                           location=center, atoms=atoms[0],
                           meta={"ifndef": "FLEXIBLE"})


def _angles_fitter(initial_center, initial_sigma, atoms, precision, R, T, convert_constraints=None):
    '''
    putting constraint converter here for ease of use below, but obviously has no function
    '''
    center = f'{np.round(initial_center, precision):.{precision}f}'
    sin_term = np.sin(np.deg2rad(float(center))) ** 2
    var = np.deg2rad(initial_sigma) ** 2
    sigma = np.round((R * T) / (sin_term * var), 2)
    # return center, sigma
    return Interaction(name='angles', func_type=2,
                       location=center, force_constant=sigma, atoms=atoms[0])

def _dihedrals_fitter(data, atoms,  n_lim = 10):

    signal = data.T[1]
    N = len(signal)
    # x = np.arange(N)

    inters = []
    n = n_lim
    res = fft.dct(signal, norm='ortho')
    res[n:] = 0

    # reconstructed = np.zeros_like(x, dtype=float)
    for m in range(n):
        C_m = np.sqrt(1 / N) if m == 0 else np.sqrt(2 / N)
        k_m = res[m] * C_m
        n_m = m / 2  # * (np.pi/(2*np.pi/N))
        # reconstructed += k_m * (1 + np.cos(n_m * np.deg2rad(x)))
        inters.append(Interaction(name='dihedrals',
                                  func_type=9,
                                  location=0,
                                  force_constant=-k_m,  # has to be -k_m to give _potential_ not _distribution_
                                  multiplicity=n_m, atoms=atoms[0]))
        print(inters[-1].parameters)

    return inters

def interaction_fitter(data, interaction, atoms, precision=3, convert_constraints=10000, T=310):

    R = constants['Boltzmann_constant']

    x = data.T[0]
    y = data.T[1]

    if interaction != 'dihedrals':

        mod = GaussianModel()

        if interaction in ['angles']:
            pars = create_params(amplitude=dict(value=y.mean(), min=0),
                                 center=dict(value=x[y.argmax()], min=x[y.argmax()] - 20, max=x[y.argmax()] + 20),
                                 sigma=dict(value=x.std(), min=x.std() / 4, max=x.std() * 1.5)
                                 )
        else:
            pars = mod.guess(y, x=x)

        fit_result = mod.fit(y, pars, x=x)

        center = fit_result.params["center"].value
        sigma = fit_result.params["sigma"].value

        func_dict = {'bonds': _bonds_fitter,
                     'constraints': _bonds_fitter,
                     'angles': _angles_fitter,
                     }

        inter = func_dict[interaction](center, sigma, atoms, precision, R, T, convert_constraints)

    else:
        inter = _dihedrals_fitter(data, atoms)

    return inter

