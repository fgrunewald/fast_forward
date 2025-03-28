

import numpy as np
from lmfit.models import GaussianModel
import matplotlib.pyplot as plt
from MDAnalysis.units import constants
import lmfit
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

def _bonds_fitter(initial_center, initial_sigma, atoms, precision, R, T, group_name, convert_constraints=10000,
                  ):
    # need this here because mdanalysis read gromacs coords in angstroms but need in nm.
    center = f'{np.round(initial_center / 10, precision):.{precision}f}'
    sigma = np.round((R * T) / ((initial_sigma / 10) ** 2), -1)
    # return center, sigma

    if sigma < convert_constraints:
        return Interaction(name='bonds', func_type=1,
                           location=center, force_constant=sigma, atoms=atoms[0],
                           meta={"comment": group_name})
    else:
        return Interaction(name='constraints', func_type=1,
                           location=center, atoms=atoms[0],
                           meta={"ifndef": "FLEXIBLE", "comment": group_name})


def _angles_fitter(initial_center, initial_sigma, atoms, precision, R, T, group_name, convert_constraints=None,
                   ):
    '''
    putting constraint converter here for ease of use below, but obviously has no function
    '''
    center = f'{np.round(initial_center, precision):.{precision}f}'
    sin_term = np.sin(np.deg2rad(float(center))) ** 2
    var = np.deg2rad(initial_sigma) ** 2
    sigma = np.round((R * T) / (sin_term * var), 2)
    # return center, sigma
    return Interaction(name='angles', func_type=2,
                       location=center, force_constant=sigma, atoms=atoms[0],
                       meta={"comment": group_name})

# Fitting function for proper dihedrals
def model_function(params, x):
    """Computes the sum of cosines with the given parameters."""
    y = np.zeros_like(x)
    num_terms = len(params) // 3  # Each term has k, n, x0
    for i in range(num_terms):
        k = params[f'k{i}']
        n = int(params[f'n{i}'].value)  # Force n to be an integer
        x0 = params[f'x0_{i}']
        y += k * (1 + np.cos(n * x - x0))
    return y

# Residual function for lmfit
def residuals(params, x, data):
    return model_function(params, x) - data

def _dihedrals_fitter(data, atoms,  group_name, max_terms = 10):

    x = data.T[0]
    y = data.T[1]

    # Iterate over different numbers of terms to find the optimal one
    best_aic = np.inf
    # best_model = None
    best_params = None
    aic_values = []

    for num_terms in range(1, max_terms + 1):
        params = lmfit.Parameters()
        for i in range(num_terms):
            params.add(f'k{i}', value=1.0)  # Initial guess
            params.add(f'n{i}', value=i, vary=False)  # Fixed integer frequency
            params.add(f'x0_{i}', value=0.0, min=-np.pi, max=np.pi)  # Phase shift

        # Perform fitting
        minimizer = lmfit.Minimizer(residuals, params, fcn_args=(x, y))
        result = minimizer.minimize()

        # Compute AIC (lower is better)
        aic = result.aic
        aic_values.append(aic)

        # Keep track of the best model
        if aic < best_aic:
            best_aic = aic
            # best_model = result
            best_params = result.params

    num_terms = len(best_params) // 3  # Each term has k, n, and x0

    factor = 1e3 # useful to scale the potential slightly
    pars_out = []
    for i in range(num_terms):
        k = -best_params[f'k{i}'].value * factor# make k negative to convert from distribution to potential
        x0 = best_params[f'x0_{i}'].value
        x0_deg = np.degrees(x0)  # Convert x0 from radians to degrees
        n = int(best_params[f'n{i}'].value)  # Ensure n is integer

        pars_out.append(Interaction(name='dihedrals',
                                    atoms=atoms[0], # contained in a list for some reason
                                    func_type=9,
                                    location=np.round(x0_deg,2),
                                    force_constant=np.round(k, 3),
                                    multiplicity=int(n),
                                    meta={"comment": group_name}
                                    ))
    return pars_out

def interaction_fitter(data, interaction, atoms, group_name,
                       precision=3, convert_constraints=10000, T=310):

    R = constants['Boltzmann_constant']

    x = data.T[0]
    y = data.T[1]

    if interaction != 'dihedrals':

        mod = GaussianModel()

        if interaction in ['angles']:
            pars = lmfit.create_params(amplitude=dict(value=y.mean(), min=0),
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

        inter = func_dict[interaction](center, sigma, atoms, precision, R, T, group_name, convert_constraints)

    else:
        inter = _dihedrals_fitter(data, atoms, group_name)

    return inter

