

'''
functions for fitting interaction distributions
'''

import numpy as np
from lmfit.models import GaussianModel
from MDAnalysis.units import constants
import lmfit
from .interaction import Interaction


def _bonds_fitter(initial_center, initial_sigma, atoms, precision, R, T, group_name, convert_constraints=10000,
                  ):
    # need this here because mdanalysis read gromacs coords in angstroms but need in nm.
    center = f'{np.round(initial_center / 10, precision):.{precision}f}'
    sigma = np.round((R * T) / ((initial_sigma / 10) ** 2), -1)

    if sigma < convert_constraints:
        return Interaction(name='bonds', func_type=1,
                           location=center, force_constant=sigma, atoms=atoms[0],
                           meta={"comment": group_name}, fit_data=[float(center), initial_sigma])
    else:
        return [Interaction(name='constraints', func_type=1,
                           location=center, atoms=atoms[0],
                           meta={"ifndef": "FLEXIBLE", "comment": group_name},
                           fit_data=[float(center), initial_sigma]),
                Interaction(name='bonds', func_type=1,
                            location=center, atoms=atoms[0], force_constant=10000,
                            meta={"ifdef": "FLEXIBLE", "comment": group_name},
                            fit_data=[float(center), initial_sigma])
            ]


def _angles_fitter(initial_center, initial_sigma, atoms, precision, R, T, group_name, convert_constraints=None,
                   ):
    '''
    putting constraint converter here for ease of use below, but obviously has no function
    '''
    center = f'{np.round(initial_center, precision):.{precision}f}'
    sin_term = np.sin(np.deg2rad(float(center))) ** 2
    var = np.deg2rad(initial_sigma) ** 2
    sigma = np.round((R * T) / (sin_term * var), 2)

    return Interaction(name='angles', func_type=2,
                       location=center, force_constant=sigma, atoms=atoms[0],
                       meta={"comment": group_name}, fit_data=[float(center), initial_sigma])

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

def _dihedrals_fitter(data, atoms,  group_name, R, T, max_terms = 10):

    x = np.linspace(-np.pi, np.pi, 360)
    y = data.T[1]

    # make a first try at fitting a gaussian to the data in case we have an improper dihedral
    gaussian = GaussianModel()
    gaussian.guess(y, x=x)
    gaussian_result = gaussian.fit(y, x=x)

    # Iterate over different numbers of terms to find the optimal one
    best_aic = np.inf
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
            best_params = result.params

    num_terms = len(best_params) // 3  # Each term has k, n, and x0

    # compare the aic values to determine which type of dihedral we have
    if best_aic < gaussian_result.aic:
        factor = 1#e2 # useful to scale the potential slightly
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
                                        force_constant=np.round(k, 6),
                                        multiplicity=int(n),
                                        meta={"comment": group_name},
                                        fit_data=[-k/factor, x0_deg, n]
                                        ))
    else:
        pars_out = Interaction(name='dihedrals',
                               atoms=atoms[0],
                               func_type=2,
                               location=np.round(np.degrees(gaussian_result.params['center']), 0),
                               force_constant=np.round((R * T) / ((gaussian_result.params['sigma']) ** 2), 0),
                               meta={"comment": group_name},
                               fit_data=[gaussian_result.params['amplitude'],
                                         gaussian_result.params['center'],
                                         gaussian_result.params['sigma']]
                               )

    return pars_out

def interaction_fitter(data, interaction, atoms, group_name,
                       precision=3, convert_constraints=10000, T=310):
    '''
    Entry function for interaction fitting
    Parameters
    ----------
    data: np.array
        nx2 array of a histogram of an interaction distribution
    interaction: str
        name of interaction
    atoms: np.array (I think?)
        indices of atoms involved in interaction
    group_name: str
        name of atoms involved in interaction
    precision: int
        number of deci
    convert_constraints: int
        force constant above which to convert bond to constraint
    T: int
        temperature of boltzmann inversion

    Returns
    -------

    '''
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
        inter = _dihedrals_fitter(data, atoms, group_name, R, T)

    return inter

