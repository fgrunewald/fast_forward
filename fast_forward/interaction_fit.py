
'''
functions for fitting interaction distributions
'''

import numpy as np
from lmfit.models import GaussianModel
from MDAnalysis.units import constants
import lmfit
from collections import defaultdict
from vermouth.molecule import Interaction

def _gaussian_fitter(x, y, initial_center, initial_sigma, initial_amplitude):
    """
    Fit a Gaussian function to an input distribution

    Parameters
    ----------
    x: np.array
        x of data
    y: np.array
        y of data
    initial_center: dict
        dictionary of values for lmfit to use as starting parameters
    initial_sigma: dict
        dictionary of values for lmfit to use as starting parameters
    initial_amplitude: dict
        dictionary of values for lmfit to use as starting parameters

    Returns
    -------
    gaussian_result: lmfit.ModelResult
        lmfit fit result
    """

    gaussian = GaussianModel()
    params = gaussian.make_params(center=initial_center,
                                  sigma=initial_sigma,
                                  amplitude=initial_amplitude
                                  )
    gaussian_result = gaussian.fit(y, params, x=x)

    return gaussian_result


class InteractionFitter:
    """
    Class to fit interactions
    """
    def __init__(self, precision, temperature, constraint_converter,
                 max_dihedrals, dihedral_scaling):
        '''

        Parameters
        ----------
        precision: int
            precision to round values to for writing
        temperature: int
            temperature of interaction distribution in boltzmann inversion
        constraint_converter: int
            threshold above which to convert bonds to constraints
        max_dihedrals: int
            maximum number of dihedrals to fit proper dihedrals with
        '''
        self.precision = precision
        self.temperature = temperature
        self.kb = constants["Boltzmann_constant"]
        self.constraint_converter = constraint_converter
        self.max_dihedrals = max_dihedrals
        self.dihedral_scaling = dihedral_scaling
        # this will store the interactions
        self.interactions_dict = defaultdict(list)
        self.fit_parameters = defaultdict(dict)

    def _bonds_fitter(self, data, group_name):
        """
        Fit bonds
        Parameters
        ----------
        data: np.array
            histogram of bond data
        group_name: str
            names of the atoms involved in the interaction joined by a "_"

        """

        x, y = data.T
        gaussian_fit = _gaussian_fitter(x, y,
                                        initial_center=dict(value=x[y.argmax()]),
                                        initial_sigma=dict(value=x.std()),
                                        initial_amplitude=dict(value=y.max())
                                        )

        initial_center = gaussian_fit.params["center"].value
        initial_sigma = gaussian_fit.params["sigma"].value

        # need this here because mdanalysis read gromacs coords in angstroms but need in nm.
        center = np.round(initial_center / 10, self.precision)
        sigma = np.round((self.kb * self.temperature) / ((initial_sigma / 10) ** 2), self.precision)

        self.fit_parameters['bonds'][group_name] = {'data': data,
                                                    'distribution_params': [center, sigma],
                                                    'fit_params': [center, initial_sigma]}

    def _angles_fitter(self, data, group_name):
        """
        Fit angles
        Parameters
        ----------
        data: np.array
            histogram of bond data
        group_name: str
            names of the atoms involved in the interaction joined by a "_"

        """
        x, y = data.T
        gaussian_fit = _gaussian_fitter(x, y,
                                        initial_center=dict(value=x[y.argmax()], min=x[y.argmax()] - 20, max=x[y.argmax()] + 20),
                                        initial_sigma=dict(value=x.std(), min=x.std() / 4, max=x.std() * 1.5),
                                        initial_amplitude=dict(value=y.max(), min=0)
                                        )

        initial_center = gaussian_fit.params["center"].value
        initial_sigma = gaussian_fit.params["sigma"].value

        center = np.round(initial_center, self.precision)
        sin_term = np.sin(np.deg2rad(float(center))) ** 2
        var = np.deg2rad(initial_sigma) ** 2
        sigma = np.round((self.kb * self.temperature) / (sin_term * var), self.precision)

        self.fit_parameters['angles'][group_name] = {'data': data,
                                                     'distribution_params': [center, sigma],
                                                     'fit_params': [center, initial_sigma]}

    # Fitting function for proper dihedrals
    def _proper_dihedral_model_function(self, params, x):
        """Computes the sum of cosines with the given parameters."""
        y = np.zeros_like(x)
        num_terms = len(params) // 3  # Each term has k, n, x0
        for i in range(1, num_terms):
            k = params[f'k{i}']
            n = int(params[f'n{i}'].value)  # Force n to be an integer
            x0 = params[f'x0_{i}']
            y += k * (1 + np.cos(n * x - x0))
        return y

    # Residual function for lmfit
    def _residuals(self, params, x, data):
        '''
        Residual function for fitting proper dihedrals
        Parameters
        ----------
        params: lmfit.Parameters
            Parameters object for fit
        x: np.array
            x variable for dihedral data
        data: np.array
            probability data for dihedral distribution

        Returns
        -------
        model - data: residuals for fitting function to optimise

        '''
        return self._proper_dihedral_model_function(params, x) - data

    def _dihedrals_fitter(self, data, group_name):
        '''
        Fitter for dihedrals.
        Will try to fit both proper and improper dihedrals, deciding which to return based on
        the Akaike information criterion of the two fits

        Parameters
        ----------
        data: np.array
            histogram of bond data
        group_name: str
            names of the atoms involved in the interaction joined by a "_"


        '''

        x = np.linspace(-np.pi, np.pi, 360)
        y = data.T[1]

        # take care of periodic effects for improper dihedrals
        x_gauss = np.concatenate((x, x + 2 * np.pi)) - np.pi
        y_gauss = np.tile(y, 2)

        # first try fitting a gaussian to the data in case we have an improper dihedral
        gaussian_result = _gaussian_fitter(x_gauss[150:-150],
                                           y_gauss[150:-150],
                                           initial_center=dict(value=x[y.argmax()],
                                                               min=-np.pi,
                                                               max=np.pi),
                                           initial_sigma=dict(value=1,
                                                              min=0,
                                                              max=np.pi/3),
                                           initial_amplitude=dict(value=y.max())
                                           )

        # now do fitting for proper dihedrals
        # Iterate over different numbers of terms to find the optimal one
        best_aic = np.inf
        best_params = None

        for num_terms in range(1, self.max_dihedrals + 1):
            params = lmfit.Parameters()
            for i in range(num_terms):
                params.add(f'k{i}', value=1.0)  # Initial guess
                params.add(f'n{i}', value=i, vary=False)  # Fixed integer frequency
                params.add(f'x0_{i}', value=0.0, min=-np.pi, max=np.pi)  # Phase shift

            # Perform fitting
            minimizer = lmfit.Minimizer(self._residuals, params, fcn_args=(x, y))
            result = minimizer.minimize()

            # Compute AIC (lower is better)
            aic = result.aic

            # Keep track of the best model
            if aic < best_aic:
                best_aic = aic
                best_params = result.params

        num_terms = len(best_params) // 3  # Each term has k, n, and x0

        condition0 = best_aic < gaussian_result.aic
        condition1 = np.isclose(gaussian_result.params['sigma'].value, gaussian_result.params['sigma'].max)

        # compare the aic values to determine which type of dihedral we have
        # also make sure we don't have a very wide gaussian, where a single periodic function will suffice
        if condition0 or condition1:
            if not condition0 and condition1:
                num_terms = 2
            pars_out = []
            for i in range(1, num_terms):
                x0 = best_params[f'x0_{i}'].value
                n = int(best_params[f'n{i}'].value)  # Ensure n is integer
                pars_out.append([best_params[f'k{i}'].value, x0, n])
            self.fit_parameters['dihedrals'][group_name] = {'data': data,
                                                            'distribution_params': {i: j for i, j in enumerate(pars_out)},
                                                            'fit_params': pars_out}

        else:
            # transform the centre back into the correct domain after fitting to account for periodicity.
            c0 = (gaussian_result.params['center'] + (2*np.pi)) % (2*np.pi) - np.pi

            center = np.round(c0, self.precision)
            sigma = np.round((self.kb * self.temperature) / ((gaussian_result.params['sigma']) ** 2), self.precision)

            self.fit_parameters['dihedrals'][group_name] = {'data': data,
                                                            'distribution_params': [center, sigma],
                                                            'fit_params': {'amp': gaussian_result.params['amplitude'],
                                                                           'center': gaussian_result.params['center'],
                                                                           'sigma': gaussian_result.params['sigma']}}

    def fit_to_gmx(self, inter_type, group_name, atoms):

        if inter_type == 'bonds':
            parameters = self.fit_parameters['bonds'][group_name]['distribution_params']
            center, sigma = parameters

            if sigma < self.constraint_converter:
                self.interactions_dict['bonds'].append(Interaction(atoms=atoms[0],
                                                                   parameters=[1, center, sigma],
                                                                   meta={"comment": group_name}))
            else:
                self.interactions_dict['bonds'].append(Interaction(atoms=atoms[0],
                                                                   parameters=[1, center, 10000],
                                                                   meta={"ifdef": "FLEXIBLE", "comment": group_name}))
                self.interactions_dict['constraints'].append(Interaction(atoms=atoms[0],
                                                                         parameters=[1, center],
                                                                         meta={"ifndef": "FLEXIBLE",
                                                                               "comment": group_name}))
        elif inter_type == 'angles':
            parameters = self.fit_parameters['angles'][group_name]['distribution_params']
            center, sigma = parameters

            # empirically derived. if sigma too big, angles get very unstable.
            if sigma > 150:
                sigma = 150

            # empirically derived. For theta_0 > 160, significant ptl energy for type 10 at equilibrium, so enforce type 1.
            if float(center) < 160:
                func_type_out = 10
            else:
                func_type_out = 1

            self.interactions_dict['angles'].append(Interaction(atoms=atoms[0],
                                                                parameters=[func_type_out, center, sigma],
                                                                meta={"comment": group_name}))

        elif inter_type == 'dihedrals':

            parameters = self.fit_parameters['dihedrals'][group_name]['distribution_params']
            if isinstance(parameters, list):
                center, sigma = parameters
                center = np.round(np.degrees(center), self.precision)

                self.interactions_dict['dihedrals'].append(Interaction(atoms=atoms[0],
                                                                       parameters=[2, center, sigma],
                                                                       meta={"comment": group_name}))
            else:
                for i in parameters.values():
                    # factors derived from the fitting directly have negligible effects (~10^-3/4),
                    # scaling them helps increase the strength of dihedral in the final interaction
                    k = - i[0] * self.dihedral_scaling
                    x0_deg = np.degrees(i[1])
                    n = i[2]
                    self.interactions_dict['dihedrals'].append(Interaction(atoms=atoms[0],
                                                                           parameters=[9,  # function type
                                                                                       np.round(x0_deg, self.precision), # center
                                                                                       np.round(k, self.precision), # force constant
                                                                                       int(n) # multiplicity
                                                                                       ],
                                                                           meta={"comment": group_name,
                                                                                 "group": group_name}))




    def fit_interaction(self, data, atoms, group_name, inter_type):
        func_dict = {'bonds': self._bonds_fitter,
                     'angles': self._angles_fitter,
                     'dihedrals': self._dihedrals_fitter
                     }
        func_dict[inter_type](data, group_name)
        self.fit_to_gmx(inter_type, group_name, atoms)
