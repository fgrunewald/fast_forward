

'''
functions for fitting interaction distributions
'''

import numpy as np
from lmfit.models import GaussianModel
from MDAnalysis.units import constants
import lmfit
from collections import defaultdict
# from .interaction import Interaction
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

    def _bonds_fitter(self, data, atoms, group_name):
        """
        Fit bonds
        Parameters
        ----------
        data: np.array
            histogram of bond data
        atoms: list
            indices of the atoms involved in the interaction
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
        sigma = np.round((self.kb * self.temperature) / ((initial_sigma / 10) ** 2), -1)

        self.fit_parameters['bonds'][group_name] = [center, initial_sigma]

        if sigma < self.constraint_converter:
            self.interactions_dict['bonds'].append(Interaction(atoms=atoms[0],
                                                               parameters=[1, center, sigma],
                                                               meta={"comment": group_name}))
        else:
            self.interactions_dict['bonds'].append(Interaction(atoms=atoms[0],
                                                               parameters=[1, center, 10000],
                                                               meta={"ifndef": "FLEXIBLE", "comment": group_name}))
            self.interactions_dict['constraints'].append(Interaction(atoms=atoms[0],
                                                                     parameters=[1, center],
                                                                     meta={"ifdef": "FLEXIBLE", "comment": group_name},))

        #     return Interaction(name='bonds',
        #                        parameters=[1, center, sigma],
        #                        atoms=atoms[0],
        #                        meta={"comment": group_name},
        #                        fit_data=[float(center), initial_sigma])
        # else:
        #     return [Interaction(name='constraints',
        #                         parameters=[1, center],
        #                         atoms=atoms[0],
        #                         meta={"ifndef": "FLEXIBLE", "comment": group_name},
        #                         fit_data=[float(center), initial_sigma]),
        #             Interaction(name='bonds',
        #                         parameters=[1, center, 10000],
        #                         atoms=atoms[0],
        #                         meta={"ifdef": "FLEXIBLE", "comment": group_name},
        #                         fit_data=[float(center), initial_sigma])
        #             ]
        # return center, sigma


    def _angles_fitter(self, data, atoms, group_name,
                       ):
        """
        Fit angles
        Parameters
        ----------
        data: np.array
            histogram of bond data
        atoms: list
            indices of the atoms involved in the interaction
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

        # if sigma too big, angles get very unstable.
        if sigma > 150:
            sigma = 150

        #empiricaly derived. For theta_0 > 160, significant ptl energy for type 10 at equilibrium, so enforce type 1.
        if float(center) < 160:
            func_type_out = 10
        else:
            func_type_out = 1

        self.interactions_dict['angles'].append(Interaction(atoms=atoms[0],
                                                            parameters=[func_type_out, center, sigma],
                                                            meta={"comment": group_name}))
        self.fit_parameters['angles'][group_name] = [center, initial_sigma]
        # return Interaction(name='angles',
        #                    atoms=atoms[0],
        #                    parameters=[func_type_out, center, sigma],
        #                    meta={"comment": group_name}, fit_data=[float(center), initial_sigma])

    # Fitting function for proper dihedrals
    def _proper_dihedral_model_function(self, params, x):
        """Computes the sum of cosines with the given parameters."""
        y = np.zeros_like(x)
        num_terms = len(params) // 3  # Each term has k, n, x0
        for i in range(1,num_terms):
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

    def _dihedrals_fitter(self, data, atoms,  group_name):
        '''
        Fitter for dihedrals.
        Will try to fit both proper and improper dihedrals, deciding which to return based on
        the Akaike information criterion of the two fits

        Parameters
        ----------
        data: np.array
            histogram of bond data
        atoms: list
            indices of the atoms involved in the interaction
        group_name: str
            names of the atoms involved in the interaction joined by a "_"


        '''

        x = np.linspace(-np.pi, np.pi, 360)
        y = data.T[1]

        # first try fitting a gaussian to the data in case we have an improper dihedral
        gaussian_result = _gaussian_fitter(x, y,
                                                initial_center=dict(value=x[y.argmax()]),
                                                initial_sigma=dict(value=15, min=0, max=np.pi),
                                                initial_amplitude=dict(value=y.max())
                                                )

        # now do fitting for proper dihedrals
        # Iterate over different numbers of terms to find the optimal one
        best_aic = np.inf
        best_params = None
        aic_values = []

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
            aic_values.append(aic)

            # Keep track of the best model
            if aic < best_aic:
                best_aic = aic
                best_params = result.params

        num_terms = len(best_params) // 3  # Each term has k, n, and x0

        # compare the aic values to determine which type of dihedral we have
        if best_aic < gaussian_result.aic:
            pars_out = []
            for i in range(1,num_terms):
                k = -best_params[f'k{i}'].value * self.dihedral_scaling # make k negative to convert from distribution to potential
                x0 = best_params[f'x0_{i}'].value
                x0_deg = np.degrees(x0)  # Convert x0 from radians to degrees
                n = int(best_params[f'n{i}'].value)  # Ensure n is integer
                # pars_out.append(Interaction(name='dihedrals',
                #                             atoms=atoms[0], # contained in a list for some reason
                #                             parameters=[9, # function type
                #                                         np.round(x0_deg, self.precision), # center
                #                                         np.round(k, self.precision), # force constant
                #                                         int(n) # multiplicity
                #                                         ],
                #                             meta={"comment": group_name, "group": group_name},
                #                             fit_data=[best_params[f'k{i}'].value, x0_deg, n]
                #                             ))
                self.interactions_dict['dihedrals'].append(Interaction(atoms=atoms[0],
                                                                       parameters=[9,  # function type
                                                                                   np.round(x0_deg, self.precision), # center
                                                                                   np.round(k, self.precision), # force constant
                                                                                   int(n) # multiplicity
                                                                                   ],
                                                                       meta={"comment": group_name,
                                                                             "group": group_name}))
                pars_out.append([best_params[f'k{i}'].value, x0_deg, n])
            self.fit_parameters['dihedrals'][group_name] = pars_out

        else:
            center = np.round(np.degrees(gaussian_result.params['center']), self.precision)
            fc = np.round((self.kb * self.temperature) / ((gaussian_result.params['sigma']) ** 2), self.precision)
            self.interactions_dict['dihedrals'].append(Interaction(atoms=atoms[0],
                                                                   parameters=[2, center, fc],
                                                                   meta={"comment": group_name}))
            self.fit_parameters['dihedrals'][group_name] = [gaussian_result.params['amplitude'],
                                                            gaussian_result.params['center'],
                                                            gaussian_result.params['sigma']]
            # pars_out = Interaction(name='dihedrals',
            #                        atoms=atoms[0],
            #                        parameters=[2, center, fc],
            #                        meta={"comment": group_name},
            #                        fit_data=[gaussian_result.params['amplitude'],
            #                                  gaussian_result.params['center'],
            #                                  gaussian_result.params['sigma']]
            #                        )

        # return pars_out

    def fit_interaction(self, data, atoms, group_name, inter_type):
        func_dict = {'bonds': self._bonds_fitter,
                     'angles': self._angles_fitter,
                     'dihedrals': self._dihedrals_fitter
                     }
        func_dict[inter_type](data, atoms, group_name)



