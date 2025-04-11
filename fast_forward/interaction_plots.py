


'''
Functions for plotting fitted interaction distributions
'''

import matplotlib.pyplot as plt
from lmfit.models import GaussianModel
from lmfit import Parameters
import numpy as np
import pickle

def _bond_plot(data, fit_params, atom_list, ax):
    '''
    plot a bonded distribution and the fit
    '''

    x = data.T[0]
    y = data.T[1]

    mod = GaussianModel(x=x)
    pars = Parameters()
    pars.add("center", fit_params[0]*10) #convert back to A
    pars.add("sigma", fit_params[1])
    fitted_distribution = mod.eval(pars,x=x)

    ax.plot(x, y, c='#6970E0', label='Distribution')
    ax.plot(x, fitted_distribution, c='#E06B69', label='Fit')

    curr_lims = ax.get_ylim()
    ax.set_ylim(0, curr_lims[1] + (curr_lims[1] * 0.1))

    ax.axvline(fit_params[0]*10, c='#506155', label=f"Center = {fit_params[0]*10: .2f}")
    ax.fill_between(np.linspace(fit_params[0]*10 - fit_params[1],
                                fit_params[0]*10 + fit_params[1],
                                100),
                    curr_lims[1] + (curr_lims[1] * 0.1),
                    color='#506155',
                    alpha=0.35,
                    label=f"Sigma = {fit_params[1]: .2f}")

    ax.legend()
    ax.set_title(f'{atom_list} bond')
    ax.set_xlabel('Distance')

    data_out = {ax.get_xlabel(): x,
                'simulated_distribution': y,
                'fitted_distribution': fitted_distribution,
                'distribution_center': fit_params[0],
                'distribution_sigma': fit_params[1]
                }
    return data_out

def _angles_plot(data, fit_params, atom_list, ax, ):
    '''
    plot an angle distribution and the fit
    '''

    x = data.T[0]
    y = data.T[1]

    mod = GaussianModel(x=x)
    pars = Parameters()
    pars.add("center", fit_params[0])
    pars.add("sigma", fit_params[1])
    fitted_distribution = mod.eval(pars, x=x)

    ax.plot(x, y, c='#6970E0', label='Distribution')
    ax.plot(x, fitted_distribution, c='#E06B69', label='Fit')

    curr_lims = ax.get_ylim()
    ax.set_ylim(0, curr_lims[1] + (curr_lims[1] * 0.1))
    ax.set_xlim(-180, 180)

    ax.axvline(fit_params[0], c='#506155', label=f"Center = {fit_params[0]: .2f}")
    ax.fill_between(np.linspace(fit_params[0] - fit_params[1],
                                fit_params[0] + fit_params[1],
                                100),
                    curr_lims[1] + (curr_lims[1] * 0.1),
                    color='#506155',
                    alpha=0.35,
                    label=f"Sigma = {fit_params[1]: .2f}")

    ax.legend()
    ax.set_title(f'{atom_list} angle')
    ax.set_xlabel('Angle')

    data_out = {ax.get_xlabel(): x,
                'simulated_distribution': y,
                'fitted_distribution': fitted_distribution,
                'distribution_center': fit_params[0],
                'distribution_sigma': fit_params[1]
                }
    return data_out

def _proper_dihedrals_plot(data, fit_params, atom_list, ax):
    '''
    plot a proper dihedral distribution and the fit
    '''
    x = np.linspace(-np.pi, np.pi, 360)
    y = data.T[1]
    y_fitted = np.zeros_like(x)
    for inter in fit_params:
        k, theta, n  = inter
        angle_rad = np.deg2rad(theta)
        y_fitted += k * (1 + np.cos(n*x - angle_rad))

    ax.plot(np.degrees(x), y, c='#6970E0', label='Distribution')
    ax.plot(np.degrees(x), y_fitted, c='#E06B69', label='Fit')

    ax.set_xlim(-180, 180)

    ax.legend()
    ax.set_title(f'{atom_list} dihedral')
    ax.set_xlabel('Angle')

    data_out = {ax.get_xlabel(): x,
                'simulated_distribution': y,
                'fitted_distribution': y_fitted,
                }
    return data_out


def _improper_dihedrals_plot(data, fit_params, atom_list, ax):
    '''
    plot an improper dihedral distribution and the fit
    '''
    x = np.linspace(-np.pi, np.pi, 360)
    y = data.T[1]

    mod = GaussianModel(x=x)
    pars = Parameters()

    pars.add("amplitude", fit_params[0])
    pars.add("center", fit_params[1])
    pars.add("sigma", fit_params[2])
    fitted_distribution = mod.eval(pars, x=x)

    ax.plot(np.degrees(x), y, c='#6970E0', label='Distribution')
    ax.plot(np.degrees(x), fitted_distribution, c='#E06B69', label='Fit')

    ax.set_xlim(-180, 180)

    ax.legend()
    ax.set_title(f'{atom_list} dihedral')
    ax.set_xlabel('Angle')


    data_out = {ax.get_xlabel(): x,
                'simulated_distribution': y,
                'fitted_distribution': fitted_distribution,
                }
    return data_out

def make_distribution_plot(fit_data, axarr=None, save_plot_data=False):
    '''

    Parameters
    ----------
    fit_data: dict
        Dictionary containing distributions and fitting parameters for interactions.
        Nested as {interaction_type: {group_name: {'data': distribution, 'fitted_params': list(params)}}
    axarr: matplotlib.pyplot.Figure.axes
        array of axes to plot the fitted distributions on
    save_plot_data: bool
        if True, save the underlying data for plots as a pickle file
    '''

    total_interactions = sum([len(fit_data[i]) for i in fit_data.keys()])

    if not axarr:
        ncols = 5
        nrows = -(total_interactions // -5) # upside-down floor division
        fig, axarr = plt.subplots(nrows=nrows, ncols=ncols,
                               figsize=(ncols*4, nrows*4))

    count = 0
    all_plot_data = {}
    for interaction_type in fit_data.keys():
        all_plot_data[interaction_type] = {}
        for atom_list in fit_data[interaction_type].keys():
            data = fit_data[interaction_type][atom_list]['data']
            fitted_parameters = fit_data[interaction_type][atom_list]['fit_params']
            ax = axarr.flatten()[count]
            if interaction_type == 'bonds':
                plot_data = _bond_plot(data, fitted_parameters, atom_list, ax)
            elif interaction_type == 'angles':
                plot_data = _angles_plot(data, fitted_parameters, atom_list, ax)
            elif interaction_type == 'dihedrals':
                if len(fitted_parameters) > 1:
                    plot_data = _proper_dihedrals_plot(data, fitted_parameters, atom_list, ax)
                else:
                    plot_data = _improper_dihedrals_plot(data, fitted_parameters, atom_list, ax)
            count += 1
            all_plot_data[interaction_type][atom_list] = plot_data

    if save_plot_data:
        pickle.dump(all_plot_data, open('plot_data.p', 'wb'))

    # remove unused axes
    for ax in axarr.flatten()[count:]:
        fig.delaxes(ax)

    # need to make room for the title
    fig.subplots_adjust(hspace = 0.3)

    fig.savefig(f'distribution_plots.png')
