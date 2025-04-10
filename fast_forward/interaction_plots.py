


'''
Functions for plotting fitted interaction distributions
'''

import matplotlib.pyplot as plt
from lmfit.models import GaussianModel
from lmfit import Parameters
import numpy as np
import pickle

def _bond_plot(data, interaction, atom_list, plot_data=False):
    #if we have a constraint
    if type(interaction) == list:
        interaction = interaction[0]
    x = data.T[0]
    y = data.T[1]

    mod = GaussianModel(x=x)
    pars = Parameters()
    pars.add("center", interaction.fit_data[0]*10) #convert back to A
    pars.add("sigma", interaction.fit_data[1])
    fitted_distribution = mod.eval(pars,x=x)

    fig, ax = plt.subplots()

    ax.plot(x, y, c='#6970E0', label='Distribution')
    ax.plot(x, fitted_distribution, c='#E06B69', label='Fit')

    curr_lims = ax.get_ylim()
    ax.set_ylim(0, curr_lims[1] + (curr_lims[1] * 0.1))

    ax.axvline(interaction.fit_data[0]*10, c='#506155', label=f"Center = {interaction.fit_data[0]*10: .2f}")
    ax.fill_between(np.linspace(interaction.fit_data[0]*10 - interaction.fit_data[1],
                                interaction.fit_data[0]*10 + interaction.fit_data[1],
                                100),
                    curr_lims[1] + (curr_lims[1] * 0.1),
                    color='#506155',
                    alpha=0.35,
                    label=f"Sigma = {interaction.fit_data[1]: .2f}")

    ax.legend()
    ax.set_title(atom_list)
    ax.set_xlabel('Distance')

    if plot_data:
        data_out = {ax.get_xlabel(): x,
                    'simulated_distribution': y,
                    'fitted_distribution': fitted_distribution,
                    'distribution_center': interaction.fit_data[0],
                    'distribution_sigma': interaction.fit_data[1]
                    }
        pickle.dump(data_out, open(f'bond_{atom_list}.p', 'wb'))

    fig.savefig(f'bond_{atom_list}.png')
    plt.close(fig)

def _angles_plot(data, interaction, atom_list, plot_data=False):

    x = data.T[0]
    y = data.T[1]

    mod = GaussianModel(x=x)
    pars = Parameters()
    pars.add("center", interaction.fit_data[0])
    pars.add("sigma", interaction.fit_data[1])
    fitted_distribution = mod.eval(pars, x=x)

    fig, ax = plt.subplots()

    ax.plot(x, y, c='#6970E0', label='Distribution')
    ax.plot(x, fitted_distribution, c='#E06B69', label='Fit')

    curr_lims = ax.get_ylim()
    ax.set_ylim(0, curr_lims[1] + (curr_lims[1] * 0.1))
    ax.set_xlim(-180, 180)

    ax.axvline(interaction.fit_data[0], c='#506155', label=f"Center = {interaction.fit_data[0]: .2f}")
    ax.fill_between(np.linspace(interaction.fit_data[0] - interaction.fit_data[1],
                                interaction.fit_data[0] + interaction.fit_data[1],
                                100),
                    curr_lims[1] + (curr_lims[1] * 0.1),
                    color='#506155',
                    alpha=0.35,
                    label=f"Sigma = {interaction.fit_data[1]: .2f}")

    ax.legend()
    ax.set_title(atom_list)
    ax.set_xlabel('Angle')

    if plot_data:
        data_out = {ax.get_xlabel(): x,
                    'simulated_distribution': y,
                    'fitted_distribution': fitted_distribution,
                    'distribution_center': interaction.fit_data[0],
                    'distribution_sigma': interaction.fit_data[1]
                    }
        pickle.dump(data_out, open(f'angle_{atom_list}.p', 'wb'))

    fig.savefig(f'angle_{atom_list}.png')
    plt.close(fig)

def _proper_dihedrals_plot(data, interaction, atom_list, plot_data=False):
    x = np.linspace(-np.pi, np.pi, 360)
    y = data.T[1]
    y_fitted = np.zeros_like(x)
    for inter in interaction:
        k, theta, n  = inter.fit_data
        angle_rad = np.deg2rad(theta)
        # k = -k/1e3 # the scaling factor is hard coded, and we need to convert back from potential to fitted angle

        y_fitted += k * (1 + np.cos(n*x - angle_rad))

    fig, ax = plt.subplots()

    ax.plot(np.degrees(x), y, c='#6970E0', label='Distribution')
    ax.plot(np.degrees(x), y_fitted, c='#E06B69', label='Fit')

    curr_lims = ax.get_ylim()
    # ax.set_ylim(0, curr_lims[1] + (curr_lims[1] * 0.1))
    ax.set_xlim(-180, 180)

    ax.legend()
    ax.set_title(atom_list)
    ax.set_xlabel('Angle')

    if plot_data:
        data_out = {ax.get_xlabel(): x,
                    'simulated_distribution': y,
                    'fitted_distribution': y_fitted,
                    }
        pickle.dump(data_out, open(f'dihedral_{atom_list}.p', 'wb'))

    fig.savefig(f'dihedral_{atom_list}.png')
    plt.close(fig)

def _improper_dihedrals_plot(data, interaction, atom_list, plot_data=False):
    x = np.linspace(-np.pi, np.pi, 360)
    y = data.T[1]

    mod = GaussianModel(x=x)
    pars = Parameters()

    pars.add("amplitude", interaction.fit_data[0])
    pars.add("center", interaction.fit_data[1])
    pars.add("sigma", interaction.fit_data[2])
    fitted_distribution = mod.eval(pars, x=x)

    fig, ax = plt.subplots()

    ax.plot(np.degrees(x), y, c='#6970E0', label='Distribution')
    ax.plot(np.degrees(x), fitted_distribution, c='#E06B69', label='Fit')

    curr_lims = ax.get_ylim()
    # ax.set_ylim(0, curr_lims[1] + (curr_lims[1] * 0.1))
    ax.set_xlim(-180, 180)

    ax.legend()
    ax.set_title(atom_list)
    ax.set_xlabel('Angle')

    if plot_data:
        data_out = {ax.get_xlabel(): x,
                    'simulated_distribution': y,
                    'fitted_distribution': fitted_distribution,
                    }
        pickle.dump(data_out, open(f'dihedral_{atom_list}.p', 'wb'))

    fig.savefig(f'dihedral_{atom_list}.png')
    plt.close(fig)


def make_distribution_plot(data, interaction, atom_list, interaction_type, plot_data=False):

    if interaction_type in ['bonds', 'constraints']:
        _bond_plot(data, interaction, atom_list)
    elif interaction_type == 'angles':
        _angles_plot(data, interaction, atom_list)
    elif interaction_type == 'dihedrals':
        if type(interaction) == list:
            _proper_dihedrals_plot(data, interaction, atom_list)
        else:
            _improper_dihedrals_plot(data, interaction, atom_list)

