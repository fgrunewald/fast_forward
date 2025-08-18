
'''
Functions for plotting fitted interaction distributions
'''

import matplotlib.pyplot as plt
import pickle
from .interaction_distribution import BINS_DICT

X_LABELS={'bonds': 'Distance',
          'angles': 'Angle',
          'dihedrals': 'Angle'}

def _plotter(data, atom_list, inter_type, ax):

    cols = ['#6970E0', '#E06B69']
    needed_keys = [key for key in list(data.keys()) if key != 'x']
    for idx, key in enumerate(needed_keys):
        ax.plot(data['x'],
                data[key],
                c = cols[idx],
                label=key)
    ax.legend()
    ax.set_title(f'{atom_list} {inter_type}')
    ax.set_xlim(BINS_DICT[inter_type][0], BINS_DICT[inter_type][-1])
    ax.set_xlabel(X_LABELS[inter_type])

def make_distribution_plot(fit_data, save_plot_data=False, axarr=None):
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
    for interaction_type in fit_data.keys():
        for atom_list in fit_data[interaction_type].keys():
            _plotter(fit_data[interaction_type][atom_list], atom_list, interaction_type, axarr.flatten()[count])
            count += 1

    if save_plot_data:
        pickle.dump(fit_data, open('plot_data.p', 'wb'))

    # remove unused axes
    for ax in axarr.flatten()[count:]:
        fig.delaxes(ax)

    # need to make room for the title
    fig.subplots_adjust(hspace = 0.3)

    fig.savefig(f'distribution_plots.png')
