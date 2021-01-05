# Copyright 2020 University of Groningen
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
import inspect
import numpy as np
import scipy.stats
import symfit as sf
from .linalg_functions import vector_angle_degrees
import matplotlib.pyplot as plt

def harmonic(x, ref, k):
    return 0.5*k*(x-ref)**2.0

def harmonic_cosine(x, ref, k):
    return 0.5*k*(np.cos(x) - np.cos(ref))**2.0


GMX_POTENTIALS = {("bonds", 1): harmonic,
                  ("bonds", 2): harmonic_cosine,
                  ("angles", 1): harmonic}

class InteractionBase:

    def __init__(self, name, func, potential, **kwargs):
        self.potential = potential
        self.parameters = kwargs
        self.name = name
        self.func = func

    def update_parameters(self, new_parameters):
        self.parameters.update(new_parameters)

    def energy(self, x):
        return self.potential(x, **self.parameters)

  # def force(self, x):
  #     return self.derivative(x, **self.parameters)

    def symfit_model(self):
        sf_parameters = sf.parameters(", ".join(self.parameters.keys()))
        x, y = sf.variables('x, y')
        model = {y: self.potential(x, *sf_parameters)}
        return model

    def vermouth_parameters(self):
        parameters = [str(self.func)] + list(map(str, self.parameters.values()))
        return parameters

    @classmethod
    def from_gmx_lib(cls, inter_type, functype, init_params=None):
        potential = GMX_POTENTIALS[(inter_type, functype)]
        if init_params:
            param_dict = dict(zip(inspect.signature(potential).parameters.keys()[1:], init_params))
        else:
            parameters = list(inspect.signature(potential).parameters.keys())[1:]
            param_dict = dict(zip(parameters, np.ones(len(parameters))))
        return cls(inter_type, functype, potential, **param_dict)


def _vector_angle_matrix(matrix_a, matrix_b):
    for v1, v2 in zip(matrix_a, matrix_b):
        yield vector_angle_degrees(v1, v2)


def _bond_distance_matrix(matrix_a):
    for v1 in matrix_a:
        yield np.linalg.norm(v1)


NORMAL_FUNCS = {"angles": _vector_angle_matrix,
                "bonds": _bond_distance_matrix,
                "constraints": _bond_distance_matrix}

def trimm_histogramm(counts, bins):
    # following answer from https://stackoverflow.com/questions/38161606/find-the-start-position-of-the-longest-sequence-of-1s
    non_zero = counts != 0
    # Get start, stop index pairs for islands/seq. of 1s
    idx_pairs = np.where(np.diff(np.hstack(([False], non_zero.astype(int) == 1,[False]))))[0].reshape(-1,2)
    # Find start stop value of longest island
    start_stop = idx_pairs[np.diff(idx_pairs,axis=1).argmax(), :]
    # trimm counts
    trimmed_counts = counts[start_stop[0]:start_stop[1]]
    # trimm bins
    trimmed_bins = bins[start_stop[0]:start_stop[1]+1]
    return trimmed_counts, trimmed_bins

def symfit_interactions(inter_type, interaction, distances, temp=298.15, gas_const=8.314):
    # compute RT
    const = temp * gas_const

    # get pair vectors
    if inter_type == 'bonds':
        vectors = [distances[:, 0, :]]
        tau = 0.005
    elif inter_type == "angles":
        vectors = [distances[:, 0, :], distances[:, 1, :]]
        tau = 2.5

    # make a histogram of the time-series
    timeseries = np.fromiter(NORMAL_FUNCS[inter_type](*vectors), dtype=float)
    counts, bins = np.histogram(timeseries, bins=np.arange(min(timeseries), max(timeseries), tau))
    # we need to clean the histogramm
    counts, bins = trimm_histogramm(counts, bins)
    # supsample the log pdf to get data for fitting
    xdata = bins[1:-1]
    hist_dist = scipy.stats.rv_histogram((counts, bins))
    log_pdf = const * -1 * hist_dist.logpdf(xdata)
    # shift minimum to zero
    log_pdf = log_pdf - min(log_pdf)
    # get the symfit model
    inter_func = InteractionBase.from_gmx_lib(inter_type, int(interaction.parameters[0]))
    fit = sf.Fit(inter_func.symfit_model(), xdata, log_pdf)
    result = fit.execute()
    np.savetxt("test.dat", log_pdf)

    # update final parameters
    inter_func.update_parameters(result.params)
    interaction.parameters[:] = inter_func.vermouth_parameters()[:]
    return interaction
