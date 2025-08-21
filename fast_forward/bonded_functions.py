import numpy as np
from numba import njit, prange
from fast_forward.linalg_functions import vector_angle_degrees, dihedral_angle
from lmfit import create_params, minimize


@njit(parallel=True)
def _fast_pair_dists(arr1, arr2):
    frames = arr1.shape[0]
    pair_dists = np.zeros((frames))
    for fdx in prange(frames):
        diff = arr1[fdx] - arr2[fdx]
        norm = (diff[0] ** 2.0 + diff[1] ** 2.0 + diff[2] ** 2.0) ** 0.5
        pair_dists[fdx] = norm
    return pair_dists


@njit(parallel=True)
def _fast_angle(arr1, arr2, arr3):
    frames = arr1.shape[0]
    angle = np.zeros((frames))
    for fdx in prange(frames):
        diff1 = arr2[fdx] - arr1[fdx]
        diff2 = arr2[fdx] - arr3[fdx]
        angle[fdx] = vector_angle_degrees(diff1, diff2)
    return angle


@njit(parallel=True)
def _fast_dih(arr1, arr2, arr3, arr4):
    frames = arr1.shape[0]
    angle = np.zeros((frames))
    for fdx in prange(frames):
        diff1 = arr1[fdx] - arr2[fdx]
        diff2 = arr3[fdx] - arr2[fdx]
        diff3 = arr3[fdx] - arr4[fdx]
        angle[fdx] = dihedral_angle(diff1, diff2, diff3)
    return angle


def _vs3fd_func(pars, x, positions):
    vals = pars.valuesdict()
    a = vals['a']
    b = vals['b']

    i, j, k = positions

    r_ij = j - i
    r_jk = k - j
    r_ijk = r_ij + a * r_jk

    r_ijk_norm = (r_ijk.T / np.linalg.norm(r_ijk, axis=1)).T
    r_predicted_vs = i + (b * r_ijk_norm)

    return x - r_predicted_vs


def _vs3fd(arr1, arr2, arr3, arr4):
    fit_params = create_params(a=1, b=1)
    pos_list = [arr2 / 10, arr3 / 10, arr4 / 10]  # convert to nm
    fit = minimize(_vs3fd_func, fit_params, args=(arr1 / 10, pos_list,))
    ab = np.array([fit.params["a"].value, fit.params["b"].value])
    return ab


def _vs3out_func(pars, x, positions):
    vals = pars.valuesdict()
    a = vals['a']
    b = vals['b']
    c = vals['c']

    i, j, k = positions

    r_ij = j - i
    r_ik = k - i
    r_predicted_vs = i + (a * r_ij) + (b * r_ik) + (c * np.cross(r_ij, r_ik))

    return x - r_predicted_vs


def _vs3out(arr1, arr2, arr3, arr4):
    fit_params = create_params(a=1, b=1, c=1)
    pos_list = [arr2 / 10, arr3 / 10, arr4 / 10]  # convert to nm
    fit = minimize(_vs3out_func, fit_params, args=(arr1 / 10, pos_list,))
    abc = np.array([fit.params["a"].value, fit.params["b"].value, fit.params["c"].value])

    return abc


def _virtual_sitesn(*args):
    return None


NORMAL_FUNCS = {"angles": _fast_angle,
                "bonds": _fast_pair_dists,
                "constraints": _fast_pair_dists,
                "dihedrals": _fast_dih,
                "virtual_sites3fd": _vs3fd,
                "virtual_sites3out": _vs3out,
                "virtual_sitesn": _virtual_sitesn,
                "distances": _fast_pair_dists
                }
