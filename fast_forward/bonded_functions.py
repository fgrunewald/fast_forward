import numpy as np
from numba import njit, prange
from fast_forward.linalg_functions import vector_angle_degrees, dihedral_angle

@njit(parallel=True)
def _fast_pair_dists(arr1, arr2):
    frames = arr1.shape[0]
    pair_dists = np.zeros((frames))
    for fdx in prange(frames):
        diff = arr1[fdx] - arr2[fdx]
        norm = (diff[0]**2.0 + diff[1]**2.0 + diff[2]**2.)**0.5
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
    angle= np.zeros((frames))
    for fdx in prange(frames):
        diff1 = arr1[fdx] - arr2[fdx]
        diff2 = arr3[fdx] - arr2[fdx]
        diff3 = arr3[fdx] - arr4[fdx]
        angle[fdx] = dihedral_angle(diff1, diff2, diff3)
    return angle

NORMAL_FUNCS = {"angles": _fast_angle,
                "bonds": _fast_pair_dists,
                "constraints": _fast_pair_dists,
                "dihedrals": _fast_dih
               }
