import argparse
import numpy as np
from numba import jit, njit, prange
import MDAnalysis as mda
from numpy import cross, dot

def _u_vect(vect):
    """
    Compute unit vector of vector.

    Parameters
    ----------
    vect: np.array

    Returns
    -----------
    np.array
    """
    unit_vect = vect/np.linalg.norm(vect)
    return unit_vect

# this is the numba implementation
u_vect = jit(_u_vect)

def _vector_angle_degrees(v1, v2):
    """
    Compute the angle between two vectors
    in degrees and between 0 180 degrees.

    Parameters
    ----------
    v1: np.ndarray[n, 3]
    v2: np.ndarray[n, 3]

    Returns
    ---------
    float
    """
    angle = np.degrees(np.arccos(np.dot(u_vect(v1), u_vect(v2))))
    return angle

# this is the numba implementation
vector_angle_degrees = jit(_vector_angle_degrees)

def _dih(r1, r2, r3):
    cross1 = cross(r1, r2)
    cross2 = cross(r2, r3)
    n1 = u_vect(cross1)
    n2 = u_vect(cross2)
    dih = vector_angle_degrees(n1, n2)
    # GROMACS specific definition of the sign of the
    # dihedral.
    if dot(r1, cross2) < 0:
        dih = - dih
    return dih


dihedral_angle = jit(_dih)

#@njit(parallel=True)
def _fast_pair_dists(arr1, arr2):
    frames = arr1.shape[0]
    pair_dists = np.zeros((frames))
    for fdx in prange(frames):
        diff = arr1[fdx] - arr2[fdx]
        norm = (diff[0]**2.0 + diff[1]**2.0 + diff[2]**2.)**0.5
        pair_dists[fdx] = norm
    return pair_dists

#@njit(parallel=True)
def _fast_angle(arr1, arr2, arr3):
    frames = arr1.shape[0]
    angle= np.zeros((frames))
    for fdx in prange(frames):
        diff1 = arr2[fdx] - arr1[fdx]
        diff2 = arr2[fdx] - arr3[fdx]
        angle[fdx] = vector_angle_degrees(diff1, diff2)
    return angle

#@njit(parallel=True)
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
