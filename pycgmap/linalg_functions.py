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
import numpy as np
from . import jit

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
    """
    Compute dihedral angle between four points,
    where `B` and `C` are in the center.
    Paramters
    ---------
    A, B, C, D:  numpy.array
    Returns
    ---------
    float
         angle in degrees
    """
    n1 = u_vect(np.cross(r1, r2))
    n2 = u_vect(np.cross(r2, r3))
    dih = vector_angle_degrees(n1, n2)
    return dih

dihedral_angle = jit(_dih)
