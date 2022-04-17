"""
Collection of linear algebra functions mostly related to
geometry operations.
"""
import numpy as np
from numba import jit, njit

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

def _vector_angle(v1, v2):
    """
    Compute the angle between two vectors.

    Parameters
    ----------
    v1: np.ndarray[n, 3]
    v2: np.ndarray[n, 3]

    Returns
    ---------
    float
    """
    angle = np.arccos(np.dot(u_vect(v1), u_vect(v2)))
    return angle

# this is the numba implementation
vector_angle = jit(_vector_angle)

def vector_angle_degrees(v1, v2):
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
    angle = np.degrees(_vector_angle(v1, v2))
    return angle

@njit
def vec2quaternion(vec, theta):
    """
    Translate a vector of 3 elements and angle theta to a quaternion.

    Parameters
    ----------
    vec : numpy 1D-array
        Vector of the quaternion.
    theta : float
        Angle of the quaternion in radian.

    Returns
    -------
    numpy 1D-array
        The full quaternion (4 elements).
    """
    w = np.cos(theta/2)
    x, y, z = np.sin(theta/2) * u_vect(vec)
    return np.array([w, x, y, z])


@njit
def calc_rotation_matrix(quaternion):
    """
    Translate a quaternion to a rotation matrix.

    Parameters
    ----------
    quaternion : numpy 1D-array of 4 elements.

    Returns
    -------
    numpy 2D-array (dimension `[3, 3]`)
        The rotation matrix.
    """
    # Initialize rotation matrix.
    matrix = np.zeros((3, 3), dtype=np.float32)
    # Get quaternion elements.
    w, x, y, z = quaternion
    # Compute rotation matrix.
    matrix[0, 0] = w**2 + x**2 - y**2 - z**2
    matrix[1, 0] = 2 * (x*y + w*z)
    matrix[2, 0] = 2 * (x*z - w*y)
    matrix[0, 1] = 2 * (x*y - w*z)
    matrix[1, 1] = w**2 - x**2 + y**2 - z**2
    matrix[2, 1] = 2 * (y*z + w*x)
    matrix[0, 2] = 2 * (x*z + w*y)
    matrix[1, 2] = 2 * (y*z - w*x)
    matrix[2, 2] = w**2 - x**2 - y**2 + z**2
    return matrix


@njit
def apply_rotation(vec_to_rotate, rotation_axis, rad_angle):
    """Rotate a vector around an axis by a given angle.
    Note
    ----
    The rotation axis is a vector of 3 elements.
    Parameters
    ----------
    vec_to_rotate : numpy 1D-array
    rotation_axis : numpy 1D-array
    rad_angle : float
    Returns
    -------
    numpy 1D-array
        The final rotated (normalized) vector.
    """
    # Generate a quaternion of the given angle (in radian).
    quaternion = vec2quaternion(rotation_axis, rad_angle)
    # Generate the rotation matrix.
    rotation_matrix = calc_rotation_matrix(quaternion)
    # Apply the rotation matrix on the vector to rotate.
    vec_rotated = np.dot(rotation_matrix, vec_to_rotate)
    return u_vect(vec_rotated)


@njit
def cross_product(A, B):
    """
    Return the cross product between vectors A & B.
    Source: http://hyperphysics.phy-astr.gsu.edu/hbase/vvec.html.

    Note
    ----
    On small vectors (i.e. of 3 elements), computing cross products
    with this functions is faster than `np.cross()`.

    Parameters
    ----------
    A : numpy 1D-array
        A vector of 3 elements.
    B : numpy 1D-array
        Another vector of 3 elements.

    Returns
    -------
    numpy 1D-array
        Cross product of A^B.
    """
    x = (A[1]*B[2]) - (A[2]*B[1])
    y = (A[0]*B[2]) - (A[2]*B[0])
    z = (A[0]*B[1]) - (A[1]*B[0])
    return np.array((x, -y, z))
