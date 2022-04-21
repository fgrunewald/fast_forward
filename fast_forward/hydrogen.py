"""
Module to reconstruct hydogens from a group of atoms.
Most of this comes from  buildh.readthedocs.io and
is proted to be compatible with the vermouth way of
handling molecules.
"""
import fnmatch
import numpy as np
from fast_forward.linalg_functions import u_vect, apply_rotation, cross_product, vector_angle

# Constants.
# From https://en.wikipedia.org/wiki/Carbon%E2%80%93hydrogen_bond
LENGTH_CH_BOND = 1.09  # in Angst
# From https://en.wikipedia.org/wiki/Tetrahedron, tetrahedral angle equals
# arccos(-1/3) ~ 1.9106 rad ~ 109.47 deg.
TETRAHEDRAL_ANGLE = np.arccos(-1/3)


def get_CH(atom, helper1, helper2, helper3):
    """
    Reconstruct the unique hydrogen of a sp3 carbon.

    Parameters
    ----------
    atom : numpy 1D-array
        Central atom on which we want to reconstruct the hydrogen.
    helper1 : numpy 1D-array
        First neighbor of central atom.
    helper2 : numpy 1D-array
        Second neighbor of central atom.
    helper3 : numpy 1D-array
        Third neighbor of central atom.

    Returns
    -------
    numpy 1D-array
        Coordinates of the rebuilt hydrogen: `([x_H, y_H, z_H])`.
    """
    # Calculate vector along tetrahedron median.
    # !!! Important !!! Use unit vectors (some bonds may have different length).
    v2 = u_vect(helper1 - atom) + u_vect(helper2 - atom) \
         + u_vect(helper3 - atom)
    # CH bond is on the opposite direction.
    unit_vect_H = u_vect(-v2)
    coor_H = LENGTH_CH_BOND * unit_vect_H + atom
    return coor_H


def get_CH2(atom, helper1, helper2):
    """
    Reconstruct the 2 hydrogens of a sp3 carbon (methylene group).

    Parameters
    ----------
    atom : numpy 1D-array
        Central atom on which we want to reconstruct hydrogens.
    helper1 : numpy 1D-array
        Heavy atom before central atom.
    helper2 : numpy 1D-array
        Heavy atom after central atom.

    Returns
    -------
    tuple of numpy 1D-arrays
        Coordinates of the two hydrogens:
        `([x_H1, y_H1, z_H1], [x_H2, y_H2, z_H2])`.
    """
    # atom->helper1 vector.
    v2 = u_vect(helper1 - atom)
    # atom->helper2 vector.
    v3 = u_vect(helper2 - atom)
    # Vector orthogonal to the helpers/atom plane.
    #v4 = normalize(np.cross(v3, v2))
    v4 = u_vect(cross_product(v3, v2))
    # Rotation axis is atom->helper1 vec minus atom->helper2 vec.
    rotation_axis = u_vect(v2 - v3)
    # Vector to be rotated by theta/2, perpendicular to rotation axis and v4.
    vec_to_rotate = u_vect(cross_product(v4, rotation_axis))
    # Reconstruct the two hydrogens.
    unit_vect_H1 = apply_rotation(vec_to_rotate, rotation_axis,
                                      -TETRAHEDRAL_ANGLE/2)

    hcoor_H1 = LENGTH_CH_BOND * unit_vect_H1 + atom

    unit_vect_H2 = apply_rotation(vec_to_rotate, rotation_axis,
                                      TETRAHEDRAL_ANGLE/2)
    hcoor_H2 = LENGTH_CH_BOND * unit_vect_H2 + atom
    return (hcoor_H1, hcoor_H2)


def get_CH3(atom, helper1, helper2):
    """
    Reconstruct the 3 hydrogens of a sp3 carbon (methyl group).

    Parameters
    ----------
    atom : numpy 1D-array
        Central atom on which we want to reconstruct hydrogens.
    helper1 : numpy 1D-array
        Heavy atom before central atom.
    helper2 : numpy 1D-array
        Heavy atom before helper1 (two atoms away from central atom).

    Returns
    -------
    tuple of numpy 1D-arrays
        Coordinates of the 3 hydrogens:
        `([x_H1, y_H1, z_H1], [x_H2, y_H2, z_H2], [x_H3, y_H3, z_H3])`.
    """
    ### Build CH3e.
    theta = TETRAHEDRAL_ANGLE
    # atom->helper1 vector.
    v2 = helper1 - atom
    # atom->helper2 vector.
    v3 = helper2 - atom
    # Rotation axis is perpendicular to the atom/helpers plane.
    #rotation_axis = normalize(np.cross(v3, v2))
    rotation_axis = u_vect(cross_product(v3, v2))
    # Rotate v2 by tetrahedral angle. New He will be in the same plane
    # as atom and helpers.
    unit_vect_He = apply_rotation(v2, rotation_axis, theta)
    coor_He = LENGTH_CH_BOND * unit_vect_He + atom
    ### Build CH3r.
    theta = (2/3) * np.pi
    rotation_axis = u_vect(helper1 - atom)
    v4 = u_vect(coor_He - atom)
    # Now we rotate atom->He bond around atom->helper1 bond by 2pi/3.
    unit_vect_Hr = apply_rotation(v4, rotation_axis, theta)
    coor_Hr = LENGTH_CH_BOND * unit_vect_Hr + atom
    ### Build CH3s.
    theta = -(2/3) * np.pi
    rotation_axis = u_vect(helper1 - atom)
    v5 = u_vect(coor_He - atom)
    # Last we rotate atom->He bond around atom->helper1 bond by -2pi/3.
    unit_vect_Hs = apply_rotation(v5, rotation_axis, theta)
    coor_Hs = LENGTH_CH_BOND * unit_vect_Hs + atom
    return coor_He, coor_Hr, coor_Hs


def get_CH_double_bond(atom, helper1, helper2):
    """Reconstruct the hydrogen of a sp2 carbon.
    Parameters
    ----------
    atom : numpy 1D-array
        Central atom on which we want to reconstruct the hydrogen.
    helper1 : numpy 1D-array
        Heavy atom before central atom.
    helper2 : numpy 1D-array
        Heavy atom after central atom.
    Returns
    -------
    tuple of numpy 1D-arrays
        Coordinates of the rebuilt hydrogen: `([x_H, y_H, z_H])`.
    """
    # calc CCC_angle helper1-atom-helper2 (in rad).
    v1 = helper1 - atom
    v2 = helper2 - atom
    CCC_angle = vector_angle(v1, v2)
    # We want to bisect the C-C-C angle ==> we take half of (2pi-CCC_angle).
    # Factorizing yields: pi - CCC_angle/2.
    theta = np.pi - (CCC_angle / 2)
    # atom->helper1 vector.
    v2 = helper1 - atom
    # atom->helper2 vector.
    v3 = helper2 - atom
    # The rotation axis is orthogonal to the atom/helpers plane.
    #rotation_axis = normalize(np.cross(v2, v3))
    rotation_axis = u_vect(cross_product(v2, v3))
    # Reconstruct H by rotating v3 by theta.
    unit_vect_H = apply_rotation(v3, rotation_axis, theta)
    coor_H = LENGTH_CH_BOND * unit_vect_H + atom
    return coor_H

def match_attributes(atom, attributes):
    for attr, value in attributes.items():
        try:
            assert fnmatch.fnmatch(getattr(atom, attr), value)
        except AttributeError:
            return False
    return True

def find_any_bonded_neighbor(atomgroup, attributes):
    """
    Find any bonded neighbor of atomgroup.

    Parameters
    ----------
    atomgroup: :class:`MDAnalysis.atomgroup`
        atom-group of single atom
    attributes: dict
        a dict of attributes the atom has to full-fill
        only string attributes are allowed

    Returns
    ------
    atom
        atomgroup of bonded neighbor
    """
    bond_list = atomgroup.get_connections("bonds")
    for atom1, atom2 in bond_list:
        if atom1 != atomgroup and match_attributes(atom1, attributes):
            return atom1

        if match_attributes(atom2, attributes):
            return atom2

REF_ATOMS =  {"CH3": [0, 1],
              "CH2": [0, 0],
              "CH": [0, 0, 0],
              "CH2d": [0, 0],
             }

BUILD_HYDRO = {"CH3": get_CH2,
               "CH2": get_CH2,
               "CH": get_CH,
               "CH2d": get_CH_double_bond,
              }

def find_helper_atoms(atom, carbon_type):
    """
    Given a carbon-type find all helper atoms
    needed in the reconstruction of the coordinates.
    """
    construction_atoms = [atom]
    for help_idx in REF_ATOMS[carbon_type]:
        anchor = construction_atoms[help_idx]
        neighbor = find_any_bonded_neighbor(anchor, {})
        construction_atoms.append(neighbor)
    return dict(zip(["atom", "helper1", "helper2", "helper3"], construction_atoms))
