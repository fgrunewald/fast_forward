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
from collections import defaultdict
import numpy as np
from tqdm import tqdm
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
    u_vect = vect/np.linalg.norm(vect)
    return u_vect

# this is the numba implementation
u_vect = jit(_u_vect)

def _vector_angle_degrees(v1, v2):
    """
    Compute the angle between two vectors
    in degrees and between 0 180 degrees.

    Parameters
    ----------
    v1: np.array
    v2: np.array

    Returns
    ---------
    float
      angle in degrees
    """
    angle = np.degrees(np.arccos(np.dot(u_vect(v1), u_vect(v2))))
    return angle

# this is the numba implementation
vector_angle_degrees = jit(_vector_angle_degrees)

def vector_angle_matrix(matrix_a, matrix_b):
    for v1, v2 in zip(matrix_a, matrix_b):
        yield vector_angle_degrees(v1, v2)

def compute_bond_base(pair_distances, pairs):
    diff = pair_distances[:, pairs[0][0], :] - pair_distances[:, pairs[0][1], :]
    dist = np.apply_along_axis(np.linalg.norm, 1, diff)/10.
    sig = np.std(dist)
    avg = np.average(dist)
    return avg, sig

def compute_angle_base(pair_distances, pairs):
    diff_a = pair_distances[:, pairs[0][0], :] - pair_distances[:, pairs[0][1], :]
    diff_b = pair_distances[:, pairs[1][0], :] - pair_distances[:, pairs[1][1], :]
    dist = list(vector_angle_matrix(diff_a, diff_b))
    sig = np.std(dist)
    avg = np.average(dist)
    return avg, sig

def compute_bond_harm(pair_distances, pairs, const, params):
    avg, sig = compute_bond_base(pair_distances, pairs)
    return [params[0], avg, const/sig**2.0]

def compute_bond_cos_harm(pair_distances, pairs, const, params):
    avg, sig = compute_bond_base(pair_distances, pairs)
    return [params[0], avg, const/(sig*2.*np.sin2(avg))]

def compute_angle_harm(pair_distances, pairs, const, params):
    avg, sig = compute_angle_base(pair_distances, pairs)
    return [params[0], avg, const/sig**2.0]

def compute_angle_cos_harm(pair_distances, pairs, const, params):
    avg, sig = compute_angle_base(pair_distances, pairs)
    return [params[0], avg, const/(sig**2.0*np.sin2(avg))]

def compute_constr(pair_distances, pairs, const, params):
    avg, _ = compute_bond_base(pair_distances, pairs)
    return [params[0], avg]

COMPUTE_INTERACTION = {
    ('bonds', '1'): compute_bond_harm,
    ('bonds', '2'): compute_bond_cos_harm,
    ('angles', '1'): compute_angle_harm,
    ('angles', '2'): compute_angle_cos_harm,
    ('constraints', '1'): compute_constr,
    ('constraints', '2'): compute_constr,
    # ('dihedrals', 1): compute_dih_prop,
    # ('dihedrals', 2): compute_dih_impor,
    # ('dihedrals', 11): compute_dih_impor,
}

def _extract_pairs_and_interactions(molecule):
    pairs = set()
    interactions = []
    for inter_type in ["bonds", "constraints", "dihedrals", "angles"]:
        for inter in molecule.interactions[inter_type]:
            pair_list = list(zip(inter.atoms[:-1], inter.atoms[1:]))
            interactions.append([inter_type, pair_list, inter])
            for pair in pair_list:
                pairs.update({pair})
    return pairs, interactions

def _update_interactions(molecule, interactions):
    molecule.interactions = defaultdict(list)
    for inter_type, _, inter in interactions:
        molecule.interactions[inter_type].append(inter)
    return molecule

def _compute_pair_distances(coordinates, pairs):
    pair_dists = np.zeros((len(pairs), 3))
    idx = 0
    for atom_1, atom_2 in pairs:
        pair_dists[idx][:] = coordinates[atom_1] - coordinates[atom_2]
        idx += 1
    return pair_dists

compute_pairs = _compute_pair_distances

def compute_interaction_parameters(universe, molecule, const):
    """
    Given positions in a `universe` and `molecule` as well
    as RxT (i.e. `const`) perform a Boltzmanntype inversion
    to give a best estimate of potentials in molecule.

    Parameters:
    -----------
    universe: :class:`mda.base.core.Universe`
    molecule: :class:`vermouth.molecule.Molecule`
    const: float

    Returns:
    --------
    :class:`vermouth.molecule.Molecule`
        molecule with best fit parameters
    """
    pairs, interactions = _extract_pairs_and_interactions(molecule)
    pair_distances = np.zeros((len(universe.trajectory), len(pairs), 3))

    fdx = 0
    for _ in tqdm(universe.trajectory):
        pair_distances[fdx, :] = compute_pairs(universe.atoms.positions, pairs)[:]
        fdx += 1

    idx = 0
    for inter_type, pairs, inter in interactions:
        functype = inter.parameters[0]
        try:
            new_params = COMPUTE_INTERACTION[(inter_type, functype)](pair_distances,
                                                             pairs,
                                                             const,
                                                             inter.parameters)
            inter.parameters[:] = new_params[:]
        except KeyError:
            msg=("WARNING - Interaction type {} function {} not implemented. "
                 "Skipping evaluation.".format(inter_type, functype))
            print(msg)
        idx+=1
    # now update the molecule
    molecule = _update_interactions(molecule, interactions)
    return molecule
