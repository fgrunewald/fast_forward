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
from tqdm import tqdm
from . import jit

def get_pairs(inter_type, inter):
    """
    Find and get all distances over trajectory
    for `pair`. This function is not order
    sensetive to the pair indices.

    Parameters:
    -----------
    pair: tuple[int, int]

    Returns:
    --------
    np.ndarray.slice
    """
    if inter_type in ["bonds", "constraints", "angles", "dihedrals"]:
        return list(zip(inter.atoms[:-1], inter.atoms[1:]))

    return []

def extract_pairs(molecule):
    """
    Lookup all interactions in the molecule attribute and write all
    pairs to a set. Note that pairs can only be generated for known
    interactions, which currently include bonds, angles, constraints,
    dihedrals of proper type and imporper type.

    Parameters:
    -----------
     :class:`MDAnalysis.Universe`

    Returns:
    --------
    set[tuple(int, int)]
        set of all pair interactions involved in interactions
    """
    pairs = set()
    for inter_type in ["bonds", "constraints", "dihedrals", "angles"]:
        for inter in molecule.interactions[inter_type]:
            pair_list = get_pairs(inter_type, inter)
            for pair in pair_list:
                pairs.update({pair})
    return pairs


class InteractionUniverse:
    """
    This calss links interactions of a :class:`vermouth.molecule.Molecule`
    to a :class:`MDAnalysis.Universe` class and provides a framework to
    lookup coordinates of the atoms involved in an interaction within
    the MDAnalysis Universe in an efficent manner. To that end all pairs
    are extracted into a set of pairs. For the set of pairs all distances
    are computed when the method `compute_pair_distances` is called. All
    pairs are linked to these distances by the following conventions:
    """

    def __init__(self, universe, molecule):
        """
        Extract the interactions from a `molecule` and
        prepare a list of indixes referring to the appropiate
        coordinates in the universe. Both molecule and universe
        become instance variables that are modifed as needed.

        Parameters:
        -----------
        molecule: :class:`vermouth.molecule.Molecule`
        universe:  :class:`MDAnalysis.Universe`
        """
        self.molecule = molecule
        self.universe = universe
        self.pairs = list(extract_pairs(molecule))
        self.dist_to_pair = dict()
        self.pair_distances = np.zeros((universe.trajectory.n_frames, len(self.pairs), 3))
        self.pair_to_idx = dict()
        # set some base arguments for the iterator
        self.type_idx = 0
        self.inter_idx = 0
        self.inter_types = list(self.molecule.interactions.keys())

        #Numba accelerated functions
        self.compute_pair_distances = self._compute_pair_distances

        total = 0
        for inter_type in self.molecule.interactions:
            total += len(self.molecule.interactions[inter_type])

        self.n_interactions = total

    def __iter__(self):
        self.type_idx = 0
        self.inter_idx = 0
        return self

    def __next__(self):
            try:
                inter_type = self.inter_types[self.type_idx]
                interaction = self.molecule.interactions[inter_type][self.inter_idx]
                idx = self.inter_idx
                self.inter_idx += 1
            except IndexError:
                self.type_idx  += 1
                idx = 0
                self.inter_idx = 0
                try:
                    inter_type = self.inter_types[self.type_idx]
                    interaction = self.molecule.interactions[inter_type][self.inter_idx]
                except IndexError:
                    raise StopIteration

            pairs_idxs = [ self.get_pair_index(pair) for pair in zip(interaction.atoms[:-1], interaction.atoms[1:])]
            dists = self.pair_distances[:, pairs_idxs, :]

            return inter_type, interaction, dists, idx

    def _compute_pair_distances(self, verbose=False):
        """
        For all pairs in self.pairs compute all distances
        between the pairs and fill in self.pairs_dists.
        """
        if verbose:
            pbar = tqdm(total=len(self.pairs))

        idx = 0
        for atom_1, atom_2 in self.pairs:
            key = frozenset([atom_1, atom_2])
            self.pair_to_idx[key] = idx
            fdx = 0
            for _ in self.universe.trajectory:
                coordinates = self.universe.atoms.positions
                self.pair_distances[fdx, idx, :] = coordinates[atom_1] - coordinates[atom_2]
                fdx += 1
            idx += 1
            if verbose:
                pbar.update(1)

        if verbose:
            pbar.close()

    def get_pair_index(self, pair):
        """
        Find and get all distances over trajectory
        for `pair`. This function is not order
        sensetive to the pair indices.

        Parameters:
        -----------
        pair: tuple[int, int]

        Returns:
        --------
        np.ndarray.slice
        """
        pair = frozenset(list(pair))
        return self.pair_to_idx[pair]
