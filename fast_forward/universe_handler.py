# Copyright 2021 University of Groningen
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
import networkx as nx
import MDAnalysis as mda
from MDAnalysis import transformations
from fast_forward.hydrogen import BUILD_HYDRO, find_helper_atoms

class UniverseHandler(mda.Universe):
    """
    Wrapper around mda.Universe which allows for
    some smart selections and exposes some useful
    properties.
    """
    def __init__(self, mol_names, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # select which molecules to treat
        self.mol_names = mol_names
        self.molecules = self.select_atoms("moltype " + " ".join(mol_names))

        # set some useful attributes
        self.mol_idxs_by_name = defaultdict(list)
        for mol_num, mol_name in dict(zip(self.atoms.molnums, self.atoms.moltypes)).items():
            if mol_name in mol_names:
                self.mol_idxs_by_name[mol_name].append(mol_num)

        mol_idxs_names = [(idx, mol_name) for idx, mol_name in self.mol_idxs_by_name.items()]
        self.normal_order = mol_idxs_names.sort()

        # Hidden attribute labels
        self.__n_atoms = None
        self.__n_residues = None
        self.__pbc_completed = False
        self.__resids = None

    def pbc_complete(self):
        if not self.__pbc_completed:
            self.trajectory.add_transformations(transformations.unwrap(self.atoms))
            self.__pbc_completed = True
        return self.__pbc_completed

    def shift_united_atom_carbons(self, association_dict):
        """
        Given an atomgroup shift it's coordiantes
        to where the center of geomtry would be if
        hydrogens were included.
        """
        # much faster to make a bonded graph even of the entire system
        bonded_graph = nx.Graph()
        bonded_graph.add_edges_from(list(self.atoms.get_connections("bonds").indices))
        # select the atoms to be treated
        select_string = "type " + " ".join(association_dict.keys())
        atoms_to_treat = self.select_atoms(select_string)
        for atom in atoms_to_treat:
            carbon_type = association_dict[atom.type]
            helper_atoms = find_helper_atoms(self, atom, carbon_type, bonded_graph)
            for ts in self.trajectory:
                hydrogen_coords = BUILD_HYDRO[carbon_type](**helper_atoms)
                new_pos = atom.position
                for hydro_coord in hydrogen_coords:
                    new_pos += hydro_coord
                atom.position = new_pos / (len(hydrogen_coords) + 1)
    @property
    def n_atoms(self):
        if self.__n_atoms is None:
            self.__n_atoms = len(self.molecules)
        return self.__n_atoms

    @property
    def resids(self):
        if self.__resids is None:
            self.__resids = self.molecules.resids
        return self.__resids

    @property
    def n_residues(self):
        if self.__n_residues is None:
            self.__n_residues = len(self.molecules.residues)
        return self.__n_residues

    def res_iter(self):
        for residue in self.molecules.residues:
            yield residue
