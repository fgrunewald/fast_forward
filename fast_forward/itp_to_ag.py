"""
Extract interactions pairs from itp and convert to atom
groups of an MDAnalysis Universe.
"""
from collections import defaultdict
import numpy as np
import networkx as nx

def find_indices(universe,
                 atoms,
                 match_attr,
                 match_values,
                 natoms):
    """
    Given a universe select all atoms that match 
    with the value of `match_attr` one or more entries in
    `match_values`. Subsequently, return indices of all
    multiples of the indices defined in atoms under the
    assumption that a single molecule has `natoms`.

    Parameters
    ----------
    universe: mda.Universe
    atoms: list[int]
    match_attr: abc.hashable
    match_values: list[abc.hashable]
    natoms: int

    Returns
    -------
    list[int]
    """
    indices = []
    atoms = np.array(atoms)
    mol_atoms = universe.atoms[np.isin(getattr(universe.atoms, match_attr), match_values)]
    if len(mol_atoms) % natoms != 0 and len(mol_atoms) != 0:
        msg = ("The number of atoms of the target molecule"
               "does not match a integer multiple of atoms"
               "selected from the universe.")
        raise IndexError(msg)
    n_mols = len(mol_atoms) // natoms
    for idx in range(0, n_mols):
        pairs = mol_atoms.indices[atoms + idx * natoms]
        indices.append(pairs)
    return indices

class ITPInteractionMapper:

    def __init__(self, universe, blocks, molnames):
        self.universe = universe
        self.blocks = dict(zip(molnames, blocks))
        # by default we try to match the molecule types
        has_molnums = hasattr(universe.atoms, "moltypes")
        self.match_attr = "moltypes"
        self.match_values = {}
        for molname in molnames:
            self.match_values[molname] = [molname]
        # if we don't have molecule types we go by residues
        # this requires there to be no overlap between the
        # target and other molecules
        if not has_molnums:
            for block, molname in self.blocks.items():
                resnames = nx.get_node_attributes(block, "resname")
                self.match_values[molname] = list(set(resnames.values()))
            self.match_attr = "resnames"

    def get_interactions_group(self, molname):
        """
        Iterate over interactions in itp file and return dict of
        grouped indices corresponding to the atoms in universe.
        """
        block = self.blocks[molname]

        indices_dict = defaultdict(dict)
        initial_parameters = defaultdict(dict)
        block_indices = defaultdict(dict)
        for inter_type in block.interactions:
            for inter in block.interactions[inter_type]:
                atoms = inter.atoms
                group = inter.meta.get("comment", None)
                if group:
                    indices = find_indices(self.universe,
                                        atoms,
                                        self.match_attr,
                                        self.match_values[molname],
                                        natoms=len(block.nodes))
                    old_indices = indices_dict[inter_type].get(group, [])
                    old_block_indices = block_indices[inter_type].get(group, [])
                    block_indices[inter_type][group] = [atoms] + old_block_indices


                    indices_dict[inter_type][group] = indices + old_indices
                    initial_parameters[inter_type][group] = inter.parameters

        return indices_dict, initial_parameters, block_indices

    def get_pairwise_interaction(self, molname):
        """
        Iterate over all atoms in itp file and return dict of
        paired indices corresponding to all pairwise distances in universe.
        group_names are given by the atomnames
        """
        block = self.blocks[molname]

        indices_dict = defaultdict(dict)
        
        for node1, name1 in block.nodes(data='atomname'):
            for node2, name2 in list(block.nodes(data='atomname'))[node1+1:]:
                atoms = np.array([node1, node2])
                group = f'{name1}_{name2}' # naming convention with node1 < node2
                indices = find_indices(self.universe,
                                        atoms,
                                        self.match_attr,
                                        self.match_values[molname],
                                        natoms=len(block.nodes))
                old_indices = indices_dict['distances'].get(group, [])
                indices_dict['distances'][group] = indices + old_indices

        return indices_dict
