"""
Extract interactions pairs from itp and convert to atom
groups of an MDAnalysis Universe.
"""
from collections import defaultdict
import numpy as np
import networkx as nx
from fast_forward.universe_handler import res_as_mol

def find_mol_indices(universe, atoms, moltype):
    """
    Given a universe select all atoms that belong to molecules 
    of the given `moltype`. Subsequently, return indices of all
    multiples of the indices defined in atoms under the assumption
    that all considered molecules have the same number of atoms.

    Parameters
    ----------
    universe: mda.Universe
    atoms: list[int]
    moltype: str

    Returns
    -------
    list[numpy.array(dtype=int)]
    """
    mol_atoms = universe.select_atoms(f'moltype {moltype}')
    n_mols = len(np.unique(mol_atoms.molnums))
    try:
        mol_atom_indices = mol_atoms.indices.reshape(n_mols, -1)
    except ValueError:
        msg = ("The target molecules passed to find_mol_indices "
               "do not seem to all have the same number of atoms.")
        raise IndexError(msg) from None
    return list(mol_atom_indices[:, atoms])

class ITPInteractionMapper:
    """
    Iterate over interactions in itp file and return dict of
    grouped indices corresponding to the atoms in universe.
    """
     # we ensure we either have molecule types, or we promote res info as such
    res_as_mol(universe)

    indices_dict = defaultdict(dict)
    initial_parameters = defaultdict(dict)
    block_indices = defaultdict(dict)
    for inter_type, block_inter in block.interactions.items():
        for inter in block_inter:
            atoms = inter.atoms
            group = inter.meta.get("comment")
            if group is not None:
                indices = find_mol_indices(universe, atoms, mol_name)
                old_indices = indices_dict[inter_type].get(group, [])
                old_block_indices = block_indices[inter_type].get(group, [])
>>>>>>> main

        indices_dict = defaultdict(dict)
        initial_parameters = defaultdict(dict)
        block_indices = defaultdict(dict)
        for inter_type in block.interactions:
            for inter in block.interactions[inter_type]:
                atoms = inter.atoms
                if itp_mode == "all":
                    atomnames=[]
                    for atom in atoms:
                        atomnames.append(block.nodes[atom]['atomname'])
                    group = "_".join(atomnames)
                    inter.meta["comment"] = group
                else:
                    group = inter.meta.get("comment", None)
                if group:
                    indices = find_indices(self.universe, atoms, molname)
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
