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

def itp_to_ag(block, mol_name, universe):
    """
    Iterate over interactions in itp file and return dict of
    grouped indices corresponding to the atoms in universe.
    """
    # by default we try to match the molecule types
    has_molnums = hasattr(universe.atoms, "moltypes")
    match_attr = "moltypes"
    match_values = [mol_name]
    # if we don't have molecule types we go by residues
    # this requires there to be no overlap between the
    # target and other molecules
    if not has_molnums:
        resnames = nx.get_node_attributes(block, "resname")
        match_values = list(set(resnames.values()))
        match_attr = "resnames"

    indices_dict = defaultdict(dict)
    initial_parameters = defaultdict(dict)
    for inter_type in block.interactions:
        for inter in block.interactions[inter_type]:
            atoms = inter.atoms
            group = inter.meta.get("comment", None)
            if group:
                indices = find_indices(universe,
                                       atoms,
                                       match_attr,
                                       match_values,
                                       natoms=len(block.nodes))
                old_indices = indices_dict[inter_type].get(group, [])
                indices_dict[inter_type][group] = indices + old_indices
                initial_parameters[inter_type][group] = inter.parameters

    return indices_dict, initial_parameters

def itp_to_ag_pairs(block, mol_name, universe):
    """
    Iterate over all atoms in itp file and return dict of
    paired indices corresponding to all pairwise distances in universe.
    group_names are given by the atomnames
    """
    # by default we try to match the molecule types
    has_molnums = hasattr(universe.atoms, "moltypes")
    match_attr = "moltypes"
    match_values = [mol_name]
    # if we don't have molecule types we go by residues
    # this requires there to be no overlap between the
    # target and other molecules
    if not has_molnums:
        resnames = nx.get_node_attributes(block, "resname")
        match_values = list(set(resnames.values()))
        match_attr = "resnames"

    indices_dict = defaultdict(dict)
    # loop through all nodes twice, with the second loop starting with the next node
    for node1, name1 in block.nodes(data='atomname'):
        for node2, name2 in list(block.nodes(data='atomname'))[node1+1:]:
            atoms = np.array([node1, node2])
            group = f'{name1}_{name2}'
            indices = find_indices(universe,
                                    atoms,
                                    match_attr,
                                    match_values,
                                    natoms=len(block.nodes))
            old_indices = indices_dict['distances'].get(group, [])
            indices_dict['distances'][group] = indices + old_indices
    return indices_dict