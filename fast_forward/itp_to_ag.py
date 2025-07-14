"""
Extract interactions pairs from itp and convert to atom
groups of an MDAnalysis Universe.
"""
import numpy as np
from collections import defaultdict

def find_indices(universe, atoms, molname, natoms):
    """
    Find the atoms indices for all molecules
    with moleculename molname in universe and
    return list of indices.
    """
    indices = []
    atoms = np.array(atoms)
    mol_atoms = universe.atoms[universe.atoms.moltypes == molname]
    n_mols = len(set(mol_atoms.molnums))
    for idx in range(0, n_mols):
        pairs = mol_atoms.indices[atoms + idx * natoms]
        indices.append(pairs)
    return indices

def itp_to_ag(block, mol_name, universe):
    """
    Iterate over interactions in itp file and return dict of
    grouped indices corresponding to the atoms in universe.
    """
    indices_dict = defaultdict(dict)
    for inter_type in block.interactions:
        for inter in block.interactions[inter_type]:
            atoms = inter.atoms
            group = inter.meta.get("comment", None)
            if group:
                indices = find_indices(universe, atoms, mol_name, natoms=len(block.nodes))
                old_indices = indices_dict[inter_type].get(group, [])
                indices_dict[inter_type][group] = indices + old_indices

    return indices_dict
