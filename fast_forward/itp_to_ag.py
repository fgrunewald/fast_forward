"""
Extract interactions pairs from itp and convert to atom
groups of an MDAnalysis Universe.
"""
from MDAnalysis.core.topologyattrs import Moltypes, Molnums
from collections import defaultdict
import numpy as np
import networkx as nx

def find_mol_indices(universe, atoms, moltypes):
    """
    Given a universe select all atoms that belong to molecules 
    with moltype matching one of the entries in `moltypes`.
    Subsequently, return indices of all multiples of the
    indices defined in atoms under the assumption that all
    molecules referenced in `moltypes` have the same number of
    atoms.

    Parameters
    ----------
    universe: mda.Universe
    atoms: list[int]
    moltypes: list[abc.hashable]

    Returns
    -------
    list[numpy.array(dtype=int)]
    """
    mol_atoms = universe.select_atoms('moltype ' + ' '.join(moltypes))
    n_mols = len(np.unique(mol_atoms.molnums))
    try:
        mol_atom_indices = mol_atoms.indices.reshape(n_mols, -1)
    except ValueError:
        msg = ("The target molecules passed to find_mol_indices "
               "do not seem to all have the same number of atoms.")
        raise IndexError(msg) from None
    return list(mol_atom_indices[:, atoms])

def res_as_mol(universe):
    """
    For a universe without moltype/molnum info, promotes residues to molecules.

    Changes universe in place. Does nothing if moltype/molnum info is already
    available.
    """
    if hasattr(universe.atoms, "moltypes"):
        return

    moltypes = Moltypes(universe.residues.resnames)
    molnums = Molnums(range(len(universe.residues)))
    universe.add_TopologyAttr(moltypes)
    universe.add_TopologyAttr(molnums)

def itp_to_ag(block, mol_name, universe):
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
                indices = find_mol_indices(universe,
                                           atoms,
                                           [mol_name])
                old_indices = indices_dict[inter_type].get(group, [])
                old_block_indices = block_indices[inter_type].get(group, [])

                indices_dict[inter_type][group] = indices + old_indices
                initial_parameters[inter_type][group] = inter.parameters

                block_indices[inter_type][group] = [atoms] + old_block_indices

    return indices_dict, initial_parameters, block_indices
