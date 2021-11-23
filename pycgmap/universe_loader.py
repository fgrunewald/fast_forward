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
from collections import OrderedDict
import numpy as np
import networkx as nx
from tqdm import tqdm
import MDAnalysis as mda
from MDAnalysis.core.universe import _TOPOLOGY_ATTRS
from vermouth.graph_utils import make_residue_graph
from .contributed import ndx_to_ag
from MDAnalysis import transformations

def load_n_frames(filenames):
    """
    Generate an mda universe from serveral filenames.
    """
    u = mda.Universe(filenames[0])
    n_frames = len(filenames)
    n_atoms = u.atoms.n_atoms
    n_residues = len(u.atoms.resids)
    atom_resindex = np.array([ atom.resindex for atom in u.atoms])
    res_seg = np.array([ 0 for _ in range(0, n_atoms)])

    new_universe = mda.Universe.empty(trajectory=True,
                                      n_atoms=n_atoms,
                                      n_residues=n_residues,
                                      atom_resindex=atom_resindex,
                                      residue_segindex=res_seg,
                                     )

    for attr in ["names", "types", "resnames", "resids"]:
        values = getattr(u.atoms, attr)
        new_universe.add_TopologyAttr(attr, values=values)

    positions = np.zeros((n_frames, n_atoms, 3))
    dimensions = np.zeros((n_frames, 6))
    positions[0, :, :] = u.atoms.positions
    dimensions[0, :] = u.atoms.dimensions
    pbar = tqdm(total=n_frames-1)
    for fdx, name in enumerate(filenames[1:]):
        try:
            u = mda.Universe(name)
            positions[fdx+1, :, :] = u.atoms.positions
            dimensions[fdx+1, :] = u.atoms.dimensions
        except ValueError:
            raise IOError(name)
        pbar.update(1)
    pbar.close()
    new_universe.trajectory.coordinate_array = positions.astype(np.float32)
    new_universe.trajectory.dimensions_array = dimensions
    new_universe.trajectory.n_frames = n_frames

    return new_universe

def create_empty_mda_universe_from_topology(molecules,
                                            mol_idxs_by_name={}):
    """
    Given a list of molecules create an empty MDAnalysis universe,
    that has the number of molecules by name as defined in the dict
    mol_idxs_by_name.

    Parameters
    -----------
    topology: :class:`pycmap.topology.Topology`

    Returns
    --------
    :class:`MDAnalysis.core.universe`
    """
    # count how many atoms/particles there are in the system
    n_atoms = count_atoms(molecules, mol_idxs_by_name)
    # make a residue graph of each molecule
    res_graphs = make_residue_graphs(molecules)
    # count how many residues there are
    n_residues = count_atoms(res_graphs)
    # initialize some MDA attributes
    res_seg = np.array([1] * n_residues)
    atom_resindex = np.fromiter(composition_attribute_iter(composition, molecules, "resid"), dtype=int) - 1

    cg_universe = mda.Universe.empty(trajectory=True,
                                     n_atoms=n_atoms,
                                     n_residues=n_residues,
                                     atom_resindex=atom_resindex,
                                     residue_segindex=res_seg,
                                     )

    # assign atom based attributes
    for attr_mda, attr_mol in {"names": "atomname", "types": "atype"}.items():
        values = np.fromiter(composition_attribute_iter(composition, molecules, attr_mol), dtype='S128').astype(str)
        cg_universe.add_TopologyAttr(attr_mda, values=values)

    # assign residue based attributes
    for attr_mda, attr_mol in {"resnames": "resname", "resids": "resid"}.items():
        values = np.fromiter(composition_attribute_iter(composition, res_graphs, attr_mol), dtype='S128').astype(str)
        cg_universe.add_TopologyAttr(attr_mda, values=values)

    return cg_universe
