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
import networkx as nx
import MDAnalysis as mda
from MDAnalysis.core.universe import _TOPOLOGY_ATTRS
from vermouth.graph_utils import make_residue_graph
from .contributed import ndx_to_ag
from MDAnalysis import transformations

def create_mda_universe_from_itp(molecule):
    """
    Take a `molecule` and generate an :class:`MDAnalysis.core.universe`
    from it, setting all relevant topology attribute
    stored in molecule. The universe is initalized with a trajectory
    but no coordinates are added even if they are in the molecule.

    Parameters:
    -----------
    molecule: :class:`vermouth.molecule.Molecule`

    Returns:
    --------
    :class:`MDAnalysis.core.universe`
    """
    n_atoms = len(molecule.nodes)
    res_graph = make_residue_graph(molecule)
    node_to_attr = nx.get_node_attributes(molecule, "resid")
    atom_resindex = np.array([node_to_attr[node]-1 for node in sorted(node_to_attr.keys())])

    cg_universe = mda.Universe.empty(trajectory=True,
                                     n_atoms=n_atoms,
                                     n_residues=len(res_graph.nodes),
                                     atom_resindex=atom_resindex
                                     )

    # assign atom based attributes
    for attr_mda, attr_mol in {"names": "atomname", "types": "atomtype"}.items():
        node_to_attr = nx.get_node_attributes(molecule, attr_mol)
        if not node_to_attr:
            continue
        values = np.array([node_to_attr[node] for node in sorted(node_to_attr.keys())])
        cg_universe.add_TopologyAttr(attr_mda, values=values)

    # assign residue based attributes
    for attr_mda, attr_mol in {"resnames": "resname", "resids": "resid"}.items():
        node_to_attr = nx.get_node_attributes(res_graph, attr_mol)
        if not node_to_attr:
            continue
        values = np.array([node_to_attr[node] for node in sorted(node_to_attr.keys())])
        cg_universe.add_TopologyAttr(attr_mda, values=values)

    return cg_universe

def _center_of_geometry(atomgroup):
    return atomgroup.center_of_geometry()

def _center_of_mass(atomgroup):
    return atomgroup.center_of_mass()

MAPPING_MODES = {"COG": _center_of_geometry,
                 "COM": _center_of_mass}

def mapping_transformation(universe, molecule, cg_universe, mapping, mode="COG"):
    """
    Take a `universe` and `cg_universe` and generate the positions
    of the cg_universe by mapping the coodinates from universe
    according to the correspondance in mapping. The mode of mapping
    can be center-of-geometry (COG) or center-of-mass (COM).

    Parameters:
    -----------
    universe: :class:`MDAnalysis.core.universe`
        the atomistic universe to be mapped
    cg_universe: :class:`MDAnalysis.core.universe`
        empty cg_universe with trajectory initialized
    mapping: dict[int]
        dict mapping atoms in cg_universe to atomgroups
        in universe
    mode: str
        mapping mode; either COG or COM

    Returns:
    --------
    :class:`MDAnalysis.core.universe`
        cg_universe with trajectory of mapped positions
    """
    n_frames = len(universe.trajectory)
    n_atoms = len(cg_universe.atoms)
    new_trajectory = np.zeros((n_frames, n_atoms, 3))

    for cg_atom, atom_group in mapping.items():
        fdx = 0
        for time_step in universe.trajectory:
            pos = MAPPING_MODES[mode](atom_group)
            new_trajectory[fdx, cg_atom, :] = pos
            fdx += 1

    dimensions = np.zeros((n_frames, 6))
    for fdx, ts in enumerate(universe.trajectory):
        dimensions[fdx, :] = ts.dimensions

    cg_universe.trajectory.coordinate_array = new_trajectory.astype(np.float32)
    cg_universe.trajectory.dimensions_array = dimensions
    cg_universe.trajectory.n_frames = n_frames

    return cg_universe

def establish_mapping(universe, molecule, ndx_file=None, res_file=None):
    """
    Given a `universe` and `molecule` use the definitions of beads to
    atoms provided by an index file (`ndx_file`) or a residue mapping
    file `res_file`, establish which node in the `molecule` corresponds
    to which atomgroup part of the universe.

    Parameters:
    -----------
    universe: :class:`MDAnalysis.core.universe`
        the atomistic universe to be mapped
    molecule: :class:`vermouth.molecule.Molecule`
    ndx_file: :class:`pathlib.Path`
    red_file: :class:`pathlib.Path`

    Returns:
    --------
    dict
    """
    if ndx_file:
        with open(ndx_file) as _file:
            lines = _file.readlines()
        atom_iter = ndx_to_ag(universe, lines)
    else:
        raise IOError("Index file or residue index file needs to be specified.")

    mapping = dict(zip(molecule.nodes, atom_iter))
    return mapping