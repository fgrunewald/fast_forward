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
from collections import OrderedDict, defaultdict
import numpy as np
from tqdm import tqdm
import networkx as nx
import MDAnalysis as mda
from MDAnalysis.core.universe import _TOPOLOGY_ATTRS
from vermouth.graph_utils import make_residue_graph
from .contributed import ndx_to_ag
from MDAnalysis import transformations

def parse_mapping_file(infile):
    """
    Read a residue mapping file and return
    a dict of dicts.
    """
    mapping = defaultdict(dict)
    for line in infile:
        # skip comments
        comment_idx = line.find(';')
        if comment_idx >= 0:
            line = line[comment_idx:]

        line = line.strip()
        # set the residue name
        if line.startswith('['):
            group_name = line[1:-1].strip()
        # update the mapping dict
        else:
            bead_mapping = line.split()
            # first bead is the bead name and the rest the AA atoms
            mapping[group_name][bead_mapping[0]] = bead_mapping[1:]

    return mapping

def res_to_ag(universe, mapping):
    """
    Convert a mapping dict to atom_groups
    """
    for residue in universe.residues:
        res_map = mapping[residue.resname]
        for bead, atom_names in res_map.items():
            res_atom_ids = [ np.where(residue.atoms.names == name)[0][0] for name in atom_names]
            atom_ids = residue.atoms.indices[res_atom_ids]
            yield universe.atoms[atom_ids]

def count_nodes(composition, molecules):
    total = 0
    for mol_name, n_mols in composition:
        total += n_mols * len(molecules[mol_name].nodes)
    return total

def make_residue_graphs(molecules):
    """
    Make residues for each molecule in molecules.

    Parameters:
    -----------
    molecules: dict[vermouth.molecule.Molecule]

    Returns:
    --------
    dict[nx.graph]
        dict of molename and residue graph
    """
    res_graphs = {}
    for mol_name, mol in molecules.items():
        res_graph = make_residue_graph(mol)
        res_graphs[mol_name] = res_graph
    return res_graphs

def composition_attribute_iter(composition, molecules, attribute):
    for molname, nmols in composition:
        molecule = molecules[molname]
        nodes = molecule.nodes
        sorted(nodes)
        for _ in range(nmols):
            for node in nodes:
                yield molecule.nodes[node][attribute]

def create_mda_universe_from_itps(composition, molecules):
    """
    Take a `molecule` and generate an :class:`MDAnalysis.core.universe`
    from it, setting all relevant topology attribute
    stored in molecule. The universe is initalized with a trajectory
    but no coordinates are added even if they are in the molecule.

    Parameters:
    -----------
    composition: tuple(int, str)
        how many molecules, molecule name
    molecules: dict[:class:`vermouth.molecule.Molecule`]
        dict of molecules by molname

    Returns:
    --------
    :class:`MDAnalysis.core.universe`
    """
    n_atoms = count_nodes(composition, molecules)
    print(n_atoms)
    res_graphs = make_residue_graphs(molecules)
    n_residues = count_nodes(composition, res_graphs)
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

def load_n_frames(filenames):
    """
    Generate an mda universe from serveral filenames.
    """
    u = mda.Universe(filenames[0])
    n_frames = len(filenames)
    n_atoms = u.atoms.n_atoms
    n_residues = len(u.atoms.resids)
    atom_resindex = np.array([atom.resindex for atom in u.atoms])
    res_seg = np.array([0 for _ in range(0, n_atoms)])

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
        u = mda.Universe(name)
        positions[fdx+1, :, :] = u.atoms.positions
        dimensions[fdx+1, :] = u.atoms.dimensions
        pbar.update(1)
    pbar.close()
    new_universe.trajectory.coordinate_array = positions.astype(np.float32)
    nS128ew_universe.trajectory.dimensions_array = dimensions
    new_universe.trajectory.n_frames = n_frames

    return new_universe


def _center_of_geometry(atomgroup):
    return atomgroup.center_of_geometry()

def _center_of_mass(atomgroup):
    return atomgroup.center_of_mass()

MAPPING_MODES = {"COG": _center_of_geometry,
                 "COM": _center_of_mass}

def mapping_transformation(universe, cg_universe, mapping, mode="COG"):
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

    pbar = tqdm(total=len(mapping)*len(universe.trajectory))
    for cg_atom, atom_group in mapping.items():
        fdx = 0
        for time_step in universe.trajectory:
            pos = MAPPING_MODES[mode](atom_group)
            new_trajectory[fdx, cg_atom, :] = pos
            fdx += 1
            pbar.update(1)
    pbar.close()

    dimensions = np.zeros((n_frames, 6))
    for fdx, ts in enumerate(universe.trajectory):
        dimensions[fdx, :] = ts.dimensions

    cg_universe.trajectory.coordinate_array = new_trajectory.astype(np.float32)
    cg_universe.trajectory.dimensions_array = dimensions
    cg_universe.trajectory.n_frames = n_frames
    return cg_universe

def establish_mapping(universe, ndx_file=None, res_file=None):
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
        nodes = range(0, len(atom_iter))
        mapping = OrderedDict(zip(nodes, atom_iter))
    elif res_file:
        with open(res_file) as _file:
            lines = _file.readlines()
        res_mappings = parse_mapping_file(lines)
        atom_iter = list(res_to_ag(universe, res_mappings))
        nodes = range(0, len(atom_iter))
        mapping = OrderedDict(zip(nodes, atom_iter))
    else:
        raise IOError("Index file or residue index file needs to be specified.")

    return mapping
