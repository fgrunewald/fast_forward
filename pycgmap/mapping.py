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
import numba
from numba import njit, prange, cuda
import MDAnalysis as mda

def _selector(atomgroup, indices, names):
    #idx_string = " ".join(indices)
    names_string = " ".join(names)
   # atoms = atomgroup.select_atoms("index " + idx_string + " and name " + names_string)
    atoms = atomgroup.select_atoms("name " + names_string)
    return atoms

def create_new_universe(universe, mapped_trajectory, mappings):
    """
    Create a new universe according to the definitions in mappings
    the old universe, and the mapped trajectory.

    Parameters
    -----------
    universe :class:`pycgmap.universe_handler.UniverseHandler`
    mapped_trajectory: np.ndarray
        coordiante array with shape (n_frames, n_atoms, 3)
    mappings: dict[:class:`pycgmap.map_parser.Mapping`]
        a dict of resname, mapping object

    Returns
    --------
    :class:`MDAnalysis.core.universe`
        the new universe that ties all information
    """
    # copy the dimensions array
    n_frames = mapped_trajectory.shape[0]
    dimensions = np.zeros((n_frames, 6))
    for fdx, ts in enumerate(universe.trajectory):
        dimensions[fdx, :] = ts.dimensions

    # create the universe to be returend
    # initalize some attributes
    n_residues = universe.n_residues
    n_atoms = mapped_trajectory.shape[1]
    res_seg = np.array([1] * n_residues)
    # to read out these we have to iterate over universe again
    atom_resindex = []
    atomnames = []
    resids = []
    resnames = []
    for idx, residue in enumerate(universe.res_iter()):
        for bead in mappings[residue.resname].beads:
            atomnames.append(bead)
            resnames.append(residue.resname)
            resids.append(idx)
            atom_resindex.append(idx)

    # now create the empty universe
    cg_universe = mda.Universe.empty(trajectory=True,
                                     n_atoms=n_atoms,
                                     n_residues=n_residues,
                                     atom_resindex=atom_resindex,
                                     residue_segindex=res_seg,
                                     )
    # add the coordinates
    cg_universe.trajectory.coordinate_array = mapped_trajectory
    cg_universe.trajectory.dimensions_array = dimensions
    cg_universe.trajectory.n_frames = n_frames

    # some more labels to make the universe sane
    cg_universe.add_TopologyAttr("names", values=atomnames)
    cg_universe.add_TopologyAttr("resnames", values=resnames)
    cg_universe.add_TopologyAttr("resids", values=resids)
    return cg_universe

def forward_map_indices(universe, mappings):
    """
    Map the universe atom indices to the new universe.
    """
    mapped_atoms = []
    bead_idxs = []
    total_beads = 0
    for residue in universe.res_iter():
        resid = residue.resid
        resname = residue.resname
        mapping = mappings[resname]
        for bead_count, bead in enumerate(mapping.beads):
            idxs = mapping.bead_to_idx[bead]
            names = mapping.bead_to_atom[bead]
            atoms = _selector(residue.atoms, idxs, names)
            mapped_atoms.append(atoms.indices)
            bead_idxs.append(total_beads)
            total_beads += 1
    return mapped_atoms, bead_idxs

@njit(parallel=True)
def forward_map_positions(mapped_atoms, bead_idxs, positions, n_frames, mode):
    new_trajectory = np.zeros((n_frames, len(mapped_atoms), 3))
    for count_lv1 in prange(len(mapped_atoms)):
        bead_idx = bead_idxs[count_lv1]
        atom_idxs = mapped_atoms[count_lv1]
        for fdx in prange(0, n_frames):
            new_pos = np.array([0.0, 0.0, 0.0], dtype=np.float32)
            for atom_idx in atom_idxs:
                vector = positions[fdx, atom_idx, :]
                new_pos = new_pos + vector
            new_pos = new_pos / len(atom_idxs)
            new_trajectory[fdx, bead_idx, :] = new_pos
    return new_trajectory
