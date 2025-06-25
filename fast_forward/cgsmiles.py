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
import networkx as nx
from cgsmiles.resolve import MoleculeResolver
from .map_file_parers import Mapping

def load_cgsmiles_library(filepath):
    cgsmiles_strs = []
    with open(filepath) as _file:
        for line in _file.readlines():
            cgsmiles_strs.append(line.strip())
    return cgsmiles_strs

def find_one_graph_match(graph1, graph2):
    """
    Returns one ismags match when graphs are isomorphic
    otherwise None.
    """

    def node_match(n1, n2):
        return n1["element"] == n2['element']

    def edge_match(e1, e2):
        return e1['order'] == e2['order']

    GM = nx.isomorphism.GraphMatcher(graph1,
                                     graph2,
                                     node_match=node_match,
                                     edge_match=edge_match)

    raw_matches = GM.subgraph_isomorphisms_iter()
    try:
        mapping = next(raw_matches)
    except StopIteration:
        mapping = []
    return mapping

def cgsmiles_to_mapping(univ, cgsmiles_strs, mol_names, mol_matching=True):
    """
    Given a list of mappings described by cgsmiles strings
    this function maps the atom indices to CG beads using,
    graph isomorphism.

    Parameters
    ----------
    univ: :class:`fast_forward.UniverseHandler`
    cgsmiles_strs: list[str]
    mol_names: list[str]


    Retunrs
    -------
    list, list, dict
    """
    mapped_atoms, bead_idxs, mappings = [], [], {}
    bead_count = 0
    for idx, mol_name in enumerate(mol_names):
        mol_graph = univ.molecule_graphs[mol_name]
        # names match the order of cgs strings
        if mol_matching:
            possible_cgs = [cgsmiles_strs[idx]]
        else:
            possible_cgs = cgsmiles_strs

        for cgs_str in possible_cgs:
            # read the cgsmiles string
            resolver = MoleculeResolver.from_string(cgs_str, last_all_atom=True)
            cg, aa = resolver.resolve()
            _match = find_one_graph_match(aa, mol_graph)
            if _match:
                break
        else:
            raise SyntaxError("No matching cgsmiles string found.")

        offset=0
        for mol_idx in univ.mol_idxs_by_name[mol_name]:
            for bead in cg.nodes:
                atoms = cg.nodes[bead]['graph'].nodes
                mapped_atoms.append([_match[atom]+offset for atom in atoms])
                bead_idxs.append(bead_count)
                bead_count += 1
                # we only need to do this once per molecule type
                if mol_idx == 0:
                    resnames = [univ.atoms[_match[atom]].resname for atom in atoms]
                    resname = resnames[0]
                    if resname not in mappings:
                        mappings[resname] = Mapping(resname, resname)
                    mapping = mappings[resname]
                    for atom in atoms:
                        mapping
                        mapping.add_atom(cg.nodes[bead]["fragname"]+f"{bead}",
                                         _match[atom],
                                         atom=univ.atoms[_match[atom]].name)

            offset += len(mol_graph)
    return mapped_atoms, bead_idxs, mappings
