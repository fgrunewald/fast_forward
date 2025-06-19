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
    mapping = next(raw_matches)
    return mapping

def cgsmiles_to_mapping(cgsmiles_strs, mol_names, univ):
    """
    Convert cgsmiles str to mapping.
    """
    print(cgsmiles_strs)
    mapped_atoms, bead_idxs = [], []
    bead_count = 0
    for cgs, mol_name in zip(cgsmiles_strs, mol_names):
        print(cgs)
        resolver = MoleculeResolver.from_string(cgs, last_all_atom=True)
        cg, aa = resolver.resolve()
        mol_graph = univ.molecule_graphs[mol_name]

        try:
            _match = find_one_graph_match(mol_graph, aa)
        except StopIteration:
            msg = f"CGsmiles string {cgs} does not match molecule with name {mol_name}."
            raise SyntaxError(msg)

        for mol_idx in univ.mol_idxs_by_name[mol_name]:
            for bead in cg.nodes:
                atoms = cg.nodes[bead]['graph'].nodes
                mapped_atoms.append([_match[atom]+offset for atom in atoms])
                bead_idxs.append(bead_count)
                bead_count += 1
            offset += len(mol_graph)

    return mapped_atoms, bead_idxs
