# Copyright 2024 University of Groningen
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

from vermouth.gmx.itp import write_molecule_itp
from vermouth import Molecule
import numpy as np
import networkx as nx

def itp_writer(interactions_dict, atomtypes, molname):
    # rearrange the interactions for now until I work out how to do it properly
    defaults = {'bonds': 1, 'angles': 2, 'dihedrals': 1, 'constraints': 1,
                'virtual_sitesn': 1, 'virtual_sites3': 2}
    vermouth_interactions = {}
    for key in interactions_dict.keys():
        l = []
        for _, value in interactions_dict[key].items():
            b = value[1]
            if key == 'dihedrals':
                c = [x for xs in [[defaults[key]], value[0], [1]] for x in xs]
            else:
                c = [x for xs in [[defaults[key]], value[0]] for x in xs]
            l.append((b, c))
        vermouth_interactions[key] = l

    # get information about the atomnames
    atomnames = {}
    for key in interactions_dict.keys():
        for k, v in interactions_dict[key].items():
            names = k.split('_')
            indices = v[1]
            for i, j in zip(indices, names):
                atomnames[i] = j

    nodes = np.arange(len(atomnames))

    # make a list of nodes to construct a molecule from
    nodes_list = []
    for node in nodes:
        if node in list(atomnames.keys()):
            l = [node, {'atomname': atomnames[node],
                        'atype': atomtypes[atomnames[node]]['atype'],
                        'charge': atomtypes[atomnames[node]]['charge'],
                        'mass': atomtypes[atomnames[node]]['mass'],
                        'resid': '1',
                        'resname': molname,
                        'charge_group': node + 1}]
        else:
            l = [node, {'atomname': "TMP",
                        'atype': "TMP",
                        "charge": '0.0',
                        'resid': '1',
                        'resname': molname,
                        'charge_group': node + 1}]

        nodes_list.append(l)

    exclusions = []
    for index, node in enumerate(nodes[:-1]):
        exclusions.append([tuple(nodes[index:]), []])
    vermouth_interactions['exclusions'] = exclusions

    # make the molecule from a graph
    G = nx.Graph()
    G.add_nodes_from(nodes_list)

    mol = Molecule(G)
    mol.nrexcl = 1

    # add the interactions
    for interaction_type in vermouth_interactions.keys():
        if interaction_type != 'constraints':
            for interaction in vermouth_interactions[interaction_type]:
                mol.add_interaction(interaction_type, interaction[0], interaction[1])
        else:
            # add to both bonds and constraints for minimization purposes
            for interaction in vermouth_interactions[interaction_type]:
                mol.add_interaction(interaction_type, interaction[0], interaction[1][:2],
                                    meta={"ifndef": "FLEXIBLE"})
                mol.add_interaction('bonds', interaction[0], interaction[1][:2] + ['10000'],
                                    meta={"ifdef": "FLEXIBLE"})

    # write the molecule out
    with open(f'{molname}.itp', 'w') as f:
        write_molecule_itp(mol, f, moltype=molname)

