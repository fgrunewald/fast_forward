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

from vermouth.file_writer import deferred_open

def itp_writer(molname, block, interactions_dict, command_used):
    # rearrange the interactions for now until I work out how to do it properly
    defaults = {'bonds': 1, 'angles': 2, 'dihedrals': 1, 'constraints': 1,
                'virtual_sitesn': 1, 'virtual_sites3': 2}
    vermouth_interactions = {}
    for key in interactions_dict.keys():
        l = []
        for _, value in interactions_dict[key].items():
            b = [i for j in value[1] for i in j]

            # need this to add multiplicity default to dihedrals
            if key == 'dihedrals':
                c = [x for xs in [[defaults[key]], value[0], [1]] for x in xs]
            else:
                c = [x for xs in [[defaults[key]], value[0]] for x in xs]
            l.append((b, c))
        vermouth_interactions[key] = l

    # add the interactions
    for interaction_type in vermouth_interactions.keys():
        if interaction_type != 'constraints':
            for interaction in vermouth_interactions[interaction_type]:
                comment = '_'.join([block.nodes[i]['atomname'] for i in interaction[0]])

                block.remove_interaction(interaction_type, interaction[0])

                block.add_interaction(interaction_type, interaction[0], interaction[1], meta={"comment": comment})
        else:
            # add to both bonds and constraints for minimization purposes
            for interaction in vermouth_interactions[interaction_type]:
                comment = '_'.join([block.nodes[i]['atomname'] for i in interaction[0]])

                block.remove_interaction('bonds', interaction[0])

                block.add_interaction(interaction_type, interaction[0], interaction[1][:2],
                                    meta={"ifndef": "FLEXIBLE", "comment": comment})
                block.add_interaction('bonds', interaction[0], interaction[1][:2] + ['10000'],
                                    meta={"ifdef": "FLEXIBLE", "comment": comment})

    header = ['This file was generated using the following command:',
              command_used, '\n',
              'initial itp generation done by Fast-Forward. Please cite:',
              'https://zenodo.org/badge/latestdoi/327071500']

    # make the block a molecule for writing
    mol_out = block.to_molecule()
    mol_out.meta['molname'] = molname

    with deferred_open(f'{molname}.itp', 'w') as fout:
        write_molecule_itp(mol_out, fout, moltype=molname, header=header)
