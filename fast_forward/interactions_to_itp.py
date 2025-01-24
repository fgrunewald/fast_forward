
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
    for interaction_type in interactions_dict.keys():
        finalised_interactions = []
        for _, value in interactions_dict[interaction_type].items():
            interaction_params, interaction_atoms = value
            b = [i for j in interaction_atoms for i in j]

            # need this to add multiplicity default to dihedrals
            if interaction_type == 'dihedrals':
                c = [x for xs in [[defaults[interaction_type]], interaction_params, [1]] for x in xs]
            else:
                c = [x for xs in [[defaults[interaction_type]], interaction_params] for x in xs]
            finalised_interactions.append((b, c))
        vermouth_interactions[interaction_type] = finalised_interactions

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
