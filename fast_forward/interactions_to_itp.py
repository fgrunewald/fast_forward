
from vermouth.gmx.itp import write_molecule_itp
from vermouth.file_writer import deferred_open

def itp_writer(molname, block, interactions, command_used):

    # remove all interactions present
    for interaction_type in list(block.interactions):
        del block.interactions[interaction_type]

    for interaction in interactions.interactions:
        block.add_interaction(type_=interaction.name,
                              atoms=interaction.atoms,
                              parameters=interaction.parameters,
                              meta=interaction.meta
                              )

    header = ['This file was generated using the following command:',
              command_used, '\n',
              'itp generation done by Fast-Forward. Please cite:',
              'https://zenodo.org/badge/latestdoi/327071500']

    # make the block a molecule for writing
    mol_out = block.to_molecule()
    mol_out.meta['molname'] = molname

    with deferred_open(f'{molname}.itp', 'w') as fout:
        write_molecule_itp(mol_out, fout, moltype=molname, header=header)
