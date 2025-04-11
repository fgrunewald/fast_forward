
from vermouth.gmx.itp import write_molecule_itp
from vermouth.file_writer import deferred_open

def itp_writer(molname, interactions, block, command_used):
    '''
    Write an itp for a block with new interactions. All existing interactions will be removed and new ones
    written in their place

    Parameters
    ----------
    molname: str
        molname for block
    interactions: dict
        The dictionary of new interactions to overwrite old ones with. Formatted as {interaction_type: [Interaction]}
        where Interaction is fast_forward.interaction.Interaction.
        Should come from fast_forward.interaction_fit.InteractionFitter.interactions_dict.
    block: vermouth.molecule.Block
        The block for which to write the itp
    command_used: str
        Command used to run the program

    '''

    # remove all interactions present
    for interaction_type in list(block.interactions):
        del block.interactions[interaction_type]

    for interaction_type in interactions:
        for interaction in interactions[interaction_type]:
            if type(interaction) == list:
                for i in interaction:
                    block.add_interaction(type_=i.name,
                                          atoms=i.atoms,
                                          parameters=i.parameters,
                                          meta=i.meta
                                          )
            else:
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
