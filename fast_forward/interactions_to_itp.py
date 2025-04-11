
from vermouth.gmx.itp import write_molecule_itp
from vermouth.file_writer import deferred_open

def itp_writer(molname, block, command_used):
    '''
    Write an itp for a block with new interactions. All existing interactions will be removed and new ones
    written in their place

    Parameters
    ----------
    molname: str
        molname for block
    block: vermouth.molecule.Block
        The block for which to write the itp
    command_used: str
        Command used to run the program

    '''

    header = ['This file was generated using the following command:',
              command_used, '\n',
              'itp generation done by Fast-Forward. Please cite:',
              'https://zenodo.org/badge/latestdoi/327071500']

    # make the block a molecule for writing
    mol_out = block.to_molecule()
    mol_out.meta['molname'] = molname

    with deferred_open(f'{molname}.itp', 'w') as fout:
        write_molecule_itp(mol_out, fout, moltype=molname, header=header)
