
import networkx as nx
from vermouth.gmx.itp import write_molecule_itp
from vermouth import Molecule

def generate_template_itp(residue_dict, links_dict, molname):
    # generate the molecule
    atomcount = 0
    rescount = 1
    nodes_list = []
    for residue in residue_dict.keys():
        atoms = residue_dict[residue]['atomtypes']
        for atom in atoms.keys():
            node_data = [atomcount, {'atomname': atom,
                                     'atype': atoms[atom]['atype'],
                                     'charge': atoms[atom]['charge'],
                                     'mass': atoms[atom]['mass'],
                                     'resid': rescount,
                                     'resname': residue,
                                     'charge_group': atomcount + 1}]
            atomcount += 1

            nodes_list.append(node_data)
        rescount += 1

    G = nx.Graph()
    G.add_nodes_from(nodes_list)

    mol = Molecule(G)
    mol.nrexcl = 1

    default_parameters = {'bonds': [1, 1, 1000],
                          'angles': [2, 100, 200],
                          'dihedrals': [1, 90, 10, 1],
                          'virtual_sitesn': [1]}

    # sort out intra residue interactions
    for residue in residue_dict.keys():
        interactions = residue_dict[residue]['interactions']
        for interaction in interactions:
            for atomset in interactions[interaction]:
                atoms = []
                for i in atomset:
                    for node in mol.nodes:
                        if (mol.nodes[node]['atomname'] == i) and (mol.nodes[node]['resname'] == residue):
                            atoms.append(node)
                mol.add_interaction(interaction,
                                    atoms,
                                    default_parameters[interaction],
                                    meta={'comment': '_'.join(atomset)})

    if links_dict is not None:
        # sort out the inter residue interactions
        for interaction in links_dict.keys():
            for atomset in links_dict[interaction]:
                atoms = []
                for i in atomset:
                    resname, atomname = i.split(':')
                    for node in mol.nodes:
                        if (mol.nodes[node]['atomname'] == atomname) and (mol.nodes[node]['resname'] == resname):
                            atoms.append(node)

                mol.add_interaction(interaction,
                                    atoms,
                                    default_parameters[interaction],
                                    meta={'comment': '_'.join(atomset)})


    header = ['initial itp generation done by Fast-Forward. Please cite:',
              'https://zenodo.org/badge/latestdoi/327071500']

    # write the molecule out
    with open(f'{molname}.itp', 'w') as fout:
        write_molecule_itp(mol, fout, moltype=f'{molname}', header=header)


