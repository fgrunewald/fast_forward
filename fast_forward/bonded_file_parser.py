
import re

def read_interactions(lines):
    interactions = {}
    atoms = []
    group_name = None
    atomtypes = {}
    for line in lines:
        if line.strip().startswith(';'):
            continue
        if line.startswith('['):
            # we do this after a section is completed unless we start
            if group_name is not None:
                if group_name != 'atomtypes':
                    interactions[group_name] = atoms
                else:
                    for atom in atoms:
                        atomname = atom[0]
                        atomtype = atom[1]

                        mass_dict = {'T': 36, 'S': 54}

                        mass = mass_dict.get(atomtype[0], 72)

                        if len(atom) == 3:
                            charge = atom[2]
                        else:
                            charge = 0
                        atomtypes[atomname] = {'atype': atomtype, 'charge': charge, 'mass': mass}

            group_name = re.search(r'[a-z]+_?[a-z1-9]+', line).group(0)
            atoms = []
        else:
            if len(line.split()) > 0:
                atoms.append(line.split())
    # do the last one
    interactions[group_name] = atoms

    return {'interactions': interactions, 'atomtypes': atomtypes}


def read_links(lines):
    interactions = {}
    atoms = []
    group_name = None
    for line in lines:
        if line.strip().startswith(';'):
            continue
        if line.startswith('['):
            # we do this after a section is completed unless we start
            if group_name is not None:
                interactions[group_name] = atoms

            group_name = re.search(r'[a-z]+_?[a-z1-9]+', line).group(0)
            atoms = []
        else:
            if len(line.split()) > 0:
                atoms.append(line.split())
    # do the last one
    interactions[group_name] = atoms

    return interactions

def read_template(fin):

    with open(fin) as file:
        lines = file.readlines()

    # work out where the different residue sections are
    sections = [index for index, line in enumerate(lines) if 'moleculetype' in line or 'links' in line]

    section_pairs = [[sections[i], sections[i+1]] for i in range(len(sections)-1)] + [[sections[-1],len(lines)]]

    # read in the initial data
    residue_dict = {}
    links_dict = None
    for i in section_pairs:
        if 'moleculetype' in lines[i[0]]:
            residue_dict[lines[i[0]+1].strip()] = read_interactions(lines[i[0]+1:i[1]])
        else:
            links_dict = read_links(lines[i[0]+1:i[1]])

    return residue_dict, links_dict

def read_bonds(infile):
    with open(infile) as _file:
        lines = _file.readlines()

    interactions = {}
    atoms = []
    group_name = None
    atomtypes = {}
    for line in lines:
        if line.strip().startswith(';'):
            continue
        if line.startswith('['):
            # we do this after a section is completed unless we start
            if group_name is not None:
                if group_name != 'atomtypes':
                    interactions[group_name] = atoms
                else:
                    for atom in atoms:
                        atomname = atom[0]
                        atomtype = atom[1]

                        mass_dict = {'T': 36, 'S': 54}

                        mass = mass_dict.get(atomtype[0], 72)

                        if len(atom) == 3:
                            charge = atom[2]
                        else:
                            charge = 0
                        atomtypes[atomname] = {'atype': atomtype, 'charge': charge, 'mass': mass}

            group_name = re.search(r'[a-z]+_?[a-z1-9]+', line).group(0)
            atoms = []
        else:
            if len(line.split()) > 0:
                atoms.append(line.split())
    # do the last one
    interactions[group_name] = atoms

    return interactions, atomtypes
