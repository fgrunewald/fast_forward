import numpy as np

def ndx_to_ag(atoms_or_universe, infile):
    """
    Create MDAnalysis AtomGroups from a Universe and an index file.

    Examples
    --------
    .. code-block::
    u = mda.Universe(TPR, XTC)
    with open('index.ndx') as infile:
        for group_name, atom_group in ndx_to_ag(u, infile):
            print(group_name, atom_group)

    Parameters
    ----------
    atoms_or_universe: mda.AtomGroup or mda.Universe
        The atoms to select from.
    infile: file
        The open index file to read.

    Yields
    ------
    group_name: str
        The name of the group as written in the index file.
    atoms: mda.AtomGroup
        An atom group containing the atoms for that group.
    """
    atoms = atoms_or_universe.atoms
    group_name = None
    indices = []
    for line in infile:
        comment_idx = line.find(';')
        if comment_idx >= 0:
            line = line[comment_idx:]
        line = line.strip()
        if line.startswith('['):
            if group_name is not None:
                indices = np.array(indices, dtype=int) - 1
                yield atoms[indices] #(group_name, atoms[indices])
            group_name = line[1:-1].strip()
            indices = []
        else:
            indices += line.split()
    if group_name is not None:
        indices = np.array(indices, dtype=int) - 1
        #yield (group_name, atoms[indices])
        yield atoms[indices]
