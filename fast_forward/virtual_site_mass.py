
MARTINI_MASS=  {"S": 54, "T": 36}

def mass_redistribution(block, pair_idx):
    """
    Redistribute the mass of a virtual site across
    Parameters
    ----------
    block
    pair_idx

    Returns
    -------

    """
    for pairs in pair_idx:
        vs = pairs[0]
        block.nodes[vs]['mass'] = 0.0
        vs_type_mass = MARTINI_MASS.get(block.nodes[vs]["atype"][0], 72)
        redistribution = vs_type_mass / len(pairs[1:])

        # if only 1 constructor (i.e. building vs directly on top of atom), no redistribution necessary.
        if len(pairs[1:]) > 1:
            constructors = pairs[1:]
            for node in constructors:
                current_mass = block.nodes[node].get("mass", MARTINI_MASS.get(block.nodes[node]["atype"][0], 72))
                block.nodes[node]["mass"] = current_mass + redistribution
