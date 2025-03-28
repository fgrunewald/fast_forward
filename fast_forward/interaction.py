

class Interaction():
    def __init__(self, name, func_type,
                 location,
                 meta = {},
                 force_constant=None,
                 atoms=None, multiplicity=None):
        self.name = name
        self.atoms = atoms
        self.func_type = func_type
        self.location = location
        self.force_constant = force_constant
        self.multiplicity = int(multiplicity) if multiplicity else multiplicity
        self.meta = meta
        # put parameters together coherently depending on what we have
        self.parameters = self.make_parameters()

    def make_parameters(self):
        return [i for i in [self.func_type, self.location, self.force_constant, self.multiplicity] if i is not None]


class InteractionHolder():
    def __init__(self):
        self.interactions = []

