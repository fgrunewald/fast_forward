

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
        self.multiplicity = float(multiplicity) if multiplicity else multiplicity
        self.meta = meta

        self.parameters = self.make_parameters()

    def make_parameters(self):
        return [i for i in [self.func_type, self.location, self.force_constant, self.multiplicity] if i is not None]


class InteractionHolder():
    def __init__(self):
        self.interactions = []

    # def finalise_interactions_list(self):
        # flattened_interactions_list = [i for j in self.interactions for i in j]
        # self.interactions = flattened_interactions_list
        # for i in self.interactions:
        #     print(i)
