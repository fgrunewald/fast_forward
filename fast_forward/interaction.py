
"""
Classes to store gromacs interactions to use for plotting, writing, etc
"""
class Interaction():
    """
    A single interaction
    """
    def __init__(self, name,
                 atoms,
                 parameters=[],
                 meta = {},
                 fit_data=None):
        '''
        Parameters
        ----------
        name: str
            name of interaction
        atoms: list
            atom indices in the molecule for the interaction to apply to
        parameters: list
            parameter list for this interaction
        meta: dict
            dictionary of meta information for vermouth associated with the interaction
        fit_data: list
            data underlying the parameters used during the fit
        '''
        self.name = name
        self.atoms = atoms
        self.parameters = parameters
        self.meta = meta
        self.fit_data = fit_data

