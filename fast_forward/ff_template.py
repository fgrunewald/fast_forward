import collections
from vermouth.ffinput import FFDirector, _tokenize
from vermouth.parser_utils import SectionLineParser

class FFTemplate(FFDirector):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @SectionLineParser.section_parser('moleculetype', 'atoms')
    def _block_atoms(self, line, lineno=0):
        tokens = collections.deque(_tokenize(line))
        _parse_template_atom(tokens, self.current_block)


def _parse_template_atom(tokens, context):
    if tokens[-1].startswith('{'):
        attributes = _parse_atom_attributes(tokens.pop())
    else:
        attributes = {}

    mass = None
    charge = None
    # we use the vermouth default way
    if len(tokens) == 6:
        # deque does not support slicing
        first_six = (tokens.popleft() for _ in range(6))
        _, atype, resid, resname, name, charge_group = first_six

        # charge and mass are optional, but charge has to be defined for mass to be
        if tokens:
            charge = float(tokens.popleft())
        if tokens:
            mass = float(tokens.popleft())
    # we have the reduced fromat
    else:
        name = tokens.popleft()
        atype = tokens.popleft()
        resname = context.name
        resid = "1"
        charge_group = "1"

    if name in context:
        msg = ('There is already an atom named "{}" in the block "{}". '
               'Atom names must be unique within a block.')
        raise IOError(msg.format(name, context.name))
    atom = {
        'atomname': name,
        'atype': atype,
        'resname': resname,
        'resid': int(resid),
        'charge_group': int(charge_group),
    }
    if mass:
        atom['mass'] = mass
    if charge:
        atom['charge'] = charge

    context.add_atom(dict(collections.ChainMap(attributes, atom)))

def read_ff_template(lines, force_field):
    director = FFTemplate(force_field)
    return list(director.parse(iter(lines)))

