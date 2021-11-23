# Copyright 2020 University of Groningen
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
from collections import defaultdict
from vermouth.parser_utils import SectionLineParser

class Mapping():

    def __init__(self, from_resname, to_resname):
        self.from_resname = from_resname
        self.to_resname = to_resname
        self.beads = set()
        self.n_beads = len(self.beads)
        self.bead_to_idx = defaultdict(list)
        self.bead_to_atom = defaultdict(list)

    def add_atom(self, bead, idx, atom=None):
        self.bead_to_atom[bead].append(atom)
        self.bead_to_idx[bead].append(idx)
        self.beads.add(bead)

class MapDirector(SectionLineParser):

    COMMENT_CHAR = ';'

    def __init__(self):
        super().__init__()
        self.current_mapping = None
        self.current_from = None
        self.mappings = {}
        self.header_actions = {
            ('molecule',): self._new_mapping
        }

    def _new_mapping(self):
        if self.current_mapping:
            self.mappings[self.current_from] = self.current_mapping
            self.current_from = None

        self.current_mapping = None

    @SectionLineParser.section_parser('molecule')
    def _molecules(self, line, lineno=0):
        """
        Parses the lines in the '[molecule]'
        directive and stores it. This should actually
        be residue but to comply with the bachwards
        format we keep molecule.
        """
        tokens = line.split()
        if len(tokens) == 2:
            from_resname, to_resname = tokens
        else:
            from_resname = to_resname = tokens

        self.current_mapping = Mapping(from_resname=from_resname,
                                       to_resname=to_resname)
        self.current_from = from_resname

    @SectionLineParser.section_parser('to')
    @SectionLineParser.section_parser('from')
    @SectionLineParser.section_parser('martini')
    @SectionLineParser.section_parser('mapping')
    @SectionLineParser.section_parser('chiral')
    def _skip_line(self, line, lineno=0):
        pass

    @SectionLineParser.section_parser('atoms')
    def _atoms(self, line, lineno=0):
        """
        Parse and store atomtypes section
        """
        idx, atom, bead = line.split()
        self.current_mapping.add_atom(idx=idx, atom=atom, bead=bead)

    def finalize(self, lineno=0):
        """
        Called at the end of the file
        """
        self._new_mapping()
        return self.mappings

def read_mapping(lines):
    """
    Parses `lines` of itp format and adds the
    molecule as a block to `force_field`.

    Parameters
    ----------
    lines: list
        list of lines of an itp file
    force_field: :class:`vermouth.forcefield.ForceField`
    """
    director = MapDirector()
    return list(director.parse(iter(lines)))
