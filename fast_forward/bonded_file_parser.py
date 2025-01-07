# Copyright 2024 University of Groningen
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
import re

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
                        if atomtype[0] == 'T':
                            mass = 36
                        elif atomtype[0] == 'S':
                            mass = 54
                        else:
                            mass = 72
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
