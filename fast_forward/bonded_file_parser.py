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
    for line in lines:
        if line.startswith('['):
            # we do this after a section is completed unless we start
            if group_name is not None:
                interactions[group_name] = atoms
            group_name = re.search(r'[a-z]+', line).group(0)
            atoms = []
        else:
            if len(line.split()) > 0:
                atoms.append(line.split())
    # do the last one
    interactions[group_name] = atoms

    return interactions
