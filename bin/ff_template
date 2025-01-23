#!/usr/bin/env python3

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
"""
Main exe for generating template itps
"""
import argparse
import pathlib

from fast_forward.bonded_file_parser import read_template
from fast_forward.template_itp import generate_template_itp


"""
generate a template itp from a sketch file that looks like e.g.


[ moleculetype ]
MOL1

[ atomtypes ]
AT1 SP1
AT2 SP2
AT3 SP3
AT4 SP4

[ bonds ]
AT1 AT2
AT2 AT3
AT3 AT4

[ angles ]
AT1 AT2 AT3
AT2 AT3 AT4

[ dihedrals ]
AT1 AT2 AT3 AT4

[ moleculetype ]
MOL2

[ atomtypes ]
AT1 SP1
AT2 SP2
AT3 SP3
AT4 SP4

[ bonds ]
AT1 AT2
AT2 AT3
AT3 AT4

[ angles ]
AT1 AT2 AT3
AT2 AT3 AT4

[ dihedrals ]
AT1 AT2 AT3 AT4

[ links ]
; list as above in a molname:atomname format
[ bonds ]
MOL1:AT4 MOL2:AT1

[ angles ]
MOL1:AT3 MOL1:AT4 MOL2:AT1
MOL1:AT4 MOL2:AT1 MOL2:AT2

[ dihedrals ]
MOL1:AT2 MOL1:AT3 MOL1:AT4 MOL2:AT1
MOL1:AT3 MOL1:AT4 MOL2:AT1 MOL2:AT2
MOL1:AT4 MOL2:AT1 MOL2:AT2 MOL2:AT3

"""
def __main__():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
    parser.add_argument('-sketch', type=pathlib.Path, dest="input_file", help="input file")
    parser.add_argument('-name', type=pathlib.Path, dest="molname", help="name of molecule")

    args = parser.parse_args()

    residues, links = read_template(args.input_file)

    generate_template_itp(residues, links, args.molname)




__main__()
