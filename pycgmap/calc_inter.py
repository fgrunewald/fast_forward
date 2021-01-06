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
from tqdm import tqdm
from .interaction_class import InteractionUniverse
from .modf_boltzmann_inversion import modf_blotzmann_inversion
from .symfit_boltzmann_inversion import symfit_interactions

METHODS = {"modf_boltz": modf_blotzmann_inversion,
           "symfit_boltz": symfit_interactions}

KOWN_INTER_TYPES = ["bonds", "angles", "constraints"]

def compute_interaction_parameters(universe, molecule, mode, temp=298.15, gas_const=8.314):
    """
    Parameters:
    -----------
    universe: :class:`mda.base.core.Universe`
    molecule: :class:`vermouth.molecule.Molecule`
    const: float

    Returns:
    --------
    :class:`vermouth.molecule.Molecule`
        molecule with best fit parameters
    """
    # initiate class
    interactions = InteractionUniverse(universe, molecule)
    interactions.compute_pair_distances()

    # for each interaction compute parameters
    pbar = tqdm(total=interactions.n_interactions)
    for inter_type, inter, pair_dists, idx in interactions:
        if inter_type in KOWN_INTER_TYPES:
            inter = METHODS[mode](inter_type,
                                  inter,
                                  pair_dists,
                                  temp=temp,
                                  gas_const=gas_const)
            # update the interaction in the molecule
            interactions.molecule.interactions[inter_type][idx] = inter
        pbar.update(1)
    pbar.close()
    return molecule
