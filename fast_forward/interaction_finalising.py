"""
Manage interaction blocks that have and haven't been commented in an input itp file
"""

from collections import defaultdict
from copy import deepcopy

def finalise_interaction_types(block_interactions, fitted_interactions):
    """
    Combine two keysets of interactions, merges them, and returns the combined
    list with a priority ordering to make sure bonds are written first into file

    Parameters
    ----------
    block_interactions: dictionary keys of interactions from a block
    fitted_interactions: dictionary keys of interactions which have been fitted

    Returns
    -------
    sorted_list: list
        interaction names sorted with priority to key interaction types
    """
    combined_set = set(block_interactions) | set(fitted_interactions)
    priority_order = ['bonds', 'constraints', 'angles', 'dihedrals']
    order_map = defaultdict(lambda: len(priority_order))  # unknown keys get max priority
    order_map.update({key: i for i, key in enumerate(priority_order)})

    sorted_list = sorted(combined_set, key=lambda x: order_map[x])
    return sorted_list

def interaction_finalising(block, fitted_interactions):
    """
    Function to write fitted interactions into a block. Ensures that interactions
    that haven't otherwise been analysed (e.g. a fragment already parameterised elsewhere)
    are not overwritten.

    Parameters
    ----------
    block: vermouth.molecule.block
        Input block analysed. Modified in place.
    fitted_interactions: fast_forward.InteractionFitter.interactions_dict
        Dictionary of fitted interactions to be assigned into the block

    """
    all_interactions = finalise_interaction_types(block.interactions.keys(), fitted_interactions.keys())

    # can't iterate over something we're modifying on the fly
    copied_interactions = deepcopy(block)

    for inter_type in all_interactions:
        # first remove commented interactions from the block, these will have been analysed
        for interaction in block.interactions[inter_type]:
            if interaction.meta.get('comment', None):
                copied_interactions.remove_interaction(inter_type,
                                                       interaction.atoms)
        # reset the block interactions to be whatever's left of the copied interactions
        block.interactions[inter_type] = copied_interactions.interactions[inter_type]
        # now add the fitted interactions back in
        for new_interaction in fitted_interactions[inter_type]:
            block.add_interaction(inter_type,
                                  new_interaction.atoms,
                                  new_interaction.parameters,
                                  new_interaction.meta)
