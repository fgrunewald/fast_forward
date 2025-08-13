"""
Manage interaction blocks that have and haven't been commented in an input itp file
"""

from collections import defaultdict
from fast_forward.virtual_site_mass import MARTINI_MASS

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

def check_vs(block, interactions):
    """
    Check whether there are virtual sites in a block. If there are, find their indices and return them in a list
    """
    virtual_sites = []
    if any(['virtual_site' in i for i in interactions]):
        virtual_sites = [i for i in block.nodes if float(block.nodes[i].get('mass',
                                                                            MARTINI_MASS.get(
                                                                                block.nodes[i].get('atype')[0], 72)
                                                                            ))<1]
    return virtual_sites


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
    virtual_sites = check_vs(block, all_interactions)

    for inter_type in all_interactions:
        # first remove commented interactions from the block, these will have been analysed
        # can't remove while iterating over, so first find interactions
        for_removal = []
        for interaction in block.interactions[inter_type]:
            if interaction.meta.get('comment', None):
                for_removal.append((inter_type, interaction.atoms))
        # then remove them
        for interaction in for_removal:
            block.remove_interaction(interaction[0],
                                     interaction[1])

        # now add the fitted interactions back in
        for new_interaction in fitted_interactions[inter_type]:
            # check that we haven't made a constraint to a virtual site, which is not allowed by Gromacs
            if inter_type == 'constraints':
                if any([i in virtual_sites for i in new_interaction.atoms]):
                    # can't work out how to reset the meta of the particular bond we want
                    # so just remove it and add it back in purely as a bond
                    block.remove_interaction('bonds',
                                             # something weird happening here between np and how vermouth
                                             # looks for interactions to remove
                                             # atoms in bonds get set as tuple(np.int64, np.int64) somewhere
                                             # but I can't track it down...
                                             # this works as expected though.
                                             tuple(new_interaction.atoms)
                                             )
                    # remove the ifndef from the interaction meta so we don't write it to the bond,
                    del new_interaction.meta['ifndef']
                    # then write the interaction together with the initially calculated high force constant.
                    block.add_interaction('bonds',
                                          new_interaction.atoms,
                                          new_interaction.parameters+[new_interaction.meta['fc']],
                                          new_interaction.meta)
                    # don't write the constraint as would otherwise happen below
                    continue

            # add the new interaction in to the block
            block.add_interaction(inter_type,
                                  new_interaction.atoms,
                                  new_interaction.parameters,
                                  new_interaction.meta)
