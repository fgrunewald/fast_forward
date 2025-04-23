ff_map
********

Background
==========

To parameterise a molecule for the martini force field, effective interactions between
beads need to be mapped from the underlying mapped interactions between groups of atoms.
This can be achieved in fast forward using the `ff_inter` program, which takes a
mapped trajectory together with a topology file describing the system to calculate
the distribution of the interactions.

The command looks as follows:

.. code-block::

    ff_inter -f <.trr/.xtc/.gro/.pdb> -i <.itp> -s <.tpr> -pref <prefix>

File format
===========
Interactions can be computed from one or multiple itp files. The molecule name
is simply matched to those found in the tpr and trajectory file and then the
interactions are computed for all commented interactions. For example, from the
following entry the first two bonds are put in a group and computed together and
the last bond is skipped.

.. code-block::

    [ bonds ]
    1 2 1 2000 ; group1
    2 3 1 2000 ; group1
    4 5 1 1000
