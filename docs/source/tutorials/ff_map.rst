ff_map
********

Background
==========

ff_map is the mapping program of fast forward. Using (a collection of) mapping files, it
will convert a simulation at all-atom resolution to one at a coarse grained resolution.

mapping file format
===================


The following is an example command for mapping trajectories.

.. code-block::

    ff_map -f <.trr/.xtc/.gro/.pdb> -m <.mapping/directory of .map> -s <.tpr> -o <.xtc> -mols <molecule names>


Note that the molecule names need to correspond to the `[moleculetype]` entries in the itp
files of the molecules that are to be mapped. Entries in mapping files are written per resiude and
follow the backwards style mapping file format.

.. code-block::

    [ molecule ]
    <resnameAA> <resnameCG>
    [ martini ]
    <list of Martini beads>
    [ atoms ]
         1     <atomname AA> <bead CG> <bead CG>

For example, the Martini3 .mapping entry for Alanine not including the
hydrogens would look as follows and maps from the CHARMM force-field:

.. code-block::

    [ molecule ]
    ALA ALA
    [ martini ]
    BB SC1
    [ atoms ]
        1     N    BB
        2    HN    BB
        3    CA    BB
        5    CB    SC1
        9     C    BB
       10     O    BB

Alternatively, a directory of .map files can be passed to the  `-m` argument, such as those used
in [Vermouth](https://github.com/marrink-lab/vermouth-martinize)-style mappings. In this format, we
use a separate map file for each residue. For example, the CHARMM map file for Alanine is:

.. code-block::

    [ molecule ]
    ALA

    [from]
    charmm

    [to]
    martini3001

    [ martini ]
    BB SC1

    [ mapping ]
    charmm27 charmm36

    [ atoms ]
        1     N    BB
        2    HN    BB
        3    CA    BB
        4    HA    !BB
        5    CB    SC1
        6   HB1    !SC1
        7   HB2    !SC1
        8   HB3    !SC1
        9     C    BB
       10     O    BB

    [ chiral ]
      CB     CA    N    C
      HB1    CA    N    C
      HB2    CA    N    C
      HB3    CA    N    C

    [ chiral ]
      HA     CA    N    CB    C ; L-Ala
    ; HA     CA    N    C    CB ; D-Ala

mapping trajectories + hydrogen reconstruction
==============================================

The Martini3 model is based on center of geometry mappings, for which
it is important to include the hydrogen atoms. However, certain atomistic
united-atom models have no explicit hydrogen atoms. In oder to still obtain
good bonded interactions fast_forward can reconstruct the hydrogen atoms
using geometric rules. The only input required is a list of atom-types and
weather they correspond to CH, CH2, CH3 or double bonded CH group (CHd).

The algorithm is directly taken from the [buildH](https://github.com/patrickfuchs/buildH) package.
If you use this feauture please cite https://joss.theoj.org/papers/10.21105/joss.03521.

An example for POPC can be found in the example folder, where a CHARMM36m
POPC was stripped of all hydrogen and mapped with reconstrucion.

Note the reconstruction only affects the mapping, atomistic coordinates with
hydrogen are not written out.


