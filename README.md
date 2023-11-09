[![DOI](https://zenodo.org/badge/327071500.svg)](https://zenodo.org/badge/latestdoi/327071500)

# Fast-Forward

A simple and fast tool for forward mapping trajectories and
computing interactions. Atomistic trajectories can be forward
mapped from mapping files that correspond to backwards mapping
[style format](http://cgmartini.nl/index.php/2021-martini-online-workshop/tutorials/566-4-backward). Interactions are computed from itp files
automatically creating the time-series and distribution. Note
that we use the numba library for paralell acceleration.

## features
- residue based mapping description
- residues can be renamed
- partial residues can be mapped
- hydrogen of united-atom models can be reconstructed
- interactions are computed from itp files
- interactions can be computed for multiple molecules
- interactions are grouped by comments

## installation
```
git clone https://github.com/fgrunewald/fast_forward.git

cd fast_forward

pip install ./
```
## mapping trajectories
The following is an examplatory command for mapping trajectories.
```
ff_map -f <.trr/.xtc/.gro/.pdb> -m <.mapping/directory of .map> -s <.tpr> -o <.xtc> -mols <molecule names>
```
Note that the molecule names need to correspond to the `[moleculetype]` entries in the itp
files of the molecules that are to be mapped. Entries in mapping files are written per resiude and 
follow the backwards style mapping file format.
```
[ molecule ]
<resnameAA> <resnameCG>
[ martini ]
<list of Martini beads>
[ atoms ]
     1     <atomname AA> <bead CG> <bead CG>
```
For example, the Martini3 .mapping entry for Alanine not including the
hydrogens would look as follows and maps from the CHARMM force-field:
```
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
```

Alternatively, a directory of .map files can be passed to the  `-m` argument, such as those used 
in [Vermouth](https://github.com/marrink-lab/vermouth-martinize)-style mappings. In this format, we
use a separate map file for each residue. For example, the CHARMM map file for Alanine is:

```
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
```
## mapping trajectories + hydrogen reconstruction
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

## computing interactions
Interactions can be computed from one or multiple itp files. The molecule name
is simply matched to those found in the tpr and trajectory file and then the 
interactions are computed for all commented interactions. For example, from the
following entry the first two bonds are put in a group and computed together and
the last bond is skipped.
```
[ bonds ]
1 2 1 2000 ; group1
2 3 1 2000 ; group1
4 5 1 1000
```
The command looks as follows:
```
ff_inter -f <.trr/.xtc/.gro/.pdb> -i <.itp> -s <.tpr> -pref <prefix>
```
