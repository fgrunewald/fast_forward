# Fast-Forward

A simple and fast tool for forward mapping trajectories and
computing interactions. Atomistic trajectories can be forward
mapped from mapping files that correspond to backwards mapping
style fromat. Interactions are computed from itp files
automatically creating the time-series and distribution. Note
that we use the numba library for paralell acceleration.

## features
- residue based mapping description
- residues can be renamed
- partial residues can be mapped
- interactions are computed from itp files
- interactions can be computed for multiple molecules
- interactions are grouped by comments

## installation
```
git clone https://github.com/fgrunewald/pycgmap.git

cd fast-forward

pip install ./
```
## mapping trajectories
The following is an examplatory command for mapping trajectories.
```
ff_map -f <.trr/.xtc/.gro/.pdb> -m <.mapping> -s <.tpr> -o <.xtc>
```
Mapping files are written per resiude and follow the backwards
style mapping file format.
```
[ molecule ]
<resnameAA> <resnameCG>
[ martini ]
<list of Martini beads>
[ atoms ]
     1     <atomname AA> <bead CG> <bead CG>
```
For example, the Martini3 mapping file for Alanine not including the
hydrogens would look as follows and maps from the CHARMM force-field:
```
[ molecule ]
ALA ALA
[ martini ]
BB
[ atoms ]
    1     N    BB
    2    HN    BB
    3    CA    BB
    5    CB    BB
    9     C    BB
   10     O    BB
```
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
