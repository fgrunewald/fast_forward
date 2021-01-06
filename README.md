Simple tool for making COG bonded distributions and some rudementary fitting. 

# installation
```
git clone https://github.com/fgrunewald/pycgmap.git

cd pycgmap 

pip install ./
```
# usage
```
pycgmap -f <.trr/.xtc/.gro/.pdb> -n <.ndx> -s <.tpr> -o <.xtc> -itp <.itp>
```
the index file needs to be a mapping of the molecule to atomisitic coordinates and the itp file a full itp file with or without parameters. To generate both have a look at the follwing tool: https://github.com/marrink-lab/pycgbuilder
