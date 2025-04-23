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


