

from fast_forward.tests.datafiles import GSH_AA_TPR, GSH_AA_TRAJ, GSH_CG_GRO, GSH_CG_TRAJ, GSH_MAP

import subprocess
import numpy as np
from MDAnalysis import Universe
from vermouth.tests.helper_functions import find_in_path

def test_ff_map(tmp_path, monkeypatch):

    monkeypatch.chdir(tmp_path)
    ff_map = find_in_path(names=('ff_map', ))

    command = [ff_map,
        '-f', GSH_AA_TRAJ,
        '-s', GSH_AA_TPR,
        '-m', GSH_MAP,
        '-mols', 'LIG', '-o', 'mapped.xtc',
    ]

    proc = subprocess.run(command, cwd='.', timeout=60, check=False,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          universal_newlines=True)

    exit_code = proc.returncode
    if exit_code:
        print(proc.stdout)
        print(proc.stderr)
        assert not exit_code

    # check we generated some files
    files = list(tmp_path.iterdir())
    assert files
    # expect a .gro and .xtc from this command
    assert len(files) == 2

    reference_universe = Universe(GSH_CG_GRO, GSH_CG_TRAJ)
    new_universe = Universe([i for i in files if '.gro' in i.name][0],
                            [i for i in files if '.xtc' in i.name][0],
                            )

    # assert that the coordinates we map to are the same as in the reference
    for ts, ts0 in zip(reference_universe.trajectory, new_universe.trajectory):
        assert np.all(np.isclose(reference_universe.atoms.positions, new_universe.atoms.positions))

