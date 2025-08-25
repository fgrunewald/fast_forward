
from fast_forward.tests.datafiles import  GSH_CG_TPR, GSH_CG_TRAJ, GSH_ITP_INTIIAL, GSH_ITP_OUTPUT, GSH_DISTS

import subprocess
import numpy as np
from vermouth.tests.helper_functions import find_in_path
from glob import glob
from vermouth.tests.integration_tests.test_integration import compare_itp
from pathlib import Path
def test_ff_inter(tmp_path, monkeypatch):

    monkeypatch.chdir(tmp_path)
    ff_inter = find_in_path(names=('ff_inter', ))

    command = [ff_inter,
        '-f', GSH_CG_TRAJ,
        '-s', GSH_CG_TPR,
        '-i', GSH_ITP_INTIIAL,
        '-max-dihedral', '5',
        '-dists',
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

    reference_dats = glob(str(GSH_DISTS))
    output_dats = [i for i in files if i.suffix == '.dat']

    # check we have the same set of output distributions as expected
    assert set([Path(i).name for i in reference_dats]) == set([i.name for i in output_dats])
    # ensure identical distributions for each
    for f0, f1 in zip(sorted(list(set([Path(i).name for i in reference_dats]))),
                      sorted(list(set([i.name for i in output_dats])))):
        data0 = np.loadtxt(f0)
        data1 = np.loadtxt(f1)
        assert np.all(np.isclose(data0, data1))

    # compare the output itps
    output_itp = [i for i in files if i.suffix == '.itp'][0]
    compare_itp(GSH_ITP_OUTPUT, output_itp)
