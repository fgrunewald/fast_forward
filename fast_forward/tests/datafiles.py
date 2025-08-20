
from pathlib import Path

try:
    import pkg_resources
except ImportError:
    import os
    TEST_DATA = os.path.join(os.path.dirname(__file__), 'data')
    del os
else:
    TEST_DATA = Path(pkg_resources.resource_filename('fast_forward.tests', 'data'))
    del pkg_resources


GSH_AA_TRAJ = TEST_DATA / 'GSH/AA/atomistic.xtc'
GSH_AA_TPR = TEST_DATA / 'GSH/AA/atomistic.tpr'

GSH_CG_TRAJ = TEST_DATA / 'GSH/CG/mapped.xtc'
GSH_CG_GRO = TEST_DATA / 'GSH/CG/mapped.gro'
GSH_CG_TPR = TEST_DATA / 'GSH/CG/mapped.tpr'

GSH_MAP = TEST_DATA / 'GSH/GSH.map'



del Path