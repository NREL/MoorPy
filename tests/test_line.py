# tests MoorPy Line functionality and results (work in progres)

import pytest

from numpy.testing import assert_allclose

import numpy as np
import moorpy as mp
#from moorpy.MoorProps import getLineProps
from moorpy.helpers import getLineProps

import matplotlib.pyplot as plt


inCBs = [0, 1.0, 10.0]  # friction coefficients as inputs for test_seabed



def test_line_stiffness():
    '''Checks stiffness of mooring lines.'''


def test_getLineProps_diameter_limits():
    '''getLineProps must apply the d_min / d_max range checks correctly.

    Regression test for two bugs in the diameter validation:
      1. ``d`` (the diameter in metres) was referenced in the range checks
         before it was assigned, so any material with ``d_min >= 0`` or
         ``d_max >= 0`` raised ``UnboundLocalError`` instead of the intended
         range Exception.
      2. ``loadLineProps`` read ``d_max`` from the key ``'d_dmax'`` (a typo),
         so a user-specified ``d_max`` was silently dropped to the -1
         "disabled" default and the upper-bound check never fired.
    '''
    base = {'mass_d2': 100.0, 'MBL_0': 0.0, 'MBL_d': 1e6}

    # within the valid range -> succeeds and returns a line type dict
    src = {'lineProps': {'rope': {**base, 'd_min': 0.01, 'd_max': 0.20}}}
    lt = getLineProps(50.0, 'rope', source=src)  # 50 mm = 0.05 m, in range
    assert isinstance(lt, dict) and lt['material'] == 'rope'

    # below d_min -> the intended Exception (not UnboundLocalError)
    with pytest.raises(Exception, match='less than the min'):
        getLineProps(5.0, 'rope', source=src)    # 5 mm = 0.005 m < 0.01

    # above d_max -> the intended Exception; this only works if the d_max
    # typo is fixed (otherwise d_max stays -1 and the check is disabled)
    with pytest.raises(Exception, match='greater than the max'):
        getLineProps(300.0, 'rope', source=src)  # 300 mm = 0.3 m > 0.20


if __name__ == '__main__':
    
    import moorpy as mp
    import matplotlib.pyplot as plt
    
    ms = mp.System(depth=100)
    ms.setLineType(100, 'chain', name='chain')
    
    ms.addPoint(1,  [1, 0, -100]) # anchor point
    ms.addPoint(-1, [0, 0, 0]) # moving point
    
    ms.addLine(99, 'chain', pointA=1, pointB=2)
    
    ms.initialize()
    
    fig, ax = ms.plot()
    
    ms.solveEquilibrium()
    f0 = ms.pointList[1].getForces()
    print(f0)
    print(ms.lineList[0].KA[1,1])
    
    ms.pointList[1].setPosition([0,0.1,0])
    ms.solveEquilibrium()
    f1 = ms.pointList[1].getForces()
    print(f1)
    print(ms.lineList[0].KA[1,1])
    
    ms.plot(ax=ax, color='red')
    
    plt.show()
    