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
    