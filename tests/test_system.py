# tests MoorPy System functionality and results

import pytest

from numpy.testing import assert_allclose

import numpy as np
import moorpy as mp
from moorpy.MoorProps import getLineProps

import matplotlib.pyplot as plt
    
def test_basic():


    depth     = 600
    angle     = np.arange(3)*np.pi*2/3  # line headings list
    anchorR   = 1600                    # anchor radius/spacing
    fair_depth= 21 
    fairR     = 20
    LineLength= 1800
    typeName  = "chain"                 # identifier string for line type

    # --------------- set up mooring system ---------------------

    # Create blank system object
    ms = mp.System()

    # Set the depth of the system to the depth of the input value
    ms.depth = depth

    # add a line type
    ms.lineTypes[typeName] = getLineProps(120, name=typeName)

    # Add a free, body at [0,0,0] to the system (including some properties to make it hydrostatically stiff)
    ms.addBody(0, np.zeros(6), m=1e6, v=1e3, rM=100, AWP=1e6)

    # Set the anchor points of the system
    anchors = []
    for i in range(len(angle)):
        ms.addPoint(1, np.array([anchorR*np.cos(angle[i]), anchorR*np.sin(angle[i]), -ms.depth], dtype=float))
        anchors.append(len(ms.pointList))

    # Set the points that are attached to the body to the system
    bodypts = []
    for i in range(len(angle)):
        ms.addPoint(1, np.array([fairR*np.cos(angle[i]), fairR*np.sin(angle[i]), -fair_depth], dtype=float))
        bodypts.append(len(ms.pointList))
        ms.bodyList[0].attachPoint(ms.pointList[bodypts[i]-1].number, ms.pointList[bodypts[i]-1].r - ms.bodyList[0].r6[:3])

    # Add and attach lines to go from the anchor points to the body points
    for i in range(len(angle)):    
        ms.addLine(LineLength, typeName)
        line = len(ms.lineList)
        ms.pointList[anchors[i]-1].attachLine(ms.lineList[line-1].number, 0)
        ms.pointList[bodypts[i]-1].attachLine(ms.lineList[line-1].number, 1)

    # ------- simulate it ----------
    ms.initialize()                                             # make sure everything's connected
    ms.bodyList[0].setPosition([10,10,1,0,0,0])                 # apply an offset
    ms.solveEquilibrium3()                                      # equilibrate - see if it goes back to zero!
    
    
    # check
    assert_allclose(ms.bodyList[0].r6, np.zeros(6), rtol=0, atol=0.01, verbose=True)
    

def test_basicU():
    '''Compares a system with a U-shape line with seabed contact with an equivalent case 
       that has a node along the U.'''
       
    # a seabed contact case with 2 lines to form a U
    ms1 = mp.System()

    ms1.depth = 100
    
    ms1.lineTypes['chain'] = getLineProps(120, name='chain')

    ms1.addPoint(1, [-200, 0 , -100])                 # anchor point
    ms1.addPoint(0, [-100, 0 ,  -50], m=0, v=50)      # float
    ms1.addPoint(0, [   0, 0 , -100])                 # midpoint
    ms1.addPoint(0, [ 100, 0 ,  -40], m=0, v=50)      # float
    ms1.addPoint(1, [ 200, 0 , -100])                 # anchor point
        
    ms1.addLine(120, 'chain', pointA=1, pointB=2)
    ms1.addLine(125, 'chain', pointA=2, pointB=3)
    ms1.addLine(125, 'chain', pointA=3, pointB=4)
    ms1.addLine(120, 'chain', pointA=4, pointB=5)
    
    
    # a seabed contact case with single U line
    msU = mp.System()

    msU.depth = 100
    
    msU.lineTypes['chain'] = getLineProps(120, name='chain')

    msU.addPoint(1, [-200, 0 , -100])                 # anchor point
    msU.addPoint(0, [-100, 0 ,  -50], m=0, v=50)      # float
    msU.addPoint(0, [ 100, 0 ,  -40], m=0, v=50)      # float
    msU.addPoint(1, [ 200, 0 , -100])                 # anchor point
        
    msU.addLine(120, 'chain', pointA=1, pointB=2)
    msU.addLine(250, 'chain', pointA=2, pointB=3)
    msU.addLine(120, 'chain', pointA=3, pointB=4)
        

    # ------- simulate it ----------
    ms1.initialize()                                             # make sure everything's connected
    msU.initialize()                                             # make sure everything's connected
    
    ms1.solveEquilibrium3(tol=0.0001)                            # equilibrate - see if it goes back to zero!
    msU.solveEquilibrium3(tol=0.0001)                            # equilibrate - see if it goes back to zero!
    
    
    # compare floating point positions
    assert_allclose(np.hstack([ms1.pointList[1].r, ms1.pointList[3].r]),
                    np.hstack([msU.pointList[1].r, msU.pointList[2].r]), rtol=0, atol=0.001, verbose=True)


if __name__ == '__main__':
    
    #test_basic()
    
    '''
      # a seabed contact case with 2 lines to form a U
    ms1 = mp.System()

    ms1.depth = 100
    
    ms1.lineTypes['chain'] = getLineProps(120, name='chain')

    ms1.addPoint(1, [-200, 0 , -100])                 # anchor point
    ms1.addPoint(1, [-100, 0 ,  -50], m=0, v=50)      # float
    ms1.addPoint(0, [   0, 0 , -100])                 # midpoint
    ms1.addPoint(1, [ 100, 0 ,  -40], m=0, v=50)      # float
    ms1.addPoint(1, [ 200, 0 , -100])                 # anchor point
        
    ms1.addLine(120, 'chain', pointA=1, pointB=2)
    ms1.addLine(125, 'chain', pointA=2, pointB=3)
    ms1.addLine(125, 'chain', pointA=3, pointB=4)
    ms1.addLine(120, 'chain', pointA=4, pointB=5)
    
    
    # a seabed contact case with single U line
    ms = mp.System()

    ms.depth = 100
    
    ms.lineTypes['chain'] = getLineProps(120, name='chain')

    ms.addPoint(1, [-200, 0 , -100])                 # anchor point
    ms.addPoint(0, [-100, 0 ,  -50], m=0, v=50)      # float
    ms.addPoint(0, [ 100, 0 ,  -40], m=0, v=50)      # float
    ms.addPoint(1, [ 200, 0 , -100])                 # anchor point
        
    ms.addLine(120, 'chain', pointA=1, pointB=2)
    ms.addLine(250, 'chain', pointA=2, pointB=3)
    ms.addLine(120, 'chain', pointA=3, pointB=4)
        

    # ------- simulate it ----------
    ms1.initialize()                                             # make sure everything's connected
    ms.initialize()                                             # make sure everything's connected
    
    fig, ax = ms1.plot(color='b')
    ms.plot(ax=ax, color='k')
    
    #ms.display=2
    
    ms1.solveEquilibrium3(maxIter=20)                                      # equilibrate - see if it goes back to zero!
    #ms.solveEquilibrium3(maxIter=20)                                      # equilibrate - see if it goes back to zero!
    
    ms1.plot(ax=ax, color='g')
    #ms.plot(ax=ax, color='r')
    
    # compare
    print(ms1.pointList[1].getForces())
    print(ms.pointList[1].getForces())
    print(ms1.pointList[3].getForces())
    print(ms.pointList[2].getForces())
    
    print(ms1.pointList[1].getStiffnessA())
    print(ms.pointList[1].getStiffnessA())
    print(ms1.pointList[3].getStiffnessA())
    print(ms.pointList[2].getStiffnessA())
    '''
    
    # a seabed contact case with 2 lines to form a U
    ms1 = mp.System()

    ms1.depth = 100
    
    ms1.lineTypes['chain'] = getLineProps(120, name='chain')

    ms1.addPoint(1, [-200, 0 , -100])                 # anchor point
    ms1.addPoint(0, [-100, 0 ,  -50], m=0, v=50)      # float
    ms1.addPoint(0, [   0, 0 , -100])                 # midpoint
    ms1.addPoint(0, [ 100, 0 ,  -40], m=0, v=50)      # float
    ms1.addPoint(1, [ 200, 0 , -100])                 # anchor point
        
    ms1.addLine(120, 'chain', pointA=1, pointB=2)
    ms1.addLine(125, 'chain', pointA=2, pointB=3)
    ms1.addLine(125, 'chain', pointA=3, pointB=4)
    ms1.addLine(120, 'chain', pointA=4, pointB=5)
    
    
    # a seabed contact case with single U line
    msU = mp.System()

    msU.depth = 100
    
    msU.lineTypes['chain'] = getLineProps(120, name='chain')

    msU.addPoint(1, [-200, 0 , -100])                 # anchor point
    msU.addPoint(0, [-100, 0 ,  -50], m=0, v=50)      # float
    msU.addPoint(0, [ 100, 0 ,  -40], m=0, v=50)      # float
    msU.addPoint(1, [ 200, 0 , -100])                 # anchor point
        
    msU.addLine(120, 'chain', pointA=1, pointB=2)
    msU.addLine(250, 'chain', pointA=2, pointB=3)
    msU.addLine(120, 'chain', pointA=3, pointB=4)
        

    # ------- simulate it ----------
    ms1.initialize()                                             # make sure everything's connected
    msU.initialize()                                             # make sure everything's connected
    
    ms1.solveEquilibrium3(tol=0.0001)                            # equilibrate - see if it goes back to zero!
    msU.solveEquilibrium3(tol=0.0001)                            # equilibrate - see if it goes back to zero!
       
    fig,ax = ms1.plot(color='g')
    msU.plot(color=[1,0,0,0.5], ax=ax)
    
    plt.show()
    