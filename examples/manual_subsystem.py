# MoorPy Example Script:
# Example of manually setting up a mooring system in MoorPy and solving equilibrium.
# This example features a subsystem to make sure it works. This will become a test after.

import numpy as np
import matplotlib.pyplot as plt
import moorpy as mp
from moorpy.MoorProps import getLineProps
from moorpy.subsystem import Subsystem


# ----- choose some system geometry parameters -----

depth     = 200                             # water depth [m]
angles    = np.radians([60, 300])      # line headings list [rad]
rAnchor   = 600                            # anchor radius/spacing [m]
zFair     = -21                             # fairlead z elevation [m]
rFair     = 20                              # fairlead radius [m]
lineLength= 650                            # line unstretched length [m]
typeName  = "chain1"                        # identifier string for the line type


# ===== First a simple Subsystem by itself =====

# create a subsystem
ss = Subsystem(depth=depth, spacing=rAnchor, rBfair=[10,0,-20])

# set up the line types
ss.setLineType(180, 'chain', name='one')
ss.setLineType( 50, 'chain', name='two')

# set up the lines and points and stuff
lengths = [350, 300]
types = ['one', 'two']
ss.makeGeneric(lengths, types)

# plotting examples
ss.setEndPosition([0  ,-40,-200], endB=0)
ss.setEndPosition([300,400,  -10], endB=1)
ss.staticSolve()
# 3D plot in the Subsystem's local frame
fig, ax = ss.plot()
# 2D view of the same
ss.plot2d()
# Line-like plot (in global reference frame)
fig = plt.figure()
ax = plt.axes(projection='3d')
ss.drawLine(0, ax, color='r')


# ===== Now a Subsystem in a larger System with a floating body =====

# Create new MoorPy System and set its depth
ms = mp.System(depth=depth)

# add a line type
ms.setLineType(dnommm=120, material='chain', name=typeName)  # this would be 120 mm chain

# Add a free, body at [0,0,0] to the system (including some properties to make it hydrostatically stiff)
ms.addBody(0, np.zeros(6), m=1e6, v=1e3, rM=100, AWP=1e3)

# For each line heading, set the anchor point, the fairlead point, and the line itself
for i, angle in enumerate(angles):

    # create end Points for the line
    ms.addPoint(1, [rAnchor*np.cos(angle), rAnchor*np.sin(angle), -depth])   # create anchor point (type 0, fixed)
    ms.addPoint(1, [  rFair*np.cos(angle),   rFair*np.sin(angle),  zFair])   # create fairlead point (type 0, fixed)
    
    # attach the fairlead Point to the Body (so it's fixed to the Body rather than the ground)
    ms.bodyList[0].attachPoint(2*i+2, [rFair*np.cos(angle), rFair*np.sin(angle), zFair]) 

    # add a Line going between the anchor and fairlead Points
    ms.addLine(lineLength, typeName, pointA=2*i+1, pointB=2*i+2)

# ----- Now add a SubSystem line! -----
ss = Subsystem(mooringSys=ms, depth=depth, spacing=rAnchor, rBfair=[10,0,-20])

# set up the line types
ms.setLineType(180, 'chain', name='one')
ms.setLineType( 50, 'chain', name='two')

# set up the lines and points and stuff
ls = [350, 300]
ts = ['one', 'two']
ss.makeGeneric(lengths=ls, types=ts)

# add points that the subSystem will attach to...
ms.addPoint(1, [-rAnchor, 100, -depth])   # Point 5 - create anchor point (type 0, fixed)
ms.addPoint(1, [ -rFair ,   0,  zFair])   # Point 6 - create fairlead point (type 0, fixed)
ms.bodyList[0].attachPoint(6, [-rFair, 0, zFair])  # attach the fairlead Point to the Body 

# string the Subsystem between the two points!
ms.lineList.append(ss)  # add the SubSystem to the System's lineList
ss.number = 3
ms.pointList[4].attachLine(3, 0)  # attach it to the respective points
ms.pointList[5].attachLine(3, 1)  # attach it to the respective points


# ----- run the model to demonstrate -----

ms.initialize()                                             # make sure everything's connected

ms.solveEquilibrium()                                       # equilibrate
fig, ax = ms.plot()                                         # plot the system in original configuration
#ms.unload("sample_from_manual.txt")                         # export to MD input file

ms.bodyList[0].f6Ext = np.array([3e6, 0, 0, 0, 0, 0])       # apply an external force on the body 
ms.solveEquilibrium3()                                      # equilibrate
fig, ax = ms.plot(ax=ax, color='red')                       # plot the system in displaced configuration (on the same plot, in red)

print(f"Body offset position is {ms.bodyList[0].r6}")

plt.show()



