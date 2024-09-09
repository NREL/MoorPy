# Hasty example with a Subsystem and the new dynamic features

import numpy as np
import matplotlib.pyplot as plt
import moorpy as mp


# ===== Create a MoorPy system to test with =====

# ----- choose some system geometry parameters -----

depth     = 200                             # water depth [m]
angles    = np.radians([60, 300])      # line headings list [rad]
rAnchor   = 600                            # anchor radius/spacing [m]
zFair     = -21                             # fairlead z elevation [m]
rFair     = 20                              # fairlead radius [m]
lineLength= 650                            # line unstretched length [m]
typeName  = "chain1"                        # identifier string for the line type


# ----- Now a Subsystem in a larger System with a floating body -----

# Create new MoorPy System and set its depth
ms = mp.System(depth=depth)

# add a line type
ms.setLineType(dnommm=120, material='chain', name=typeName, Cd=1.4, Ca=1, CdAx=0.4, CaAx=0)  # this would be 120 mm chain

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
ss = mp.Subsystem(mooringSys=ms, depth=depth, spacing=rAnchor, rBfair=[10,0,-20])

# set up the line types
ms.setLineType(180, 'chain', name='one', Cd=1.4, Ca=1, CdAx=0.4, CaAx=0)
ms.setLineType( 50, 'chain', name='two', Cd=1.4, Ca=1, CdAx=0.4, CaAx=0)


# set up the lines and points and stuff
ls = [350, 300]
ts = ['one', 'two']
ss.makeGeneric(lengths=ls, types=ts)
ss.lineList[1].nNodes = 10  # reduce nodes on rope line for easier viewing
ss.initialize()  # update things after changing node number

# add points that the subSystem will attach to...
ms.addPoint(1, [-rAnchor, 100, -depth])   # Point 5 - create anchor point (type 0, fixed)
ms.addPoint(1, [ -rFair ,   0,  zFair])   # Point 6 - create fairlead point (type 0, fixed)
ms.bodyList[0].attachPoint(6, [-rFair, 0, zFair])  # attach the fairlead Point to the Body 

# string the Subsystem between the two points!
ms.lineList.append(ss)  # add the SubSystem to the System's lineList
ss.number = 3
ms.pointList[4].attachLine(3, 0)  # attach it to the respective points
ms.pointList[5].attachLine(3, 1)  # attach it to the respective points


# ----- run the model to check that the Subsystem is working -----

ms.initialize()                                             # make sure everything's connected

ms.solveEquilibrium()                                       # equilibrate
fig, ax = ms.plot()                                         # plot the system in original configuration

ms.bodyList[0].f6Ext = np.array([3e6, 0, 0, 0, 0, 0])       # apply an external force on the body 
ms.solveEquilibrium()                                      # equilibrate
fig, ax = ms.plot(ax=ax, color='red')                       # plot the system in displaced configuration (on the same plot, in red)

print(f"Body offset position is {ms.bodyList[0].r6}")



# ===== Now try out Serag's dynamic tension functions with it =====

# dynamic solve of some lines
kbot = 3E+06
cbot = 3E+05

moorpy_freqs = []
moorpy_fairten_psds = []
moorpy_ten_stds = []

omegas = np.array([ 0.02,  0.04,  0.06,  0.08 ])
S_zeta = np.array([ 10.0 ,  10.0 ,  10.0 ,  10.0  ])
RAO_fl = np.array([[ 2.0 ,  0.0 ,  0.0 ],
                   [ 2.0 ,  0.0 ,  0.0 ],
                   [ 2.0 ,  0.0 ,  0.0 ],
                   [ 2.0 ,  0.0 ,  0.0 ]])

T_nodes_psd_fd,T_nodes_std_fd,s,r_static,r_dynamic,r_total,X = ms.lineList[1].dynamicSolve(
                                omegas,S_zeta,RAO_A=0,RAO_B=RAO_fl,depth=-ms.depth,
                                kbot=kbot,cbot=cbot, seabed_tol=1e-4, tol=0.01, iters=100, w=0.8)
    
fig2,ax2 = plt.subplots(1,1)

plt.plot(T_nodes_std_fd,'-r',label = 'MoorPy (FD)')
plt.xlabel('Node number (anchor to fairlead)')
plt.ylabel('Fairlead tension std. dev [$N$]')

# now try a Subsystem
T_nodes_psd_fd,T_nodes_std_fd,s,r_static,r_dynamic,r_total,X = ms.lineList[2].dynamicSolve(
                                omegas,S_zeta,RAO_A=0,RAO_B=RAO_fl,depth=-ms.depth,
                                kbot=kbot,cbot=cbot, seabed_tol=1e-4, tol=0.01, iters=100, w=0.8)
    
fig2,ax2 = plt.subplots(1,1)
plt.plot(T_nodes_std_fd,'-r',label = 'MoorPy (FD)')
plt.xlabel('Node number (anchor to fairlead)')
plt.ylabel('Fairlead tension std. dev [$N$]')


# Some mode shape plots
ms.lineList[1].getModes(plot_modes=7, amp_factor=1000)  # with a Line
ms.lineList[2].getModes(plot_modes=7, amp_factor=1000)  # with a Subsystem

plt.show()
