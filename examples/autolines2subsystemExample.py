"""
Created on Mon Jul 1, 2024

@author: ralkarem
"""
import moorpy as mp
import numpy as np
from moorpy import helpers
from copy import deepcopy
import matplotlib.pyplot as plt
from moorpy.helpers import lines2subsystem
# ms = mp.System(file='MoorDynSemiTautUpd_simple.dat')  # test 0 - pass
# ms = mp.System(file='MoorDynSemiTautUpd.dat')  # test 1 - pass
# ms = mp.System(file='MoorDynSemiTautUpd_flipped_lines.dat') # test 2 - pass
# ms = mp.System(file='SharedMooringR.dat') # test 4 - pass
# ms = mp.System(file='SharedMooring2 (1) 1.dat') # test 5 - pass

# activate this section if test 5 is to be tested.
# for body in ms.bodyList:
#     body.m = 19911423.956678286
#     body.v = 19480.104108645974
#     body.rCG = np.array([ 1.49820657e-15,  1.49820657e-15, -2.54122031e+00])
#     body.AWP = 446.69520543229874
#     body.rM = np.array([2.24104273e-15, 1.49402849e-15, 1.19971829e+01])
#     body.type = -1
# ms.bodyList[1].setPosition([1600,0,0,0,0,0])

ms.initialize()
ms.solveEquilibrium()
settings = {}
settings["linelabels"] = True
settings["pointlabels"] = True  
ms.plot(**settings)
plt.show()
ms1 = deepcopy(ms) # copy the moorpy system
ms1 = helpers.lines2ss(ms1)
ms1.plot(**settings)
plt.show()

