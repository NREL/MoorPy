"""
Created on Mon Jul 1, 2024

@author: ralkarem
"""
import moorpy as mp
import numpy as np
from moorpy import helpers
from copy import deepcopy
import matplotlib.pyplot as plt
ms = mp.System(file='MoorDynSemiTautUpd.dat')

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

