# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 09:35:50 2024

@author: lsirkis
"""
import moorpy as mp
import numpy as np
from moorpy import helpers

# original line, unchanged
ms = mp.System(file='C:/Users/LSIRKIS/Documents/OpenFAST/IEA-15-240-RWT-UMaineSemi/Mooring/MoorDynSemiTautUpd.dat')
ms.initialize()
ms.solveEquilibrium()
settings = {}
settings["linelabels"] = True
settings["pointlabels"] = True  
ms.plot(**settings)

from copy import deepcopy
ms1 = deepcopy(ms) # copy the moorpy system
############## Create one subsystem ###############
ms1 = helpers.lines2subsystem([0,1],ms1)
# lines 0 and 1 are now replaced by subsystem, can be found in ms.lineList[4]

# ############ Create all subsystems ##############
print('\nReplacing all lines with subsystems\n')
nsecs = [2,2,2] # number of sections in each mooring line
nlines = len(nsecs) # number of lines

# replace all lines with subsystems
for i in range(0,nlines):
    sdelList = []
    for j in range(0,nsecs[i]): # always starts at 0 since old lines deleted and new subsystem added at end
        sdelList.append(j)
    ms = helpers.lines2subsystem(sdelList,ms)

ms.initialize()
ms.solveEquilibrium()
  
ms.plot(**settings)