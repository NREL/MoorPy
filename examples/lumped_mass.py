import moorpy as mp
import os
import numpy as np
from moorpy.helpers import lines2subsystem, lines2ss
try:
    import pickle5 as pickle
except:
    import pickle
import matplotlib.pyplot as plt

def JONSWAP(ws, Hs, Tp, Gamma=None):
    '''JONSWAP function copied from RAFT
    '''
    # If peak shape parameter gamma is not specified, use the recommendation 
    # from IEC 61400-3 as a function of Hs and Tp. For PM spectrum, use 1.
    if not Gamma:
        TpOvrSqrtHs = Tp/np.sqrt(Hs)
        if TpOvrSqrtHs <= 3.6:
            Gamma = 5.0
        elif TpOvrSqrtHs >= 5.0:
            Gamma = 1.0
        else:
            Gamma = np.exp( 5.75 - 1.15*TpOvrSqrtHs )
    
    # handle both scalar and array inputs
    if isinstance(ws, (list, tuple, np.ndarray)):
        ws = np.array(ws)
    else:
        ws = np.array([ws])

    # initialize output array
    S = np.zeros(len(ws))


    # the calculations
    f        = 0.5/np.pi * ws                         # wave frequencies in Hz
    fpOvrf4  = pow((Tp*f), -4.0)                      # a common term, (fp/f)^4 = (Tp*f)^(-4)
    C        = 1.0 - ( 0.287*np.log(Gamma) )          # normalizing factor
    Sigma = 0.07*(f <= 1.0/Tp) + 0.09*(f > 1.0/Tp)    # scaling factor

    Alpha = np.exp( -0.5*((f*Tp - 1.0)/Sigma)**2 )

    return  0.5/np.pi *C* 0.3125*Hs*Hs*fpOvrf4/f *np.exp( -1.25*fpOvrf4 )* Gamma**Alpha



current_dir = os.path.dirname(os.path.abspath(__file__))
ms = mp.System(os.path.join(current_dir, 'volturn_chain.dat'))
ms.initialize()
ms.solveEquilibrium()
# ms = lines2ss(ms) # For cases with multisegment lines, need to convert each of them to a subsystem

# Updates the dynamic matrices of all the lines in the system.
# This function can only properly update the inertia, added mass, and stiffness matrices of each line.        
# Though it updates the damping matrix, this is done considering unitary amplitude motions of the nodes and
# no fluid kinematics, so it is not correct.
ms.updateSystemDynamicMatrices() 
M, A, B, K_dyn = ms.getCoupledDynamicMatrices() # Get the dynamic matrices of the system
K_qsA  = ms.getCoupledStiffnessA(lines_only=True) # We can also compute the stiffness matrix without the lumped mass model. K_dyn should be similar to K_qsa


# Get the dynamic tension along the line
line = ms.lineList[0] # Get the line object
RAO_data = pickle.load(open(os.path.join(current_dir, 'RAO_fl.pkl'), 'rb')) # Read the nFreq x 3 RAO matrix (nFreq x 4 complex numpy array, first column are the frequencies in rad/s)
RAO_fl = RAO_data[:, 1:] # Motion RAOs of the fairlead
w = RAO_data[:, 0] # Frequencies of the RAO data
Sw = JONSWAP(ws = w, Hs = 6, Tp = 8) # Arbitrary wave spectrum
T_nodes_amp, T_nodes_psd,T_nodes_std,s,r_static,r_dynamic,r_total,X = line.dynamicSolve(w, Sw, RAO_A=0,RAO_B=RAO_fl, depth=np.abs(line.rA[2]))

fig, ax = plt.subplots(1, 1)
ax.plot(w, T_nodes_psd[:,-1], '-k')
ax.set_xlabel('Frequency (rad/s)')
ax.set_ylabel('PSD fairlead tension (N^2.s/rad)')
plt.show()