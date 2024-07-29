
import numpy as np
import time
import yaml
import os
import re

 
# base class for MoorPy exceptions
class Error(Exception):
    ''' Base class for MoorPy exceptions'''
    pass

# Catenary error class
class CatenaryError(Error):
    '''Derived error class for catenary function errors. Contains an error message.'''
    def __init__(self, message):
        self.message = message

# Line Object error class
class LineError(Error):
    '''Derived error class for Line object errors. Contains an error message and the line number with the error.'''
    def __init__(self, num, message):
        self.line_num = num
        self.message = message

# Solve error class for any solver process
class SolveError(Error):
    '''Derived error class for various solver errors. Contains an error message'''
    def __init__(self, message):
        self.message = message
    

# Generic MoorPy error
class MoorPyError(Error):
    '''Derived error class for MoorPy. Contains an error message'''
    def __init__(self, message):
        self.message = str(message)




def printMat(mat):
    '''Prints a matrix to a format that is specified

    Parameters
    ----------
    mat : array
        Any matrix that is to be printed.

    Returns
    -------
    None.

    '''
    for i in range(mat.shape[0]):
        print( "\t".join(["{:+8.3e}"]*mat.shape[1]).format( *mat[i,:] ))
        
def printVec(vec):
    '''Prints a vector to a format that is specified

    Parameters
    ----------
    vec : array
        Any vector that is to be printed.

    Returns
    -------
    None.

    '''
    print( "\t".join(["{:+9.4e}"]*len(vec)).format( *vec ))



def unitVector(r):
    '''Returns the unit vector along the direction of input vector r.'''

    L = np.linalg.norm(r)

    return r/L



def getInterpNums(xlist, xin, istart=0):  # should turn into function in helpers
    '''
    Paramaters
    ----------
    xlist : array
        list of x values
    xin : float
        x value to be interpolated
    istart : int
        first lower index to try
    
    Returns
    -------
    i : int
        lower index to interpolate from
    fout : float
        fraction to return   such that y* = y[i] + fout*(y[i+1]-y[i])
    '''
    
    nx = len(xlist)
  
    if xin <= xlist[0]:  #  below lowest data point
        i = 0
        fout = 0.0
  
    elif xlist[-1] <= xin:  # above highest data point
        i = nx-1
        fout = 0.0
  
    else:  # within the data range
 
        # if istart is below the actual value, start with it instead of 
        # starting at 0 to save time, but make sure it doesn't overstep the array
        if xlist[min(istart,nx)] < xin:
            i1 = istart
        else:
            i1 = 0

        for i in range(i1, nx-1):
            if xlist[i+1] > xin:
                fout = (xin - xlist[i] )/( xlist[i+1] - xlist[i] )
                break
    
    return i, fout
    
        

def getH(r):
    '''function gets the alternator matrix, H, that when multiplied with a vector,
    returns the cross product of r and that vector

    Parameters
    ----------
    r : array
        the position vector that another vector is from a point of interest.

    Returns
    -------
    H : matrix
        the alternator matrix for the size-3 vector, r.

    '''
    
    H = np.array([[ 0   , r[2],-r[1]],
                  [-r[2], 0   , r[0]],
                  [ r[1],-r[0], 0   ]])
    return H
    
    
def rotationMatrix(x3,x2,x1):
    '''Calculates a rotation matrix based on order-z,y,x instrinsic (tait-bryan?) angles, meaning
    they are about the ROTATED axes. (rotation about z-axis would be (0,0,theta) )
    
    Parameters
    ----------
    x3, x2, x1: floats
        The angles that the rotated axes are from the nonrotated axes. Normally roll,pitch,yaw respectively. [rad]

    Returns
    -------
    R : matrix
        The rotation matrix
    '''
    # initialize the sines and cosines
    s1 = np.sin(x1) 
    c1 = np.cos(x1)
    s2 = np.sin(x2) 
    c2 = np.cos(x2)
    s3 = np.sin(x3) 
    c3 = np.cos(x3)
    
    # create the rotation matrix
    R = np.array([[ c1*c2,  c1*s2*s3-c3*s1,  s1*s3+c1*c3*s2],
                  [ c2*s1,  c1*c3+s1*s2*s3,  c3*s1*s2-c1*s3],
                  [   -s2,           c2*s3,           c2*c3]])
    
    return R    


def rotatePosition(rRelPoint, rot3):
    '''Calculates the new position of a point by applying a rotation (rotates a vector by three angles)
    
    Parameters
    ----------
    rRelPoint : array
        x,y,z coordinates of a point relative to a local frame [m]
    rot3 : array
        Three angles that describe the difference between the local frame and the global frame/ Normally roll,pitch,yaw. [rad]

    Returns
    -------
    rRel : array
        The relative rotated position of the point about the local frame [m]
    '''
    
    # get rotation matrix from three provided angles
    RotMat = rotationMatrix(rot3[0], rot3[1], rot3[2])     
    
    # find location of point in unrotated reference frame about reference point
    rRel = np.matmul(RotMat,rRelPoint)    
    
    return rRel
    

def transformPosition(rRelPoint, r6):
    '''Calculates the position of a point based on its position relative to translated and rotated 6DOF body
    
    Parameters
    ----------
    rRelPoint : array
        x,y,z coordinates of a point relative to a local frame [m]
    r6 : array
        6DOF position vector of the origin of the local frame, in the global frame coorindates [m, rad]

    Returns
    -------
    rAbs : array
        The absolute position of the point about the global frame [m]
    '''
    # note: r6 should be in global orientation frame
    
    # absolute location = rotation of relative position + absolute position of reference point
    rAbs = rotatePosition(rRelPoint, r6[3:]) + r6[:3]
        
    return rAbs
    
    
def translateForce3to6DOF(r, Fin):
    '''Takes in a position vector and a force vector (applied at the positon), and calculates 
    the resulting 6-DOF force and moment vector.    
    
    Parameters
    ----------
    r : array
        x,y,z coordinates at which force is acting [m]
    Fin : array
        x,y,z components of force [N]

    Returns
    -------
    Fout : array
        The resulting force and moment vector [N, Nm]
    '''
    
    # initialize output vector as same dtype as input vector (to support both real and complex inputs)
    Fout = np.zeros(6, dtype=Fin.dtype) 
    
    # set the first three elements of the output vector the same as the input vector
    Fout[:3] = Fin
    
    # set the last three elements of the output vector as the cross product of r and Fin
    Fout[3:] = np.cross(r, Fin)
    
    return Fout


def set_plot_center(ax, x=None, y=None, z=None):
    '''Sets the center point in x and y of a 3d plot'''
    
    # adjust the center point of the figure if requested, by moving out one of the bounds
    if not x is None:
        xlims = ax.get_xlim3d()
        if   x > np.mean(xlims): ax.set_xlim([xlims[0], x + (x - xlims[0])])
        elif x < np.mean(xlims): ax.set_xlim([x - (xlims[1] - x), xlims[1]])
        
    if not y is None:
        ylims = ax.get_ylim3d()
        if   y > np.mean(ylims): ax.set_ylim([ylims[0], y + (y - ylims[0])])
        elif y < np.mean(ylims): ax.set_ylim([y - (ylims[1] - y), ylims[1]])
    
    if not z is None:
        zlims = ax.get_zlim3d()
        if   z > np.mean(zlims): ax.set_zlim([zlims[0], z + (z - zlims[0])])
        elif z < np.mean(zlims): ax.set_zlim([z - (zlims[1] - z), zlims[1]])
        
    # make sure the aspect ratio stays equal
    set_axes_equal(ax)
        
    '''    
        # set the AXIS bounds on the axis (changing these bounds can change the perspective of the matplotlib figure)
        if xbounds != None:
            ax.set_xbound(xbounds[0], xbounds[1])
            ax.autoscale(enable=False, axis='x')
        if ybounds != None:
            ax.set_ybound(ybounds[0], ybounds[1])
            ax.autoscale(enable=False, axis='y')
        if zbounds != None:
            ax.set_zbound(zbounds[0], zbounds[1])
            ax.autoscale(enable=False, axis='x')
    '''

def set_axes_equal(ax):
    '''Sets 3D plot axes to equal scale

    Parameters
    ----------
    ax : matplotlib.pyplot axes
        the axes that are to be set equal in scale to each other.

    Returns
    -------
    None.

    '''
    
    rangex = np.diff(ax.get_xlim3d())[0]
    rangey = np.diff(ax.get_ylim3d())[0]
    rangez = np.diff(ax.get_zlim3d())[0]
    
    ax.set_box_aspect([rangex, rangey, rangez])  # note: this may require a matplotlib update
    

def quiver_data_to_segments(X, Y, Z, u, v, w, scale=1):
    '''function to help with animation of 3d quivers'''
    
    if scale < 0.0:  # negative scale input will be treated as setting the desired RMS quiver length
        scale = -scale/np.sqrt(np.mean(u**2 + v**2 + w**2))

    segments = (X, Y, Z, X+u*scale, Y+v*scale, Z+w*scale)
    segments = np.array(segments).reshape(6,-1)
    return [[[x1, y1, z1], [x2, y2, z2]] for x1, y1, z1, x2, y2, z2 in zip(*list(segments))]
    

def dsolve2(eval_func, X0, Ytarget=[], step_func=None, args=[], tol=0.0001, ytol=0, maxIter=20, 
           Xmin=[], Xmax=[], a_max=2.0, dX_last=[], stepfac=4, display=0, dodamping=False):
    '''
    PARAMETERS
    ----------    
    eval_func : function
        function to solve (will be passed array X, and must return array Y of same size)
    X0 : array
        initial guess of X
    Ytarget : array (optional)
        target function results (Y), assumed zero if not provided
    stp_func : function (optional)
        function use for adjusting the variables (computing dX) each step. 
        If not provided, Netwon's method with finite differencing is used.
    args : list
        A list of variables (e.g. the system object) to be passed to both the eval_func and step_func
    tol : float or array
        If scalar, the*relative* convergence tolerance (applied to step size components, dX).
        If an array, must be same size as X, and specifies an absolute convergence threshold for each variable.
    ytol: float, optional
        If specified, this is the absolute error tolerance that must be satisfied. This overrides the tol setting which otherwise works based on x values.
    Xmin, Xmax 
        Bounds. by default start bounds at infinity
    a_max
        maximum step size acceleration allowed
    dX_last
        Used if you want to dictate the initial step size/direction based on a previous attempt
    '''
    success = False
    start_time = time.time()
    # process inputs and format as arrays in case they aren't already
    
    X = np.array(np.atleast_1d(X0), dtype=float)         # start off design variable
    N = len(X)
    
    Xs = np.zeros([maxIter,N]) # make arrays to store X and error results of the solve
    Es = np.zeros([maxIter,N])
    dXlist = np.zeros([maxIter,N])
    dXlist2 = np.zeros([maxIter,N])
    
    damper = 1.0   # used to add a relaxation/damping factor to reduce the step size and combat instability
    
    
    # check the target Y value input
    if len(Ytarget)==N:
        Ytarget = np.array(Ytarget, dtype=float)
    elif len(Ytarget)==0:
        Ytarget = np.zeros(N, dtype=float)
    else:
        raise TypeError("Ytarget must be of same length as X0")
        
    # ensure all tolerances are positive
    if ytol==0:  # if not using ytol
        if np.isscalar(tol) and tol <= 0.0:
            raise ValueError('tol value passed to dsovle2 must be positive')
        elif not np.isscalar(tol) and any(np.array(tol) <= 0):
            raise ValueError('every tol entry passed to dsovle2 must be positive')
        
    # if a step function wasn't provided, provide a default one
    if step_func==None:
        if display>1:
            print("Using default finite difference step func")
        
        def step_func(X, args, Y, oths, Ytarget, err, tols, iter, maxIter):
            ''' this now assumes tols passed in is a vector and are absolute quantities'''
            J = np.zeros([N,N])       # Initialize the Jacobian matrix that has to be a square matrix with nRows = len(X)
            
            for i in range(N):             # Newton's method: perturb each element of the X variable by a little, calculate the outputs from the
                X2 = np.array(X)                # minimizing function, find the difference and divide by the perturbation (finding dForce/d change in design variable)
                deltaX = stepfac*tols[i]                  # note: this function uses the tols variable that is computed in dsolve based on the tol input
                X2[i] += deltaX
                Y2, _, _ = eval_func(X2, args)    # here we use the provided eval_func
                
                J[:,i] = (Y2-Y)/deltaX             # and append that column to each respective column of the Jacobian matrix
               
            if N > 1:
                dX = -np.matmul(np.linalg.inv(J), Y-Ytarget)   # Take this nth output from the minimizing function and divide it by the jacobian (derivative)
            else:
                if J[0,0] == 0.0:
                    raise ValueError('dsolve2 found a zero gradient')
                    
                dX = np.array([-(Y[0]-Ytarget[0])/J[0,0]])
                
                if display > 1:
                    print(f" step_func iter {iter} X={X[0]:9.2e}, error={Y[0]-Ytarget[0]:9.2e}, slope={J[0,0]:9.2e}, dX={dX[0]:9.2e}")

            return dX                              # returns dX (step to make)

    
    
    # handle bounds
    if len(Xmin)==0:
        Xmin = np.zeros(N)-np.inf
    elif len(Xmin)==N:
        Xmin = np.array(Xmin, dtype=float)
    else:
        raise TypeError("Xmin must be of same length as X0")
        
    if len(Xmax)==0:
        Xmax = np.zeros(N)+np.inf
    elif len(Xmax)==N:
        Xmax = np.array(Xmax, dtype=float)
    else:
        raise TypeError("Xmax must be of same length as X0")
    
    
    
    if len(dX_last)==0:
        dX_last = np.zeros(N)
    else:
        dX_last = np.array(dX_last, dtype=float)

    if display>0:
        print(f"Starting dsolve iterations>>>   aiming for Y={Ytarget}")

    
    for iter in range(maxIter):

        
        # call evaluation function
        Y, oths, stop = eval_func(X, args)
        
        # compute error
        err = Y - Ytarget
        
        if display==2:
            print(f"  new iteration #{iter} with RMS error {np.linalg.norm(err):8.3e}")
        if display>2:
            print(f"  new iteration #{iter} with X={X} and Y={Y}")

        Xs[iter,:] = X
        Es[iter,:] = err

        # stop if commanded by objective function
        if stop:
            break
        
        # handle tolerances input
        if np.isscalar(tol):
            tols = tol*(np.abs(X)+tol)
        else:
            tols = np.array(tol)
        
        # check maximum iteration
        if iter==maxIter-1:
            if display>0:
                print("Failed to find solution after "+str(iter)+" iterations, with error of "+str(err))
                
            # looks like things didn't converge, so if N=1 do a linear fit on the last 30% of points to estimate the soln
            if N==1:
            
                m,b = np.polyfit(Es[int(0.7*iter):iter,0], Xs[int(0.7*iter):iter,0], 1)            
                X = np.array([b])
                Y = np.array([0.0]) 
                if display>0: print(f"dsolve is using linear fit to estimate solution at X={b}")
                
            break

        #>>>> COULD ALSO HAVE AN ITERATION RESTART FUNCTION? >>> 
        #  that returns a restart boolean, as well as what values to use to restart things if true. How?
        
        else: 
            dX = step_func(X, args, Y, oths, Ytarget, err, tols, iter, maxIter)
        

        #if display>2:
        #    breakpoint()

        # Make sure we're not diverging by keeping things from reversing too much.
        # Track the previous step (dX_last) and if the current step reverses too much, stop it part way.
        # Stop it at a plane part way between the current X value and the previous X value (using golden ratio, why not).  
        
        # get the point along the previous step vector where we'll draw the bounding hyperplane (could be a line, plane, or more in higher dimensions)
        Xlim = X - 0.62*dX_last
        
        # the equation for the plane we don't want to recross is then sum(X*dX_last) = sum(Xlim*dX_last)
        if np.sum((X+dX)*dX_last) < np.sum(Xlim*dX_last):         # if we cross are going to cross it
            
            alpha = np.sum((Xlim-X)*dX_last)/np.sum(dX*dX_last)    # this is how much we need to scale down dX to land on it rather than cross it
               
            if display > 2:
                print("  limiting oscillation with alpha="+str(alpha))
                print(f"   dX_last was {dX_last}, dX was going to be {dX}, now it'll be {alpha*dX}")
            
            dX = alpha*dX  # scale down dX
            
        # also avoid extreme accelerations in the same direction        
        for i in range(N):
            
            # should update the following for ytol >>>
            if abs(dX_last[i]) > tols[i]:                           # only worry about accelerations if the last step was non-negligible
        
                dX_max = a_max*dX_last[i]                           # set the maximum permissible dx in each direction based an an acceleration limit
                
                if dX_max == 0.0:                                   # avoid a divide-by-zero case (if dX[i] was zero to start with)
                    breakpoint()
                    dX[i] = 0.0                     
                else:    
                    a_i = dX[i]/dX_max                              # calculate ratio of desired dx to max dx
              
                    if a_i > 1.0:
                    
                        if display > 2:
                            print(f"    limiting acceleration ({1.0/a_i:6.4f}) for axis {i}")
                            print(f"     dX_last was {dX_last}, dX was going to be {dX}")
                        
                        #dX = dX*a_max/a_i  # scale it down to the maximum value
                        dX[i] = dX[i]/a_i  # scale it down to the maximum value (treat each DOF individually)
                        
                        if display > 2:
                            print(f"     now dX will be {dX}")
        
        dXlist[iter,:] = dX
        #if iter==196:
            #breakpoint() 
        
        
        # add damping if cyclic behavior is detected at the halfway point
        if dodamping and iter == int(0.5*maxIter):
            if display > 2:   print(f"dsolve2 is at iteration {iter} (50% of maxIter)")
                    
            for j in range(2,iter-1):
                iterc = iter - j
                if all(np.abs(X - Xs[iterc,:]) < tols):
                    print(f"dsolve2 is going in circles detected at iteration {iter}")
                    print(f"last similar point was at iteration {iterc}")
                    damper = damper * 0.9
                    break
                    
        dX = damper*dX
            
            
        # enforce bounds
        for i in range(N):
            
            if X[i] + dX[i] < Xmin[i]:
                dX[i] = Xmin[i] - X[i]
                
            elif X[i] + dX[i] > Xmax[i]:
                dX[i] = Xmax[i] - X[i]

        dXlist2[iter,:] = dX
        # check for convergence
        if (ytol==0 and all(np.abs(dX) < tols)) or (ytol > 0 and all(np.abs(err) < ytol)):
        
            if display>0:
                print("Iteration converged after "+str(iter)+" iterations with error of "+str(err)+" and dX of "+str(dX))
                print("Solution X is "+str(X))
            
                #if abs(err) > 10:
                #    breakpoint()
                
                if display > 0:
                    print("Total run time: {:8.2f} seconds = {:8.2f} minutes".format((time.time() - start_time),((time.time() - start_time)/60)))

            
            if any(X == Xmin) or any(X == Xmax):
                success = False
                print("Warning: dsolve ended on a bound.")
            else:
                success = True
                
            break

        dX_last = 1.0*dX # remember this current value
        
           
        X = X + dX
         
    # truncate empty parts of these arrays
    Xs      = Xs     [:iter+1]
    Es      = Es     [:iter+1]
    dXlist  = dXlist [:iter+1]
    dXlist2 = dXlist2[:iter+1]

    return X, Y, dict(iter=iter, err=err, dX=dX_last, oths=oths, Xs=Xs, Es=Es, success=success, dXlist=dXlist, dXlist2=dXlist2)


def dsolvePlot(info):
    '''Plots dsolve or dsolve solution process based on based dict of dsolve output data'''

    import matplotlib.pyplot as plt

    n = info['Xs'].shape[1]  # number of variables

    if n < 8:
        fig, ax = plt.subplots(2*n, 1, sharex=True)
        for i in range(n):
            ax[  i].plot(info['Xs'][:info['iter']+1,i])
            ax[n+i].plot(info['Es'][:info['iter']+1,i])
        ax[-1].set_xlabel("iteration")
    else:
        fig, ax = plt.subplots(n, 2, sharex=True)
        for i in range(n):
            ax[i,0].plot(info['Xs'][:info['iter']+1,i])
            ax[i,1].plot(info['Es'][:info['iter']+1,i])
        ax[-1,0].set_xlabel("iteration, X")
        ax[-1,1].set_xlabel("iteration, Error")
    plt.show()



def getLineProps(dnommm, material, lineProps=None, source=None, name="", rho=1025.0, g=9.81, **kwargs):
    '''Sets up a dictionary that represents a mooring line type based on the 
    specified diameter and material type. The returned dictionary can serve as
    a MoorPy line type. Data used for determining these properties is a MoorPy
    lineTypes dictionary data structure, created by loadLineProps. This data
    can be passed in via the lineProps parameter, or a new data set can be
    generated based on a YAML filename or dictionary passed in via the source 
    parameter. The lineProps dictionary should be error-checked at creation,
    so it is not error check in this function for efficiency.
        
    Parameters
    ----------
    dnommm : float
        nominal diameter [mm].
    material : string
        string identifier of the material type be used.
    lineProps : dictionary
        A MoorPy lineProps dictionary data structure containing the property scaling coefficients.
    source : dict or filename (optional)
        YAML file name or dictionary containing line property scaling coefficients
    name : any dict index (optional)
        Identifier for the line type (otherwise will be generated automatically).
    rho : float (optional)
        Water density used for computing apparent (wet) weight [kg/m^3].
    g : float (optional)
        Gravitational constant used for computing weight [m/s^2].
    '''
    
    if lineProps==None and source==None:
        raise Exception("Either lineProps or source keyword arguments must be provided")
    
    # deal with the source (is it a dictionary, or reading in a new yaml?)
    if not source==None:
        if not lineProps==None:
            print('Warning: both lineProps and source arguments were passed to getLineProps. lineProps will be ignored.')
        lineProps = loadLineProps(source)
        
    # raise an error if the material isn't in the source dictionary
    if not material in lineProps:
        raise ValueError(f'Specified mooring line material, {material}, is not in the database.')
    
    # calculate the relevant properties for this specific line type
    mat = lineProps[material]       # shorthand for the sub-dictionary of properties for the material in question    
    d = dnommm*0.001                # convert nominal diameter from mm to m      
    mass = mat['mass_d2']*d**2
    MBL  = mat[ 'MBL_0'] + mat[ 'MBL_d']*d + mat[ 'MBL_d2']*d**2 + mat[ 'MBL_d3']*d**3 
    EA   = mat[  'EA_0'] + mat[  'EA_d']*d + mat[  'EA_d2']*d**2 + mat[  'EA_d3']*d**3 + mat['EA_MBL']*MBL 
    cost =(mat['cost_0'] + mat['cost_d']*d + mat['cost_d2']*d**2 + mat['cost_d3']*d**3 
                         + mat['cost_mass']*mass + mat['cost_EA']*EA + mat['cost_MBL']*MBL)
    # add in drag and added mass coefficients if available, if not, use defaults
    if 'Cd' in mat:
        Cd   = mat['Cd']
    else:
        Cd = 1.2
    if 'Cd_ax' in mat:
        CdAx = mat['Cd_ax']
    else:
        CdAx = 0.2
    if 'Ca' in mat:
        Ca = mat['Ca']
    else:
        Ca = 1.0
    if 'Ca_ax' in mat:
        CaAx = mat['Ca_ax']
    else:
        CaAx = 0.0
        
    # internally calculate the volumetric diameter using a ratio
    d_vol = mat['dvol_dnom']*d  # [m]

    # use the volumetric diameter to calculate the apparent weight per unit length 
    w = (mass - np.pi/4*d_vol**2 *rho)*g
    
    # stiffness values for viscoelastic approach 
    EAd = mat['EAd_MBL']*MBL     # dynamic stiffness constant: Krd alpha term x MBL [N]
    EAd_Lm = mat['EAd_MBL_Lm']   # dynamic stiffness Lm slope: Krd beta term (to be multiplied by mean load) [-]
    
    # Set up a main identifier for the linetype unless one is provided
    if name=="":
        typestring = f"{material}{dnommm:.0f}"  # note: previously was type instead of material, undefined
    else:
        typestring = name
    
    notes = f"made with getLineProps"

    lineType = dict(name=typestring, d_vol=d_vol, m=mass, EA=EA, w=w,
                    MBL=MBL, EAd=EAd, EAd_Lm=EAd_Lm, input_d=d,
                    cost=cost, notes=notes, material=material, 
                    Cd=Cd, CdAx=CdAx, Ca=Ca, CaAx=CaAx)
    
    lineType.update(kwargs)   # add any custom arguments provided in the call to the lineType's dictionary
          
    return lineType


def loadLineProps(source):
    '''Loads a set of MoorPy mooring line property scaling coefficients from
    a specified YAML file or passed dictionary. Any coefficients not included
    will take a default value (zero for everything except diameter ratio, 
    which is 1). It returns a dictionary containing the complete mooring line
    property scaling coefficient set to use for any provided mooring line types.
    
    Parameters
    ----------
    source : dict or filename
        YAML file name or dictionary containing line property scaling coefficients
    
    Returns
    -------
    dictionary
        LineProps dictionary listing each supported mooring line type and 
        subdictionaries of scaling coefficients for each.
    '''

    if type(source) is dict:
        source = source
        
    elif source is None or source=="default":
        import os
        mpdir = os.path.dirname(os.path.realpath(__file__))
        with open(os.path.join(mpdir,"MoorProps_default.yaml")) as file:
            source = yaml.load(file, Loader=yaml.FullLoader)
        
    elif type(source) is str:
        with open(source) as file:
            source = yaml.load(file, Loader=yaml.FullLoader)

    else:
        raise Exception("loadLineProps supplied with invalid source")

    if 'lineProps' in source:
        lineProps = source['lineProps']
    else:
        raise Exception("YAML file or dictionary must have a 'lineProps' field containing the data")

    
    output = dict()  # output dictionary combining default values with loaded coefficients
    
    # combine loaded coefficients and default values into dictionary that will be saved for each material
    for mat, props in lineProps.items():  
        output[mat] = {}
        output[mat]['mass_d2'  ] = getFromDict(props, 'mass_d2')  # mass must scale with d^2
        output[mat]['EA_0'     ] = getFromDict(props, 'EA_0'     , default=0.0)
        output[mat]['EA_d'     ] = getFromDict(props, 'EA_d'     , default=0.0)
        output[mat]['EA_d2'    ] = getFromDict(props, 'EA_d2'    , default=0.0)
        output[mat]['EA_d3'    ] = getFromDict(props, 'EA_d3'    , default=0.0)
        output[mat]['EA_MBL'   ] = getFromDict(props, 'EA_MBL'   , default=0.0)
        output[mat]['EAd_MBL'  ] = getFromDict(props, 'EAd_MBL'  , default=0.0)
        output[mat]['EAd_MBL_Lm']= getFromDict(props, 'EAd_MBL_Lm',default=0.0)
        output[mat]['Cd'       ] = getFromDict(props, 'Cd'       , default=0.0)
        output[mat]['Cd_ax'    ] = getFromDict(props, 'Cd_ax'    , default=0.0)
        output[mat]['Ca'       ] = getFromDict(props, 'Ca'       , default=0.0)
        output[mat]['Ca_ax'    ] = getFromDict(props, 'Ca_ax'    , default=0.0)
        
        output[mat]['MBL_0'    ] = getFromDict(props, 'MBL_0'    , default=0.0)
        output[mat]['MBL_d'    ] = getFromDict(props, 'MBL_d'    , default=0.0)
        output[mat]['MBL_d2'   ] = getFromDict(props, 'MBL_d2'   , default=0.0)
        output[mat]['MBL_d3'   ] = getFromDict(props, 'MBL_d3'   , default=0.0)
        output[mat]['dvol_dnom'] = getFromDict(props, 'dvol_dnom', default=1.0)

        # special handling if material density is provided
        if 'density' in props:
            if 'dvol_dnom' in props:
                raise ValueError("Only one parameter can be specified to calculate the volumetric diameter. Choose either 'dvol_dnom' or 'density'.")
            else:
                mass_d2 = output[mat]['mass_d2']
                material_density = getFromDict(props, 'density')
                output[mat]['dvol_dnom'] = np.sqrt((mass_d2/material_density)*(4/np.pi))
        
        # cost coefficients
        output[mat]['cost_0'   ] = getFromDict(props, 'cost_0'   , default=0.0)
        output[mat]['cost_d'   ] = getFromDict(props, 'cost_d'   , default=0.0)
        output[mat]['cost_d2'  ] = getFromDict(props, 'cost_d2'  , default=0.0)
        output[mat]['cost_d3'  ] = getFromDict(props, 'cost_d3'  , default=0.0)
        output[mat]['cost_mass'] = getFromDict(props, 'cost_mass', default=0.0)
        output[mat]['cost_EA'  ] = getFromDict(props, 'cost_EA'  , default=0.0)
        output[mat]['cost_MBL' ] = getFromDict(props, 'cost_MBL' , default=0.0)

    return output




def getPointProps(weight, rho=1025.0, g=9.81, **kwargs):
    '''for now this is just getClumpMV put in a place where it could grow 
    into a fully versatile equivalent to getMoorProps.
    '''
    
    '''A function to provide a consistent scheme for converting a clump weight/float magnitude to the 
    mass and volume to use in a MoorPy Point.'''
    
    if weight >= 0:                          # if the top point of the intermediate line has a clump weight
        pointvol = 0.0
        pointmass = weight*1000.0           # input variables are in units of tons (1000 kg), convert to kg
    else:
        pointvol = -weight*1200.0/rho  # input variables are still in tons. Assume additional 20% of BM mass
        pointmass = -weight*200.0

    return dict(m=pointmass, v=pointvol)


def loadPointProps(source):
    '''Loads a set of MoorPy point property scaling coefficients from
    a specified YAML file or passed dictionary. 
    
    Parameters
    ----------
    source : dict or filename
        YAML file name or dictionary containing line property scaling coefficients
    
    Returns
    -------
    dictionary
        PointProps dictionary listing each supported mooring line type and 
        subdictionaries of scaling coefficients for each.
    '''
    
    '''
    if type(source) is dict:
        source = source
        
    elif source is None or source=="default":
        import os
        mpdir = os.path.dirname(os.path.realpath(__file__))
        with open(os.path.join(mpdir,"PointProps_default.yaml")) as file:
            source = yaml.load(file, Loader=yaml.FullLoader)
        
    elif type(source) is str:
        with open(source) as file:
            source = yaml.load(file, Loader=yaml.FullLoader)

    else:
        raise Exception("loadLineProps supplied with invalid source")
    
    if 'lineProps' in source:
        lineProps = source['lineProps']
    else:
        raise Exception("YAML file or dictionary must have a 'lineProps' field containing the data")
    '''
    
    output = dict()  # output dictionary combining default values with loaded coefficients
    
    #output['generic'] = dict(rho = , m_v = , )
    
    return output



def getFromDict(dict, key, shape=0, dtype=float, default=None):
    '''
    Function to streamline getting values from design dictionary from YAML file, including error checking.

    Parameters
    ----------
    dict : dict
        the dictionary
    key : string
        the key in the dictionary
    shape : list, optional
        The desired shape of the output. If not provided, assuming scalar output. If -1, any input shape is used.
    dtype : type
        Must be a python type than can serve as a function to format the input value to the right type.
    default : number, optional
        The default value to fill in if the item isn't in the dictionary. Otherwise will raise error if the key doesn't exist.
    '''
    # in future could support nested keys   if type(key)==list: ...

    if key in dict:
        val = dict[key]                                      # get the value from the dictionary
        if shape==0:                                         # scalar input expected
            if np.isscalar(val):
                return dtype(val)
            else:
                raise ValueError(f"Value for key '{key}' is expected to be a scalar but instead is: {val}")
        elif shape==-1:                                      # any input shape accepted
            if np.isscalar(val):
                return dtype(val)
            else:
                return np.array(val, dtype=dtype)
        else:
            if np.isscalar(val):                             # if a scalar value is provided and we need to produce an array (of any shape)
                return np.tile(dtype(val), shape)

            elif np.isscalar(shape):                         # if expecting a 1D array
                if len(val) == shape:
                    return np.array([dtype(v) for v in val])
                else:
                    raise ValueError(f"Value for key '{key}' is not the expected size of {shape} and is instead: {val}")

            else:                                            # must be expecting a multi-D array
                vala = np.array(val, dtype=dtype)            # make array

                if list(vala.shape) == shape:                      # if provided with the right shape
                    return vala
                elif len(shape) > 2:
                    raise ValueError("Function getFromDict isn't set up for shapes larger than 2 dimensions")
                elif vala.ndim==1 and len(vala)==shape[1]:   # if we expect an MxN array, and an array of size N is provided, tile it M times
                    return np.tile(vala, [shape[0], 1] )
                else:
                    raise ValueError(f"Value for key '{key}' is not a compatible size for target size of {shape} and is instead: {val}")

    else:
        if default == None:
            raise ValueError(f"Key '{key}' not found in input file...")
        else:
            if shape==0 or shape==-1:
                return default
            else:
                return np.tile(default, shape)


def addToDict(dict1, dict2, key1, key2, default=None):
    '''
    Function to streamline getting values from one dictionary and 
    putting them in another dictionary (potentially under a different key),
    including error checking.

    Parameters
    ----------
    dict1 : dict
        the input dictionary
    dict2 : dict
        the output dictionary
    key1 : string
        the key in the input dictionary
    key2 : string
        the key in the output dictionary
    default : number, optional
        The default value to fill in if the item isn't in the input dictionary.
        Otherwise will raise error if the key doesn't exist.
    '''
    
    if key1 in dict1:
        val = dict1[key1]  
    else:
        if default == None:
            raise ValueError(f"Key '{key1}' not found in input dictionary...")
        else:
            val = default
    
    dict2[key2] = val


def drawBox(ax, r1, r2, color=[0,0,0,0.2]):
    '''Draw a box along the x-y-z axes between two provided corner points.'''
    
    
    ax.plot([r1[0], r2[0]], [r1[1], r1[1]], [r1[2], r1[2]], color=color) # along x
    ax.plot([r1[0], r2[0]], [r2[1], r2[1]], [r1[2], r1[2]], color=color)
    ax.plot([r1[0], r2[0]], [r1[1], r1[1]], [r2[2], r2[2]], color=color)
    ax.plot([r1[0], r2[0]], [r2[1], r2[1]], [r2[2], r2[2]], color=color)
    ax.plot([r1[0], r1[0]], [r1[1], r2[1]], [r1[2], r1[2]], color=color) # along y
    ax.plot([r2[0], r2[0]], [r1[1], r2[1]], [r1[2], r1[2]], color=color)
    ax.plot([r1[0], r1[0]], [r1[1], r2[1]], [r2[2], r2[2]], color=color)
    ax.plot([r2[0], r2[0]], [r1[1], r2[1]], [r2[2], r2[2]], color=color)
    ax.plot([r1[0], r1[0]], [r1[1], r1[1]], [r1[2], r2[2]], color=color) # along z
    ax.plot([r1[0], r1[0]], [r2[1], r2[1]], [r1[2], r2[2]], color=color)
    ax.plot([r2[0], r2[0]], [r1[1], r1[1]], [r1[2], r2[2]], color=color)
    ax.plot([r2[0], r2[0]], [r2[1], r2[1]], [r1[2], r2[2]], color=color)


def makeTower(twrH, twrRad):
    '''Sets up mesh points for visualizing a cylindrical structure (should align with RAFT eventually.'''
    
    n = 8
    X = []
    Y = []
    Z = []
    ax=np.zeros(n+1)
    ay=np.zeros(n+1)
    for jj in range(n+1):
        ax[jj] = np.cos(float(jj)/float(n)*2.0*np.pi)
        ay[jj] = np.sin(float(jj)/float(n)*2.0*np.pi)
        
    for ii in range(int(len(twrRad)-1)):
        z0 = twrH*float(ii)/float(len(twrRad)-1)
        z1 = twrH*float(ii+1)/float(len(twrRad)-1)
        for jj in range(n+1):
            X.append(twrRad[ii]*ax[jj])
            Y.append(twrRad[ii]*ay[jj])
            Z.append(z0)            
            X.append(twrRad[ii+1]*ax[jj])
            Y.append(twrRad[ii+1]*ay[jj])
            Z.append(z1)
    
    Xs = np.array(X)
    Ys = np.array(Y)
    Zs = np.array(Z)    
    
    return Xs, Ys, Zs

def lines2ss(ms):
    '''
    This function automatically detects multi-segmented 
    mooring lines in the MoorPy system object, convert 
    them into subsystems, and updates the MoorPy system. 
    It detects whether the line is suspended or anchored,
    as well as whether it is in the right order (if not
    it will re-oreder).

    Parameters
    ----------
    ms (object): 
        MoorPy system object

    Returns
    ----------
    ms : object
        an updated MoorPy system object with the replaced 
        multi-segmented mooring lines with subsystems.

    '''
    
    i = 0 
    while True:
        subsys_line_id = []
        subsys_point_id = []
        line_ID_of_interest = []
        point_ID_of_interest = []
        pointi = ms.pointList[i]
        if len(pointi.attached) > 2:
            raise ValueError("f point number {pointi.number} branches out.")
        # 1) define the connected lines if any
        subsys_line_id.append(pointi.attached[0])
        subsys_point_id.append(pointi.number)
        # 2) check where the line with line_ID has been repeated in other points
        while True:
            for line_id in subsys_line_id:
                for pointj in ms.pointList:
                    if line_id in pointj.attached:
                        line_ID_of_interest.append(pointj.attached)
                        point_ID_of_interest.append(pointj.number)
                        # if len(pointj.attached) > 2:  # this is the case where we end the subsystem chain if the subsystem line is branching
                            # continue
            old_subsys_line = subsys_line_id
            old_subsys_point = subsys_point_id
            # 3) get the unique values
            subsys_line_id = np.unique(np.concatenate(line_ID_of_interest))
            subsys_point_id = np.unique(point_ID_of_interest)
            if len(subsys_line_id) == len(old_subsys_line) and len(subsys_point_id) == len(old_subsys_point):
                break

        # 4) check if the subsystem is at its lowest state (line and two points), in that case, move to the next point.
        if len(subsys_line_id) == 1 and len(subsys_point_id) == 2:
            i += 1
            if i >= len(ms.pointList):
                break
            continue
        # 5) define the case for this subsys: (case=0: anchored, case=1: suspended)
        ends_z = []
        for pointk in subsys_point_id:
            ends_z.append(ms.pointList[pointk - 1].r[-1])

        dist_to_seabed = np.abs(np.min(ends_z)) - ms.depth  # distance from the lowest point to the seabed
        # if the abs(distance) is below 20%, and it is of fixed type, we consider it anchored (seabed can be variant and anchor might not be at exact ms.seabed)
        if np.abs(dist_to_seabed) < 0.2 * ms.depth and ms.pointList[np.argmin(ends_z)].type == 1:  
            case = 0
        else:
            case = 1  # suspended

        # 6) rearrange lines.
        if case==0:
            anchored_point = ms.pointList[subsys_point_id[np.argmin(ends_z)]-1].r
            # find the distance between the anchored point and the middle of each line:
            dist_to_anchor = []

            for line_id in subsys_line_id:
                rA = ms.lineList[line_id-1].rA
                rB = ms.lineList[line_id-1].rB
                rAdist_to_anchor = np.linalg.norm(rA - anchored_point)
                rBdist_to_anchor = np.linalg.norm(rB - anchored_point)
                dist_to_anchor.append(np.mean([rAdist_to_anchor, rBdist_to_anchor]))
                if rAdist_to_anchor > rBdist_to_anchor:  # Find a way to switch rA and rB so that rA is the point closer to the anchor
                    pass

            subsys_lines = subsys_line_id[np.argsort(dist_to_anchor)] - 1

            # find the distance between the anchored point and the mooring points
            pdist_to_anchor = []
            for point_id in subsys_point_id:
                pdist_to_anchor.append(np.linalg.norm(ms.pointList[point_id - 1].r - anchored_point))

            subsys_points = subsys_point_id[np.argsort(pdist_to_anchor)] - 1
        else:
            subsys_lines = subsys_line_id - 1 # no need to rearrange because it could work from either end
            subsys_points = subsys_point_id - 1

        
        lines = list(subsys_lines)
        points = list(subsys_points)
        ms = lines2subsystem(lines, ms, span=None, case=case)
        ms.initialize()
        ms.solveEquilibrium()
        i += 1
        if i >= len(ms.pointList):
            break

    return ms

def lines2subsystem(lines,ms,span=None,case=0):
    '''Takes a set of connected lines (in order from rA to rB) in a moorpy system and creates a subsystem equivalent.
    The original set of lines are then removed from the moorpy system and replaced with the 
    subsystem.

    Parameters
    ----------
    lines : list
        List of indices in the ms.lineList to replace.
    ms : object
        MoorPy system object the lines are part of
    span : float (optional)
        Span of the total line (from start to end of subsystem)
    case : int (optional)
        0 = end A on seabed
        1 = suspended line with end A at another floater
        2 = suspended line is symmetric, end A is assumed the midpoint

    Returns
    -------
    ms : object
        MoorPy system object with new subsystem line

    '''
    from moorpy.subsystem import Subsystem
    from copy import deepcopy
    # save a deepcopy of the line list to delete
    originalList = deepcopy(lines)
    
    # # check that all lines connect (all are sections of one full mooring line)
    # for i in range(0,len(lines)):
    #     if i>0:
    #         if not all(b==0 for b in ms.lineList[lines[i]].rB == ms.lineList[lines[i-1]].rA):
    #             raise Exception('Lines indices must be provided in order from rA to rB.')
    # get the span of the subsystem line
    if not span:
        span = np.sqrt((ms.lineList[lines[0]].rA[0]-ms.lineList[lines[-1]].rB[0])**2+(ms.lineList[lines[0]].rA[1]-ms.lineList[lines[-1]].rB[1])**2)
    # make a subsystem object
    ss = Subsystem(depth=ms.depth, span=span,rBFair=ms.lineList[lines[-1]].rB)
    lengths = []
    types = []
    pt = [] # list of points that 
    # ptB = []

    # go through each line
    for i in lines:
        # see which point is connected to end A and end B
        for j in range(0,len(ms.pointList)):
            for k in range(0,len(ms.pointList[j].attached)):
                if ms.pointList[j].attached[k] == i+1: #and ms.pointList[j].attachedEndB[k] == 0:
                    if not j in pt:
                        pt.append(j)
                # elif ms.pointList[j].attached[k] == i+1 and ms.pointList[j].attachedEndB[k] == 1:
                #     ptB.append(j)

        # collect line lengths and types                                          
        lengths.append(ms.lineList[i].L)
        types.append(ms.lineList[i].type['name'])
        ss.lineTypes[types[-1]] = ms.lineTypes[types[-1]]
        
    # use makeGeneric to build the subsystem line
    ss.makeGeneric(lengths,types,suspended=case)
    ss.setEndPosition(ms.lineList[lines[0]].rA,endB=0)
    ss.setEndPosition(ms.lineList[lines[-1]].rB,endB=1)
    
    # add in any info on the points connected to the lines
    # currently, mass, volume, and diameter but others could be added
    for i in range(0,len(pt)):
        ss.pointList[i].m = ms.pointList[pt[i]].m
        ss.pointList[i].v = ms.pointList[pt[i]].v
        ss.pointList[i].d = ms.pointList[pt[i]].d
        # ss.pointList[i].m = ms.pointList[ptB[i]].m
        # ss.pointList[i].v = ms.pointList[ptB[i]].v
        # ss.pointList[i].d = ms.pointList[ptB[i]].d
    
    from moorpy import helpers
    # delete old line
    for i in range(0,len(lines)):
        decB = 0 # boolean to check if ptB has been decreased already for this line
        decA = 0 # boolean to check if ptA has been decreased already for this line
        if i == 0 and i < len(lines) - 1:
            # first line of multiple (keep only point A)
            delpts = 2
            for j in range(0,len(ms.pointList)):
                if lines[i]+1 in ms.pointList[j].attached:
                    if pt[-1]>j and decB == 0:                    
                        pt[-1] -= 1
                        decB = 1
                    if pt[0]>j and decA == 0:
                        pt[0] -= 1 
                        decA = 1
        elif i == 0 and i == len(lines) - 1:
            # first line, only line (keep point A and B)
            delpts = 0 
        elif i == len(lines) - 1:
            # last line, keep point A because already been deleted, and keep point B (fairlead)
            delpts = 0
        else:
            # not beginning or end line, point A (previous line pointB) will have already been deleted so don't delete point A
            delpts = 2
            # reduce index of last point B in ptB list and first point A in ptA list (only care about last ptB and first ptA now) by one
            for j in range(0,len(ms.pointList)):
                if lines[i]+1 in ms.pointList[j].attached:
                    if pt[-1]>j and decB == 0:                    
                        pt[-1] -= 1
                        decB = 1
                    if pt[0]>j and decA == 0:
                        pt[0] -= 1 
                        decA = 1
        # adjust index of any lines that have a higher index than the line to delete
        for j in range(0,len(lines)):
            if lines[i]<lines[j]:
                lines[j] -= 1
        # delete old line and any necessary points
        helpers.deleteLine(ms,lines[i],delpts)
    
    # print('Replacing lines ',originalList,' with a subsystem appended to the end of the lineList ')
    # append subsystem to ms
    ms.lineList.append(ss)
    ssNum = len(ms.lineList)
    # attach subystem line to the end points
    ms.pointList[pt[0]].attachLine(ssNum,0) # rA
    ms.pointList[pt[-1]].attachLine(ssNum,1) # rB     
        
    return(ms)

def deleteLine(ms,ln,delpts=0):
    '''
    Deletes a line from the linelist, and updates the points to have the correct line
    index listed as attached. If delpts=True, also deletes all points associated with 
    that line.

    Parameters
    ----------
    ms : system object
        MoorPy system object the line is within
    ln : int
        Line index number to delete from ms.lineList
    delpts : int
        Set to 0 to keep all points associated with deleted line, 1 to delete only the point A of the line,
        2 to delete only point B of the line, and 3 to delete all points on the line

    Returns
    -------
    None.

    '''
    # delete line
    ms.lineList.pop(ln)
    # reduce line.number for any line after deleted line
    for i in range(0,len(ms.lineList)):
        if ms.lineList[i].number > ln + 1:
            ms.lineList[i].number = ms.lineList[i].number - 1
    
    # adjust attached line index number in any point that is attached to a line index after the deleted line index
    numpts = len(ms.pointList)
    i=0
    reset = 0
    while i < numpts:
        j = 0
        while j < len(ms.pointList[i].attached):
            Bbool = ms.pointList[i].attachedEndB[j]
            reset = 0 # turn off boolean to reset i
            if ms.pointList[i].attached[j] > ln+1 :
                ms.pointList[i].attached[j] = ms.pointList[i].attached[j] - 1 # new index will be one less
            # delete points if wanted
            elif ms.pointList[i].attached[j] == ln+1:
                # remove line number from attached list
                ms.pointList[i].attached.pop(j)
                if delpts == 0:
                    # keep the point but delete the line from the attachedEndB list
                    ms.pointList[i].attachedEndB.pop(j)
                if delpts == 1 or delpts == 3:
                    if Bbool == 0:
                    # if ms.pointList[i].attachedEndB[j] == 0:  
                        # reduce number of times through the loop
                        numpts = numpts-1
                        # reduce point.number for each point after deleted point
                        for k in range(0,len(ms.pointList)):
                            if ms.pointList[k].number > i + 1:
                                ms.pointList[k].number = ms.pointList[k].number - 1
                        # lower index of any body attached points after deleted point, remove deleted point from body attached points
                        for k in range(0,len(ms.bodyList)):
                            ii = 0
                            while ii < len(ms.bodyList[k].attachedP):
                                if ms.bodyList[k].attachedP[ii] == i+1 :
                                    # remove point
                                    ms.bodyList[k].attachedP.pop(ii)
                                    # remove relative points from list
                                    ms.bodyList[k].rPointRel.pop(ii)
                                    ii = ii - 1 # reduce iter because attachedP[1] is now attachedP[0] so need to get that one
                                elif ms.bodyList[k].attachedP[ii] > i+1 :
                                    # reduce index by one
                                    ms.bodyList[k].attachedP[ii] = ms.bodyList[k].attachedP[ii] - 1
                                ii += 1
                        # remove point
                        ms.pointList.pop(i)
                        # trigger boolean to reset i and j back one (since now point x+1 will be point x)
                        reset = 1
                        
                if delpts == 2 or delpts == 3:
                    if Bbool == 1:
                    # if ms.pointList[i].attachedEndB[j] == 1:
                        # reduce number of times through the loop
                        numpts = numpts-1
                        # reduce point.number for each point after deleted point
                        for k in range(0,len(ms.pointList)):
                            if ms.pointList[k].number > i + 1:
                                ms.pointList[k].number = ms.pointList[k].number - 1
                        # lower index of any body attached points after deleted point, remove deleted point from body attached points
                        for k in range(0,len(ms.bodyList)):
                            ii = 0
                            while ii < len(ms.bodyList[k].attachedP):
                                if ms.bodyList[k].attachedP[ii] == i+1 :
                                    # remove relative points from list
                                    ms.bodyList[k].rPointRel.pop(ii)
                                    # remove point
                                    ms.bodyList[k].attachedP.pop(ii)
                                    ii = ii - 1 # reduce iter because attachedP[1] is now attachedP[0] so need to get that one
                                elif ms.bodyList[k].attachedP[ii] > i+1 :
                                    # reduce index by one
                                    ms.bodyList[k].attachedP[ii] = ms.bodyList[k].attachedP[ii] - 1
                                ii += 1
                        # remove point
                        ms.pointList.pop(i)
                        # trigger boolean to reset i and j back one (since now point x+1 will be point x)
                        reset = 1
                        
            j += 1
            if reset:
                j -= 1
            # need to get out of inner loop if the last point in the list was deleted (for statement for inner loop would throw an error otherwise)
            if i >= len(ms.pointList):
                break
        # reset i if any points were removed
        if reset:
            i -= 1
        # increment i
        i += 1
    return(ms)
                    
def subsystem2Line(ms,ssNum,nsegs=10):
    '''Replace a subsystem with equivalent set of lines

    Parameters
    ----------
    ms : system object
        MoorPy system object that contains the subsystem
    ssNum : int
        index in the lineList of ms which points to the subsystem object to replace
    nsegs : list OR int, optional
        Number of segments per line for each line. Can be an integer (all line sections have the same # of segments)
        OR can be a list (# of segments for each section of line in order from A to B)

    Returns
    -------
    None.

    '''
    # get subsystem object
    ss = ms.lineList[ssNum]   
    types = []
    lengths = []
    points = []
    # record line types, lines, and points in the subsystem
    for i in range(0,len(ss.lineList)):
        types.append(ss.lineList[i].type)
        lengths.append(ss.lineList[i].L)
        if not types[-1]['name'] in ms.lineTypes:
            # add type to lineTypes list
            ms.lineTypes[types[-1]['name']] = types[-1]
    for i,spt in enumerate(ss.pointList):
        # gather all info about the points in the subsystem
        points.append({'r':spt.r,'m':spt.m,'v':spt.v,'CdA':spt.CdA,'d':spt.d,'type':spt.type,'Ca':spt.Ca})
    # points[0]['r'] = ss.rA
    # points[-1]['r'] = ss.rB
        if spt.attachedEndB[-1]:
            endB = i
            points[endB]['r'] = ss.rB
        if spt.attachedEndB[0] == 0:
            endA = i
            points[endA]['r'] = ss.rA
    # get actual r of end points (r in subsystem is not true location)
    for i in range(0,len(ms.pointList)):
        # check if point is attached to the subsystem line
        for j in range(0,len(ms.pointList[i].attached)):
            if ms.pointList[i].attached[j] == ssNum+1:
                if ms.pointList[i].attachedEndB[j]:
                    # for k in range(0,len(ms.bodyList)):
                    #     if i+1 in ms.bodyList[k].attachedP:
                    #         points[-1]['body'] = k    
                    # update end B r
                    points[endB]['r'] = ms.pointList[i].r
                    points[endB]['type'] = ms.pointList[i].type
                    # check if end points are attached to a body
                    for k in range(0,len(ms.bodyList)):
                        if i+1 in ms.bodyList[k].attachedP:
                            points[endB]['body'] = k
                else:
                    # update end A r
                    points[endA]['r'] = ms.pointList[i].r
                    points[endA]['type'] = ms.pointList[i].type 
                    # check if end points are attached to a body
                    for k in range(0,len(ms.bodyList)):
                        if i+1 in ms.bodyList[k].attachedP:
                            points[endA]['body'] = k
    # approximate midpoint r with depth of subsystem point r and angle from two end points
    aang = np.arctan2(points[0]['r'][1] - points[-1]['r'][1],points[0]['r'][0] - points[-1]['r'][0])
    # update x-y location of any midpoints if they exist
    if len(points)>2:
        for i in range(1,len(points)-1):
            ll = np.sqrt(points[i]['r'][0]**2+points[i]['r'][1]**2)
            points[i]['r'][0] = (ll)*np.cos(aang)+points[-1]['r'][0]# poits[-1]['r][0]
            points[i]['r'][1] = (ll)*np.sin(aang)+points[-1]['r'][1]
            #points[i]['r'][0] = (ll+np.abs(points[-1]['r'][0]))*np.cos(aang)
            #points[i]['r'][1] = (ll+np.abs(points[-1]['r'][1]))*np.sin(aang)

    from moorpy import helpers
    # remove subsystem line, delete all associated points
    helpers.deleteLine(ms,ssNum,delpts=3)
    # add in new lines to replace subsystem
    for i in range(0,len(types)):
        # determine # of segments for this line section
        if isinstance(nsegs,list):
            NSegs = nsegs[i]
        elif isinstance(nsegs,int):
            NSegs = nsegs
        else:
            raise Exception('Input nsegs must be either a list or an integer')
        # add point A
        ms.addPoint(points[i]['type'], points[i]['r'], points[i]['m'], points[i]['v'], d=points[i]['d'])
        ms.pointList[-1].CdA = points[i]['CdA']
        ms.pointList[-1].Ca = points[i]['Ca']
        # add line
        ms.addLine(lengths[i], types[i]['name'], nSegs=NSegs)
        # attach new line to its point A
        ms.pointList[-1].attachLine(len(ms.lineList),endB=0) # end A for line just created
        # attach to any bodies the point was originally attached to
        if 'body' in points[i]:
            ms.bodyList[points[i]['body']].attachPoint(len(ms.pointList),[ms.pointList[-1].r[0]-ms.bodyList[points[i]['body']].r6[0],ms.pointList[-1].r[1]-ms.bodyList[points[i]['body']].r6[1],ms.pointList[-1].r[2]])
        if i>0:
            # new point is end B for previous line, attach as point B
            ms.pointList[-1].attachLine(len(ms.lineList)-1,endB=1)
    
    # add last point (should be the fairlead)
    ms.addPoint(points[-1]['type'], points[-1]['r'], points[-1]['m'], points[-1]['v'], d=points[-1]['d'])
    ms.pointList[-1].CdA = points[-1]['CdA']
    ms.pointList[-1].Ca = points[-1]['Ca']
    # attach to last line as point B
    ms.pointList[-1].attachLine(len(ms.lineList),endB=1)
    # attach to a body if applicable
    if 'body' in points[-1]:
        ms.bodyList[points[-1]['body']].attachPoint(len(ms.pointList),[ms.pointList[-1].r[0]-ms.bodyList[points[-1]['body']].r6[0],ms.pointList[-1].r[1]-ms.bodyList[points[-1]['body']].r6[1],ms.pointList[-1].r[2]])
   
            
    
                

def readBathymetryFile(filename):
    '''Read a MoorDyn-style bathymetry input file (rectangular grid of depths)
    and return the lists of x and y coordinates and the matrix of depths.
    '''
    f = open(filename, 'r')

    # skip the header
    line = next(f)
    # collect the number of grid values in the x and y directions from the second and third lines
    line = next(f)
    nGridX = int(line.split()[1])
    line = next(f)
    nGridY = int(line.split()[1])
    # allocate the Xs, Ys, and main bathymetry grid arrays
    bathGrid_Xs = np.zeros(nGridX)
    bathGrid_Ys = np.zeros(nGridY)
    bathGrid = np.zeros([nGridY, nGridX])  # MH swapped order June 30
    # read in the fourth line to the Xs array
    line = next(f)
    bathGrid_Xs = [float(line.split()[i]) for i in range(nGridX)]
    # read in the remaining lines in the file into the Ys array (first entry) and the main bathymetry grid
    for i in range(nGridY):
        line = next(f)
        entries = line.split()
        bathGrid_Ys[i] = entries[0]
        bathGrid[i,:] = entries[1:]
    
    return bathGrid_Xs, bathGrid_Ys, bathGrid


def read_mooring_file(dirName,fileName):
    # Taken from line system.... maybe should be a helper function?
    # load data from time series for single mooring line
    
    print('attempting to load '+dirName+fileName)
    
    f = open(dirName+fileName, 'r')
    
    channels = []
    units = []
    data = []
    i=0
    
    for line in f:          # loop through lines in file
    
        if (i == 0):
            for entry in line.split():      # loop over the elemets, split by whitespace
                channels.append(entry)      # append to the last element of the list
                
        elif (i == 1):
            for entry in line.split():      # loop over the elemets, split by whitespace
                units.append(entry)         # append to the last element of the list
        
        elif len(line.split()) > 0:
            data.append([])  # add a new sublist to the data matrix
            import re
            r = re.compile(r"(?<=\d)\-(?=\d)")  # catch any instances where a large negative exponent has been written with the "E"
            line2 = r.sub("E-",line)            # and add in the E
            
            
            for entry in line2.split():      # loop over the elemets, split by whitespace
                data[-1].append(entry)      # append to the last element of the list
            
        else:
            break
    
        i+=1
    
    f.close()  # close data file
    
    # use a dictionary for convenient access of channel columns (eg. data[t][ch['PtfmPitch'] )
    ch = dict(zip(channels, range(len(channels))))
    
    data2 = np.array(data)
    
    data3 = data2.astype(float)
    
    return data3, ch, channels, units    

def read_output_file(dirName,fileName, skiplines=-1, hasunits=1, chanlim=999, dictionary=True):

    # load data from FAST output file
    # looks for channel names, then units (if hasunits==1), then data lines after first skipping [skiplines] lines.
    # skiplines == -1 signals to search for first channel names line based on starting channel "Time".
    
#   print('attempting to load '+dirName+fileName)
    f = open(dirName+fileName, 'r')
    
    channels = []
    units = []
    data = []
    i=0
    
    for line in f:          # loop through lines in file
    
        if (skiplines == -1):               # special case signalling to search for "Time" at start of channel line
            entries = line.split()          # split elements by whitespace
            print(entries)
            if entries[0].count('Time') > 0 or entries[0].count('time') > 0:  # if we find the time keyword
                skiplines = i
                print("got skiplines="+str(i))
            else:
                pass
    
        if (i < skiplines or skiplines < 0):        # if we haven't gotten to the first channel line or we're in search mode, skip
            pass
            
        elif (i == skiplines):
            for entry in line.split():      # loop over the elemets, split by whitespace
                channels.append(entry)      # append to the last element of the list
                
        elif (i == skiplines+1 and hasunits == 1):
            for entry in line.split():      # loop over the elemets, split by whitespace
                if entry.count('kN') > 0 and entry.count('m') > 0:  # correct for a possible weird character
                    entry = '(kN-m)'
                    
                units.append(entry)         # append to the last element of the list
        
        elif len(line.split()) > 0:
            data.append([])  # add a new sublist to the data matrix
            
            r = re.compile(r"(?<=\d)\-(?=\d)")  # catch any instances where a large negative exponent has been written with the "E"
            line2 = r.sub("E-",line)           # and add in the E
            
            j=0
            for entry in line2.split():      # loop over the elements, split by whitespace
                if j > chanlim:
                    break
                j+=1    
                data[-1].append(entry)      # append to the last element of the list
    
        else:
            break
    
        i+=1
    
    f.close()  # close data file
    
    
    # use a dictionary for convenient access of channel columns (eg. data[t][ch['PtfmPitch'] )
    ch = dict(zip(channels, range(len(channels))))
    
    #print ch['WindVxi']

    data2 = np.array(data)
    
    data3 = data2.astype(float)
    
    if dictionary:
        dataDict = {}
        unitDict = {}
        for i in range(len(channels)):
            dataDict[channels[i]] = data3[:,i]
            unitDict[channels[i]] = units[i]
        return dataDict, unitDict
    else:
        return data3, ch, channels, units