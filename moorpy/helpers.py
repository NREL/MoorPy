
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
    '''Prints the provided matrix with elements separated by tabs.

    Parameters
    ----------
    mat : array
        Any matrix that is to be printed.
    '''
    
    if mat.dtype in [int, bool]:
        for i in range(mat.shape[0]):
            print( "\t".join(["{}"]*mat.shape[1]).format( *mat[i,:] ))
    else:
        for i in range(mat.shape[0]):
            print( "\t".join(["{:+8.3e}"]*mat.shape[1]).format( *mat[i,:] ))


def printVec(vec):
    '''Prints a vector with 4 decimal scientific notation and consisten spaces

    Parameters
    ----------
    vec : array
        Any vector that is to be printed.
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
    
    Returns
    -------
    lineType : dictionary
        A lineType dictionary
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

    # Checking valid diameter ranges for MBL
    if mat['MBL_dmax'] >= 0 and mat['MBL_dmin'] >= 0: # if min and max values given and the diameter outside of the range
        if d > mat['MBL_dmax'] or d < mat['MBL_dmin']:
            raise Exception(f"Input diameter {d} m is outside of the valid range for the MBL curve of {mat['MBL_dmin']}-{mat['MBL_dmax']} m")
    elif mat['MBL_dmax'] >= 0 and d > mat['MBL_dmax']: # if a max value is given and the diameter is greater than the max
        raise Exception(f"Input diameter {d} m is greater than the max valid value for the MBL curve of {mat['MBL_dmax']} m")
    elif mat['MBL_dmin'] >= 0 and d < mat['MBL_dmin']: # if a min value is given and the diameter is less than the min
        raise Exception(f"Input diameter {d} m is less than the min valid value for the MBL curve of {mat['MBL_dmin']} m")
    MBL  = mat[ 'MBL_0'] + mat[ 'MBL_d']*d + mat[ 'MBL_d2']*d**2 + mat[ 'MBL_d3']*d**3 

    # Checking valid diameter ranges fo EA
    if mat['EA_dmax'] >= 0 and mat['EA_dmin'] >= 0: # if min and max values given and the diameter outside of the range
        if d > mat['EA_dmax'] or d < mat['EA_dmin']:
            raise Exception(f"Input diameter {d} m is outside of the valid range for the EA curve of {mat['EA_dmin']}-{mat['EA_dmax']} m")
    elif mat['EA_dmax'] >= 0 and d > mat['EA_dmax']: # if a max value is given and the diameter is greater than the max
        raise Exception(f"Input diameter {d} m is greater than the max valid value for the EA curve of {mat['EA_dmax']} m")
    elif mat['EA_dmin'] >= 0 and d < mat['EA_dmin']: # if a min value is given and the diameter is less than the min
        raise Exception(f"Input diameter {d} m is less than the min valid value for the EA curve of {mat['EA_dmin']} m")
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
                    MBL=MBL, EAd=EAd, EAd_Lm=EAd_Lm, d_nom=d,
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
    output : dictionary
        LineProps dictionary listing each supported mooring line type and 
        subdictionaries of scaling coefficients for each.
    '''

    if type(source) is dict:
        source = source
        
    elif source is None or source=="default":
        import os
        mpdir = os.path.dirname(os.path.realpath(__file__))
        with open(os.path.join(mpdir,"library/MoorProps_default.yaml")) as file:
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
        output[mat]['EA_dmin'  ] = getFromDict(props, 'EA_dmin'  , default=-1.0) # -1 to disable checking
        output[mat]['EA_dmax'  ] = getFromDict(props, 'EA_dmax'  , default=-1.0) # -1 to disable checking
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
        output[mat]['MBL_dmin' ] = getFromDict(props, 'MBL_dmin' , default=-1.0) # -1 to disable checking
        output[mat]['MBL_dmax' ] = getFromDict(props, 'MBL_dmax' , default=-1.0) # -1 to disable checking
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


def getPointProps(design, Props = None, source = None, name="", **kwargs):
    '''Sets up a dictionary that represents a mooring point type based on the 
    specified design. Data used for determining these properties is a MoorPy
    pointTypes dictionary data structure, created by loadPointProps. This data
    can be passed in via the Props parameter, or a new data set can be
    generated based on a YAML filename or dictionary passed in via the source 
    parameter. The pointProps dictionary should be error-checked at creation,
    so it is not error check in this function for efficiency.
    
    Parameters
    ----------
    design : string or dict
        design keyword from DesignProps or dictionary with num_a_<anchor key>, num_b_<buoy key>, num_c_<connect key> entries
    Props : dictionary (optional)
        A MoorPy PointProps dictionary data structure containing designs and property scaling coefficients.
    source : dict or filename (optional)
        YAML file name or dictionary containing designs and property scaling coefficients.
    name : string (optional)
        Identifier for the point type (otherwise wll be generated automatically).

    Returns
    -------
    dictionary
        A pointType dictionary 
    '''

    pointinfo = {} # Dictionary to hold any information statements

    if Props==None  and source==None:
        raise MoorPyError('Either Props or source keyword arguments must be provided')

    # deal with the source (is it a dictionary, or reading in a new yaml?)
    if not source==None:
        if not Props==None :
            print('Warning: both Props and source arguments were passed to getLineProps. Props will be ignored.')
        Props = loadPointProps(source)

    # if design is a keyword, load it from design props, otherwise just use the passed dict
    if type(design) == str: 
        # raise an error if the material isn't in the source dictionary
        if not design in Props["DesignProps"]:
            raise MoorPyError(f'Specified point design, {design}, is not in the database.')
        dtype = Props["DesignProps"][design]       # shorthand for the sub-dictionary of properties for the design in question
    else:
        dtype = design

    # Check that fractions of anchor size and buoy size are realistic (and create them as num * 1/total num object if they dont exist)
    num_aTypes = 0
    num_bTypes = 0
    for key in dtype.keys(): # count total num anchors and buoys
        if 'num_a_' in key:
            num_aTypes += 1
        if 'num_b_' in key:
            num_bTypes += 1

    aFrac_sum = 0.0
    bFrac_sum = 0.0
    key_list = list(dtype.keys()) # becasue we are changing the dict in the loop, we only want to iterate over the original keys
    for key in key_list:
        if 'num_a_' in key:
            if not (f"frac_a_{key[6:]}" in dtype.keys()): # if we dont have frac_a_<key>. This only can happen if design is passed in as dict and not keyword
                # set frac as 1/total num anchs
                dtype[f"frac_a_{key[6:]}"] = 1/num_aTypes

            aFrac_sum += dtype[f"frac_a_{key[6:]}"]

        if 'num_b_' in key:
            if not (f"frac_b_{key[6:]}" in dtype.keys()): # if we dont have frac_a_<key>. This only can happen if design is passed in as dict and not keyword
                # set frac as 1/total num buoys
                dtype[f"frac_b_{key[6:]}"] = 1/num_bTypes

            bFrac_sum += dtype[f"frac_b_{key[6:]}"]

    if aFrac_sum > 1:
        raise MoorPyError(f"Anchor fractions sum to {aFrac_sum}, which is greater than 1") 
    if bFrac_sum > 1:
        raise MoorPyError(f"Buoy fractions sum to {bFrac_sum}, which is greater than 1") 

    # Set up a name identifier for the pointtype unless one is provided
    if name=="":
        if type(design) == str: # if design is a keyword, use that
            typestring = design 
        else: # otherwise :(
            typestring = "no_name"
    else:
        typestring = name

    # initialize list to hold FOS. 0 is default, but this shouldn't ever be used because c_bool will be false if nothing is added to the list
    pointfos_list = [0.0] # FOS

    # initialize pointType outputs
    a_bool = False
    b_bool = False
    c_bool = False
    anchorList = []
    aprops_out = {}
    Bcost = dict(cost_b0 = 0.0, cost_b1 = 0.0, cost_b2 = 0.0, cost_b3 = 0.0)
    Ccost = dict(cost_load0 = 0.0, cost_load1 = 0.0, cost_load2 = 0.0, cost_load3 = 0.0)

    # load ABC props and assign values to the point.

    for key in dtype.keys():
        
        # anchors
        if 'num_a_' in key:
            '''
            anchor stuff is handled in MoorProps and point.getCost_and_MBL, so just agregate general things here
            '''
            a_bool = True
            anchorList.append({"name" : key[6:], "num" : dtype[key], "frac" : dtype[f"frac_a_{key[6:]}"]})
            aprops_out[key[6:]] = Props["AnchorProps"][key[6:]] # Error checking of valid keys is already done in loadPointProps

        # buoys
        if 'num_b_' in key:
            '''
            Number of this buoy type * the cost coefficient to get our general cost coefficients for the point
            Each buoy has a buoyancy of (total buoyancy of point * % of buoyancy). 

            Note: 
            This assumes if there are multiple buoys of the same type they share the same % of buoyancy. I.e. frac_b_key / num_b_key

            Note:
            Num of type * (unit cost) = Num of type * (B1 * (buoyancy * % of buoyancy ) + B2 * (buoyancy * % of buoyancy )^2 + ...) = Num of type * B1 * % of buoyancy * buoyancy + Num of type * B2 * % of buoyancy^2 * buoyancy^2 + ...
            '''

            b_bool = True
            Bcost["cost_b0"] += dtype[key] * Props["BuoyProps"][key[6:]]['cost_b0']                                                # Error checking of valid keys is already done in loadPointProps
            Bcost["cost_b1"] += dtype[key] * Props["BuoyProps"][key[6:]]['cost_b1'] * (dtype[f"frac_b_{key[6:]}"] / dtype[key])    # Error checking of valid keys is already done in loadPointProps
            Bcost["cost_b2"] += dtype[key] * Props["BuoyProps"][key[6:]]['cost_b2'] * (dtype[f"frac_b_{key[6:]}"] / dtype[key])**2 # Error checking of valid keys is already done in loadPointProps
            Bcost["cost_b3"] += dtype[key] * Props["BuoyProps"][key[6:]]['cost_b3'] * (dtype[f"frac_b_{key[6:]}"] / dtype[key])**3 # Error checking of valid keys is already done in loadPointProps

        # connections
        if 'num_c_' in key:
            '''
            Number of this connect type * the cost coefficient to get our general cost coefficients for the point

            Note:
            Num * (C1 * (Load * FOS) + C2 * (Load * FOS)**2 + ...) = Num * C1 * FOS * Load + Num * C2 * FOS^2 * Load^2 + ...
            '''
            c_bool = True
            FOS = Props["ConnectProps"][key[6:]]["FOS"]
            Ccost["cost_load0"] += dtype[key] * Props["ConnectProps"][key[6:]]['cost_MBL0']          # Error checking of valid keys is already done in loadPointProps
            Ccost["cost_load1"] += dtype[key] * Props["ConnectProps"][key[6:]]['cost_MBL1'] * FOS    # Error checking of valid keys is already done in loadPointProps
            Ccost["cost_load2"] += dtype[key] * Props["ConnectProps"][key[6:]]['cost_MBL2'] * FOS**2 # Error checking of valid keys is already done in loadPointProps
            Ccost["cost_load3"] += dtype[key] * Props["ConnectProps"][key[6:]]['cost_MBL3'] * FOS**3 # Error checking of valid keys is already done in loadPointProps
            pointfos_list.append(FOS)

    pointType = dict(name=typestring, Anchors = a_bool, Buoys = b_bool, Connections = c_bool, 
                     anchor_list = anchorList, aprops = aprops_out, buoy_cost = Bcost, 
                     connector_cost = Ccost, FOS = np.max(pointfos_list), info=pointinfo)

    pointType.update(kwargs)   # add any custom arguments provided in the call to the pointType's dictionary

    return pointType


def loadPointProps(source):
    '''Loads a set of MoorPy point property scaling coefficients from
    a specified YAML file or passed dictionary. Any coefficients not included will take a 
    default value. It returns two dictionaries containing the 
    complete point design descriptions and the connection component scaling 
    coefficient set to use for any provided point data types.
    
    Parameters
    ----------
    source : dict or filename
        YAML file name or dictionary containing point property dictionaries of designs and connection hardware
    
    Returns
    -------
    output : dictionary
        PointProps dictionary listing each supported connection design type and 
        subdictionaries of anchor, buoy, and connection properties.
    '''

    # Note that the assumption that connection hardware matches line MBL is not necessarily correct. Typical hardware FOS are 5:1 (tri-plates) or 6:1 (shackles). 
    # A more accurate approach would find nominal size to MBL coefficients for each type, and find a relation between attached line diam and nominal size.

    if type(source) is dict:
        source = source
        
    elif source is None or source=="default":
        import os
        mpdir = os.path.dirname(os.path.realpath(__file__))
        with open(os.path.join(mpdir,"library/PointProps_default.yaml")) as file:
            source = yaml.load(file, Loader=yaml.FullLoader)
        
    elif type(source) is str:
        with open(source) as file:
            source = yaml.load(file, Loader=yaml.FullLoader)

    else:
        raise MoorPyError("loadPointProps supplied with invalid source")
    
    # read the Anchor dict
    if 'AnchorProps' in source:
        AnchorProps = source['AnchorProps']
    else:
        raise MoorPyError("YAML file or dictionary must have a 'AnchorProps' field containing the data") # We might want to remove this requirement becasue if a user just gives point mass volume and cost they wont need a anchor dict
    
    anchor = dict()  # output dictionary combining default values with loaded coefficients
    
     # combine loaded coefficients and default values into dictionary that will be saved for each material
    for atype, props in AnchorProps.items():  
        anchor[atype] = {}

        anchor[atype]['matcost_m'   ] = getFromDict(props, 'matcost_m'   , default=0.0)
        anchor[atype]['matcost_m2'  ] = getFromDict(props, 'matcost_m2'  , default=0.0)
        anchor[atype]['matcost_m3'  ] = getFromDict(props, 'matcost_m3'  , default=0.0)
        anchor[atype]['matcost_a'   ] = getFromDict(props, 'matcost_a'   , default=0.0)
        anchor[atype]['matcost_a2'  ] = getFromDict(props, 'matcost_a2'  , default=0.0)
        anchor[atype]['matcost_a3'  ] = getFromDict(props, 'matcost_a3'  , default=0.0)
        anchor[atype]['matcost'     ] = getFromDict(props, 'matcost'     , default=0.0)
        anchor[atype]['instcost_m'  ] = getFromDict(props, 'instcost_m'  , default=0.0)
        anchor[atype]['instcost_m2' ] = getFromDict(props, 'instcost_m2' , default=0.0)
        anchor[atype]['instcost_m3' ] = getFromDict(props, 'instcost_m3' , default=0.0)
        anchor[atype]['instcost_a'  ] = getFromDict(props, 'instcost_a'  , default=0.0)
        anchor[atype]['instcost_a2' ] = getFromDict(props, 'instcost_a2' , default=0.0)
        anchor[atype]['instcost_a3' ] = getFromDict(props, 'instcost_a3' , default=0.0)
        anchor[atype]['instcost'    ] = getFromDict(props, 'instcost'    , default=0.0)
        anchor[atype]['decomcost_m' ] = getFromDict(props, 'decomcost_m' , default=0.0)
        anchor[atype]['decomcost_m2'] = getFromDict(props, 'decomcost_m2', default=0.0)
        anchor[atype]['decomcost_m3'] = getFromDict(props, 'decomcost_m3', default=0.0)
        anchor[atype]['decomcost_a' ] = getFromDict(props, 'decomcost_a' , default=0.0)
        anchor[atype]['decomcost_a2'] = getFromDict(props, 'decomcost_a2', default=0.0)
        anchor[atype]['decomcost_a3'] = getFromDict(props, 'decomcost_a3', default=0.0)
        anchor[atype]['decomcost'   ] = getFromDict(props, 'decomcost'   , default=0.0)

    # read the Buoy dict
    if 'BuoyProps' in source:
        BuoyProps = source['BuoyProps']
    else:
        raise MoorPyError("YAML file or dictionary must have a 'BuoyProps' field containing the data") # We might want to remove this requirement becasue if a user just gives point mass volume and cost they wont need a buoy dict
    
    buoy = dict()  # output dictionary combining default values with loaded coefficients
    
     # combine loaded coefficients and default values into dictionary that will be saved for each material
    for btype, props in BuoyProps.items():  
        buoy[btype] = {}
        # buoy[btype]['mass_b'   ] = getFromDict(props, 'mass_b'   , default=0.0)  
        # buoy[btype]['v_b'      ] = getFromDict(props, 'v_b'      , default=0.0)
        buoy[btype]['cost_b0'  ] = getFromDict(props, 'cost_b0'  , default=0.0)
        buoy[btype]['cost_b1'  ] = getFromDict(props, 'cost_b1'  , default=0.0)
        buoy[btype]['cost_b2'  ] = getFromDict(props, 'cost_b2'  , default=0.0)
        buoy[btype]['cost_b3'  ] = getFromDict(props, 'cost_b3'  , default=0.0)

    # read the connection hardware dict
    if 'ConnectProps' in source:
        ConnectProps = source['ConnectProps']
    else:
        raise MoorPyError("YAML file or dictionary must have a 'ConnectProps' field containing the data") # We might want to remove this requirement becasue if a user just gives point mass volume and cost they wont need a connect dict
    
    connect = dict()  # output dictionary combining default values with loaded coefficients
    
     # combine loaded coefficients and default values into dictionary that will be saved for each material
    for ctype, props in ConnectProps.items():  
        connect[ctype] = {}
        # connect[ctype]['mass_mbl' ] = getFromDict(props, 'mass_mbl' , default=0.0)  
        # connect[ctype]['v_mbl'    ] = getFromDict(props, 'v_mbl'    , default=0.0)
        connect[ctype]['FOS'      ] = getFromDict(props, 'FOS'      , default=1.0) # defaults to 1, to return design load as MBL
        connect[ctype]['cost_MBL0'] = getFromDict(props, 'cost_MBL0', default=0.0)
        connect[ctype]['cost_MBL1'] = getFromDict(props, 'cost_MBL1', default=0.0)
        connect[ctype]['cost_MBL2'] = getFromDict(props, 'cost_MBL2', default=0.0)
        connect[ctype]['cost_MBL3'] = getFromDict(props, 'cost_MBL3', default=0.0)

    # read the point design dict
    if 'DesignProps' in source:
        DesignProps = source['DesignProps']
    else:
        raise MoorPyError("YAML file or dictionary must have a 'DesignProps' field containing the data")

    point = dict()  # output dictionary combining default values with loaded designs
    
    # combine loaded coefficients and default values into dictionary that will be saved for each material
    for dtype, props in DesignProps.items():  
        point[dtype] = {}
        num_anchs = 0 # counter for allocating fraction defaults
        num_buoys = 0 # counter for allocating fraction defaults

        # load the number of each component in the design
        for key in props.keys():
            # handle the numbers of each Anchor type
            if 'num_a_' in key:
                atype = key[6:]
                if atype in AnchorProps.keys():
                    point[dtype][key] = getFromDict(props, key, default=0) 
                else:
                    raise MoorPyError(f'Anchor type {atype} not found in AnchorProps dictionary')
                
                num_anchs += 1
                    
            # handle the numbers of each Buoy type
            if 'num_b_' in key:
                btype = key[6:]
                if btype in BuoyProps.keys():
                    point[dtype][key] = getFromDict(props, key, default=0) 
                else:
                    raise MoorPyError(f'Buoy type {btype} not found in BuoyProps dictionary')
                
                num_buoys += 1
                
            # handle the numbers of each Connection type
            if 'num_c_' in key:
                ctype = key[6:]
                if ctype in ConnectProps.keys():
                    point[dtype][key] = getFromDict(props, key, default=0) 
                else:
                    raise MoorPyError(f'Connection type {ctype} not found in ConnectProps dictionary')

        # load the fractions, which default as 1/total number of type. This checks for num_<key> because frac_<key> should only be loaded if num_<key> is provided
        if num_anchs > 0 or num_buoys > 0:
            for key in props.keys():
                if 'num_a_' in key:
                    # handle the fractions for each Anchor type
                    point[dtype][f'frac_a_{key[6:]}'] = getFromDict(props, f'frac_a_{key[6:]}', default = 1/num_anchs)    
                if 'num_b_' in key:
                    # handle the fractions for each Buoy type
                    point[dtype][f'frac_b_{key[6:]}'] = getFromDict(props, f'frac_b_{key[6:]}', default = 1/num_buoys)

    # return ABCD dict
    return dict(AnchorProps = anchor, BuoyProps = buoy, ConnectProps = connect, DesignProps = point)


def getAnchorCost(fx, fz, type='drag-embedment'):
    '''Simple interface function for getting the cost of an anchor given the
    anchor load. Points to other functions in MoorPy that are in a state of 
    flux.
    Inputs: fx [N], fy [N], anchor type ('drag-embedment' or 'suction').
    Returns: anchor material cost [$]
    '''
    
    # Currently using the original cost function
    from moorpy.MoorProps import getAnchorCostOld
    
    anchorMatCost, anchorInstCost, anchorDecomCost, info = getAnchorCostOld(fx, fz, type=type)
    
    return anchorMatCost


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
        ms = lines2subsystem(lines, points, ms, span=None, case=case)
        ms.initialize()
        ms.solveEquilibrium()
        i += 1
        if i >= len(ms.pointList):
            break

    return ms

def lines2subsystem(lines,points, ms,span=None,case=0):
    '''Takes a set of connected lines (in order from rA to rB) in a moorpy system and creates a subsystem equivalent.
    The original set of lines are then removed from the moorpy system and replaced with the 
    subsystem.

    Parameters
    ----------
    lines : list
        List of indices in the ms.lineList to replace.
    points : list
        List of indices in the ms.pointList to replace
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

    # go through each line
    for i in lines:
        # collect line lengths and types                                          
        lengths.append(ms.lineList[i].L)
        types.append(ms.lineList[i].type['name'])
        ss.lineTypes[types[-1]] = ms.lineTypes[types[-1]]
        
    # use makeGeneric to build the subsystem line
    nSegs = [ms.lineList[i].nNodes-1 for i in lines]
    ss.makeGeneric(lengths,types,suspended=case, nSegs=nSegs)
    ss.setEndPosition(ms.lineList[lines[0]].rA,endB=0)
    ss.setEndPosition(ms.lineList[lines[-1]].rB,endB=1)
    
    # add in any info on the points connected to the lines
    # currently, mass, volume, and diameter but others could be added
    for i in range(0,len(points)):
        ss.pointList[i].m = ms.pointList[points[i]].m
        ss.pointList[i].v = ms.pointList[points[i]].v
        ss.pointList[i].d = ms.pointList[points[i]].d
    
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
                    if points[-1]>j and decB == 0:                    
                        points[-1] -= 1
                        decB = 1
                    if points[0]>j and decA == 0:
                        points[0] -= 1 
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
                    if points[-1]>j and decB == 0:                    
                        points[-1] -= 1
                        decB = 1
                    if points[0]>j and decA == 0:
                        points[0] -= 1 
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
    ms.pointList[points[0]].attachLine(ssNum,0) # rA
    ms.pointList[points[-1]].attachLine(ssNum,1) # rB     
        
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
   
def ss2lines(ms, nsegs=10):
    '''Replace all subsystem in ms with equivalent set of lines

    Parameters
    ----------
    ms : system object
        MoorPy system object that contains the subsystem
    nsegs : list OR int, optional
        Number of segments per line for each line. Can be an integer (all line sections have the same # of segments)
        OR can be a list (# of segments for each section of line in order from A to B)

    Returns
    -------
    new mooring system object with all subsystems replaced with lines

    '''
    from copy import deepcopy
    original_ms = deepcopy(ms)
    subsystemCount = len(original_ms.lineList)
    sub_idx = 0
    newly_created_lines = []
    shared_point_map = {}
    # Use a stable key per original point (coordinates are fine if stable)
    def _pt_key(P):
        # round to avoid tiny float diffs; include type to be safer
        return (round(P.r[0], 6), round(P.r[1], 6), round(P.r[2], 6), int(P.type))    
    for _ in range(subsystemCount):
        # get subsystem object
        ss = ms.lineList[sub_idx]   
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
                if ms.pointList[i].attached[j] == sub_idx+1:
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
                        key = _pt_key(ms.pointList[i])
                        if len(ms.pointList[i].attached) > 1:
                            # This is a shared anchor
                            other_attached_ss = []
                            for att in ms.pointList[i].attached:
                                if att != sub_idx + 1 and att not in newly_created_lines:
                                    other_attached_ss.append(att - 1)

                            # --- mark this endpoint as shared & key it for reuse later
                            points[endA]['shared'] = True
                            # build a stable key based on the ORIGINAL system's point data if possible
                            key = _pt_key(ms.pointList[i])
                            points[endA]['share_key'] = key
                            points[endA]['shared_with'] = other_attached_ss
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
        for i in range(0, len(ms.pointList)):
            key = _pt_key(ms.pointList[i])
            if key in shared_point_map:
                # This point has already been created for another subsystem - detach from current subsystem
                if sub_idx + 1 in ms.pointList[i].attached:
                    ms.pointList[i].attached.remove(sub_idx + 1)
        # remove subsystem line, delete all associated points
        helpers.deleteLine(ms,sub_idx,delpts=3)  # for some reason, in the delete line, it attached 19 to 13 again. Please investigate
        # add in new lines to replace subsystem
        for i in range(0,len(types)):
            # determine # of segments for this line section
            if isinstance(nsegs,list):
                NSegs = nsegs[i]
            elif isinstance(nsegs,int):
                NSegs = nsegs
            else:
                raise Exception('Input nsegs must be either a list or an integer')
            # add point A (if it wasn't added already by other subsystems since shared anchors are a possibility)
            # --- Reuse shared anchors (do not add again)
            use_existing = False
            curr_pt_idx = None  # 1-based index in ms.pointList            
            if points[i].get('shared', False):
                key = points[i]['share_key']
                if key in shared_point_map:
                    # already created for another subsystem – find it then reuse it
                    # find it
                    for j, point in enumerate(ms.pointList):
                        k_j = _pt_key(point)
                        if k_j == key:
                            curr_pt_idx = j + 1
                            break
                    use_existing = True            

            if not use_existing:
                # Not shared, or shared but not created yet -> create it now
                ms.addPoint(points[i]['type'], points[i]['r'], points[i]['m'], points[i]['v'], d=points[i]['d'])
                ms.pointList[-1].CdA = points[i]['CdA']
                ms.pointList[-1].Ca  = points[i]['Ca']
                curr_pt_idx = len(ms.pointList)
                if points[i].get('shared', False):
                    # This is a shared anchor, and has just been created so let's attach it to the other subsystems now
                    for other_ss in points[i]['shared_with']:
                        ms.pointList[curr_pt_idx - 1].attachLine(other_ss, endB=0)
                    # remember it if this is a shared anchor
                    shared_point_map[points[i]['share_key']] = curr_pt_idx

            # add line
            ms.addLine(lengths[i], types[i]['name'], nSegs=NSegs)
            newly_created_lines.append(len(ms.lineList))
            # attach point A to the just-created line
            ms.pointList[curr_pt_idx - 1].attachLine(len(ms.lineList), endB=0)

            # attach to any bodies the point was originally attached to
            if (not use_existing) and ('body' in points[i]):
                ms.bodyList[points[i]['body']].attachPoint(
                    curr_pt_idx,
                    [ms.pointList[curr_pt_idx - 1].r[0] - ms.bodyList[points[i]['body']].r6[0],
                    ms.pointList[curr_pt_idx - 1].r[1] - ms.bodyList[points[i]['body']].r6[1],
                    ms.pointList[curr_pt_idx - 1].r[2]]
                )
            if i > 0:
                # this same point is end B for previous line
                ms.pointList[curr_pt_idx - 1].attachLine(len(ms.lineList) - 1, endB=1)
        
        # add last point (should be the fairlead)
        ms.addPoint(points[-1]['type'], points[-1]['r'], points[-1]['m'], points[-1]['v'], d=points[-1]['d'])
        ms.pointList[-1].CdA = points[-1]['CdA']
        ms.pointList[-1].Ca = points[-1]['Ca']
        # attach to last line as point B
        ms.pointList[-1].attachLine(len(ms.lineList),endB=1)
        # attach to a body if applicable
        if 'body' in points[-1]:
            ms.bodyList[points[-1]['body']].attachPoint(len(ms.pointList),[ms.pointList[-1].r[0]-ms.bodyList[points[-1]['body']].r6[0],ms.pointList[-1].r[1]-ms.bodyList[points[-1]['body']].r6[1],ms.pointList[-1].r[2]])        
    return ms
            
def duplicateSyntheticLines(ms):
    '''reads in a MoorPy system and duplicates linetypes with nonzero EAd. needed for system unload to work with 
    multiple mean load values'''
    import copy
    
    # list of line types with nonzero EAd
    types = []
    for key, lineType in ms.lineTypes.items():
        if 'EAd' in lineType.keys() and lineType['EAd'] > 0:
            types.append(lineType['name'])
    

    if len(types) > 0:
        
        #iterate through types with nonzero EAd
        for t in types:

            #store indexes of lines that use that line type
            inds = []
            for i, line in enumerate(ms.lineList):
                if line.type['name'] == t:
                    inds.append(i)
            
            names = [t]
            
            # make copies of lineType (so that each segment with nonzero EAd has unique LineType)
            for i in inds[1:]:
                
                # insert the copies right below the existing linetype to make ordering more logical
                
                pos = list(ms.lineTypes.keys()).index(t) + 1
                items = list(ms.lineTypes.items())     
                items.insert(pos, (t+str(i), copy.deepcopy(ms.lineTypes[t])))
                ms.lineTypes = dict(items)
                ms.lineTypes[t + str(i)]['name'] = t + str(i)
                
                names.append(t + str(i))
            
            #make sure each line points to the correct lineType
            for j, i in enumerate(inds):
                ms.lineList[i].type = ms.lineTypes[names[j]]
        return(ms)
    else:
        print('No lines have nonzero EAd')                

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
    
    # Check if *** characters in data2 and if so replace with 0
    for i in range(len(data2)):
        for j in range(len(data2[i])):
            if data2[i][j].count('***') > 0:
                data2[i][j] = '0.0'

    data3 = data2.astype(float)
    
    return data3, ch, channels, units    

def read_file(dirName,fileName, skiplines=-1, hasunits=1, chanlim=999, dictionary=True):

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


def get_horizontal_oop_vec(p1,p2):
    """
    Evaluates the horizontal out of plane vector given the coordinates of two points.
    """
    hor_vec = p2 - np.array([p1[0],p1[1],p2[2]])
    ver_vec = p1 - np.array([p1[0],p1[1],p2[2]])

    if np.isclose(np.linalg.norm(hor_vec),0): # vertical line
        n_op = np.array([1,0,0]) 
    elif np.isclose(np.linalg.norm(ver_vec),0): # horizontal line
        oop_vec = np.cross(hor_vec,np.array([0,0,1])) 
        n_op = oop_vec/np.linalg.norm(oop_vec)
    else:
        oop_vec = np.cross(hor_vec,ver_vec)
        n_op = oop_vec/np.linalg.norm(oop_vec)
    return n_op


def guyan_reduce(K, n_coupled):
    """
    Perform Guyan reduction on stiffness matrix K, keeping the first n_coupled DOFs as the loaded dofs.
    Returns the reduced stiffness matrix with dimensions (n_coupled, n_coupled).
    See https://hal.science/hal-01711552v1/document
    """
 
    # Partitions of the stiffness matrix
    Kcc = K[:n_coupled, :n_coupled]
    Kcu = K[:n_coupled, n_coupled:]
    Kuc = K[n_coupled:, :n_coupled]
    Kuu = K[n_coupled:, n_coupled:]

    # Using the partitions to compute the reduced matrix
    try:
        Kuu_inv = np.linalg.inv(Kuu)
    except:
        Kuu_inv = np.linalg.pinv(Kuu)
    K_reduced = Kcc - Kcu @ Kuu_inv @ Kuc
    return K_reduced


def get_dynamic_matrices(Line, omegas, S_zeta,r_dynamic,depth,kbot,cbot,seabed_tol=1e-4):
    """
    Evaluates dynamic matrices for a Line object.

    Parameters
    ----------
    Line : Line
        An instance of MoorPy's Line class
    omegas : ndarray
        Array of frequencies in rad/s.
    S_zeta : ndarray
        Wave spectrum array in m^2/(rad/s), must be of the same length as omegas.
    r_dynamic : ndarray
        A 3d array of the frequency dependent complex amplitudes of line nodes. The array has a shape of (m,n,3) where m is the number of frequencies,
        n is the number of nodes, and the last dimension of 3 correspond to the x,y,z components.
    depth : float
        Water depth.
    kbot : float
        Vertical stiffness for points lying on the seabed.
    cbot : float
        Vertical damping for points lying on the seabed.
    seabed_tol : float, optional
        Distance from seabed within which a node is considered to be lying on the seabed, by default 1e-4 m.

    Returns
    -------
    M: ndarray
        Mass matrix.
    A: ndarray
        Added mass matrix.
    B: ndarray
        Damping matrix.
    K: ndarray
       Stiffness matrix.
    M: ndarray
        Mass matrix.
    r_mean: ndarray
        Mean static positions of the nodes given as a (m,3) array where m is the number of nodes.
    EA_segs: ndarray
        Extensional stiffness of segments.
    """
    # extract line properties
    N = Line.nNodes
    mden = Line.type['m'] # line mass density function
    deq = Line.type['d_vol'] # line volume equivalent diameter
    EA = Line.type['EA'] # extensional stiffness
    Ca = Line.type['Ca'] # normal added mass coeff
    CaAx = Line.type['CaAx'] # tangential added mass coeff
    Cd = Line.type['Cd'] # normal drag coeff
    CdAx = Line.type['CdAx'] # tangential drag coeff

    # extract mean node coordinates
    X_mean,Y_mean,Z_mean,T_mean = Line.getLineCoords(0.0,n=N) # coordinates of line nodes and tension values
    r_mean = np.vstack((X_mean,Y_mean,Z_mean)).T # coordinates of line nodes

    # evaluate node velocities
    r_dynamic = np.ones((len(omegas),N,3),dtype='float')*r_dynamic
    v_dynamic = 1j*omegas[:,None,None]*r_dynamic

    # define out of plane normal
    h_op = get_horizontal_oop_vec(r_mean[0],r_mean[-1]) # horizontal out-of-plane vector
    hh_op = np.outer(h_op,h_op)

    # intialize line matrices
    M = np.zeros([3*N, 3*N], dtype='float') # mass matrix
    A = np.zeros([3*N, 3*N], dtype='float') # added mass matrix
    B = np.zeros([3*N, 3*N], dtype='float') # linearized viscous damping matrix
    K = np.zeros([3*N, 3*N], dtype='float') # stiffness matrix

    # Node 0
    dr_e1 = r_mean[1] - r_mean[0]
    L_e1 = np.linalg.norm(dr_e1) # element 1 length
    t_e1 = (dr_e1)/L_e1 # tangential unit vector
    p_e1 = np.cross(t_e1,h_op) # in plane normal unit vector


    ut_e1 = np.einsum('ij,j->i',v_dynamic[:,0,:],t_e1) # tangential velocity
    uh_e1 = np.einsum('ij,j->i',v_dynamic[:,0,:],h_op) # normal horizontal out of plane velocity
    up_e1 = np.einsum('ij,j->i',v_dynamic[:,0,:],p_e1) # normal in plane velocity

    sigma_ut_e1 = np.sqrt(np.trapezoid(np.abs(ut_e1)**2*S_zeta,omegas))
    sigma_uh_e1 = np.sqrt(np.trapezoid(np.abs(uh_e1)**2*S_zeta,omegas))
    sigma_up_e1 = np.sqrt(np.trapezoid(np.abs(up_e1)**2*S_zeta,omegas))

    tt_e1 = np.outer(t_e1,t_e1) # local tangential to global components transformation matrix
    pp_e1 = np.outer(p_e1,p_e1) # local normal inplane to global components transformation matrix

    M[0:3,0:3] += mden*L_e1/2*np.eye(3) # element 1 mass contribution

    A_e1 = 1025*np.pi/4*deq**2*L_e1/2*(Ca*(hh_op+pp_e1) + CaAx*tt_e1) # element 1 added mass contribution

    B_e1 = 0.5*1025*deq*L_e1/2*np.sqrt(8/np.pi)*(Cd*(sigma_uh_e1*hh_op + sigma_up_e1*pp_e1) +
                                                CdAx*sigma_ut_e1*tt_e1) # element 1 damping contribution 

    K_e1 = EA/L_e1*tt_e1 + (T_mean[0]/L_e1)*(hh_op+pp_e1) # element 1 stiffness (axial + geometric)

    ## assembling element 1 contributions (rows corresponding to node 0)
    A[0:3,0:3] += A_e1 
    B[0:3,0:3] += B_e1
    K[0:3,0:3] += K_e1
    K[0:3,3:6] += -K_e1

    ## add seabed contribution to node 0
    if np.isclose(r_mean[0,2],-depth,seabed_tol):
        K[2,2] += kbot
        B[2,2] += cbot 

    # Internal nodes loop (each internal node has contributions from two elements n-1/2 and n+1/2)
    for n in range(1, N-1):
        
        ## backward element (n-1/2) contributions
        dr_bw = r_mean[n-1] - r_mean[n]
        L_bw = np.linalg.norm(dr_bw) # element 1 length
        t_bw = (dr_bw)/L_bw # tangential unit vector
        p_bw = np.cross(t_bw,h_op) # in plane normal unit vector

        ut_bw = np.einsum('ij,j->i',v_dynamic[:,n,:],t_bw) # tangential velocity
        uh_bw = np.einsum('ij,j->i',v_dynamic[:,n,:],h_op) # normal horizontal out of plane velocity
        up_bw = np.einsum('ij,j->i',v_dynamic[:,n,:],p_bw) # normal in plane velocity

        sigma_ut_bw = np.sqrt(np.trapezoid(np.abs(ut_bw)**2*S_zeta,omegas))
        sigma_uh_bw = np.sqrt(np.trapezoid(np.abs(uh_bw)**2*S_zeta,omegas))
        sigma_up_bw = np.sqrt(np.trapezoid(np.abs(up_bw)**2*S_zeta,omegas))

        tt_bw = np.outer(t_bw,t_bw) # local tangential to global components transformation matrix
        pp_bw = np.outer(p_bw,p_bw) # local normal inplane to global components transformation matrix

        M[3*n:3*n+3,3*n:3*n+3] += mden*L_bw/2*np.eye(3) # mass contribution from adjacent elements

        A_bw = 1025*np.pi/4*deq**2*L_bw/2*(Ca*(hh_op+pp_bw) + CaAx*tt_bw) # backward element added mass contribution

        B_bw = 0.5*1025*deq*L_bw/2*np.sqrt(8/np.pi)*(Cd*(sigma_uh_bw*hh_op + sigma_up_bw*pp_bw) +
                                                        CdAx*sigma_ut_bw*tt_bw) # backward element damping contribution 

        K_bw = EA/L_bw*tt_bw + (T_mean[n]/L_bw)*(hh_op+pp_bw) # backward element stiffness (axial + geometric)

        ## forward element (n+1/2) contributions
        dr_fw = r_mean[n+1] - r_mean[n]
        L_fw = np.linalg.norm(dr_fw) # element 1 length
        t_fw = (dr_fw)/L_fw # tangential unit vector
        p_fw = np.cross(t_fw,h_op) # in plane normal unit vector


        ut_fw = np.einsum('ij,j->i',v_dynamic[:,n,:],t_fw) # tangential velocity
        uh_fw = np.einsum('ij,j->i',v_dynamic[:,n,:],h_op) # normal horizontal out of plane velocity
        up_fw = np.einsum('ij,j->i',v_dynamic[:,n,:],p_fw) # normal in plane velocity

        sigma_ut_fw = np.sqrt(np.trapezoid(np.abs(ut_fw)**2*S_zeta,omegas))
        sigma_uh_fw = np.sqrt(np.trapezoid(np.abs(uh_fw)**2*S_zeta,omegas))
        sigma_up_fw = np.sqrt(np.trapezoid(np.abs(up_fw)**2*S_zeta,omegas))

        tt_fw = np.outer(t_fw,t_fw) # local tangential to global components transformation matrix
        pp_fw = np.outer(p_fw,p_fw) # local normal inplane to global components transformation matrix
        
        M[3*n:3*n+3,3*n:3*n+3] += mden*L_fw/2*np.eye(3) # mass contribution from adjacent elements

        A_fw = 1025*np.pi/4*deq**2*L_fw/2*(Ca*(hh_op+pp_fw) + CaAx*tt_fw) # forward element added mass contribution

        B_fw = 0.5*1025*deq*L_fw/2*np.sqrt(8/np.pi)*(Cd*(sigma_uh_fw*hh_op + sigma_up_fw*pp_fw) +
                                                    CdAx*sigma_ut_fw*tt_fw) # forward element damping contribution 

        K_fw = EA/L_fw*tt_fw + (T_mean[n]/L_fw)*(hh_op+pp_fw) # forward element stiffness (axial + geometric)

        ## assembling bwd and fwd elements contributions (rows corresponding to node n)
        A[3*n:3*n+3,3*n:3*n+3] += A_bw + A_fw
        B[3*n:3*n+3,3*n:3*n+3] += B_bw + B_fw
        K[3*n:3*n+3,3*n:3*n+3] += K_bw + K_fw 
        K[3*n:3*n+3,3*n-3:3*n] += -K_bw
        K[3*n:3*n+3,3*n+3:3*n+6] += -K_fw
        
        ## add seabed contribution to node n
        if np.isclose(r_mean[n,2],-depth,seabed_tol):
            K[3*n+2,3*n+2] += kbot
            B[3*n+2,3*n+2] += cbot

    # Node N
    dr_eN = r_mean[N-1] - r_mean[N-2]
    L_eN = np.linalg.norm(dr_eN) # element N length
    t_eN = (dr_eN)/L_eN # tangential unit vector
    p_eN = np.cross(t_eN,h_op) # in plane normal unit vector

    ut_eN = np.einsum('ij,j->i',v_dynamic[:,N-1,:],t_eN) # tangential velocity
    uh_eN = np.einsum('ij,j->i',v_dynamic[:,N-1,:],h_op) # normal horizontal out of plane velocity
    up_eN = np.einsum('ij,j->i',v_dynamic[:,N-1,:],p_eN) # normal in plane velocity

    sigma_ut_eN = np.sqrt(np.trapezoid(np.abs(ut_eN)**2*S_zeta,omegas))
    sigma_uh_eN = np.sqrt(np.trapezoid(np.abs(uh_eN)**2*S_zeta,omegas))
    sigma_up_eN = np.sqrt(np.trapezoid(np.abs(up_eN)**2*S_zeta,omegas))

    tt_eN = np.outer(t_eN,t_eN) # local tangential to global components transformation matrix
    pp_eN = np.outer(p_eN,p_eN) # local normal inplane to global components transformation matrix

    M[3*(N-1):3*(N-1)+3,3*(N-1):3*(N-1)+3] += mden*L_eN/2*np.eye(3) # element N mass contribution

    A_eN = 1025*np.pi/4*deq**2*L_eN/2*(Ca*(hh_op+pp_eN) + CaAx*tt_eN) # element N added mass contribution

    B_eN = 0.5*1025*deq*L_eN/2*np.sqrt(8/np.pi)*(Cd*(sigma_uh_eN*hh_op + sigma_up_eN*pp_eN) +
                                                CdAx*sigma_ut_eN*tt_eN) # element N damping contribution 

    K_eN = EA/L_eN*tt_eN + (T_mean[N-1]/L_eN)*(hh_op+pp_eN) # element N stiffness (axial + geometric)

    ## assembling element N contributions (rows corresponding to node N)
    A[3*(N-1):3*(N-1)+3,3*(N-1):3*(N-1)+3] += A_eN 
    B[3*(N-1):3*(N-1)+3,3*(N-1):3*(N-1)+3] += B_eN 
    K[3*(N-1):3*(N-1)+3,3*(N-1):3*(N-1)+3] += K_eN
    K[3*(N-1):3*(N-1)+3,3*(N-1)-3:3*(N-1)] += -K_eN

    ## add seabed contribution to node N
    if np.isclose(r_mean[N-1,2],-depth,seabed_tol):
        K[3*(N-1)+2,3*(N-1)+2] += kbot
        B[3*(N-1)+2,3*(N-1)+2] += cbot

    EA_segs = Line.type['EA']*np.ones(Line.nNodes - 1)

    return M,A,B,K,r_mean,EA_segs


def get_dynamic_tension(line,omegas,S_zeta,RAO_A,RAO_B,depth,kbot,cbot,seabed_tol=1e-4,tol = 0.01,iters=100, w = 0.8, conv_time=False, returnMatrices=False):
    """
    Evaluates dynamic tension at all the nodes for an instance of MoorPy's line or CompositeLine classes.

    Parameters
    ----------
    line : line/CompositeLine
        An instance of MoorPy's line or CompositeLine classes.
    omegas : ndarray
        Array of frequencies in rad/s.
    S_zeta : ndarray
        Wave spectrum array in m^2/(rad/s), must be of the same length as omegas.
    RAO_A : ndarray
        Translational RAOs for end node A given as a (m,3) array where m is the number of frequencies (must be equal to the length of omegas) .
    RAO_B : ndarray
        Translational RAOs for end node B given as a (m,3) array where m is the number of frequencies (must be equal to the length of omegas) .
    depth : float
        Water depth.
    kbot : float
        Vertical stiffness for points lying on the seabed.
    cbot : float
        Vertical damping for points lying on the seabed.
    seabed_tol : float, optional
        Distance from seabed within which a node is considered to be lying on the seabed, by default 1e-4 m.
    tol : float, optional
        Relative tolerance for iteration, by default 0.01
    iters : int, optional
        Maximum number of iterations, by default 100
    w : float, optional
        Succesive relaxation coefficient, by default 0.8

    Returns (if returnMatrices = False)
    -------
    T_nodes_psd: ndarray
        Tension amplitude at nodes given as (m,n) array where m is the number of frequencies and n is the number of nodes.
    T_nodes_psd: ndarray
        Tension spectra at nodes given as (m,n) array where m is the number of frequencies and n is the number of nodes.
    T_nodes_std: ndarray
        Tension standard deviation at nodes.
    s: ndarray:
        Node locations along the line.
    r_static: ndarray
        Nodes mean static position given as (n,3) array where n the number of nodes.
    r_dynamic: ndarray
        Nodes complex dynamic amplitudes given as (m,n,3) array where m the number of frequencies, n is the number of nodes
    r_total: ndarray
        Combined static and dynamic positions.
    X: ndarray
        Solution of the dynamic problem.

    Returns (if returnMatrices = True)
    -------
    M: (n,n) inertia matrix, where n is the number of nodes
    A: (n,n) added mass matrix
    B: (n,n) damping matrix
    K: (n,n) stiffness matrix
    """
    N = line.nNodes
    n_dofs = 3*N

    if np.all(RAO_A == 0):
        RAO_A = np.zeros_like(RAO_B)
    if np.all(RAO_B == 0):
        RAO_B = np.zeros_like(RAO_A)

    # intialize iteration matrices
    r_dynamic_init = np.ones((len(omegas),N,3))
    M,A,B,K,r_static,EA_segs = line.getDynamicMatrices(omegas,S_zeta,r_dynamic_init,depth,kbot,cbot,seabed_tol=seabed_tol) # TODO: return EA_segs
    X = np.zeros((len(omegas),n_dofs),dtype = 'complex')
    r_dynamic = np.zeros(((len(omegas),int(n_dofs/3),3)),dtype = 'complex')
    S_Xd = np.zeros((len(omegas),n_dofs),dtype = 'float')
    sigma_Xd = np.zeros(n_dofs,dtype = 'float')
    sigma_Xd0 = np.zeros(n_dofs,dtype = 'float')
    X[:, :3] = RAO_A
    X[:,-3:] = RAO_B

    # solving dynamics
    start = time.time()
    for ni in range(iters):
        H = - omegas[:,None,None]**2*(M+A)[None,:,:]\
            + 1j*omegas[:,None,None]*B[None,:,:]\
            + K[None,:,:]\
        
        F_A = np.einsum('nij,njk->ni',-H[:,3:-3, :3],RAO_A[:,:,None])
        F_B = np.einsum('nij,njk->ni',-H[:,3:-3,-3:],RAO_B[:,:,None])
        F = F_A + F_B

        # Solve for each frequency
        for i in range(H.shape[0]):
            X[i, 3:-3] = np.linalg.solve(H[i, 3:-3, 3:-3], F[i])
            
        S_Xd[:] = np.abs(1j*omegas[:,None]*X)**2*S_zeta[:,None]
        sigma_Xd[:] = np.sqrt(np.trapezoid(S_Xd,omegas,axis=0)) 
        r_dynamic[:] = X.reshape(X.shape[0],int(X.shape[1]/3),3)
        if (np.abs(sigma_Xd-sigma_Xd0) <= tol*np.abs(sigma_Xd0)).all():
            break
        else:
            sigma_Xd0[:] = w * sigma_Xd + (1.-w) * sigma_Xd0
            _,_,B[:],_,_,_ = line.getDynamicMatrices(omegas,S_zeta,r_dynamic,depth,kbot,cbot,seabed_tol=seabed_tol)
    if conv_time:
        print(f'Finished {ni} dynamic tension iterations in {time.time()-start} seconds (w = {w}).')

    # evaluate tension
    dw = np.diff(omegas,
            prepend= omegas[0] - (omegas[1]-omegas[0]),
            append= omegas[-1] + (omegas[-1]-omegas[-2]))
    dw = (dw[1:]+dw[:-1])/2
    wave_amps = np.sqrt(2*S_zeta*dw) #evaluate wave amplitudes of harmonic components from wave spectrum

    r_dynamic *= wave_amps[:,None,None]
    r_total = r_static[None,:,:] + r_dynamic
    dr_static = r_static[:-1] - r_static[1:]
    dr_dynamic = r_dynamic[:,:-1,:] - r_dynamic[:,1:,:]
    tangents = dr_static/np.linalg.norm(r_static[:-1] - r_static[1:], axis=-1)[:,None]
    L_static = np.linalg.norm(dr_static, axis=-1)
    dL_dynamic = np.einsum('mni,ni->mn', dr_dynamic, tangents)
    eps_segs = dL_dynamic/L_static

    T_segs = EA_segs * eps_segs
    T_nodes_amp = np.zeros((len(omegas),N), dtype='complex')
    T_nodes_amp[:,0] = T_segs[:,0]
    T_nodes_amp[:,1:-1] = (T_segs[:,:-1] + T_segs[:,1:])/2
    T_nodes_amp[:,-1] = T_segs[:,-1]

    # S_T = np.zeros((len(omegas),N))
    # S_T[:,1:] = T_e**2/dw[:,None]
    # S_T[:,0] = S_T[:,1]

    T_nodes_psd = np.abs(T_nodes_amp)**2/(2*dw[:,None])
    T_nodes_std = np.sqrt(np.trapezoid(T_nodes_psd,omegas,axis=0))


    dr = np.diff(r_static,axis=0)
    ds = np.linalg.norm(dr,axis=1)
    s = np.zeros_like(T_nodes_std)
    s = np.cumsum(ds)

    if returnMatrices:
        return M, A, B, K
    else:
        return T_nodes_amp, T_nodes_psd, T_nodes_std, s, r_static, r_dynamic, r_total, X


def get_modes(line,fix_A=True,fix_B=True,plot_modes=False,amp_factor=1,adj_view = False,kbot=3E+06,cbot=3E+05,seabed_tol=1E-04):
    import matplotlib.pyplot as plt
    from collections.abc import Iterable

    M,A,_,K,r_nodes,_ = line.getDynamicMatrices(np.ones(1), np.ones(1),0.,line.sys.depth,kbot,cbot,seabed_tol=seabed_tol)

    if fix_A:
        n1 = 1
    else:
        n1 = 0
    
    if fix_B:
        n2 = -1
    else:
        n2 = r_nodes.shape[0]
    
    eigvals,eigvecs = la.eig(K[3*n1:3*n2,3*n1:3*n2],M[3*n1:3*n2,3*n1:3*n2]+A[3*n1:3*n2,3*n1:3*n2])
    stable_eigvals = eigvals[np.real(eigvals)>0]
    stable_eigvecs = eigvecs[:,np.real(eigvals)>0]
    
    idx = stable_eigvals.argsort()[::-1]   
    stable_eigvals = stable_eigvals[idx]
    stable_eigvecs = stable_eigvecs[:,idx]
   
    freqs = np.sqrt(np.real(stable_eigvals))/2/np.pi
    mode_shapes = np.zeros(stable_eigvecs.shape,dtype='float')
    
    for i in range(stable_eigvecs.shape[1]):
        mode_shapes[:,i] = r_nodes[n1:n2].flatten('C') + stable_eigvecs[:,i]

    freqs = np.flip(freqs)
    mode_shapes = np.flip(mode_shapes,axis=1)

    if plot_modes:
        cols = 4
        rows = plot_modes//cols + bool(plot_modes%cols)
        fig,ax = plt.subplots(rows,cols,subplot_kw={"projection": "3d"},figsize=(5*cols,5*rows))

        i = 0
        for axes in ax:
            if not isinstance(axes,Iterable):
                axes = [axes]

            for axis in axes:
                if i >= plot_modes:
                    break

                r = r_nodes.copy()
                r[n1:n2] = mode_shapes[:,i].reshape([int(len(mode_shapes[:,i])/3),3])
                r = (r-r_nodes)*amp_factor
                x = r[:,0]
                y = r[:,1]
                z = r[:,2]
                
                x_0 = r_nodes[:,0]
                y_0 = r_nodes[:,1]
                z_0 = r_nodes[:,2]
                    
                axis.plot(x_0,y_0,z_0,'-ko',label='initial')
                axis.plot(x+x_0,y+y_0,z+z_0,'--ro',label='mode shape')
                
                # R_0 = np.sqrt(x_0**2 + y_0**2)
                if adj_view:
                    # h_min = np.min((x_0,y_0))
                    # h_max = np.max((x_0,y_0))
                    # axis.set_xlim(h_min,h_max)
                    # axis.set_ylim(h_min,h_max)
                    # axis.set_zlim(z_0.min(),z_0.max())
                    sigma_x = x.std() 
                    sigma_y = y.std()
                    sigma_z = z.std()
                    azim = np.arctan2(sigma_x,sigma_y)*180/np.pi
                    elev = np.arctan2(np.hypot(sigma_x,sigma_y),sigma_z)*180/np.pi
                    axis.view_init(elev=elev,azim=azim)

                # axis.set_box_aspect([np.ptp(coord) for coord in [x, y, z]])
                axis.set_xlabel('X (m)')
                axis.set_ylabel('Y (m)')
                axis.set_zlabel('Z (m)')
                axis.set_title(f'f = {freqs[i]:.3f} Hz, T = {1/freqs[i]:.3f} s')

                i+=1

        fig.tight_layout()
        return freqs,mode_shapes,r_nodes,M,A,K,fig,ax        
    else:
        return freqs,mode_shapes,r_nodes,M,A,K
