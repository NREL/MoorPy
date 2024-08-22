import numpy as np
import yaml

from moorpy.system import System
from moorpy.body import Body
from moorpy.point import Point
from moorpy.line import Line, from2Dto3Drotated
from moorpy.lineType import LineType
from moorpy.helpers import (rotationMatrix, rotatePosition, getH, printVec, 
                            set_axes_equal, dsolve2, SolveError, MoorPyError, 
                            loadLineProps, getLineProps, read_mooring_file, 
                            printMat, printVec, getInterpNums, unitVector,
                            getFromDict, addToDict)



class Subsystem(System, Line):
    '''A class for a mooring line or dynamic cable subsystem.
    It includes Line sections but also can fit into a larger System
    the same way a Line can.
    
    A subsystem tracks its own objects in its local coordinate system.
    This local coordinate system puts the x axis along the line heading,
    from anchor to fairlead, with the anchor at x=0.
    
    For a multi-section line or cable (the main reason to use a SubSystem),
    the free DOFs of the points are contained in the SubSystem and are not 
    seen by the parent System. Solving their equilibrium is an internal
    solve nested within the larger System's equilibrium solve.
    
    The first and last points in the subsystem are set as coupled so that
    their stiffness matrices can be computed and passed back the same way as
    those of a line's two ends.
    
    Line types should be referenced from the overall System... [how?]
    
    Key adjustments: replace staticSolve from Line to do an internal 
    equilibrium solve using System.solveEquilibrium.
    '''
    
    
    def __init__(self, mooringSys=None, num=0, depth=0, rho=1025, g=9.81, 
                 lineProps=None, **kwargs):
        '''Shortened initializer for just the SubSystem aspects.
        
        How should this be initialized???
        '''
        
        # some basics for compatibility with Line objects when used in a System lineList
        self.sys = mooringSys  # store a reference to the overall mooring system (instance of System class)
        self.number = num
        self.isRod = False
        self.type = {}  # because Line objects have a type dict...
        
        # lists to hold mooring system objects
        self.bodyList = []
        self.rodList = []
        self.pointList = []
        self.lineList = []
        self.lineTypes = {}  # dict indexed as a list [0:n-1] corresponding to lineList. Contains *references* to lineType dicts.
        self.rodTypes = {}
        
        # some basic subsystem things
        self.rA = np.zeros(3)
        self.rB = np.zeros(3)
        self.LH = 0
        self.th = 0
        self.sin_th = 0
        self.cos_th = 1
        
        # load mooring line property scaling coefficients for easy use when creating line types
        self.lineProps = loadLineProps(lineProps)
        
        # the ground body (number 0, type 1[fixed]) never moves but is the parent of all anchored things
        self.groundBody = Body(self, 0, 1, np.zeros(6))
        
        # constants used in the analysis
        self.depth = depth  # water depth [m]
        self.rho   = rho    # water density [kg/m^3]
        self.g     = g      # gravitational acceleration [m/s^2]
        
        # Set position tolerance to use in equilibrium solves [m]
        self.eqtol = getFromDict(kwargs, 'eqtol', default=0.01)
        
        # water current - currentMod 0 = no current; 1 = steady uniform current
        self.currentMod = 0         # flag for current model to use
        self.current = np.zeros(3)  # current velocity vector [m/s]
        if 'current' in kwargs:
            self.currentMod = 1
            self.current = getFromDict(kwargs, 'current', shape=3)
            
        # flag to indicate shared line or suspended cable being modeled as symmetric
        self.shared = getFromDict(kwargs, 'shared', dtype=bool, default=False)
        
        self.span   = getFromDict(kwargs, 'span', default=0)  # horizontal end-end distance [m]
        self.rad_fair = getFromDict(kwargs, 'rad_fair', default=0)  # [m] fairlead radius [m]
        self.z_fair   = getFromDict(kwargs, 'z_fair'  , default=0)  # [m] fairlead z coord [m]
        
        # the old rAFair input for the position of a second fairlead of a shaerd mooring line is gone
        # currently we assume rad_fair and z_fair are same for both ends of a shared line...
        
        # seabed bathymetry - seabedMod 0 = flat; 1 = uniform slope, 2 = grid
        self.seabedMod = 0
        
        if 'xSlope' in kwargs or 'ySlope' in kwargs:
            self.seabedMod = 1
            self.xSlope = getFromDict(kwargs, 'xSlope', default=0)
            self.ySlope = getFromDict(kwargs, 'ySlope', default=0)
        
        if 'bathymetry' in kwargs:
            self.seabedMod = 2
            self.bathGrid_Xs, self.bathGrid_Ys, self.bathGrid = self.readBathymetryFile(kwargs['bathymetry'])
        # Note, System and Subsystem bathymetry can be set after creation 
        # by setBathymetry, which uses link to existing data, for efficiency.
        
        # initializing variables and lists        
        self.nDOF = 0       # number of (free) degrees of freedom of the mooring system (needs to be set elsewhere)        
        self.freeDOFs = []  # array of the values of the free DOFs of the system at different instants (2D list)
        
        self.nCpldDOF = 0   # number of (coupled) degrees of freedom of the mooring system (needs to be set elsewhere)        
        self.cpldDOFs = []  # array of the values of the coupled DOFs of the system at different instants (2D list)
        
        self.display = 0    # a flag that controls how much printing occurs in methods within the System (Set manually. Values > 0 cause increasing output.)
        
        self.MDoptions = {} # dictionary that can hold any MoorDyn options read in from an input file, so they can be saved in a new MD file if need be
        self.qs = 1         # flag that it's a MoorPy analysis, so System methods don't complain. One day should replace this <<<
        
    
    def initialize(self, daf_dict={}):
        '''Initialize the Subsystem including DAFs.'''
        
        # Count number of line sections and total number of nodes
        self.nLines = len(self.lineList)
        self.nNodes = np.sum([l.nNodes for l in self.lineList]) - self.nLines + 1
        
        # Use the System initialize method
        System.initialize(self)
        
        # dynamic amplication factor for each line section, and anchor forces 
        # (DAFS[-2] is for vertical load, DAFS[-1] is for horizontal load)
        self.DAFs = getFromDict(daf_dict, 'DAFs', shape=self.nLines+2, default=1.0)  
        
        # mean tension [N] of each line section end [section #, end A/B] for any given mean offset
        self.Te0 = np.zeros([self.nLines,2])  # undispalced values
        self.TeM = np.zeros([self.nLines,2])  # mean offset values
        self.TeD = np.zeros([self.nLines,2])  # dynamic values
        
        # adjustment on laylength... positive means that the dynamic lay length is greater than linedesign laylength 
        self.LayLen_adj = getFromDict(daf_dict, 'LayLen_adj', shape=0, default=0.0) 
        
        # Fatigue damage over all nodes (placeholder for now)
        self.damage = np.zeros(self.nNodes)
        
    
    def makeGeneric(self, lengths, types, connectors=[], suspended=0):

        '''Creates a cable of n components going between an anchor point and
        a floating body (or a bridle point). If shared, it goes between two
        floating bodies.

        Parameters
        ----------
        lengths : list of floats
            List of each section's length. This also implies the number of
            sections.
        types : list of strings or dicts
            List of lineType names or dicts for each section. If strings, 
            these names must match keys in the parent system lineTypes 
            dictionary or the subsystem's lineTypes dictionary. If dicts,
            these dicts are referred to for each lineType (by reference).
        connectors : list of dicts
            List of length nLines-1 with dicts of optional properties for
            any interior points (between sections).
        suspended : int
            Selector shared/suspended cases: 
            - 0 (default): end A is on the seabed,
            - 1: the assembly is suspended and end A is at another floating system,
            - 2: the assembly is suspended and assumed symmetric, end A is the midpoint.
        '''
        
        # some initialization steps.
        self.nLines = len(lengths)
        if len(connectors) == 0:
            connectors = [{}]*(self.nLines - 1)
        elif not len(connectors) == self.nLines - 1:
            raise Exception('Length of connectors must be nLines - 1')
        
        if not len(types)==self.nLines:
            raise Exception("The specified number of lengths and types is inconsistent.")
        
        # get cumulative sum of line lengths, starting from anchor segment
        Lcsum = np.cumsum(np.array(lengths))
        
        # set end A location depending on whether configuration is suspended/symmetrical
        if suspended==2:  # symmetrical suspended case
            rA = np.array([-0.5*self.span-self.rad_fair, 0, -1])  # shared line midpoint coordinates
            self.shared = True  # flag that it's being modeled as symmetric
        elif suspended==1:  # general suspended case
            rA = np.array([-self.span-self.rad_fair, 0, self.z_fair])  # other suspended end
        else:  # normal anchored line case
            rA = np.array([-self.span-self.rad_fair, 0, -self.depth])  # anchor coordinates
        rB = np.array([-self.rad_fair, 0, self.z_fair])     # fairlead coordinates

        self.rA = rA
        self.rB = rB

        if suspended==2:
            self.addPoint(0, rA, DOFs=[2]) # add shared line point, free only to move in z
        else:
            self.addPoint(-1, rA, DOFs=[0,2])                # add anchor point
        
        # Go through each line segment and add its upper point, add the line, and connect the line to the points
        for i in range(self.nLines):

            # find the specified lineType dict and save a reference to it
            if type(types[i]) == dict:  # if it's a dictionary, just point to it
                self.lineTypes[i] = types[i]
            # otherwise we're assuming it's a string of the lineType name
            elif types[i] in self.lineTypes:  # first look for the name in the subsystem
                self.lineTypes[i] = self.lineTypes[types[i]]
            elif self.sys: # otherwise look in the parent system, if there is one
                if types[i] in self.sys.lineTypes:  # first look for the name in the subsystem
                    self.lineTypes[i] = self.sys.lineTypes[types[i]]
                else:
                    raise Exception(f"Can't find lineType '{types[i]}' in the SubSystem or parent System.")
            else:
                raise Exception(f"Can't find lineType '{types[i]}' in the SubSystem.")
            
            # add the line segment using the reference to its lineType dict
            self.addLine(lengths[i],self.lineTypes[i])

            # add the upper end point of the segment
            if i==self.nLines-1:                            # if this is the upper-most line
                self.addPoint(-1, rB, DOFs=[0,2])  # add the fairlead point (make it coupled)
                #self.bodyList[0].attachPoint(i+2, rB)       # attach the fairlead point to the body (two points already created)
            else:                                           # if this is an intermediate line
                m = connectors[i].get('m', 0)
                v = connectors[i].get('v', 0)
                # add the point, initializing linearly between anchor and fairlead/midpoint
                self.addPoint(0, rA + (rB-rA)*Lcsum[i]/Lcsum[-1], m=m, v=v, DOFs=[0,2])

            # attach the line to the points
            self.pointList[-2].attachLine(i+1, 0)       # attach end A of the line
            self.pointList[-1].attachLine(i+1, 1)       # attach end B of the line
        
    
    def setEndPosition(self, r, endB, sink=False):
        '''Sets either end position of the subsystem in the global/system
        reference frame. This is included mainly to mimic the Line method.

        Parameters
        ----------
        r : array
            x,y,z coorindate position vector of the line end [m].
        endB : boolean
            An indicator of whether the r array is at the end or beginning of the line
        sink : bool
            If true, and if there is a subsystem, the z position will be on the seabed.

        Raises
        ------
        LineError
            If the given endB value is not a 1 or 0
        '''
        
        if sink: # set z coordinate on seabed
            z, _ = self.getDepthFromBathymetry(r[0], r[1])
            r = np.array([r[0], r[1], z])  # (making a copy of r to not overwrite it)
        
        # set end coordinates in global frame just like for a Line
        if endB == 1:
            self.rB = np.array(r, dtype=float)
        elif endB == 0:
            self.rA = np.array(r, dtype=float)
        else:
            raise LineError("setEndPosition: endB value has to be either 1 or 0")
    
    
    def staticSolve(self, reset=False, tol=0, profiles=0):
        '''Solve internal equilibrium of the Subsystem and saves the forces
        and stiffnesses at the ends in the global reference frame. All the 
        This method mimics the behavior of the Line.staticSolve method and 
        the end coordinates rA and rB must be set beforehand. The solve 
        equilibrium happens in the local 2D plane. Values in this local 
        frame are also saved. 
        '''
        
        if tol==0:
            tol=self.eqtol
        
        # transform end positions to SubSystem internal coordinate system
        # inputs are self.rA, rB in global frame
        # outputs should be pointList[0] and [N] .r
        dr =  self.rB - self.rA
        LH = np.hypot(dr[0], dr[1])         # horizontal spacing of line ends
        LV = dr[2]                          # vertical rise from end A to B
        if self.shared:
            self.pointList[ 0].setPosition([  -self.span/2   , 0, self.rA[2]])
            self.pointList[-1].setPosition([ -self.span/2+LH, 0, self.rB[2]])
        else:
            self.pointList[ 0].setPosition([ -self.span   , 0, self.rA[2]])
            self.pointList[-1].setPosition([ -self.span+LH, 0, self.rB[2]])
            
        # get equilibrium
        self.solveEquilibrium(tol=tol)
        
        # get 2D stiffness matrices of end points
        K = self.getCoupledStiffnessA()
        
        # transform coordinates and forces back into global frame
        if LH > 0:
            self.cos_th = dr[0]/LH          # cos of line heading
            self.sin_th = dr[1]/LH          # sin of line heading
            th = np.arctan2(dr[1],dr[0])    # heading (from A to B) [rad]
            self.R = rotationMatrix(0,0,th)  # rotation matrix about z that goes from +x direction to heading
            
        else:   # special case of vertical line: line heading is undefined - use zero as default
            self.cos_th = 1.0
            self.sin_th = 0.0
            self.R = np.eye(3)
        
        # save end forces and stiffness matrices (first in local frame)
        self.fA_L = self.pointList[ 0].getForces(xyz=True) # force at end A
        self.fB_L = self.pointList[-1].getForces(xyz=True) # force at end B
        
        # Compute transverse (out-of-plane) stiffness term
        if LH < 0.01*abs(LV):  # if line is nearly vertical (note: this theshold is unverified)
            Kt = 0.5*(self.fA_L[2] - self.fB_L[2])/LV  # compute Kt based on vertical tension/span
        else:  # otherwise use the classic horizontal approach
            Kt = -self.fB_L[0]/LH  
        
        # expand to get 3D stiffness matrices
        '''
        R = np.eye(3)
        self.KA_L  = from2Dto3Drotated(K[:2,:2],  Kt, R.T)  # reaction at A due to motion of A
        self.KB_L  = from2Dto3Drotated(K[2:,2:],  Kt, R.T)  # reaction at B due to motion of B
        self.KBA_L = from2Dto3Drotated(K[2:,:2], -Kt, R.T)  # reaction at B due to motion of A
        '''
        self.KA_L = np.array([[K[0,0], 0 , K[0,1]],
                              [  0   , Kt,   0   ],
                              [K[1,0], 0 , K[1,1]]])
        
        # If symmetrical model, ignore midpoint stiffness and force
        if self.shared:  
            self.KB_L = np.array(self.KA_L)  # same stiffness as A
            self.fB_L = np.array([-self.fA_L[0], 0, self.fA_L[1]]) # mirror of fA
        else:
            self.KB_L = np.array([[K[2,2], 0 , K[2,3]],
                                  [  0   , Kt,   0   ],
                                  [K[3,2], 0 , K[3,3]]])
        
        self.KBA_L = -self.KB_L
        
        # Save tension magnitudes
        self.TA = np.linalg.norm(self.fA_L)  # tensions [N]
        self.TB = np.linalg.norm(self.fB_L)
        
        # Save rotated quantities for larger system in global reference frame
        self.fA = np.matmul(self.R, self.fA_L)
        self.fB = np.matmul(self.R, self.fB_L)
        
        self.KA  = np.matmul(np.matmul(self.R, self.KA_L ), self.R.T)
        self.KB  = np.matmul(np.matmul(self.R, self.KB_L ), self.R.T)
        self.KBA = np.matmul(np.matmul(self.R, self.KBA_L), self.R.T)
        
        #self.LBot = info["LBot"]  <<< this should be calculated considering all lines
        
        # ----- plot the profile if requested -----
        if profiles > 1:
            import matplotlib.pyplot as plt
            self.plot2d()
            plt.show()    
    
    
    def drawLine2d(self, Time, ax, color="k", endpoints=False, Xuvec=[1,0,0], Yuvec=[0,0,1], Xoff=0, Yoff=0, colortension=False, plotnodes=[], plotnodesline=[],label="",cmap='rainbow', alpha=1.0):
        '''wrapper to System.plot2d with some transformation applied'''
        
        for i, line in enumerate(self.lineList):
            
            # color and width settings
            if color == 'self':
                colorplot = line.color  # attempt to allow custom colors
                lw = line.lw
            elif color == None:
                colorplot = [0.3, 0.3, 0.3]  # if no color, default to grey
                lw = 1
            else:
                colorplot = color
                lw = 1
            
            # get the Line's local coordinates
            Xs0, Ys0, Zs, tensions = line.getLineCoords(Time)
            
            # transform to global coordinates (Ys0 should be zero for now)
            Xs = self.rA[0] + (Xs0 + self.span)*self.cos_th - Ys0*self.sin_th
            Ys = self.rA[1] + (Xs0 + self.span)*self.sin_th + Ys0*self.cos_th
            
            # apply 3D to 2D transformation to provide desired viewing angle
            Xs2d = Xs*Xuvec[0] + Ys*Xuvec[1] + Zs*Xuvec[2] + Xoff
            Ys2d = Xs*Yuvec[0] + Ys*Yuvec[1] + Zs*Yuvec[2] + Yoff
            
            
            if colortension:    # if the mooring lines want to be plotted with colors based on node tensions
                maxT = np.max(tensions); minT = np.min(tensions)
                for i in range(len(Xs)-1):          # for each node in the line
                    color_ratio = ((tensions[i] + tensions[i+1])/2 - minT)/(maxT - minT)  # ratio of the node tension in relation to the max and min tension
                    cmap_obj = cm.get_cmap(cmap_tension)    # create a cmap object based on the desired colormap
                    rgba = cmap_obj(color_ratio)    # return the rbga values of the colormap of where the node tension is
                    ax.plot(Xs2d[i:i+2], Ys2d[i:i+2], color=rgba, zorder=100)
            else:
                ax.plot(Xs2d, Ys2d, color=colorplot, lw=lw, zorder=100, label=label, alpha=alpha)
            
            if len(plotnodes) > 0:
                for i,node in enumerate(plotnodes):
                    if self.number==plotnodesline[i]:
                        ax.plot(Xs2d[node], Ys2d[node], 'o', color=colorplot, markersize=5)
            
            if endpoints == True:
                ax.scatter([Xs2d[0], Xs2d[-1]], [Ys2d[0], Ys2d[-1]], color = colorplot)


    def drawLine(self, Time, ax, color="k", endpoints=False, shadow=True, colortension=False, cmap_tension='rainbow'):
        '''wrapper to System.plot with some transformation applied'''
        
        for i, line in enumerate(self.lineList):
            
            # color and width settings
            if color == 'self':
                color = line.color  # attempt to allow custom colors
                lw = line.lw
            elif color == None:
                color = [0.3, 0.3, 0.3]  # if no color, default to grey
                lw = 1
            else:
                lw = 1
            
            # get the Line's local coordinates
            Xs0, Ys0, Zs, tensions = line.getLineCoords(Time)
            
            # transform to global coordinates (Ys0 should be zero for now)
            Xs = self.rA[0] + (Xs0 + self.span)*self.cos_th - Ys0*self.sin_th
            Ys = self.rA[1] + (Xs0 + self.span)*self.sin_th + Ys0*self.cos_th
            
            if colortension:    # if the mooring lines want to be plotted with colors based on node tensions
                maxT = np.max(tensions); minT = np.min(tensions)
                for i in range(len(Xs)-1):          # for each node in the line
                    color_ratio = ((tensions[i] + tensions[i+1])/2 - minT)/(maxT - minT)  # ratio of the node tension in relation to the max and min tension
                    cmap_obj = cm.get_cmap(cmap_tension)    # create a cmap object based on the desired colormap
                    rgba = cmap_obj(color_ratio)    # return the rbga values of the colormap of where the node tension is
                    #linebit.append(ax.plot(Xs[i:i+2], Ys[i:i+2], Zs[i:i+2], color=rgba, zorder=100))
                    ax.plot(Xs[i:i+2], Ys[i:i+2], Zs[i:i+2], color=rgba, zorder=100)
            else:
                #linebit.append(ax.plot(Xs, Ys, Zs, color=color, lw=lw, zorder=100))
                ax.plot(Xs, Ys, Zs, color=color, lw=lw, zorder=100)
            
            if shadow:
                ax.plot(Xs, Ys, np.zeros_like(Xs)-self.depth, color=[0.5, 0.5, 0.5, 0.2], lw=lw, zorder = 1.5) # draw shadow
            
            if endpoints == True:
                #linebit.append(ax.scatter([Xs[0], Xs[-1]], [Ys[0], Ys[-1]], [Zs[0], Zs[-1]], color = color))
                ax.scatter([Xs[0], Xs[-1]], [Ys[0], Ys[-1]], [Zs[0], Zs[-1]], color = color)
    
    
    def setOffset(self, offset, z=0):
        '''Moves end B of the Subsystem to represent an offset from the 
        undisplaced position of the endpoint. End A is set based on the 
        'span' (shouldn't change), and B is set based on offset and the
        rad_fair/z_fair setting. Optional argument z can be added for a z offset.
        '''
        
        # Use static EA values and unstretched lengths
        self.revertToStaticStiffness()

        # Ensure end A position and set end B position to offset values
        if self.shared:
            self.rA = np.array([-self.span/2-self.rad_fair, 0, self.rA[2]])
            self.rB = np.array([-self.rad_fair + offset/2, 0, self.z_fair+z]) 
            
        else:
            self.rA = np.array([-self.span-self.rad_fair, 0, self.rA[2]])
            self.rB = np.array([-self.rad_fair + offset, 0, self.z_fair+z]) 
            
        self.staticSolve(tol=self.eqtol)  # solve the subsystem
        
        # Store some values at this offset position that may be used later
        for i, line in enumerate(self.lineList):
            self.TeM[i,0] = np.linalg.norm(line.fA)
            self.TeM[i,1] = np.linalg.norm(line.fB)
                
        self.anchorFx0 = self.lineList[0].fA[0]
        self.anchorFz0 = self.lineList[0].fA[2]
    
        self.TeD = np.copy(self.TeM)  # set the dynamic values as well, in case they get queried right away
        
    
    def setDynamicOffset(self, offset, z=0):
        '''Moves end B of the Subsystem to represent a dynamic offset from
        the previous offset location. Uses dynamic stiffness values and also
        applies dynamic amplification factors (DAFs) on the difference from
        the mean tensions (which would have been calculated in getOffset).
        End A is set based on the 'span' (shouldn't change), 
        and B is set based on offset and the rad_fair/z_fair setting. 
        Optional argument z can be added for a z offset.
        '''
        
        if not self.dynamic_stiffness_activated: # if not already using them,
            System.activateDynamicStiffness(self)  # switch to dynamic EA values
        
        # adjust end B to the absolute offsets specified
        self.rB = np.array([-self.rad_fair + offset, 0, self.z_fair+z]) 
        
        self.staticSolve(tol=self.eqtol)  # solve the subsystem
        
        # Store dynamic values at this offset position that may be used later
        for i, line in enumerate(self.lineList):
            self.TeD[i,:] = self.TeM[i,:] + self.DAFs[i]*( np.array([line.TA, line.TB]) - self.TeM[i,:] )
        
        
    
    def activateDynamicStiffness(self, display=0):
        '''Calls the dynamic stiffness method from System rather than from Line.'''
        System.activateDynamicStiffness(self, display=display)
    
    
    def revertToStaticStiffness(self):
        '''Calls the static stiffness method from System rather than from Line.'''
        System.revertToStaticStiffness(self)
    
    
    # ---- Extra convenience functions (subsystem should be in equilibrium) -----


    def getLayLength(self, iLine=0):
        '''Computes how much of a line is on the seabed.
           Finds the total LBot of the whole Mooring looking at all lines in the Mooring
        '''

        uplift_ang = np.degrees(np.arctan2(self.lineList[iLine].fA[2], self.lineList[iLine].fA[0]))

        LBot = sum([line.LBot for line in self.lineList if line.number <= self.nLines])

        if LBot == 0.0:  # if there is zero lay length, enter a negative value equal to the negative sign of the anchor inclination
            #LBot = -self.lineList[iLine].fA[2]/np.linalg.norm(self.lineList[iLine].fA) # this is to give the optimizer a gradient to help GET out of the constraint violation

            LBot = -uplift_ang
            # UPDATED: calculate the angle of uplift when LBot is negative to keep constraint 'continuous'
            # If input constraint value is negative, that's the highest angle (degrees) of vertical uplift allowed
        
        return LBot


    def getPointHeight(self, iPoint):
        '''Computes how high a point is off the seabed.
           Default is assuming point 1.
        '''

        if iPoint==0:   # if the input index of the point that needs to be constrained is the anchor point (0)
            Xs, Ys, Zs, Ts = self.lineList[iPoint].getLineCoords(0) # get the current coordinates of the nodes in the bottom line
            height_off_seabed = Zs[1] + self.depth              # the rope point that needs to be constrained is the first node of the line away from the anchor
        else:
            height_off_seabed = self.pointList[iPoint].r[2] + self.depth # measure the height of the rope-chain connecting point

        # if it's on the seabed, include some other factor to avoid a zero jacobian
        if height_off_seabed <= 0:
            fudge_factor = -self.lineList[iPoint].LBot       # add the negative of the contact length of the next line
        else:
            fudge_factor = 0.0

        if height_off_seabed < 0:
            print("Something is wrong")
            breakpoint()

        return height_off_seabed + fudge_factor

    
    def getHog(self, iLine):
        '''Compute the z elevation of the highest point along a line section.'''
        line = self.lineList[iLine]
        return max([line.z_extreme, line.rA[2], line.rB[2]])
    
    def getSag(self, iLine):
        '''Compute the z elevation of the lowest point along a line section.'''
        line = self.lineList[iLine]
        return min([line.z_extreme, line.rA[2], line.rB[2]])
    
    def getMinSag(self):
        '''Compute the z elevation of the lowest point of the whole subsystem.'''
        return min([ min([line.z_extreme, line.rA[2], line.rB[2]]) for line in self.lineList ])

    
    def getTen(self, iLine):
        '''Compute the end (maximum) tension for a specific line, including
        a dynamic amplification factor.'''
        line = self.lineList[iLine]
        '''
        dynamicTe = max([ self.Te0[iLine,0] + self.DAFs[iLine]*(np.linalg.norm(line.fA) - self.Te0[iLine,0]),
                          self.Te0[iLine,1] + self.DAFs[iLine]*(np.linalg.norm(line.fB) - self.Te0[iLine,1])])
        '''
        dynamicTe = max( self.TeD[iLine,:] )
        
        return dynamicTe 
    
    
    def getTenSF(self, iLine):
        '''Compute MBL/tension for a specific line.'''
        return self.lineList[iLine].type['MBL'] / self.getTen(iLine) 
        

    def getMinTenSF(self, display=0):
        '''Compute minimum MBL/tension across all line sections.'''
        return min([self.getTenSF(i) for i in range(len(self.lineList))])
        
        
    def getCurv(self, iLine):
        '''Compute the maximum curvature (minimum bending radius) of a line
        section. In the future this could include a DAF.'''        
        
        # For now, this gets the max of the nodal curvatures for this line
        # which should already be calculated by calcCurvature (Is==iLine is
        # a mask). In future this could be replaced with an analytic calc.
        return max(self.Ks[self.Is==iLine]) 
    
    def getCurvSF(self, iLine):
        '''Compute minimum allowable bending radius / actual bend radius for
        a line/cable section.'''
        return 1 /( self.lineList[iLine].type['MBR'] * self.getCurv(iLine) )
    
    def getMinCurvSF(self, display=0):
        '''Compute minimum MBL/tension across all line sections.'''
        return min([self.getCurvSF(i) for i in range(len(self.lineList))])
        
    
    def getYawStiffness(self):
        '''Compute the contribution to platform yaw stiffness. Should already 
        be in equilibrium.'''
        
        tau0 = -self.fB[0]  # horizontal tension component [N]

        yaw_stiff = (tau0/l)*self.rad_fair**2 + tau0*self.rad_fair  # [N-m]
        
        return yaw_stiff


    def calcCurvature(self):
        '''Get the curvature [1/m] at each node of a Subsystem. Should already 
        be in equilibrium. Also populates nodal values for S, X, Y, Z, T, and K
        along the full length.'''
        
        n = self.nNodes
        
        self.Ss = np.zeros(n)  # arc length along line assembly (unstretched) [m]
        self.Xs = np.zeros(n)
        self.Ys = np.zeros(n)
        self.Zs = np.zeros(n)
        self.Ts = np.zeros(n)  # tensions [N]
        self.Ks = np.zeros(n)  # curvatures [1/m]
        self.Is = np.zeros(n, dtype=int)  # what line ID this node corresponds to (end A inclusive)
        
        ia = 0  # index of a Line's first node, as counted along the Subsystem
        La = 0  # length from end A of the Subsystem to end A of the Line
        
        for i, line in enumerate(self.lineList):
            ni = line.nNodes # number of nodes in this Line
            Ss = np.linspace(0, line.L, ni) + La
            Xs, Ys, Zs, Ts = line.getLineCoords(0)
            
            # paste in values for this line into the full array
            self.Ss[ia:ia+ni] = Ss
            self.Xs[ia:ia+ni] = Xs
            self.Ys[ia:ia+ni] = Ys
            self.Zs[ia:ia+ni] = Zs
            self.Ts[ia:ia+ni] = Ts
            self.Is[ia:ia+ni] = i
            
            # increment the starting index and length for the next line 
            ia += ni-1  # one point will be overwritten since lines share end points
            La += line.L
            
        
        #ds = np.diff(self.Ss)  # unstretched length of each segment
        ds = np.zeros(n-1)  # stretched length of each segment
        
        q = np.zeros([n-1, 3])  # unit vector of each segment
        
        
        # get length and direction of each segment
        for i in range(n-1):
            
            dr = np.array([self.Xs[i+1]-self.Xs[i], self.Ys[i+1]-self.Ys[i], self.Zs[i+1]-self.Zs[i]])
            
            ds[i] = np.linalg.norm(dr) # stretched length of each segment [m]
            
            q[i,:] = dr/ds[i]  # unit direction vector of the segment
            
            
        # now go through and get curvature at each node
        curvature = np.zeros(n)  # curvature at each node [1/m]
        for i in range(1,n-1):
            dl = (ds[i-1] + ds[i])/2  # average of two segment lengths
            
            curvature[i] = (2 / dl) * np.sqrt((1 - np.dot(q[i-1,:], q[i,:])) / 2)
        
        # assume the curvature of the end nodes is the same as the adjacent nodes
        curvature[ 0] = curvature[ 1]
        curvature[-1] = curvature[-2]
        
        # save values
        self.Ks = curvature

        return curvature
    
    
    def loadData(self, dirname, rootname, line_IDs):
        '''Load data from MoorDyn into the Subsystem lineList.
        
        dirname 
            folder where output files are located
        rootname
            root name of the OpenFAST primary file (.fst)
        line_IDs : list
            list of the line numbers (of MoorDyn) that correspond to the 
            lineList. These are the numbers of the line.out file to load.
        
        '''
        
        self.qs = 0  # indicate this Subsystem is reading in MoorDyn data
        
        for i, line in enumerate(self.lineList):
            line.loadData(dirname, rootname, line_IDs[i])

