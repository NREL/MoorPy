
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
        #self.rodTypes = {}
        
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
        
        # water current - currentMod 0 = no current; 1 = steady uniform current
        self.currentMod = 0         # flag for current model to use
        self.current = np.zeros(3)  # current velocity vector [m/s]
        if 'current' in kwargs:
            self.currentMod = 1
            self.current = getFromDict(kwargs, 'current', shape=3)
            
        self.shared     = getFromDict(kwargs, 'shared', dtype=bool, default=False)  # flag to indicate shared line
        self.span    = getFromDict(kwargs, 'span', default=0)                 # spacing (to rename as span<<<)
        self.rBFair     = getFromDict(kwargs, 'rBFair', shape=-1, default=[0,0,0])  # [m] end coordinates relative to attached body's ref point

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
    
    
    def makeGeneric(self, lengths, types, points=[], shared=False):
        '''Creates a cable of n components going between an anchor point and
        a floating body (or a bridle point). If shared, it goes between two
        floating bodies.

        Parameters
        ----------
        lengths
            List of each section's length. This also implies the number of
            sections.
        types
            List of lineType names for each section. These names must match
            keys in the parent system lineTypes dictionary or the subsystem's lineTypes dictionary...
        i_buoy
            list of which sections have buoyancy
        cable_type_name
            name of the cable type of the cable
        '''
        
        # some initialization steps.
        self.nLines = len(lengths)
        
        if not len(types)==self.nLines:
            raise Exception("The specified number of lengths and types is inconsistent.")
        
        # get cumulative sum of line lengths, starting from anchor segment
        Lcsum = np.cumsum(np.array(lengths))

        if shared:
            rA = np.array([-0.5*self.span, 0, -1])           # shared line midpiont coordinates
        else:
            rA = np.array([-self.span, 0, -self.depth])    # anchor coordinates
        rB = np.array([-self.rBFair[0], 0, self.rBFair[2]])     # fairlead coordinates

        # add the imaginary body and anchor/mid point
        #self.addBody(-1, np.zeros(6)) <<<<<<<<<<<<<<<  we'll probaly NOT have a body, just do calcs as if there were one
        self.rA = rA
        self.rB = rB

        if shared:
            self.addPoint(-1, rA, DOFs=[2]) # add shared line point, free only to move in z
        else:
            self.addPoint(-1, rA, DOFs=[0,2])  # add anchor point

        # Go through each line segment and add its upper point, add the line, and connect the line to the points
        for i in range(self.nLines):

            # find the specified lineType dict and save a reference to it
            if types[i] in self.lineTypes:  # first look for the name in the subsystem
                self.lineTypes[i] = self.lineTypes[types[i]]
            elif self.sys: # otherwise look in the parent system, if there is one
                if types[i] in self.sys.lineTypes:  # first look for the name in the subsystem
                    self.lineTypes[i] = self.sys.lineTypes[types[i]]
                else:
                    raise Exception(f"Can't find lineType '{types[i]}' in the SubSystem or parent System.")
            else:
                raise Exception(f"Can't find lineType '{types[i]}' in the SubSystem.")
            
            # add the line segment using the reference to its lineType dict
            self.addLine(lengths[i], self.lineTypes[i])

            # add the upper end point of the segment
            if i==self.nLines-1:                            # if this is the upper-most line
                self.addPoint(-1, rB, DOFs=[0,2])                        # add the fairlead point (make it coupled)
                #self.bodyList[0].attachPoint(i+2, rB)       # attach the fairlead point to the body (two points already created)
            else:                                           # if this is an intermediate line
                # add the point, initializing linearly between anchor and fairlead/midpoint
                self.addPoint(0, rA + (rB-rA)*Lcsum[i]/Lcsum[-1], DOFs=[0,2])


            # attach the line to the points
            self.pointList[-2].attachLine(i+1, 0)       # attach end A of the line
            self.pointList[-1].attachLine(i+1, 1)       # attach end B of the line

        # if a points list is provided, apply any mass or other properties it contains?
    
        self.nLines = len(self.lineList)
        self.nNodes = np.sum([line.nNodes for line in self.lineList]) - self.nLines + 1
    
    
    def makeFromSystem(self, sys, point_id):
        '''Populate a Subsystem based on a series of mooring lines in an
        existing system (sys), starting at an end point (point_id).
        Taken from Serag's CompositeLine class.
        In this version, the subsystem still needs to be added to a system
        afterward.
        '''
        
        if len(self.pointList)+len(self.lineList) > 0:
            raise Exception('makeFromSystem can only be called for an empty Subsystem.')
        
        point = sys.pointList[point_id-1] # Starting point id

        # check starting point to make sure it is either coupled or fixed and that it is not connected to more than 1 line.
        if point.type == 0:
            raise Exception(f'Starting point ({point.number}) must not be free (change point type to -1 or 1).')
        elif len(point.attached)>1:
            raise Exception(f'Starting point ({point.number}) cannot be connected to more than 1 line')

        # Save point A info
        self.rA = np.array(point.r)
        self.addPoint(1, self.rA)  # add anchor point

        # Move through the points along the composite
        while True:
            # make sure that the point is free
            if len(point.attached) > 2:
                raise Exception(f'Point {point.number} is attached to more than two lines.')
            
            # get the line id and line object
            line_id = point.attached[-1] # get next line's id
            line = sys.lineList[line_id - 1] # get line object
            self.lineTypes[line.type['name']] = dict(line.type)  # copy the lineType dict
            self.addLine(line.L0, self.lineTypes[line.type['name']]) # add the line
            

            # get the next point
            attached_points = line.attached.copy() # get the IDs of the points attached to the line
            pointA_id = point.number # get first point ID
            attached_points.remove(pointA_id) # remove first point from attached point list
            pointB_id = attached_points[0] # get second point id
            point = sys.pointList[pointB_id-1] # get second point object

            if line(point.attached) == 1:  # must be the endpoint
                self.addPoint(-1, point.r)
            else:  # intermediate point along line
                self.addPoint(0, point.r, DOFs=[0,2]) # may need to ensure there's no y component
                
            # Attach the line ends to the points
            self.pointList[-2].attachLine(len(self.lineList), 0)
            self.pointList[-1].attachLine(len(self.lineList), 1)
             
            # break from the loop when a point with a single attachment is reached
            if len(point.attached) == 1:
                break
        
        # make sure that the last point is not a free point
        if point.type == 0:
            raise Exception(f'Last point ({point.number}) is a free point.')
        
        self.rB = np.array(point.r)
        self.nLines = len(self.lineList)
        self.nNodes = np.sum([line.nNodes for line in self.lineList]) - self.nLines + 1
    
        
    def initialize(self):
        '''Initializes the subsystem objects to their initial positions, and
        counts number of nodes.'''
        
        self.nDOF, self.nCpldDOF, _ = self.getDOFs()
        
        self.nNodes = np.sum([line.nNodes for line in self.lineList]) - self.nLines + 1
        
        for point in self.pointList:
            point.setPosition(point.r)
            
        self.staticSolve()
    
    
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
            self.rB = np.array(r, dtype=np.float_)
        elif endB == 0:
            self.rA = np.array(r, dtype=np.float_)
        else:
            raise LineError("setEndPosition: endB value has to be either 1 or 0")
    
    
    def staticSolve(self, reset=False, tol=0.01, profiles=0):
        '''Solve internal equilibrium of the Subsystem and saves the forces
        and stiffnesses at the ends in the global reference frame. All the 
        This method mimics the behavior of the Line.staticSolve method and 
        the end coordinates rA and rB must be set beforehand. The solve 
        equilibrium happens in the local 2D plane. Values in this local 
        frame are also saved. 
        '''
        
        # transform end positions to SubSystem internal coordinate system
        # inputs are self.rA, rB in global frame
        # outputs should be pointList[0] and [N] .r
        dr =  self.rB - self.rA
        LH = np.hypot(dr[0], dr[1])         # horizontal spacing of line ends
        self.pointList[ 0].setPosition([ -self.span   , 0, self.rA[2]])
        self.pointList[-1].setPosition([ -self.span+LH, 0, self.rB[2]])
        
        # get equilibrium
        self.solveEquilibrium(tol=tol)
        
        # get 2D stiffness matrices of end points
        K = self.getCoupledStiffness()
        
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
        '''
        self.KA_L  = K[:3,:3]                # reaction at A due to motion of A
        self.KB_L  = K[3:,3:]                # reaction at B due to motion of B
        self.KAB_L = K[3:,:3]                # reaction at B due to motion of A
        '''
        
        # expand to get 3D stiffness matrices
        R = np.eye(3)
        self.KA_L  = from2Dto3Drotated(K[:2,:2], -self.fB_L[0], LH, R.T)  # reaction at A due to motion of A
        self.KB_L  = from2Dto3Drotated(K[2:,2:], -self.fB_L[0], LH, R.T)  # reaction at B due to motion of B
        self.KBA_L = from2Dto3Drotated(K[2:,:2], -self.fB_L[0], LH, R.T)  # reaction at B due to motion of A
        
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
    
    
    def drawLine2d(self, Time, ax, color="k", Xuvec=[1,0,0], Yuvec=[0,0,1], Xoff=0, Yoff=0, colortension=False, cmap='rainbow'):
        '''wrapper to System.plot2d with some transformation applied'''
        
        for i, line in enumerate(self.lineList):
            
            Xs0, Ys0, Zs, tensions = self.getLineCoords(Time)
            
            # transform to global coordinates
            Xs = self.rA[0] + Xs0*self.cos_th - Ys0*self.sin_th
            Ys = self.rA[1] + Xs0*self.sin_th + Ys0*self.cos_th
            
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
                ax.plot(Xs2d, Ys2d, color=color, lw=lw, zorder=100)
            
            if endpoints == True:
                ax.scatter([Xs2d[0], Xs2d[-1]], [Ys2d[0], Ys2d[-1]], color = color)


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
        rBFair setting. Optional argument z can be added for a z offset.
        '''
        
        self.rA = np.array([-self.span, 0, self.rA[2]])
        self.rB = np.array([-self.rBFair[0] + offset, 0, self.rBFair[2]+z]) 
        
        # should we also do the static solve now?
        
    
    def activateDynamicStiffness(self, display=0):
        '''Switch mooring system model to dynamic line stiffness
        values and adjust the unstretched line lengths to maintain the
        same tensions. This only has an effect when dynamic line properties
        are used.'''
        System.activateDynamicStiffness(self, display=display)
    
    
    def revertToStaticStiffness(self):
        '''Revert mooring system model back to the static stiffness
        values and the original unstretched lenths.'''
        System.revertToStaticStiffness(self)
    
    
    # ----- Function for dynamic frequency-domain tensions -----
    
    def getDynamicMatrices(self, omegas, S_zeta, r_dynamic, depth, kbot, cbot,
                           seabed_tol=1e-4):
        '''Compute M,A,B,K matrices for the Subsystem. This calls 
        get_dynamic_matrices() for each Line in the Subsystem then combines
        the results. Note that this method overrides the Line method. Other
        Line methods used for dynamics can be used directly in Subsystem.
        '''
        self.nNodes = np.sum([line.nNodes for line in self.lineList]) - self.nLines + 1
        
        EA_segs = np.zeros(self.nNodes-1) # extensional stiffness of the segments
        n_dofs = 3*self.nNodes # number of dofs
        M = np.zeros([n_dofs,n_dofs], dtype='float')
        A = np.zeros([n_dofs,n_dofs], dtype='float')
        B = np.zeros([n_dofs,n_dofs], dtype='float')
        K = np.zeros([n_dofs,n_dofs], dtype='float')
        r_mean = np.zeros([self.nNodes,3], dtype='float')
        r_dynamic = np.ones((len(omegas),self.nNodes,3),dtype='float')*r_dynamic
        v_dynamic = 1j*omegas[:,None,None]*r_dynamic

        n = 0  # starting index of the next line's entries in the matrices
        
        for line in self.lineList:
            n1 = int(n/3) 
            n2 = n1 + line.nNodes

            # Filling matrices for line (dof n to dof 3xline_nodes+n)
            M_n,A_n,B_n,K_n,r_n,EA_segs_n = line.getDynamicMatrices(omegas, S_zeta,r_dynamic[:,n1:n2,:],depth,kbot,cbot,seabed_tol=seabed_tol)
            M[n:3*line.nNodes+n,n:3*line.nNodes+n] += M_n
            A[n:3*line.nNodes+n,n:3*line.nNodes+n] += A_n
            B[n:3*line.nNodes+n,n:3*line.nNodes+n] += B_n
            K[n:3*line.nNodes+n,n:3*line.nNodes+n] += K_n

            # Attachment point properties
            attachment = self.pointList[line.attached[-1]-1] # attachment point
            attachment_idx = n2 - 1 # last node index
            sigma_vp = np.sqrt(np.trapz(np.abs(v_dynamic[:,attachment_idx,:])**2*S_zeta[:,None],omegas,axis=0)) # standard deviations of the global components of the attachment point's velocity

            M[3*line.nNodes - 3:3*line.nNodes, 3*line.nNodes - 3:3*line.nNodes] += attachment.m*np.eye(3) 
            A[3*line.nNodes - 3:3*line.nNodes, 3*line.nNodes - 3:3*line.nNodes] += attachment.Ca* self.rho * attachment.v * np.eye(3)
            B[3*line.nNodes - 3:3*line.nNodes, 3*line.nNodes - 3:3*line.nNodes] += 0.5* self.rho * attachment.CdA * np.pi*(3*attachment.v/4/np.pi)**(2/3) \
                                                                                   * np.eye(3) * np.sqrt(8/np.pi) * np.diag(sigma_vp)

            # Static line properties
            r_mean[n1:n2,:] = r_n
            EA_segs[n1:n2-1] = EA_segs_n 

            n += 3*line.nNodes - 3 # next line starting node add the number of dofs of the current line minus 3 to get the last shared node
        
        return M,A,B,K,r_mean,EA_segs
    
    
    # ---- Extra convenience functions (subsystem should be in equilibrium) -----


    def getLayLength(self, iLine=0):
        '''Computes how much of a line is on the seabed.
           Finds the total LBot of the whole Mooring looking at all lines in the Mooring
        '''

        uplift_ang = np.degrees(np.arctan(self.lineList[iLine].fA[2]/self.lineList[iLine].fA[0]))

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
        return max([line.info['Zextreme'], line.rA[2], line.rB[2]])
    
    def getSag(self, iLine):
        '''Compute the z elevation of the lowest point along a line section.'''
        line = self.lineList[iLine]
        return min([line.info['Zextreme'], line.rA[2], line.rB[2]])
    
    def getMinSag(self):
        '''Compute the z elevation of the lowest point of the whole subsystem.'''
        return min([line.info['Zextreme'] + 
                    min(line.rA[2], line.rB[2]) for line in self.lineList])

    
    def getTen(self, iLine):
        '''Compute the end (maximum) tension for a specific line, including
        a dynamic amplification factor.'''        
        line = self.lineList[iLine]
        
        dynamicTe = max([ self.Te0[iLine,0] + self.DAFs[iLine]*(np.linalg.norm(line.fA) - self.Te0[iLine,0]),
                          self.Te0[iLine,1] + self.DAFs[iLine]*(np.linalg.norm(line.fB) - self.Te0[iLine,1])])
        
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

        yaw_stiff = (tau0/l)*self.rBFair[0]**2 + tau0*self.rBFair[0]  # [N-m]
        
        return yaw_stiff


    def calcCurvature(self):
        '''Get the curvature [1/m] at each node of a Subsystem. Should already 
        be in equilibrium. Also populates nodal values for S, X, Y, Z, T, and K
        along the full length.'''
        
        self.nNodes = np.sum([l.nNodes for l in self.lineList]) - self.nLines + 1
        
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
  