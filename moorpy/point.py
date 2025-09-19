

import numpy as np

class Point():
    '''A class for any object in the mooring system that can be described by three translational coorindates'''
    
    def __init__(self, mooringSys, num, type, r, typeData = {}, m=0, v=0, a=0, fExt=np.zeros(3), DOFs=[0,1,2], d=0, zSpan=[-1,1], CdA=0.0, Ca=0.0):
        '''Initialize Point attributes

        Parameters
        ----------
        mooringSys : system object
            The system object that contains the point object
        num : int
            indentifier number
        type : int
            the point type: 0 free to move, 1 fixed, -1 coupled externally
        r : array
            x,y,z coorindate position vector [m].
        typeData : dict, optional
            structure that holds the PointType info (see setPointType)
        m : float, optional
            mass [kg]. The default is 0.
        v : float, optional
            submerged volume [m^3]. The default is 0.
        a : float, optional
            plate area if point is a plate anchor. For use in cost calcs.
        CdA : float, optional
            Product of drag coefficient and cross sectional area in any direction [m^2]. The default is 0.
        Ca : float, optional
            Added mass coefficient in any direction.
        fExt : array, optional
            applied external force vector in global orientation (not including weight/buoyancy). The default is np.zeros(3).
        DOFs: list
            list of which coordinate directions are DOFs for this point (default 0,1,2=x,y,z). E.g. set [2] for vertical motion only.
        d : float, optional
            diameter [m]. The default is 0.
        zSpan : [float, float], optional
            The lower and upper limits of the Point's volume relative to its coordinate [m]. 
            This only affects the change in buoyancy when crossing the free surface. The 
            default is [-1,1], i.e. a 2-m tall volume.
        
        Returns
        -------
        None.

        '''
    
        self.sys    = mooringSys        # store a reference to the overall mooring system (instance of System class)
    
        self.number = num
        self.type = type                # 1: fixed/attached to something, 0 free to move, or -1 coupled externally
        self.r = np.array(r, dtype=float)
        self.entity = typeData          # dict for entity (e.g. anchor, point, buoy) info. Should be pointType structure. See return from helpers.getPointProps
        self.cost = {}                  # empty dictionary to contain cost info
        self.loads = {}                 # empty dictionary to contain load info
        
        self.m  = float(m)
        self.a  = float(a)
        self.v  = float(v)
        self.CdA= float(CdA)
        self.Ca = float(Ca)
        self.fExt = fExt                # external forces plus weight/buoyancy
        self.fBot = 10.0                # this is a seabed contact force that will be added if a point is specified below the seabed
        self.zSub = 0.0                 # this is the depth that the point is positioned below the seabed (since r[2] will be capped at the depth)
        self.zTol = 2.0                 # depth tolerance to be used when updating the point's position relative to the seabed
        
        self.DOFs = DOFs
        self.nDOF = len(DOFs)
        
        self.d = d                      # the diameter of the point, if applicable. Used for hydrostatics [m]
        
        self.attached     = []         # ID numbers of any Lines attached to the Point
        self.attachedEndB = []         # specifies which end of the line is attached (1: end B, 0: end A)
        
        self.cable = False             # specifies whether the point should be modeled as a Rod (for dynamic cables) or not
    
        if len(zSpan)==2:
            self.zSpan = np.array(zSpan, dtype=float)
        else:
            raise ValueError("Point zSpan parameter must contain two numbers.")
    
        #print("Created Point "+str(self.number))
    
    
    def attachLine(self, lineID, endB):
        '''Adds a Line end to the Point

        Parameters
        ----------
        lineID : int
            The identifier ID number of a line
        endB : boolean
            Determines which end of the line is attached to the point
        '''
    
        self.attached.append(lineID)  
        self.attachedEndB.append(endB)
        #print("attached Line "+str(lineID)+" to Point "+str(self.number))
    
    def detachLine(self, lineID, endB):
        '''Detaches a Line end from the Point

        Parameters
        ----------
        lineID : int
            The identifier ID number of a line
        endB : boolean
            Determines which end of the line is to be detached from the point
        '''
        
        # get attachment index
        i1 = self.attached.index(lineID)
        i2 = self.attachedEndB.index(endB)
        if not i1==i2:
            raise Exception("issue with the right end of the line to detach...")
        
        self.attached.pop(i1)
        self.attachedEndB.pop(i1)
        print("detached Line "+str(lineID)+" from Point "+str(self.number))
    
    
    def setPosition(self, r):
        '''Sets the position of the Point, along with that of any dependent objects.

        Parameters
        ----------
        r : array
            x,y,z coordinate position vector of the point [m]

        Raises
        ------
        ValueError
            If the length of the input r array is not of length 3

        Returns
        -------
        None.

        '''
        
        # update the position of the Point itself
        if len(r) == 3:   # original case, setting all three coordinates as normal, asuming x,y,z
            self.r = np.array(r)
        elif len(r) == self.nDOF:
            self.r[self.DOFs] = r          # this does a mapping based on self.DOFs, to support points with e.g. only a z DOF or only x and z DOFs
        else:
            raise ValueError(f"Point setPosition method requires an argument of size 3 or nDOF, but size {len(r):d} was provided")
        
        # update the point's depth and position based on relation to seabed
        depth, _ = self.sys.getDepthFromBathymetry(self.r[0], self.r[1]) 
        
        self.zSub = np.max([-self.zTol, -self.r[2] - depth])   # depth of submergence in seabed if > -zTol
        self.r = np.array([self.r[0], self.r[1], np.max([self.r[2], -depth])]) # don't let it sink below the seabed
        
        # update the position of any attached Line ends
        for LineID,endB in zip(self.attached,self.attachedEndB):
            self.sys.lineList[LineID-1].setEndPosition(self.r, endB)
            
        if len(self.r) < 3:
            print("Double check how this point's position vector is calculated")
            breakpoint()
            
            
    
    def getForces(self, lines_only=False, seabed=True, xyz=False):
        '''Sums the forces on the Point, including its own plus those of any attached Lines.

        Parameters
        ----------
        lines_only : boolean, optional
            An option for calculating forces from just the mooring lines or not. The default is False.
        seabed : bool, optional
            if False, will not include the effect of the seabed pushing the point up
        xyz : boolean, optional
            if False, returns only forces corresponding to enabled DOFs. If true, returns forces in x,y,z regardless of DOFs. 
            
        Returns
        -------
        f : array
            The force vector applied to the point in its current position [N]

        '''
    
        f = np.zeros(3)         # create empty force vector on the point
        
        if lines_only==False:
            '''
            radius = self.d/2           # can do this, or find the radius using r=(3*self.v/(4*np.pi))**(1/3)
            x = max(0, radius**2 - self.r[2]**2)
            dWP = 2*np.sqrt(x)          # diameter at the waterplane [m]
            AWP = np.pi/4 * dWP**2      # waterplane area [m]
            #v_half = (4/3)*np.pi*(np.sqrt(x)**3) * 0.5  # volume of the half sphere that is cut by the waterplane [m^3]
            #v = abs(-min(0, np.sign(self.r[2]))*self.v - v_half)    # submerged volume of the point [m^3]
            '''
            f[2] += -self.m*self.sys.g  # add weight 

            #f[2] += self.v*self.sys.rho*self.sys.g   # add buoyancy using submerged volume
            
            if self.r[2] + self.zSpan[1] < 0.0:                # add buoyancy if fully submerged
                f[2] +=  self.v*self.sys.rho*self.sys.g
            elif self.r[2] + self.zSpan[0] < 0.0:    # add some buoyancy if part-submerged (linear variation, constant Awp)
                f[2] +=  self.v*self.sys.rho*self.sys.g * (self.r[2] + self.zSpan[0])/(self.zSpan[0]-self.zSpan[1])
            # (no buoyancy force added if it's fully out of the water, which would be very exciting for the Point)
            
            f += np.array(self.fExt) # add external forces
            #f[2] -= self.sys.rho*self.sys.g*AWP*self.r[2]   # hydrostatic heave stiffness
            
            # handle case of Point resting on or below the seabed, to provide a restoring force
            # add smooth transition to fz=0 at seabed (starts at zTol above seabed)
            f[2] += max(self.m - self.v*self.sys.rho, 0)*self.sys.g * (self.zSub + self.zTol)/self.zTol

                
        # add forces from attached lines
        for LineID,endB in zip(self.attached,self.attachedEndB):
            # f += self.sys.lineList[LineID-1].getEndForce(endB)
            if endB:
                f += self.sys.lineList[LineID-1].fB
            else:
                f += self.sys.lineList[LineID-1].fA
        
        if xyz:
            return f
        else:
            return f[self.DOFs]    # return only the force(s) in the enable DOFs
        
    
    
    def getStiffness(self, X = [], tol=0.0001, dx = 0.01):
        '''Gets the stiffness matrix of the point due only to mooring lines with all other objects free to equilibrate.
        NOTE: This method currently isn't set up to worry about nDOF and DOFs settings of the Point. It only works for DOFs=[0,1,2].

        Parameters
        ----------
        X1 : array
            The position vector of the Point at which the stiffness matrix is to be calculated.
        dx : float, optional
            The change in displacement to be used for calculating the change in force. The default is 0.01.

        Returns
        -------
        K : matrix
            The stiffness matrix of the point at the given position X1.

        '''
        
        #print("Getting Point "+str(self.number)+" stiffness matrix...")
        
        if len(X) == 3:
            X1 = np.array(X)
        elif len(X)==0:
            X1 = self.r
        else:
            raise ValueError('Point.getStiffness expects the optional X parameter to be size 3')
        
        # set this Point's type to fixed so mooring system equilibrium response to its displacements can be found
        type0 = self.type                         # store original type to restore later
        self.type = 1                             # set type to 1 (not free) so that it won't be adjusted when finding equilibrium
        
        # if this Point is attached to a Body, set that Body's type to fixed so equilibrium can be found
        for body in self.sys.bodyList:            # search through all the bodies in the mooring system
            if self.number in body.attachedP:     # find the one that this Point is attached to (if at all)
                num = body.number                 # store body number to index later
                Btype0 = body.type                # store original body type to restore later
                body.type = 1                     # set body type to 1 (not free) so that it won't be adjusted when finding equilibrium 
        
        # ensure this Point is positioned at the desired linearization point
        self.setPosition(X1)                      # set position to linearization point
        self.sys.solveEquilibrium3(tol=tol)       # find equilibrium of mooring system given this Point in current position
        f = self.getForces(lines_only=True)       # get the net 6DOF forces/moments from any attached lines 

        # Build a stiffness matrix by perturbing each DOF in turn
        K = np.zeros([3,3])
        
        for i in range(len(K)):
            X2 = X1 + np.insert(np.zeros(2),i,dx) # calculate perturbed Point position by adding dx to DOF in question            
            self.setPosition(X2)                  # perturb this Point's position
            self.sys.solveEquilibrium3(tol=tol)   # find equilibrium of mooring system given this Point's new position
            f_2 =self.getForces(lines_only=True)  # get the net 3DOF forces/moments from any attached lines 

            K[:,i] = -(f_2-f)/dx                  # get stiffness in this DOF via finite difference and add to matrix column
            
        # ----------------- restore the system back to previous positions ------------------
        self.setPosition(X1)                      # set position to linearization point
        self.sys.solveEquilibrium3(tol=tol)       # find equilibrium of mooring system given this Point in current position
        self.type = type0                         # restore the Point's type to its original value
        for body in self.sys.bodyList:
            if self.number in body.attachedP:
                num = body.number
                self.sys.bodyList[num-1].type = Btype0    # restore the type of the Body that the Point is attached to back to original value

        
        return K
    
    
    
    def getStiffnessA(self, lines_only=False, xyz=False):
        '''Gets analytical stiffness matrix of Point due only to mooring lines with other objects fixed.

        Returns
        -------
        K : matrix
            3x3 analytic stiffness matrix.

        '''
        
        #print("Getting Point "+str(self.number)+" analytic stiffness matrix...")
        
        K = np.zeros([3,3])         # create an empty 3x3 stiffness matrix
        
        # append the stiffness matrix of each line attached to the point
        for lineID,endB in zip(self.attached,self.attachedEndB):
            line = self.sys.lineList[lineID-1]
            #KA, KB = line.getStiffnessMatrix()
            
            if endB == 1:                  # assuming convention of end A is attached to the point, so if not,
                #KA, KB = KB, KA            # swap matrices of ends A and B                                
                K += line.KB
            else:
                K += line.KA 
            
        # NOTE: can rotate the line's stiffness matrix in either Line.getStiffnessMatrix() or here in Point.getStiffnessA()
        
        # add seabed or hydrostatic terms if needed
        if lines_only==False:
        
            # if partially submerged, apply a hydrostatic stiffness based on buoyancy
            if self.r[2] + self.zSpan[1] > 0.0 and self.r[2] + self.zSpan[0] < 0.0: 
                K[2,2] += self.sys.rho*self.sys.g * self.v/(self.zSpan[1]-self.zSpan[0])  # assumes volume is distributed evenly across zSpan
            
            # if less than zTol above the seabed (could even be below the seabed), apply a stiffness (should bring wet weight to zero at seabed)
            if self.r[2] < self.zTol - self.sys.depth:
                K[2,2] += max(self.m - self.v*self.sys.rho, 0)*self.sys.g / self.zTol
                
            # if on seabed, apply a large stiffness to help out system equilibrium solve (if it's transitioning off, keep it a small step to start with)    
            if self.r[2] == -self.sys.depth:
                K[2,2] += 1.0e12
        if sum(np.isnan(K).ravel()) > 0:
            raise ValueError("Something is wrong")
        if xyz:                     # if asked to output all DOFs, do it
            return K
        else:                       # otherwise only return rows/columns of active DOFs
            return K[:,self.DOFs][self.DOFs,:]

    def getDynamicMatrices(self):
        '''Gets inertia, added mass, damping, and stiffness matrices of Point due only to mooring lines (with other objects fixed)
        using a lumped mass model.

        Returns
        -------
        '''
        M = np.zeros([3,3])
        A = np.zeros([3,3])
        B = np.zeros([3,3])
        K = np.zeros([3,3])

        # append the stiffness matrix of each line attached to the point
        for lineID,endB in zip(self.attached,self.attachedEndB):            
            line = self.sys.lineList[lineID-1]     

            M_all, A_all, B_all, K_all = line.getDynamicMatricesLumped()
             
            M += M_all[-3:,-3:] if endB == 1 else M_all[:3, :3]
            A += A_all[-3:,-3:] if endB == 1 else A_all[:3, :3]
            B += B_all[-3:,-3:] if endB == 1 else B_all[:3, :3]
            K += K_all[-3:,-3:] if endB == 1 else K_all[:3, :3]

        return M, A, B, K
        
    def getCost(self):
        '''Calls getCost_and_MBL, for backwards compatability
        Returns the outputs of moorprops.getAnchorCost()

        Returns
        -------
        anchorMatCost : float
            The anchor material cost
        anchorInstCost : float
            The anchor installation cost
        anchorDecomCost : float
            The anchor decomissioning cost
        info : dict
            An info dictionary
        '''
        self.getCost_and_MBL()
        
        # Returns outputs from getAnchorCost: anchorMatCost, anchorInstCost, anchorDecomCost, info 
        return self.entity['INFO']['Anchors - anchorMatCost'], self.entity['INFO']['Anchors - anchorInstCost'], self.entity['INFO']['Anchors - anchorDecomCost'], self.entity['INFO']['Anchors - getAnchorCost info']
        
    def getCost_and_MBL(self, fx = 0.0, fz = 0.0, peak_tension = None, buoyancy = None):
        '''Calculates the cost and MBL of a point defined by the point.entity structure. 

        If point.entity is a pointType structure, this gives the total cost including all 
        the components in the point design

        Parameters
        ----------
        fx : float (optional)
            The maximum horizontal force on the point [N]
        fz : float (optional)
            The maximum vertical force on the point [N]
        peak_tension : float (optional)
            The peak tension seen at the point. Used for calculating the MBL of points. If not given, vector sum of fx and fz is used. [N]
        buoyancy : float (optional) 
            The total buoyancy of the point used for buoy costs. This is divided amoung the buoys using frac_b_key. [N]

        Returns
        -------
        cost : float
            The total cost of the point [$]
        MBL : float
            The minimum MBL of the point [N]
        INFO : dict
            An info dictionary with things like the outputs of getAnchorCost
        '''
        
        self.entity['INFO'] = {} # info structure

        self.cost = 0.0 # clear any old cost numbers and start with 0

        if any((key != "INFO" and key != "type" and key != "anchor_type") for key in self.entity): # if there are keys other than "INFO" and "type" and "anchor_type" then this must be PointType Struct
            
            # ------ Anchors ------
            '''
            Unlike buoys and connections, getPointProps does not agregate anchor curves because of unique use cases and UHC data being in MoorProps. 
            This is handled by getAnchorCost and the code below, hence why there is more code here than buoy and connection sections. 
            '''

            aUHC = 0.0 # zero initial
            if self.entity["Anchors"]:
                
                from moorpy.MoorProps import getAnchorCost

                aUHC_list = []
                for anchor in self.entity["anchor_list"]:
                                    
                    if fx > 0.0 or fz > 0.0: # if forces given, find anchor size and calc cost via getAnchorCost
                        aCost = getAnchorCost(type = anchor["name"], fx = fx, fz = fz, aprops = self.entity["aprops"]) # not looking for anchor UHC

                        aUHC_list.append(aCost[-1]["UHC"]) # extract UHC from the info dict
                        
                    else: # use the point mass or area as defined by user. Mass and area multiplied by fraction of total size per anchor type. Default is frac_a_key/num_a_key where frac_a_key is 1/total number of anchor types in point

                        aCost = getAnchorCost(type = anchor["name"], mass = self.m * (anchor["frac"] / anchor["num"]), area = self.a * (anchor["frac"] / anchor["num"]), aprops = self.entity["aprops"]) # not looking for anchor UHC

                        self.entity['INFO']["Anchors - Note"] = "UHC not included in point MBL calc as no data was provided"

                    if sum(aCost[:3]) == 0.0 and (self.m > 0.0 or self.a > 0.0):
                        raise ValueError(f"{anchor['name']} anchor costs are not yet supported when mass = {self.m * (anchor['frac'] / anchor['num'])} kg and area = {self.a  * (anchor['frac'] / anchor['num'])} m^2")
                    
                    self.cost += anchor["num"]*np.sum(aCost[:3]) # number of these anchors * unit cost, where aCost[:3] is the mat, inst, and decom costs

                    if self.a > 0.0:
                        self.entity['INFO']["Anchors - Area"] = self.a

                    # store the outputs of getAnchorCost for backwards compatability with getCost function
                    self.entity['INFO']['Anchors - anchorMatCost']      = aCost[0]
                    self.entity['INFO']['Anchors - anchorInstCost']     = aCost[1]
                    self.entity['INFO']['Anchors - anchorDecomCost']    = aCost[2]
                    self.entity['INFO']['Anchors - getAnchorCost info'] = aCost[3]

                # minimum aUHC reported back, 0.0 if empty list
                if len(aUHC_list) > 0: 
                    aUHC = np.min(aUHC_list)
                    
            # ------ Buoys ------
            if self.entity["Buoys"]: 
                if buoyancy == None:
                    # TODO: can we calculate the buoyancy based on the point mass and volume and use that instead?
                    raise ValueError("Buoyancy is required to find the cost of a point with buoys")
                
                self.cost += self.entity["buoy_cost"]["cost_b0"] + self.entity["buoy_cost"]["cost_b1"] * buoyancy + self.entity["buoy_cost"]["cost_b2"] * buoyancy**2 + self.entity["buoy_cost"]["cost_b3"] * buoyancy**3

            # ------ Connections ------
            cMBL = 0.0 # intitialize connect MBL to use for MBL
            if self.entity["Connections"]: 
                if peak_tension == None:
                    if fx > 0.0 or fz > 0.0: # if forces are provided use those to find peak tension
                        peak_tension = np.sqrt(fx**2 + fz**2)
                    else:
                        raise ValueError("Peak tension or fx and fz is required to find the cost of a point with connection hardware")
                
                self.cost += self.entity["connector_cost"]["cost_load0"] + self.entity["connector_cost"]["cost_load1"] * peak_tension + self.entity["connector_cost"]["cost_load2"] * peak_tension**2 + self.entity["connector_cost"]["cost_load3"] * peak_tension**3

                cMBL = self.entity["FOS"]*peak_tension

            # MBL of point is smallest capacity between component MBL and anchor UHC
            if aUHC == 0 and cMBL > 0:
                MBL = cMBL # MBL [N]
            elif aUHC > 0 and cMBL == 0:
                MBL = aUHC # MBL [N]
            else:
                MBL = np.min([cMBL, aUHC]) # MBL [N]
            
        else: # No other keys this must be the old method
            '''Fill in and returns a cost dictionary for this Point object.
            So far it only applies for if the point is an anchor. This is an 
            old method.
            '''

            self.entity["INFO"]["Anchors - Note"] = "Costs calculated using max fx and fz loads on point if point is on seabed (old method)"
            from moorpy.MoorProps import getAnchorCost
            
            # figure out if it should be an anchor if it isn't already defined
            if self.entity['type'] == '':
                depth, _ = self.sys.getDepthFromBathymetry(self.r[0], self.r[1]) 
                if self.r[3] == depth and self.type==1:  # if it's fixed on the seabed
                    self.entity['type'] = 'anchor'       # assume it's an anchor
                    if self.FA[2] == 0:
                        self.entity['anchor_type'] = 'drag-embedment'
                    else:
                        self.entity['anchor_type'] = 'suction'
            
            # calculate costs if it's an anchor (using simple model)
            if self.entity['type'] == 'anchor':
                aCost = getAnchorCost(self.loads['fx_max'], 
                                                    self.loads['fz_max'],
                                                type=self.entity['anchor_type'])
                
            self.cost = np.sum(aCost[:3]) # unit cost, including material, install, decomissioning

            # store the outputs of getAnchorCost for backwards compatability with getCost function
            self.entity['INFO']['Anchors - anchorMatCost']      = aCost[0]
            self.entity['INFO']['Anchors - anchorInstCost']     = aCost[1]
            self.entity['INFO']['Anchors - anchorDecomCost']    = aCost[2]
            self.entity['INFO']['Anchors - getAnchorCost info'] = aCost[3]

        return self.cost, MBL, self.entity["INFO"]


