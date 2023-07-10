import pdb

import numpy as np
from matplotlib import cm
from moorpy.Catenary import catenary
from moorpy.nonlinear import nonlinear                                      
from moorpy.helpers import (unitVector, LineError, CatenaryError, 
                     rotationMatrix, makeTower, read_mooring_file, 
                     quiver_data_to_segments, printVec, printMat)
from os import path

 
class Line():
    '''A class for any mooring line that consists of a single material'''

    def __init__(self, mooringSys, num, L, lineType, nSegs=100, cb=0, isRod=0, attachments = [0,0]):
        '''Initialize Line attributes

        Parameters
        ----------
        mooringSys : system object
            The system object that contains the point object
        num : int
            indentifier number
        L : float
            line unstretched length [m]
        lineType : dict
            dictionary containing the coefficients needed to describe the line (could reference an entry of System.lineTypes).
        nSegs : int, optional
            number of segments to split the line into. Used in MoorPy just for plotting. The default is 100.
        cb : float, optional
            line seabed friction coefficient (will be set negative if line is fully suspended). The default is 0.
        isRod : boolean, optional
            determines whether the line is a rod or not. The default is 0.
        attachments : TYPE, optional
            ID numbers of any Points attached to the Line. The default is [0,0]. << consider removing

        Returns
        -------
        None.

        '''
        
        self.sys    = mooringSys       # store a reference to the overall mooring system (instance of System class)
        
        self.number = num
        self.isRod = isRod
            
        self.L = L  # line unstretched length
        self.type = lineType    # dictionary of a System.lineTypes entry
        
        self.nNodes = int(nSegs) + 1
        self.cb = float(cb)    # friction coefficient (will automatically be set negative if line is fully suspended)
        self.sbnorm = []    # Seabed Normal Vector (to be filled with a 3x1 normal vector describing seabed orientation)
        
        self.rA = np.zeros(3) # end coordinates
        self.rB = np.zeros(3)
        self.fA = np.zeros(3) # end forces
        self.fB = np.zeros(3)
        
        #Perhaps this could be made less intrusive by defining it using a line.addpoint() method instead, similar to point.attachline().
        self.attached = attachments  # ID numbers of the Points at the Line ends [a,b] >>> NOTE: not fully supported <<<<
        self.th = 0           # heading of line from end A to B
        self.HF = 0           # fairlead horizontal force saved for next solve
        self.VF = 0           # fairlead vertical force saved for next solve
        self.KA = []          # to be filled with the 2x2 end stiffness matrix from catenary
        self.KB = []          # to be filled with the 2x2 end stiffness matrix from catenary
        self.info = {}        # to hold all info provided by catenary
        
        self.qs = 1  # flag indicating quasi-static analysis (1). Set to 0 for time series data
        self.show = True      # a flag that will be set to false if we don't want to show the line (e.g. if results missing)
        #print("Created Line "+str(self.number))
        self.color = 'k'
        self.lw=0.5
        
        

    
    def loadData(self, dirname, rootname, sep='.MD.'):
        '''Loads line-specific time series data from a MoorDyn output file'''
        
        self.qs = 0 # signals time series data
        
        if self.isRod==1:
            strtype='Rod'
        elif self.isRod==0:
            strtype='Line'

        filename = dirname+rootname+sep+strtype+str(self.number)+'.out'
        
        if path.exists(filename):


        # try:
        
            # load time series data
            data, ch, channels, units = read_mooring_file("", filename) # remember number starts on 1 rather than 0

            # get time info
            if ("Time" in ch):
                self.Tdata = data[:,ch["Time"]]
                self.dt = self.Tdata[1]-self.Tdata[0]
            else:
                raise LineError("loadData: could not find Time channel for mooring line "+str(self.number))
        
            
            nT = len(self.Tdata)  # number of time steps
            
            # check for position data <<<<<<
            
            self.xp = np.zeros([nT,self.nNodes])
            self.yp = np.zeros([nT,self.nNodes])
            self.zp = np.zeros([nT,self.nNodes])
            
            
            for i in range(self.nNodes):
                self.xp[:,i] = data[:, ch['Node'+str(i)+'px']]
                self.yp[:,i] = data[:, ch['Node'+str(i)+'py']]
                self.zp[:,i] = data[:, ch['Node'+str(i)+'pz']]
            
            '''
            if self.isRod==0:
                self.Te = np.zeros([nT,self.nNodes-1])   # read in tension data if available
                if "Seg1Te" in ch:
                    for i in range(self.nNodes-1):
                        self.Te[:,i] = data[:, ch['Seg'+str(i+1)+'Te']]
                        
                self.Ku = np.zeros([nT,self.nNodes])   # read in curvature data if available
                if "Node0Ku" in ch:
                    for i in range(self.nNodes):
                        self.Ku[:,i] = data[:, ch['Node'+str(i)+'Ku']]
            else:
                # read in Rod buoyancy force data if available
                if "Node0Box" in ch:
                    self.Bx = np.zeros([nT,self.nNodes])   
                    self.By = np.zeros([nT,self.nNodes])
                    self.Bz = np.zeros([nT,self.nNodes])
                    for i in range(self.nNodes):
                        self.Bx[:,i] = data[:, ch['Node'+str(i)+'Box']]
                        self.By[:,i] = data[:, ch['Node'+str(i)+'Boy']]
                        self.Bz[:,i] = data[:, ch['Node'+str(i)+'Boz']]

            if "Node0Ux" in ch:
                self.Ux = np.zeros([nT,self.nNodes])   # read in fluid velocity data if available
                self.Uy = np.zeros([nT,self.nNodes])
                self.Uz = np.zeros([nT,self.nNodes])
                for i in range(self.nNodes):
                    self.Ux[:,i] = data[:, ch['Node'+str(i)+'Ux']]
                    self.Uy[:,i] = data[:, ch['Node'+str(i)+'Uy']]
                    self.Uz[:,i] = data[:, ch['Node'+str(i)+'Uz']]
            
            #Read in tension data if available
            if "Seg1Ten" in ch:
                self.Ten = np.zeros([nT,self.nNodes-1])   
                for i in range(self.nNodes-1):
                    self.Ten[:,i] = data[:, ch['Seg'+str(i+1)+'Ten']]
            '''
            
            
            
            # --- Read in additional data if available ---

            # segment tension  <<< to be changed to nodal tensions in future MD versions
            if "Seg1Ten" in ch:
                self.Tendata = True
                self.Te = np.zeros([nT,self.nNodes-1])
                for i in range(self.nNodes-1):
                    self.Te[:,i] = data[:, ch['Seg'+str(i+1)+'Ten']]
            elif "Seg1Te" in ch:
                self.Tendata = True
                self.Te = np.zeros([nT,self.nNodes-1])
                for i in range(self.nNodes-1):
                    self.Te[:,i] = data[:, ch['Seg'+str(i+1)+'Te']]
            else:
                self.Tendata = False
                        
            # curvature at node
            if "Node0Ku" in ch:
                self.Kudata = True
                self.Ku = np.zeros([nT,self.nNodes])   
                for i in range(self.nNodes):
                    self.Ku[:,i] = data[:, ch['Node'+str(i)+'Ku']]
            else:
                self.Kudata = False
            
            # water velocity data 
            if "Node0Ux" in ch:  
                self.Udata = True
                self.Ux = np.zeros([nT,self.nNodes])
                self.Uy = np.zeros([nT,self.nNodes])
                self.Uz = np.zeros([nT,self.nNodes])
                for i in range(self.nNodes):
                    self.Ux[:,i] = data[:, ch['Node'+str(i)+'Ux']]
                    self.Uy[:,i] = data[:, ch['Node'+str(i)+'Uy']]
                    self.Uz[:,i] = data[:, ch['Node'+str(i)+'Uz']]
            else:
                self.Udata = False
                
            # buoyancy force data
            if "Node0Box" in ch:  
                self.Bdata = True
                self.Bx = np.zeros([nT,self.nNodes])
                self.By = np.zeros([nT,self.nNodes])
                self.Bz = np.zeros([nT,self.nNodes])
                for i in range(self.nNodes):
                    self.Bx[:,i] = data[:, ch['Node'+str(i)+'Box']]
                    self.By[:,i] = data[:, ch['Node'+str(i)+'Boy']]
                    self.Bz[:,i] = data[:, ch['Node'+str(i)+'Boz']]
            else:
                self.Bdata = False
                
            # hydro drag data
            if "Node0Dx" in ch: 
                self.Ddata = True
                self.Dx = np.zeros([nT,self.nNodes])   # read in fluid velocity data if available
                self.Dy = np.zeros([nT,self.nNodes])
                self.Dz = np.zeros([nT,self.nNodes])
                for i in range(self.nNodes):
                    self.Dx[:,i] = data[:, ch['Node'+str(i)+'Dx']]
                    self.Dy[:,i] = data[:, ch['Node'+str(i)+'Dy']]
                    self.Dz[:,i] = data[:, ch['Node'+str(i)+'Dz']]
            else:
                self.Ddata = False
                
            # weight data
            if "Node0Wx" in ch: 
                self.Wdata = True
                self.Wx = np.zeros([nT,self.nNodes])   # read in fluid velocity data if available
                self.Wy = np.zeros([nT,self.nNodes])
                self.Wz = np.zeros([nT,self.nNodes])
                for i in range(self.nNodes):
                    self.Wx[:,i] = data[:, ch['Node'+str(i)+'Wx']]
                    self.Wy[:,i] = data[:, ch['Node'+str(i)+'Wy']]
                    self.Wz[:,i] = data[:, ch['Node'+str(i)+'Wz']]
            else:
                self.Wdata = False
            
            
            
            # initialize positions (is this used?)
            self.xpi= self.xp[0,:]
            self.ypi= self.yp[0,:]
            self.zpi= self.zp[0,:]
            
            # calculate the dynamic LBot !!!!!!! doesn't work for sloped bathymetry yet !!!!!!!!!!
            for i in range(len(self.zp[0])):
                if np.max(self.zp[:,i]) > self.zp[0,0]:
                    inode = i
                    break
                else:
                    inode = i
            self.LBotDyn = (inode-1)*self.L/(self.nNodes-1)
            
            # get length (constant)
            #self.L = np.sqrt( (self.xpi[-1]-self.xpi[0])**2 + (self.ypi[-1]-self.ypi[0])**2 + (self.zpi[-1]-self.zpi[0])**2 )
            # ^^^^^^^ why are we changing the self.L value to not the unstretched length specified in MoorDyn?
            # moved this below the dynamic LBot calculation because I wanted to use the original self.L
            # >>> this is probably needed for Rods - should look into using for Rods only <<<
            
            # check for tension data <<<<<<<
            
            self.show = True
            
        else:
            self.Tdata = []
            self.show = False
            print(f"Error geting data for {'Rod' if self.isRod else 'Line'} {self.number}: {filename}")
            print("dirname: {} or rootname: {} is incorrect".format(dirname, rootname))
            
         
        # >>> this was another option for handling issues - maybe no longer needed <<<
        #except Exception as e:
        #    # don't fail if there's an issue finding data, just flag that the line shouldn't be shown/plotted
        #    print(f"Error geting data for {'Rod' if self.isRod else 'Line'} {self.number}: ")
        #    print(e)
        #    self.show = False
        
        

    def getTimestep(self, Time):
        '''Get the time step to use for showing time series data'''
        
        if Time < 0: 
            ts = np.int_(-Time)  # negative value indicates passing a time step index
        else:           # otherwise it's a time in s, so find closest time step
            if len(self.Tdata) > 0:
                for index, item in enumerate(self.Tdata):                
                    ts = -1
                    if item > Time:
                        ts = index
                        break
                if ts==-1:
                    raise LineError(self.number, "getTimestep: requested time likely out of range")
            else:
                raise LineError(self.number, "getTimestep: zero time steps are stored")

        return ts
        
        

    def getLineCoords(self, Time, n=0):    # formerly UpdateLine
        '''Gets the updated line coordinates for drawing and plotting purposes.'''
        
        if n==0: n = self.nNodes
    
        # special temporary case to draw a rod for visualization. This assumes the rod end points have already been set somehow
        if self.qs==1 and self.isRod > 0:
        
            # make points for appropriately sized cylinder
            d = self.type['d_vol']
            Xs, Ys, Zs = makeTower(self.L, np.array([d/2, d/2]))   # add in makeTower method once you start using Rods
            
            # get unit vector and orientation matrix
            k = (self.rB-self.rA)/self.L
            Rmat = np.array(rotationMatrix(0, np.arctan2(np.hypot(k[0],k[1]), k[2]), np.arctan2(k[1],k[0])))
        
            # translate and rotate into proper position for Rod
            coords = np.vstack([Xs, Ys, Zs])
            newcoords = np.matmul(Rmat,coords)
            Xs = newcoords[0,:] + self.rA[0]
            Ys = newcoords[1,:] + self.rA[1]
            Zs = newcoords[2,:] + self.rA[2]
            
            return Xs, Ys, Zs, None
        
    
        # if a quasi-static analysis, just call the catenary function to return the line coordinates
        elif self.qs==1:
            
            self.staticSolve(profiles=1) # call with flag to tell Catenary to return node info
            
            #Xs = self.rA[0] + self.info["X"]*self.cosBeta 
            #Ys = self.rA[1] + self.info["X"]*self.sinBeta 
            #Zs = self.rA[2] + self.info["Z"]
            #Ts = self.info["Te"]
            Xs = self.Xs
            Ys = self.Ys
            Zs = self.Zs
            Ts = self.Ts
            return Xs, Ys, Zs, Ts
            
        # otherwise, count on read-in time-series data
        else:

            # figure out what time step to use
            ts = self.getTimestep(Time)
            
            # drawing rods
            if self.isRod > 0:
            
                k1 = np.array([ self.xp[ts,-1]-self.xp[ts,0], self.yp[ts,-1]-self.yp[ts,0], self.zp[ts,-1]-self.zp[ts,0] ]) / self.L # unit vector
                
                k = np.array(k1) # make copy
            
                Rmat = np.array(rotationMatrix(0, np.arctan2(np.hypot(k[0],k[1]), k[2]), np.arctan2(k[1],k[0])))  # <<< should fix this up at some point, MattLib func may be wrong
                
                # make points for appropriately sized cylinder
                d = self.type['d_vol']
                Xs, Ys, Zs = makeTower(self.L, np.array([d/2, d/2]))   # add in makeTower method once you start using Rods
                
                # translate and rotate into proper position for Rod
                coords = np.vstack([Xs, Ys, Zs])
                newcoords = np.matmul(Rmat,coords)
                Xs = newcoords[0,:] + self.xp[ts,0]
                Ys = newcoords[1,:] + self.yp[ts,0]
                Zs = newcoords[2,:] + self.zp[ts,0]
                
                return Xs, Ys, Zs, None
                
            # drawing lines
            else:
                
                # handle whether or not there is tension data
                try:  # use average to go from segment tension to node tensions <<< can skip this once MD is updated to output node tensions
                    Te = 0.5*(np.append(self.Te[ts,0], self.Te[ts,:]) +np.append(self.Te[ts,:], self.Te[ts,-1]))
                except: # otherwise return zeros to avoid an error (might want a warning in some cases?)
                    Te = np.zeros(self.nNodes)
                
                return self.xp[ts,:], self.yp[ts,:], self.zp[ts,:], Te
    
    
    def getCoordinate(self, s, n=100):
        '''Returns position and tension at a specific point along the line's unstretched length'''
        
        dr =  self.rB - self.rA                 
        LH = np.hypot(dr[0], dr[1])  
            
        Ss = np.linspace(0, self.L, n)
        Xs, Ys, Zs, Ts = self.getLineCoords(0.0, n=n)
        
        X = np.interp(s, Ss, Xs)*dr[0]/LH
        Y = np.interp(s, Ss, Ys)*dr[1]/LH
        Z = np.interp(s, Ss, Zs)
        T = np.interp(s, Ss, Ts)
        
        return X, Y, Z, T
        
    
    
    def drawLine2d(self, Time, ax, color="k", Xuvec=[1,0,0], Yuvec=[0,0,1], Xoff=0, Yoff=0, colortension=False, cmap='rainbow', plotnodes=[], plotnodesline=[], label="", alpha=1.0):
        '''Draw the line on 2D plot (ax must be 2D)

        Parameters
        ----------
        Time : float
            time value at which to draw the line
        ax : axis
            the axis on which the line is to be drawn
        color : string, optional
            color identifier in one letter (k=black, b=blue,...). The default is "k".
        Xuvec : list, optional
            plane at which the x-axis is desired. The default is [1,0,0].
        Yuvec : list, optional
            plane at which the y-axis is desired. The default is [0,0,1].
        colortension : bool, optional
            toggle to plot the lines in a colormap based on node tensions. The default is False
        cmap : string, optional
            colormap string type to plot tensions when colortension=True. The default is 'rainbow'

        Returns
        -------
        linebit : list
            list of axes and points on which the line can be plotted

        '''
        
        linebit = []  # make empty list to hold plotted lines, however many there are
        
        if self.isRod > 0:
            
            Xs, Ys, Zs, Te = self.getLineCoords(Time)
        
            # apply any 3D to 2D transformation here to provide desired viewing angle
            Xs2d = Xs*Xuvec[0] + Ys*Xuvec[1] + Zs*Xuvec[2] 
            Ys2d = Xs*Yuvec[0] + Ys*Yuvec[1] + Zs*Yuvec[2] 
        
            for i in range(int(len(Xs)/2-1)):
                linebit.append(ax.plot(Xs2d[2*i:2*i+2]    ,Ys2d[2*i:2*i+2]    , lw=0.5, color=color))  # side edges
                linebit.append(ax.plot(Xs2d[[2*i,2*i+2]]  ,Ys2d[[2*i,2*i+2]]  , lw=0.5, color=color))  # end A edges
                linebit.append(ax.plot(Xs2d[[2*i+1,2*i+3]],Ys2d[[2*i+1,2*i+3]], lw=0.5, color=color))  # end B edges
        
        # drawing lines...
        else:            
            # >>> can probably streamline the next bit of code a fair bit <<<
            if self.qs==1:
                Xs, Ys, Zs, tensions = self.getLineCoords(Time)
            elif self.qs==0:
                Xs, Ys, Zs, Ts = self.getLineCoords(Time)
                self.rA = np.array([Xs[0], Ys[0], Zs[0]])
                self.rB = np.array([Xs[-1], Ys[-1], Zs[-1]])
                tensions = self.getLineTens()   
            
            # apply any 3D to 2D transformation here to provide desired viewing angle
            Xs2d = Xs*Xuvec[0] + Ys*Xuvec[1] + Zs*Xuvec[2] + Xoff
            Ys2d = Xs*Yuvec[0] + Ys*Yuvec[1] + Zs*Yuvec[2] + Yoff
            
            if colortension:    # if the mooring lines want to be plotted with colors based on node tensions
                maxt = np.max(tensions); mint = np.min(tensions)
                for i in range(len(Xs)-1):          # for each node in the line
                    color_ratio = ((tensions[i] + tensions[i+1])/2 - mint)/(maxt - mint)  # ratio of the node tension in relation to the max and min tension
                    cmap_obj = cm.get_cmap(cmap)    # create a cmap object based on the desired colormap
                    rgba = cmap_obj(color_ratio)    # return the rbga values of the colormap of where the node tension is
                    linebit.append(ax.plot(Xs2d[i:i+2], Ys2d[i:i+2], color=rgba))
            else:
                linebit.append(ax.plot(Xs2d, Ys2d, lw=1, color=color, label=label, alpha=alpha)) # previously had lw=1 (linewidth)
            
            if len(plotnodes) > 0:
                for i,node in enumerate(plotnodes):
                    if self.number==plotnodesline[i]:
                        linebit.append(ax.plot(Xs2d[node], Ys2d[node], 'o', color=color, markersize=5))   
            
        self.linebit = linebit # can we store this internally?
        
        self.X = np.array([Xs, Ys, Zs])
            
        return linebit

    

    def drawLine(self, Time, ax, color="k", endpoints=False, shadow=True, colortension=False, cmap_tension='rainbow'):
        '''Draw the line in 3D
        
        Parameters
        ----------
        Time : float
            time value at which to draw the line
        ax : axis
            the axis on which the line is to be drawn
        color : string, optional
            color identifier in one letter (k=black, b=blue,...). The default is "k".
        endpoints : bool, optional
            toggle to plot the end points of the lines. The default is False
        shadow : bool, optional
            toggle to plot the mooring line shadow on the seabed. The default is True
        colortension : bool, optional
            toggle to plot the lines in a colormap based on node tensions. The default is False
        cmap : string, optional
            colormap string type to plot tensions when colortension=True. The default is 'rainbow'
            
        Returns
        -------
        linebit : list
            list of axes and points on which the line can be plotted
        '''
        
        if not self.show:  # exit if this line isn't set to be shown
            return 0
        
        if color == 'self':
            color = self.color  # attempt to allow custom colors
            lw = self.lw
        else:
            lw = 1
        
        linebit = []  # make empty list to hold plotted lines, however many there are
    
        if self.isRod > 0:
            
            if color==None:
                color = [0.3, 0.3, 0.3]  # if no color provided, default to dark grey rather than rainbow rods
                
            Xs, Ys, Zs, Ts = self.getLineCoords(Time)
            
            for i in range(int(len(Xs)/2-1)):
                linebit.append(ax.plot(Xs[2*i:2*i+2],Ys[2*i:2*i+2],Zs[2*i:2*i+2]            , color=color))  # side edges
                linebit.append(ax.plot(Xs[[2*i,2*i+2]],Ys[[2*i,2*i+2]],Zs[[2*i,2*i+2]]      , color=color))  # end A edges
                linebit.append(ax.plot(Xs[[2*i+1,2*i+3]],Ys[[2*i+1,2*i+3]],Zs[[2*i+1,2*i+3]], color=color))  # end B edges
            
            # scatter points for line ends 
            #if endpoints == True:
            #    linebit.append(ax.scatter([Xs[0], Xs[-1]], [Ys[0], Ys[-1]], [Zs[0], Zs[-1]], color = color))
        
        # drawing lines...
        else:
            # >>> can probably streamline the next bit of code a fair bit <<<
            if self.qs==1:  # returns the node positions and tensions of the line, doesn't matter what time
                Xs, Ys, Zs, tensions = self.getLineCoords(Time)
            elif self.qs==0: # returns the node positions and time data at the given time
                Xs, Ys, Zs, Ts = self.getLineCoords(Time)
                self.rA = np.array([Xs[0], Ys[0], Zs[0]])
                self.rB = np.array([Xs[-1], Ys[-1], Zs[-1]])
                tensions = self.getLineTens()
            
            if colortension:    # if the mooring lines want to be plotted with colors based on node tensions
                maxt = np.max(tensions); mint = np.min(tensions)
                for i in range(len(Xs)-1):          # for each node in the line
                    color_ratio = ((tensions[i] + tensions[i+1])/2 - mint)/(maxt - mint)  # ratio of the node tension in relation to the max and min tension
                    cmap_obj = cm.get_cmap(cmap_tension)    # create a cmap object based on the desired colormap
                    rgba = cmap_obj(color_ratio)    # return the rbga values of the colormap of where the node tension is
                    linebit.append(ax.plot(Xs[i:i+2], Ys[i:i+2], Zs[i:i+2], color=rgba, zorder=100))
            else:
                linebit.append(ax.plot(Xs, Ys, Zs, color=color, lw=lw, zorder=100))
            
            if shadow:
                ax.plot(Xs, Ys, np.zeros_like(Xs)-self.sys.depth, color=[0.5, 0.5, 0.5, 0.2], lw=lw, zorder = 1.5) # draw shadow
            
            if endpoints == True:
                linebit.append(ax.scatter([Xs[0], Xs[-1]], [Ys[0], Ys[-1]], [Zs[0], Zs[-1]], color = color))
                
                    
            # draw additional data if available (should make this for rods too eventually - drawn along their axis nodes)
            if self.qs == 0:              
                ts = self.getTimestep(Time)
                
                if self.Tendata:
                    pass
                if self.Kudata:
                    pass        
                if self.Udata:
                    self.Ubits = ax.quiver(Xs, Ys, Zs, self.Ux[ts,:], self.Uy[ts,:], self.Uz[ts,:], color="blue")  # make quiver plot and save handle to line object
                if self.Bdata:
                    self.Bbits = ax.quiver(Xs, Ys, Zs, self.Bx[ts,:], self.By[ts,:], self.Bz[ts,:], color="red")
                if self.Ddata:
                    self.Dbits = ax.quiver(Xs, Ys, Zs, self.Dx[ts,:], self.Dy[ts,:], self.Dz[ts,:], color="green")
                if self.Wdata:
                    self.Wbits = ax.quiver(Xs, Ys, Zs, self.Wx[ts,:], self.Wy[ts,:], self.Wz[ts,:], color="orange")
                
                
        self.linebit = linebit # can we store this internally?
        
        self.X = np.array([Xs, Ys, Zs])
        
            
        return linebit
    
    

        
        
    def redrawLine(self, Time, colortension=False, cmap_tension='rainbow', drawU=True):  #, linebit):
        '''Update 3D line drawing based on instantaneous position'''
        
        linebit = self.linebit
        
        if self.isRod > 0:
            
            Xs, Ys, Zs, Ts = self.getLineCoords(Time)
            
            for i in range(int(len(Xs)/2-1)):
                        
                linebit[3*i  ][0].set_data(Xs[2*i:2*i+2],Ys[2*i:2*i+2])    # side edges (x and y coordinates)
                linebit[3*i  ][0].set_3d_properties(Zs[2*i:2*i+2])         #            (z coordinates)             
                linebit[3*i+1][0].set_data(Xs[[2*i,2*i+2]],Ys[[2*i,2*i+2]])           # end A edges
                linebit[3*i+1][0].set_3d_properties(Zs[[2*i,2*i+2]])                    
                linebit[3*i+2][0].set_data(Xs[[2*i+1,2*i+3]],Ys[[2*i+1,2*i+3]])   # end B edges
                linebit[3*i+2][0].set_3d_properties(Zs[[2*i+1,2*i+3]])
        
        # drawing lines...
        else:
        
            Xs, Ys, Zs, Ts = self.getLineCoords(Time)
            
            if colortension:
                self.rA = np.array([Xs[0], Ys[0], Zs[0]])       # update the line ends based on the MoorDyn data
                self.rB = np.array([Xs[-1], Ys[-1], Zs[-1]])
                tensions = self.getLineTens()                   # get the tensions of the line calculated quasi-statically
                maxt = np.max(tensions); mint = np.min(tensions)
                cmap_obj = cm.get_cmap(cmap_tension)               # create the colormap object
                
                for i in range(len(Xs)-1):  # for each node in the line, find the relative tension of the segment based on the max and min tensions
                    color_ratio = ((tensions[i] + tensions[i+1])/2 - mint)/(maxt - mint)
                    rgba = cmap_obj(color_ratio)
                    linebit[i][0]._color = rgba         # set the color of the segment to a new color based on its updated tension
                    linebit[i][0].set_data(Xs[i:i+2],Ys[i:i+2])     # set the x and y coordinates
                    linebit[i][0].set_3d_properties(Zs[i:i+2])      # set the z coorindates
            
            else:
                linebit[0][0].set_data(Xs,Ys)    # (x and y coordinates)
                linebit[0][0].set_3d_properties(Zs)         # (z coordinates) 
                    
            
        
            # draw additional data if available (should make this for rods too eventually - drawn along their axis nodes)
            if self.qs == 0:
                ts = self.getTimestep(Time)                    
                s = 0.0002
                
                if self.Tendata:
                    pass
                if self.Kudata:
                    pass        
                if self.Udata:
                    self.Ubits.set_segments(quiver_data_to_segments(Xs, Ys, Zs, self.Ux[ts,:], self.Uy[ts,:], self.Uz[ts,:], scale=10.))
                if self.Bdata:
                    self.Bbits.set_segments(quiver_data_to_segments(Xs, Ys, Zs, self.Bx[ts,:], self.By[ts,:], self.Bz[ts,:], scale=s))
                if self.Ddata:
                    self.Dbits.set_segments(quiver_data_to_segments(Xs, Ys, Zs, self.Dx[ts,:], self.Dy[ts,:], self.Dz[ts,:], scale=s))
                if self.Wdata:
                    self.Wbits.set_segments(quiver_data_to_segments(Xs, Ys, Zs, self.Wx[ts,:], self.Wy[ts,:], self.Wz[ts,:], scale=s))
                
                
        
        return linebit
        
        
    
    
    def setEndPosition(self, r, endB):
        '''Sets the end position of the line based on the input endB value.

        Parameters
        ----------
        r : array
            x,y,z coorindate position vector of the line end [m].
        endB : boolean
            An indicator of whether the r array is at the end or beginning of the line

        Raises
        ------
        LineError
            If the given endB value is not a 1 or 0

        Returns
        -------
        None.

        '''
        
        if endB == 1:
            self.rB = np.array(r, dtype=np.float_)
        elif endB == 0:
            self.rA = np.array(r, dtype=np.float_)
        else:
            raise LineError("setEndPosition: endB value has to be either 1 or 0")
        
        
    def staticSolve(self, reset=False, tol=0.0001, profiles=0):
        '''Solves static equilibrium of line. Sets the end forces of the line based on the end points' positions.

        Parameters
        ----------
        reset : boolean, optional
            Determines if the previous fairlead force values will be used for the catenary iteration. The default is False.

        tol : float
            Convergence tolerance for catenary solver measured as absolute error of x and z values in m.
            
        profiles : int
            Values greater than 0 signal for line profile data to be saved (used for plotting, getting distributed tensions, etc).

        Raises
        ------
        LineError
            If the horizontal force at the fairlead (HF) is less than 0

        Returns
        -------
        None.

        '''
        
        # ensure line profile information is computed if needed for computing current loads
        if self.sys.currentMod == 1 and profiles == 0:
            profiles = 1

        #depth = self.sys.depth
        # get seabed depth and slope under each line end
        depthA, nvecA = self.sys.getDepthFromBathymetry(self.rA[0], self.rA[1])
        depthB, nvecB = self.sys.getDepthFromBathymetry(self.rB[0], self.rB[1])
        depth = 0.5*(depthA + depthB)  # temporary approach <<<
        
        #breakpoint()
        
        dr =  self.rB - self.rA
        LH = np.hypot(dr[0], dr[1])     # horizontal spacing of line ends
        LV = dr[2]                # vertical offset from end A to end B
        if LH >0:
            cosBeta = dr[0]/LH                 # cos of line heading
            sinBeta = dr[1]/LH                 # sin of line heading
            self.th = np.arctan2(dr[1],dr[0])  # line heading
        else:   # special case of vertical line: line heading is undefined - use zero as default
            cosBeta = 0.0
            sinBeta = 0.0
            self.th = 0.0

        if self.rA[2] < -depthA:
            self.rA[2] = -depthA
            self.cb = 0
            #raise LineError("Line {} end A is lower than the seabed.".format(self.number)) <<< temporarily adjust to seabed depth
        elif self.rB[2] < -depthB:
            raise LineError("Line {} end B is lower than the seabed.".format(self.number))
        else:
            self.cb = -depthA - self.rA[2]  # when cb < 0, -cb is defined as height of end A off seabed (in catenary)
        '''
        elif np.min([self.rA[2],self.rB[2]]) > -depth:
            self.cb = -depth - np.min([self.rA[2],self.rB[2]])   # if this line's lower end is off the seabed, set cb negative and to the distance off the seabed
        elif self.cb < 0:   # if a line end is at the seabed, but the cb is still set negative to indicate off the seabed
            self.cb = 0.0     # set to zero so that the line includes seabed interaction.
        '''
        
        if self.HF < 0:  # or self.VF < 0:  <<<<<<<<<<< it shouldn't matter if VF is negative - this could happen for buoyant lines, etc.
            raise LineError("Line HF cannot be negative") # this could be a ValueError too...
            
        if reset==True:   # Indicates not to use previous fairlead force values to start catenary 
            self.HF = 0   # iteration with, and insteady use the default values.
        
        
        # ----- Do some preprocessing to determine the seabed slope -----

        # check if there's seabed contact AND seabed slope at end A
        #if self.rA[2] + depthA < 1 and nvecA[2] < 1.0:
        if self.rA[2] + depthA < 1:
            
            self.cb = 0 # hasty override for safety
            
            #Determine the heading of the line and construct the d_hat unit vector
            d_hat = [cosBeta, sinBeta,0]   
            
            #Determine the v vector which is the d vector projectd onto the seabed
            v_hat = [np.sqrt(1-((d_hat[0]*nvecA[0]+d_hat[1]*nvecA[1])/nvecA[2])**2)*d_hat[0], 
                     np.sqrt(1-((d_hat[0]*nvecA[0]+d_hat[1]*nvecA[1])/nvecA[2])**2)*d_hat[1], 
                     -(d_hat[0]*nvecA[0]+d_hat[1]*nvecA[1])/nvecA[2]]  
                     
            #Determine the seabed slope
            if v_hat[2] == 0:
                alpha = 0 
            else: 
                cosArg = np.dot(d_hat, v_hat)/(np.linalg.norm(d_hat)*np.linalg.norm(v_hat))
                if cosArg > 1:
                    cosArg = 1
                if cosArg < -1:
                    cosArg = -1
                alpha = np.sign(v_hat[2])*(180/np.pi)*np.arccos(cosArg) 
            #Throw an error if the seabed slope is more than the arctan(Zf/Xf) then line would be fully on the slope and or would need to go through the slope which cannot be handled by our equations
            if alpha > (180/np.pi)*np.arctan((self.rB[2]-self.rA[2])/LH):    
                raise LineError(self.number, "Fairlead/Anchor Position not compatible with Positive Seabed Slope")    
        else:
            alpha = 0
            
       
        # ----- get line results for linear or nonlinear elasticity -----
        
        #If EA is found in the line properties we will run the original catenary function 
        if 'EA' in self.type:
            try:
                (fAH, fAV, fBH, fBV, info) = catenary(LH, LV, self.L, self.type['EA'], 
                                              self.type['w'], CB=self.cb, alpha=alpha, Tol=tol, 
                                              HF0=self.HF, VF0=self.VF, nNodes = self.nNodes, plots=profiles)   
                                                                                                                                    
            except CatenaryError as error:
                raise LineError(self.number, error.message)       
       #If EA isnt found then we will use the ten-str relationship defined in the input file 
        else:
             (fAH, fAV, fBH, fBV, info) = nonlinear(LH, LV, self.L, self.type['Str'], self.type['Ten'],self.type['w']) 
            # should have a profiles=1 option for this too (could be linear, or use rod style)
            
        #Postion of the line
        if profiles > 0:
            Xs = self.rA[0]+cosBeta*info["X"]
            Ys = self.rA[1]+sinBeta*info["X"]
            Zs = self.rA[2]+info["Z"]
            Ts = info["Te"]
        
            
        self.HF = info["HF"]
        self.VF = info["VF"]
        self.KA2 = info["stiffnessA"]
        self.KB2 = info["stiffnessB"]
        self.LBot = info["LBot"]
        self.info = info
            
        self.fA[0] = fAH*cosBeta
        self.fA[1] = fAH*sinBeta
        self.fA[2] = fAV
        self.fB[0] = fBH*cosBeta
        self.fB[1] = fBH*sinBeta
        self.fB[2] = fBV
        self.TA = np.sqrt(fAH*fAH + fAV*fAV) # end tensions
        self.TB = np.sqrt(fBH*fBH + fBV*fBV)
        
        # a lazy save for now
        self.sinBeta = sinBeta
        self.cosBeta = cosBeta
        
        
        # ----- compute 3d stiffness matrix for both line ends (3 DOF + 3 DOF) -----
        
        # may want to skip this next bit when just getting profiles for plotting...
        
        # create the rotation matrix based on the heading angle that the line is from the horizontal
        R = rotationMatrix(0,0,self.th)

        # line end stiffness matrices
        self.KA  = from2Dto3Drotated(self.info['stiffnessA'], -fBH, LH, R)  # reaction at A due to motion of A
        self.KB  = from2Dto3Drotated(self.info['stiffnessB'], -fBH, LH, R)  # reaction at B due to motion of B
        self.KAB = from2Dto3Drotated(self.info['stiffnessAB'], fBH, LH, R)  # reaction at B due to motion of A
        # <<< move the above into an else block for currentMod below
        
        # ----- apply current loads if applicable -----
        
        if self.sys.currentMod == 1: 
            #Calculate the current loading on the line
            
            # -------------- calculate new current drag and resultant W ----------

            U = self.sys.current
            
            FCurrent = np.zeros(3)  # total current force on line in x, y, z [N]        
            
            # Loop through each segment along the line and add up the drag forces.
            # This is in contrast to MoorDyn calculating for nodes.
            for i in range(self.nNodes-1):
                #For each segment find the tangent vector and then calculate the current loading
                dr_seg = np.array([Xs[i+1] - Xs[i], Ys[i+1] - Ys[i], Zs[i+1] - Zs[i]])  # segment vector
                ds_seg = np.linalg.norm(dr_seg)
                
                if ds_seg > 0:                   # only include if segment length > 0
                    q = dr_seg/ds_seg
                    # transverse and axial current velocity components
                    Uq = np.dot(U, q) * q
                    Up = U - Uq          
                    # transverse and axial drag forces on segment
                    dp = 0.5*self.sys.rho*self.type["Cd"]        *self.type["d_vol"]*ds_seg*np.linalg.norm(Up)*Up
                    dq = 0.5*self.sys.rho*self.type["CdAx"]*np.pi*self.type["d_vol"]*ds_seg*np.linalg.norm(Uq)*Uq
                    # add to total current force on line
                    FCurrent += dp + dq    
            
            
            #Define Equivalent Weight Vector
            w_eqv = FCurrent/self.L + np.array([0, 0, -self.type["w"]])
            w_eqv_norm = w_eqv/np.linalg.norm(w_eqv)
            
            
            # ------------ Perform rotation/transformation -------------
            
            #Be careful if vectors are the same.  Bad things will happen
            if np.array_equal(w_eqv_norm,np.array([0, 0, -1])):
                R_curr = np.eye(3,3)
            else:
                R_curr = RotFrm2Vect(w_eqv_norm, np.array([0, 0, -1]))  # rotation vector to make w_eqv vertical
            
            #rotate the seabed normal and the line headings
            rot_sbnorm = np.matmul(R_curr, nvecA)
            rot_dr = np.matmul(R_curr, dr)  # vector from A to B in rotated frame
            
            # figure out what rotation about Z' is needed to align end A onto the X'-Z' plane
            theta_z = -np.arctan2(rot_dr[1], rot_dr[0])
            
            # add this rotation to the previous one to get the complete rotation needed
            R_z = rotationMatrix(0, 0, theta_z)
            R_current = np.matmul(R_z, R_curr)
            
            # figure out slope in plane (only if contacting the seabed)
            if self.rA[2] <= -depthA or self.rB[2] <= -depthB:
                rot_sbnorm2 = np.matmul(R_current, nvecA)
            
                dz_dx = -rot_sbnorm2[0]*(1.0/rot_sbnorm2[2])  # seabed slope components
                dz_dy = -rot_sbnorm2[1]*(1.0/rot_sbnorm2[2])  # seabed slope components
                # we only care about dz_dx since the line is in the X-Z plane in this rotated situation
                alpha = np.degrees(np.arctan(dz_dx))
                cb = self.cb
            else:
                alpha = 0
                cb = min(0, rot_dr[2]) - 100  # put the seabed out of reach
            
            #Reassign the values for LH and LV 
            LH = np.linalg.norm(rot_dr[:2])
            LV = rot_dr[2]
            
            
            # ----- get line results for linear or nonlinear elasticity -----
            
            #NOTE FOR  CAT SOLVE THAT WE ARE GOING TO USE W_EQV INSTEAD OF W
            #If EA is found in the line properties we will run the original catenary function 
            if 'EA' in self.type:
                try:
                    (fAH, fAV, fBH, fBV, info) = catenary(LH, LV, self.L, self.type['EA'], 
                                                          np.linalg.norm(w_eqv),
                                                          CB=cb, alpha=alpha, HF0=self.HF, VF0=self.VF, 
                                                          nNodes=self.nNodes, plots=profiles)                                                                      
                except CatenaryError as error:
                    raise LineError(self.number, error.message)       
            #If EA isnt found then we will use the ten-str relationship defined in the input file 
            else:
                (fAH, fAV, fBH, fBV, info) = nonlinear(LH, LV, self.L, self.type['Str'], self.type['Ten'],np.linalg.norm(w_eqv)) 
        

        
            #Unrotate the system (Anti-rotation ;p)
            
            if profiles > 0:
                for i in range(0,self.nNodes):
                    temp_array = np.array([info['X'][i], 0 ,info['Z'][i]])
                    unrot_pos = np.matmul(temp_array, R_current)
                    #overwrite line positions
                    Xs[i] = unrot_pos[0]
                    Ys[i] = unrot_pos[1]
                    Zs[i] = unrot_pos[2]
                
                if self.number==3 and abs(self.rB[2] - (-1.55813956e+02)) < 0.001:
                    breakpoint()
                
                # global line coordinates
                Xs = self.rA[0] + Xs 
                Ys = self.rA[1] + Ys
                Zs = self.rA[2] + Zs
                Ts = info["Te"]
            
            # line end stiffness matrices transformed into global orientation
            self.KA  = from2Dto3Drotated(info['stiffnessA'], -fBH, LH, R_current)  # reaction at A due to motion of A
            self.KB  = from2Dto3Drotated(info['stiffnessB'], -fBH, LH, R_current)  # reaction at B due to motion of B
            self.KAB = from2Dto3Drotated(info['stiffnessAB'], fBH, LH, R_current)  # reaction at B due to motion of A

            #Reassign our unrotated forces
            self.fA = np.matmul(np.array([fAH, 0, fAV]), R_current)
            self.fB = np.matmul(np.array([fBH, 0, fBV]), R_current)
            self.TA = np.linalg.norm(self.fA) # end tensions
            self.TB = np.linalg.norm(self.fB)
        
        
        if profiles > 0:
            self.Xs = Xs
            self.Ys = Ys
            self.Zs = Zs
            self.Ts = Ts
        
        if profiles > 1:
            import matplotlib.pyplot as plt
            plt.plot(self.info['X'], self.info['Z'])
            plt.show()
        
    
        
    def getEndForce(self, endB):
        '''Returns the force of the line at the specified end based on the endB value

        Parameters
        ----------
        endB : boolean
            An indicator of which end of the line is the force wanted

        Raises
        ------
        LineError
            If the given endB value is not a 1 or 0

        Returns
        -------
        fA or fB: array
            The force vector at the end of the line

        '''
        
        if endB == 1:
            return self.fB
        elif endB == 0:
            return self.fA
        else:
            raise LineError("getEndForce: endB value has to be either 1 or 0")
            
            
            
    def getStiffnessMatrix(self):
        '''Returns the stiffness matrix of a line derived from analytic terms in the jacobian of catenary

        Raises
        ------
        LineError
            If a singluar matrix error occurs while taking the inverse of the Line's Jacobian matrix.

        Returns
        -------
        K2_rot : matrix
            the analytic stiffness matrix of the line in the rotated frame.

        '''

        # take the inverse of the Jacobian to get the starting analytic stiffness matrix
        '''
        if np.isnan(self.jacobian[0,0]): #if self.LBot >= self.L and self.HF==0. and self.VF==0.  << handle tricky cases here?
            K = np.array([[0., 0.], [0., 1.0/self.jacobian[1,1] ]])
        else:
            try:
                K = np.linalg.inv(self.jacobian)
            except:
                raise LineError(self.number, f"Check Line Length ({self.L}), it might be too long, or check catenary ProfileType")
        '''
        
        # solve for required variables to set up the perpendicular stiffness. Keep it horizontal
        L_xy = np.linalg.norm(self.rB[:2] - self.rA[:2])
        T_xy = np.linalg.norm(self.fB[:2])
        Kt = T_xy/L_xy
        
        # initialize the line's analytic stiffness matrix in the "in-line" plane
        KA = np.array([[self.KA2[0,0], 0 , self.KA2[0,1]],
                       [     0      , Kt,      0      ],
                       [self.KA2[1,0], 0 , self.KA2[1,1]]])
                       
        KB = np.array([[self.KB2[0,0], 0 , self.KB2[0,1]],
                       [     0      , Kt,      0      ],
                       [self.KB2[1,0], 0 , self.KB2[1,1]]])
        
        # create the rotation matrix based on the heading angle that the line is from the horizontal
        R = rotationMatrix(0,0,self.th)
        
        # rotate the matrix to be about the global frame [K'] = [R][K][R]^T
        KA_rot = np.matmul(np.matmul(R, KA), R.T)
        KB_rot = np.matmul(np.matmul(R, KB), R.T)
        
        return KA_rot, KB_rot

    
    def getLineTens(self):
        '''Calls the catenary function to return the tensions of the Line for a quasi-static analysis'''

        # >>> this can probably be done using data already generated by static Solve <<<
        '''
        depth = self.sys.depth
    
        dr =  self.rB - self.rA                 
        LH = np.hypot(dr[0], dr[1])     # horizontal spacing of line ends
        LV = dr[2]                      # vertical offset from end A to end B
        
        if np.min([self.rA[2],self.rB[2]]) > -depth:
            self.cb = -depth - np.min([self.rA[2],self.rB[2]])   # if this line's lower end is off the seabed, set cb negative and to the distance off the seabed
        elif self.cb < 0:   # if a line end is at the seabed, but the cb is still set negative to indicate off the seabed
            self.cb = 0.0     # set to zero so that the line includes seabed interaction.
    
        tol = 0.0001
        # --------------------
        
        #Determine the heading of the line and construct the d_hat unit vector
        x_excursion = self.rB[0]-self.rA[0]
        y_excursion = self.rB[1]-self.rA[1]
        total_xy_dist = np.sqrt( x_excursion* x_excursion + y_excursion*y_excursion)
        d_hat = [x_excursion/total_xy_dist,y_excursion/total_xy_dist,0]   
        #Determine the v vector which is the d vector projectd onto the seabed
        v_hat = [np.sqrt(1-((d_hat[0]*self.sbnorm[0]+d_hat[1]*self.sbnorm[1])/self.sbnorm[2])**2)*d_hat[0], np.sqrt(1-((d_hat[0]*self.sbnorm[0]+d_hat[1]*self.sbnorm[1])/self.sbnorm[2])**2)*d_hat[1], -(d_hat[0]*self.sbnorm[0]+d_hat[1]*self.sbnorm[1])/self.sbnorm[2]]    
        #Determine the seabed slope
        if v_hat[2] == 0:
            alpha = 0 
        else: 
            cosArg = np.dot(d_hat, v_hat)/(np.linalg.norm(d_hat)*np.linalg.norm(v_hat))
            if cosArg > 1:
                cosArg = 1
            if cosArg < -1:
                cosArg = -1
            alpha = np.sign(v_hat[2])*(180/np.pi)*np.arccos(cosArg) 
        #Throw an error if the seabed slope is more than the arctan(Zf/Xf) then line would be fully on the slope and or would need to go through the slope which cannot be handled by our equations
        if alpha > (180/np.pi)*np.arctan((self.rB[2]-self.rA[2])/total_xy_dist):    
            raise LineError(22,"Fairlead/Anchor Position not compatible with Positive Seabed Slope")    
        # --------------------                                                                                                  
    
        #If EA is found in the line properties we will run the original catenary function 
        if 'EA' in self.type:
            try:
                tol = 0.000001 #TODO figure out why tol and profiles are not defined. These values are hardcoded from defaults in other function calls
                profiles = 1
                (fAH, fAV, fBH, fBV, info) = catenary(LH, LV, self.L, self.type['EA'], self.type['w'], CB=self.cb, Tol=tol, HF0=self.HF, VF0=self.VF, plots=profiles)   # call line model
                                                                                                                  
            except CatenaryError as error:
                raise LineError(self.number, error.message) 
        #If EA isnt found then we will use the ten-str relationship defined in the input file 
        else:
             (fAH, fAV, fBH, fBV, info) = nonlinear(LH, LV, self.L, self.type['Str'], self.type['Ten'],self.type['w']) 
        '''

        self.staticSolve(profiles=1) # call with flag to tell Catenary to return node info (may be unnecessary)

        Ts = self.info["Te"]
        return Ts
    

    def getTension(self, s):
        '''Returns tension at a given point along the line
        
        Parameters
        ----------
        
        s : scalar or array-like
            Value or array of values for the arc length along the line from end A to end B at which
            the information is desired. Positive values are arc length in m, negative values are a
            relative location where 0 is end A, -1 is end B, and -0.5 is the midpoint.
        
        Returns
        -------
        
        tension value(s)
        
        '''
        #if s < 0:
        #    s = -s*self.L            
        #if s > self.L:
        #    raise ValueError('Specified arc length is larger than the line unstretched length.')
        
        Te = np.interp(s, self.info['s'], self.info['Te'])
        
        return Te


    def getPosition(self, s):
        '''Returns position at a given point along the line
        
        Parameters
        ----------
        
        s : scalar or array-like
            Value or array of values for the arc length along the line from end A to end B at which
            the information is desired. Positive values are arc length in m, negative values are a
            relative location where 0 is end A, -1 is end B, and -0.5 is the midpoint.
        
        Returns
        -------
        
        position vector(s)
        
        '''
        
        # >>> should be merged with getLineCoords and getCoordinate functionality <<<
        
        x = np.interp(s, self.info['s'], self.info['X'])
        z = np.interp(s, self.info['s'], self.info['Z'])
        
        
        dr =  self.rB - self.rA                 
        LH = np.hypot(dr[0], dr[1])
        Xs = self.rA[0] + x*dr[0]/LH
        Ys = self.rA[1] + x*dr[1]/LH
        Zs = self.rA[2] + z
        
        return np.vstack([ Xs, Ys, Zs])
    
    def attachLine(self, lineID, endB):
        pass


def from2Dto3Drotated(K2D, F, L, R): 
    '''Initialize a line end's analytic stiffness matrix in the 
    plane of the catenary then rotate the matrix to be about the 
    global frame using [K'] = [R][K][R]^T
    
    Parameters
    ----------
    K2D : 2x2 matrix
        Planar stiffness matrix of line end [N/m]
    F : float
        Line horizontal tension component [N]
    L : float
        Line horizontal distance end-to-end [m]
    R : 3x3 matrix
        Rotation matrix from global frame to plane to the local
        X-Z plane of the line
        
    Returns
    -------
    3x3 stiffness matrix in global orientation [N/m].
    '''

    if L > 0:
        Kt = F/L         # transverse stiffness term
    else:
        Kt = 0.0
    
    K2 = np.array([[K2D[0,0], 0 , K2D[0,1]],
                   [  0     , Kt,   0     ],
                   [K2D[1,0], 0 , K2D[1,1]]])
    
    return np.matmul(np.matmul(R, K2), R.T)    
    


def RotFrm2Vect( A, B):
    '''Rodriguez rotation function, which returns the rotation matrix 
    that transforms vector A into Vector B.
    '''
    
    v = np.cross(A,B)
    ssc = np.array([[0, -v[2], v[1]],
                [v[2], 0, -v[0]],
                [-v[1], v[0], 0]])
         
    R =  np.eye(3,3) + ssc + np.matmul(ssc,ssc)*(1-np.dot(A,B))/(np.linalg.norm(v)*np.linalg.norm(v))            

    return R