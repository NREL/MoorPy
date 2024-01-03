import numpy as np
from moorpy.dynamic_tension_functions import get_dynamic_tension, get_modes

class CompositeLine():
    def __init__(self, sys, point_id, rho=1025.):

        self.rho = rho
        point = sys.pointList[point_id-1] # Starting point id

        # check starting point to make sure it is either coupled or fixed and that it is not connected to more than 1 line.
        if point.type == 0:
            raise ValueError(f'Starting point ({point.number}) must not be free (change point type to -1 or 1).')
        elif len(point.attached)>1:
            raise ValueError(f'Starting point ({point.number}) cannot be connected to more than 1 line')

        self.pointList = []
        self.lineList = []

        self.pointA = point

        # Move through the points along the composite
        while True:
            # make sure that the point is free
            if len(point.attached)>2:
                raise ValueError(f'Point {point.number} is attached to more than two lines.')
            
            # get the line id and line object
            line_id = point.attached[-1] # get next line's id
            line = sys.lineList[line_id - 1] # get line object

            # append the starting point and the line object
            self.pointList.append(point)
            self.lineList.append(line)

            # get the next point
            attached_points = line.attached.copy() # get the IDs of the points attached to the line
            pointA_id = point.number # get first point ID
            attached_points.remove(pointA_id) # remove first point from attached point list
            pointB_id = attached_points[0] # get second point id
            point = sys.pointList[pointB_id-1] # get second point object

            # break from the loop when a point with a single attachment is reached
            if len(point.attached) == 1:
                self.pointList.append(point)
                break
        # make sure that the last point is not a free point
        if point.type == 0:
            raise ValueError(f'Last point ({point.number}) is a free point.') # make sure that the last point is not free
        
        self.pointB = point
        self.nNodes = (np.sum([line.nNodes for line in self.lineList]) - (len(self.lineList)-1))
        self.sys = sys

    def getDynamicMatrices(self, omegas, S_zeta,r_dynamic,depth,kbot,cbot,seabed_tol=1e-4):
        n_nodes = self.nNodes # number of nodes
        EA_segs = np.zeros(n_nodes-1) # extensional stiffness of the segments
        n_dofs = 3*n_nodes # number of dofs
        M = np.zeros([n_dofs,n_dofs], dtype='float')
        A = np.zeros([n_dofs,n_dofs], dtype='float')
        B = np.zeros([n_dofs,n_dofs], dtype='float')
        K = np.zeros([n_dofs,n_dofs], dtype='float')
        n = 0
        r_mean = np.zeros([n_nodes,3], dtype='float')
        r_dynamic = np.ones((len(omegas),n_nodes,3),dtype='float')*r_dynamic
        v_dynamic = 1j*omegas[:,None,None]*r_dynamic

        for line in self.lineList:
            line_type = line.type
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
    
    def dynamicSolve(self,omegas,S_zeta,RAO_A,RAO_B,depth,kbot,cbot,seabed_tol=1e-4,tol = 0.01,iters=100, w = 0.8):
        T_nodes_psd,T_nodes_std,s,r_static,r_dynamic,r_total,X = get_dynamic_tension(self,omegas,S_zeta,RAO_A,RAO_B,depth,kbot,cbot,
                                                                                     seabed_tol=seabed_tol,tol = tol,iters=iters, w = w)
        return T_nodes_psd,T_nodes_std,s,r_static,r_dynamic,r_total,X

    def getModes(self,fix_A=True,fix_B=True,plot_modes=False,amp_factor=1,adj_view = False,kbot=3E+06,cbot=3E+05,seabed_tol=1E-04):
        
        if plot_modes:
            freqs,mode_shapes,r_nodes,M,A,K,fig,ax = get_modes(self,fix_A=fix_A,fix_B=fix_B,plot_modes=plot_modes,amp_factor=amp_factor,adj_view = adj_view,
                                                               kbot=kbot,cbot=cbot,seabed_tol=seabed_tol)
            return freqs,mode_shapes,r_nodes,M,A,K,fig,ax
        else:
            freqs,mode_shapes,r_nodes,M,A,K = get_modes(self,fix_A=fix_A,fix_B=fix_B,plot_modes=plot_modes,amp_factor=amp_factor,adj_view = adj_view,
                                                        kbot=kbot,cbot=cbot,seabed_tol=seabed_tol)
            return freqs,mode_shapes,r_nodes,M,A,K