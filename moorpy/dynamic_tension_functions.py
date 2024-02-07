import numpy as np
from numpy.linalg import solve,norm
import scipy.linalg as la
from datetime import datetime
import matplotlib.pyplot as plt
from collections.abc import Iterable


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
    L_e1 = la.norm(dr_e1) # element 1 length
    t_e1 = (dr_e1)/L_e1 # tangential unit vector
    p_e1 = np.cross(t_e1,h_op) # in plane normal unit vector


    ut_e1 = np.einsum('ij,j->i',v_dynamic[:,0,:],t_e1) # tangential velocity
    uh_e1 = np.einsum('ij,j->i',v_dynamic[:,0,:],h_op) # normal horizontal out of plane velocity
    up_e1 = np.einsum('ij,j->i',v_dynamic[:,0,:],p_e1) # normal in plane velocity

    sigma_ut_e1 = np.sqrt(np.trapz(np.abs(ut_e1)**2*S_zeta,omegas))
    sigma_uh_e1 = np.sqrt(np.trapz(np.abs(uh_e1)**2*S_zeta,omegas))
    sigma_up_e1 = np.sqrt(np.trapz(np.abs(up_e1)**2*S_zeta,omegas))

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
        L_bw = la.norm(dr_bw) # element 1 length
        t_bw = (dr_bw)/L_bw # tangential unit vector
        p_bw = np.cross(t_bw,h_op) # in plane normal unit vector

        ut_bw = np.einsum('ij,j->i',v_dynamic[:,n,:],t_bw) # tangential velocity
        uh_bw = np.einsum('ij,j->i',v_dynamic[:,n,:],h_op) # normal horizontal out of plane velocity
        up_bw = np.einsum('ij,j->i',v_dynamic[:,n,:],p_bw) # normal in plane velocity

        sigma_ut_bw = np.sqrt(np.trapz(np.abs(ut_bw)**2*S_zeta,omegas))
        sigma_uh_bw = np.sqrt(np.trapz(np.abs(uh_bw)**2*S_zeta,omegas))
        sigma_up_bw = np.sqrt(np.trapz(np.abs(up_bw)**2*S_zeta,omegas))

        tt_bw = np.outer(t_bw,t_bw) # local tangential to global components transformation matrix
        pp_bw = np.outer(p_bw,p_bw) # local normal inplane to global components transformation matrix

        M[3*n:3*n+3,3*n:3*n+3] += mden*L_bw/2*np.eye(3) # mass contribution from adjacent elements

        A_bw = 1025*np.pi/4*deq**2*L_bw/2*(Ca*(hh_op+pp_bw) + CaAx*tt_bw) # backward element added mass contribution

        B_bw = 0.5*1025*deq*L_bw/2*np.sqrt(8/np.pi)*(Cd*(sigma_uh_bw*hh_op + sigma_up_bw*pp_bw) +
                                                        CdAx*sigma_ut_bw*tt_bw) # backward element damping contribution 

        K_bw = EA/L_bw*tt_bw + (T_mean[n]/L_bw)*(hh_op+pp_bw) # backward element stiffness (axial + geometric)

        ## forward element (n+1/2) contributions
        dr_fw = r_mean[n+1] - r_mean[n]
        L_fw = la.norm(dr_fw) # element 1 length
        t_fw = (dr_fw)/L_fw # tangential unit vector
        p_fw = np.cross(t_fw,h_op) # in plane normal unit vector


        ut_fw = np.einsum('ij,j->i',v_dynamic[:,n,:],t_fw) # tangential velocity
        uh_fw = np.einsum('ij,j->i',v_dynamic[:,n,:],h_op) # normal horizontal out of plane velocity
        up_fw = np.einsum('ij,j->i',v_dynamic[:,n,:],p_fw) # normal in plane velocity

        sigma_ut_fw = np.sqrt(np.trapz(np.abs(ut_fw)**2*S_zeta,omegas))
        sigma_uh_fw = np.sqrt(np.trapz(np.abs(uh_fw)**2*S_zeta,omegas))
        sigma_up_fw = np.sqrt(np.trapz(np.abs(up_fw)**2*S_zeta,omegas))

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
    L_eN = la.norm(dr_eN) # element N length
    t_eN = (dr_eN)/L_eN # tangential unit vector
    p_eN = np.cross(t_eN,h_op) # in plane normal unit vector

    ut_eN = np.einsum('ij,j->i',v_dynamic[:,N-1,:],t_eN) # tangential velocity
    uh_eN = np.einsum('ij,j->i',v_dynamic[:,N-1,:],h_op) # normal horizontal out of plane velocity
    up_eN = np.einsum('ij,j->i',v_dynamic[:,N-1,:],p_eN) # normal in plane velocity

    sigma_ut_eN = np.sqrt(np.trapz(np.abs(ut_eN)**2*S_zeta,omegas))
    sigma_uh_eN = np.sqrt(np.trapz(np.abs(uh_eN)**2*S_zeta,omegas))
    sigma_up_eN = np.sqrt(np.trapz(np.abs(up_eN)**2*S_zeta,omegas))

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


def get_dynamic_tension(Line,omegas,S_zeta,RAO_A,RAO_B,depth,kbot,cbot,seabed_tol=1e-4,tol = 0.01,iters=100, w = 0.8, conv_time=True):
    """
    Evaluates dynamic tension at all the nodes for an instance of MoorPy's Line or CompositeLine classes.

    Parameters
    ----------
    Line : Line/CompositeLine
        An instance of MoorPy's Line or CompositeLine classes.
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

    Returns
    -------
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
    """
    N = Line.nNodes
    n_dofs = 3*N

    if np.all(RAO_A == 0):
        RAO_A = np.zeros_like(RAO_B)
    if np.all(RAO_B == 0):
        RAO_B = np.zeros_like(RAO_A)

    # intialize iteration matrices
    r_dynamic_init = np.ones((len(omegas),N,3))
    M,A,B,K,r_static,EA_segs = Line.getDynamicMatrices(omegas,S_zeta,r_dynamic_init,depth,kbot,cbot,seabed_tol=seabed_tol) # TODO: return EA_segs
    X = np.zeros((len(omegas),n_dofs),dtype = 'complex')
    r_dynamic = np.zeros(((len(omegas),int(n_dofs/3),3)),dtype = 'complex')
    S_Xd = np.zeros((len(omegas),n_dofs),dtype = 'float')
    sigma_Xd = np.zeros(n_dofs,dtype = 'float')
    sigma_Xd0 = np.zeros(n_dofs,dtype = 'float')
    X[:, :3] = RAO_A
    X[:,-3:] = RAO_B

    # solving dynamics
    start = datetime.now()
    for ni in range(iters):
        H = - omegas[:,None,None]**2*(M+A)[None,:,:]\
            + 1j*omegas[:,None,None]*B[None,:,:]\
            + K[None,:,:]\
        
        F_A = np.einsum('nij,njk->ni',-H[:,3:-3, :3],RAO_A[:,:,None])
        F_B = np.einsum('nij,njk->ni',-H[:,3:-3,-3:],RAO_B[:,:,None])
        F = F_A + F_B

        X[:,3:-3] = solve(H[:,3:-3,3:-3],F)
        S_Xd[:] = np.abs(1j*omegas[:,None]*X)**2*S_zeta[:,None]
        sigma_Xd[:] = np.sqrt(np.trapz(S_Xd,omegas,axis=0)) 
        r_dynamic[:] = X.reshape(X.shape[0],int(X.shape[1]/3),3)
        if (np.abs(sigma_Xd-sigma_Xd0) <= tol*np.abs(sigma_Xd0)).all():
            break
        else:
            sigma_Xd0[:] = w * sigma_Xd + (1.-w) * sigma_Xd0
            _,_,B[:],_,_,_ = Line.getDynamicMatrices(omegas,S_zeta,r_dynamic,depth,kbot,cbot,seabed_tol=seabed_tol)
    if conv_time:
        print(f'Finished {ni} dynamic tension iterations in {datetime.now()-start} seconds (w = {w}).')

    # evaluate tension
    dw = np.diff(omegas,
            prepend= omegas[0] - (omegas[1]-omegas[0]),
            append= omegas[-1] + (omegas[-1]-omegas[-2]))
    dw = (dw[1:]+dw[:-1])/2
    wave_amps = np.sqrt(S_zeta*dw) #evaluate wave amplitudes of harmonic components from wave spectrum

    r_dynamic *= wave_amps[:,None,None]
    r_total = r_static[None,:,:] + r_dynamic
    dr_static = r_static[:-1] - r_static[1:]
    dr_dynamic = r_dynamic[:,:-1,:] - r_dynamic[:,1:,:]
    tangents = dr_static/la.norm(r_static[:-1] - r_static[1:], axis=-1)[:,None]
    L_static = la.norm(dr_static, axis=-1)
    dL_dynamic = np.einsum('mni,ni->mn', dr_dynamic, tangents)
    eps_segs = np.abs(dL_dynamic)/L_static

    T_segs = EA_segs * eps_segs
    T_nodes = np.zeros((len(omegas),N))
    T_nodes[:,0] = T_segs[:,0]
    T_nodes[:,1:-1] = (T_segs[:,:-1] + T_segs[:,1:])/2
    T_nodes[:,-1] = T_segs[:,-1]

    # S_T = np.zeros((len(omegas),N))
    # S_T[:,1:] = T_e**2/dw[:,None]
    # S_T[:,0] = S_T[:,1]

    T_nodes_psd = T_nodes**2/dw[:,None]
    T_nodes_std = np.sqrt(np.trapz(T_nodes_psd,omegas,axis=0))


    dr = np.diff(r_static,axis=0)
    ds = la.norm(dr,axis=1)
    s = np.zeros_like(T_nodes_std)
    s = np.cumsum(ds)

    return T_nodes_psd,T_nodes_std,s,r_static,r_dynamic,r_total,X


def get_modes(line,fix_A=True,fix_B=True,plot_modes=False,amp_factor=1,adj_view = False,kbot=3E+06,cbot=3E+05,seabed_tol=1E-04):
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
