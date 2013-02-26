try:
    import numpy as np
except:
    print 'ERROR'
    print ' Numpy cannot be imported '
    print ' To use the python module of RADMC-3D you need to install Numpy'


try:
    from matplotlib.pylab import *
except:
    print ' WARNING'
    print ' matploblib.pylab cannot be imported ' 
    print ' To used the visualization functionality of the python module of RADMC-3D you need to install matplotlib'
    print ' Without matplotlib you can use the python module to set up a model but you will not be able to plot things or'
    print ' display images'


try:
    ppdisk = __import__('model_ppdisk')
except:
    try:
        ppdisk  = __import__('radmc3d.model_ppdisk', fromlist=['']) 
    except:
        print 'ERROR'
        print ' model_ppdisk.py could not be imported (required by the selected warp_model "m97")'
        print ' The model files should either be in the current working directory or'
        print ' in the radmc3d python module directory'

from radmc3dPy.natconst import *
import sys


"""
PYTHON module for RADMC3D 
(c) Attila Juhasz 2011,2012,2013

Warped protoplanetary disk model of Mouillet et al. 1997

"""

def get_desc():
    """
    A one line description of the model
    """

    return 'Warped circumstellar disk model of Mouillet et al. 1997'

def get_default_params():
    """
    Function to provide default parameter values 

    OUTPUT:
    -------

    Returns a list whose elements are also lists with three elements:
    1) parameter name, 2) parameter value, 3) parameter description
    All three elements should be strings. The string of the parameter
    value will be directly written out to the parameter file if requested,
    and the value of the string expression will be evaluated and be put
    to radmc3dData.ppar. The third element contains the description of the
    parameter which will be written in the comment field of the line when
    a parameter file is written. 

    """

    defpar = ppdisk.get_default_params()

    defpar.append(['warp_model', "'m97'", ' Name of the warp model'])
    defpar.append(['warp_dcomp', '5.0*au', ' Distance of the companion, causing the perturbation, from the central star'])
    defpar.append(['warp_icomp', '10./180.*pi', ' Incliation of the companion with respecto the disk plane in radians'])
    defpar.append(['warp_mcomp', '1.0*mju', ' Mass of the perturber companion'])
    defpar.append(['warp_t', '1e5*year', ' Time at which the perturbation should be calculated'])

    return defpar

def get_warp_z0_m97(rcyl=None, phi=None, grid=None, mstar=None, dcomp=None, icomp=None, \
                        mcomp=None, t=None):
    """
     Function to create a twisted warp (Mouillet et al. 1997, MNRAS, 292, 896)
       This is a purely kinematic model to describe warps in the outer
       disk caused by a companion within the disk (R_DISK>D_COMP)
     USAGE
    
     Z0 = get_warp_z0_mouillet97(GRID=GRID, M_STAR=M_STAR, D_COMP=D_COMP, $
                                 I_COMP=I_COMP, M_COMP=M_COMP, T=T)
    
     OUTPUT
    
              Z0     : 3D numpy array, height of the midplane above the unperturbed value (z0)
    
     KEYWORDS 
    
              crd    : A numpy array the returning value of crd_trans.ctrans_sph2cyl containg
                       the cylindrical coordinates at each point
              MSTAR  : Mass of the central star [g] 
              DCOMP  : Distance of the companion from the central star [cm]
              ICOMP  : Inclination angle of the obit of the secondary
                       with respect to the plane of the disk [radian]
              MCOMP  : Mass of the secondary [g]
              T      : Time at which the twisted warp should be calculated 
                       NOTE : If T=0 no perturbation is present
    """
#
# Get the coordinates right;
#   Spherical coordinates in RADMC3D used [r, theta, phi] while the description of this 
#   model is in cylindrical coordinates [rcyl, phi, z]. So let's rename the coordinate
#   variables;
#   r     = grid.x     rcyl = r * sin(theta)
#   theta = grid.y  -> phi  = phi 
#   phi   = grid.z     z    = r * cos(theta)
#
    if (grid!=None):
        rr, th = np.meshgrid(grid.x, grid.y)
        rcyl   = rr * np.sin(th)
        z      = rr * np.cos(th)
        phi    = grid.z
        nr     = grid.nx
        nz     = grid.ny
        nphi   = grid.nz
    else:
        nr     = rcyl.shape[0]
        nphi   = phi.shape[0]
        nz     = 1 
        z      = array([0.0], dtype=float64)
#
# Now do the calculations 
#

    z0 = np.zeros([nr, nphi, nz], dtype=np.float64)
    omega      = sqrt(gg * mstar[0] / rcyl.T**3.)
    omega_prec = -3.0/4.0 * gg * mcomp * dcomp**2. / (omega * rcyl.T**5.)
    dum_z0     = -3./8.0 * gg * mcomp * dcomp**2. * sin(2.0*icomp) / \
                     (omega * omega_prec * rcyl.T**4.)

    for ip in range(nphi):
        z0[:,ip,:] = dum_z0 * (sin(phi[ip]) - sin(phi[ip] - omega_prec*t))
                
    return z0

# -----------------------------------------------------------------------------------------------
# 
# -----------------------------------------------------------------------------------------------
def get_warp_omega_prec(rcyl=None, phi=None, mstar=None, dcomp=None, icomp=None, \
                                mcomp=None, t=None):

    nc  = natconst()

    omega      = sqrt(gg * mstar / rcyl.T**3.)
    omega_prec = -3.0/4.0 * gg * mcomp * dcomp**2. / (omega * rcyl.T**5.)


#    nx  = rcyl.shape[0]
#    ny  = phi.shape[0]
#    
##
## Calculate the cylindrical coordinates
##
#    omega_prec = np.zeros([nx, ny], dtype=np.float64)
#    for ix in range(nx):
#        omega      = np.sqrt(nc.gg*mstar/rcyl[ix]**3.)
#        omega_prec[ix,:] = -3.0/4.0 * nc.gg * mcomp * dcomp**2. / (omega * rcyl[ix]**5.)

    return omega_prec

# -----------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------
def get_gas_density(grid=None, ppar=None):
    """
    Calculates dust density in g/cm^3
    """

#
# Apply perturbations 
#
    
    z0      = np.zeros([grid.nx, grid.nz, grid.ny], dtype=np.float64)


    if ppar.has_key('warp_model'):
        if ppar['warp_model'].lower()=='m97':
            z0 =get_warp_z0_m97(grid=grid, mstar=ppar['mstar'], \
                                    dcomp=ppar['warp_dcomp'],  icomp=ppar['warp_icomp'], mcomp=ppar['warp_mcomp'],\
                                    t=ppar['warp_t'])

    rr, th = np.meshgrid(grid.x, grid.y)
    rcyl = rr * np.sin(th)
    zz   = rr * np.cos(th)

    for ix in range(grid.nx):
        for iy in range(grid.ny):
            if (rcyl[iy,ix]<ppar['warp_dcomp']):
                z0[ix,:,iy] = 0.0


    rho = ppdisk.get_gas_density(z0=z0, grid=grid, ppar=ppar)


    return rho
# -----------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------
def get_dust_density(grid=None, ppar=None):
    """
    Calculates dust density in g/cm^3
    """

#
# Apply perturbations 
#
    
    z0      = np.zeros([grid.nx, grid.nz, grid.ny], dtype=np.float64)


    if ppar.has_key('warp_model'):
        if ppar['warp_model'].lower()=='m97':
            z0 =get_warp_z0_m97(grid=grid, mstar=ppar['mstar'], \
                                    dcomp=ppar['warp_dcomp'],  icomp=ppar['warp_icomp'], mcomp=ppar['warp_mcomp'],\
                                    t=ppar['warp_t'])

    rr, th = np.meshgrid(grid.x, grid.y)
    rcyl = rr * np.sin(th)
    zz   = rr * np.cos(th)

    for ix in range(grid.nx):
        for iy in range(grid.ny):
            if (rcyl[iy,ix]<ppar['warp_dcomp']):
                z0[ix,:,iy] = 0.0


    rho = ppdisk.get_dust_density(z0=z0, grid=grid, ppar=ppar)


    return rho


# -----------------------------------------------------------------------------------------------
# 
# -----------------------------------------------------------------------------------------------
def get_velocity(grid=None, ppar=None):

    

    # Calculate the twist angle
    # First figure out the tilt angle of that shell
    omp = -3.0/4.0 * gg * ppar['mcomp'] * ppar['dcomp']**2. / (sqrt(gg*ppar['mstar'][0]/grid.x**3) * grid.x**5.)
    alp = omp * ppar['t']

    warp_ang_twist = np.zeros([grid.nx], dtype=np.float64)
    warp_ang_tilt  = np.zeros([grid.nx], dtype=np.float64)
    for ix in range(grid.nx):
        if (alp[ix]==0.):
            warp_ang_twist[ix] = 0.
            warp_ang_tilt[ix] = 0.
        elif(alp[ix]==np.pi):
            z0  = 9./32. * grid.x[ix] * np.sin(2.*ppar['icomp']) * (1. - sin(np.pi/2. - alp[ix]))
            if (z0>0.):
                warp_ang_tilt[ix]  = arctan(z0/grid.x[ix])
            else:
                warp_ang_twist[ix] = 3.*np.pi/2.

            warp_ang_twist[ix] = np.pi/2.
        else:
            warp_ang_twist[ix] = arctan((1.0 - cos(alp[ix])) / sin(alp[ix])) 
            z0  = 9./32. * grid.x[ix] * np.sin(2.*ppar['icomp']) * (warp_ang_twist[ix] - sin(warp_ang_twist[ix] - alp[ix]))
            warp_ang_tilt[ix] = arctan(z0 / grid.x[ix])
        
        print ix, omp[ix], ppar['t'], alp[ix], warp_ang_twist[ix], warp_ang_tilt[ix]

    # Calculate the keplerian velocity field 
    rcyl = zeros([grid.nx, grid.ny, grid.nz], dtype=float64)
    for ix in range(grid.nx):
        for iz in range(grid.nz):
            rcyl[ix,:,iz] = abs(grid.x[ix]*sin(grid.y + warp_ang_tilt[ix] * cos(grid.z[iz])))

    fig = figure()
    plot(grid.x/au, warp_ang_twist/np.pi*180.)
    
    fig = figure()
    plot(grid.x/au, warp_ang_tilt/np.pi*180.)

    dddd = raw_input('')
    # Calculate the keplerian velocity 
    vkep = sqrt(gg*ppar['mstar'][0]/rcyl)

    # Now we know the amplitude of the velocities so I have to figure out what their direction
    # is. Let's take the unperturbet velocity field, where only the phi component of the velocity
    # vector is not zero and then rotate the vectors of each shell in the spherical mesh according
    # to the tilt and twist of that sphere.

    vel = zeros([grid.nx,grid.ny,grid.nz,3], dtype=float64)

    # For performance reasons let's create 3d arrays from the rotation angles and the coordinates
    # in the exact same dimensions as we wish to have the velocity field [ix, iy, iz]. If all variables
    # are in the form of arrays with the same dimensions we can use the performance advantage of numpy
    
    incl  = zeros([grid.nx, grid.ny, grid.nz], dtype=float64)
    twist = zeros([grid.nx, grid.ny, grid.nz], dtype=float64)
   
    r     = zeros([grid.nx, grid.ny, grid.nz], dtype=float64)
    phi   = zeros([grid.nx, grid.ny, grid.nz], dtype=float64)
    theta = zeros([grid.nx, grid.ny, grid.nz], dtype=float64)
   

    for ix in range(grid.nx):
        incl[ix,:,:]  = warp_ang_tilt[ix]
        twist[ix,:,:] = warp_ang_twist[ix] 
        r[ix,:,:]     = grid.x[ix]
        
    for iz in range(grid.nz):
        phi[:,:,iz]   = grid.z[iz]
    for iy in range(grid.ny):
        theta[:,iy,:] = grid.y[iy]

    # Here comes the rotation of the vector components itself 
    vphi = vkep * cos(incl * sin(phi-twist))
    vtheta = vkep * (sin(incl * sin(phi-twist)))

    # Fill up the returning array and we're done
    vel[:,:,:,1] = vtheta
    vel[:,:,:,2] = vphi


    return vel

