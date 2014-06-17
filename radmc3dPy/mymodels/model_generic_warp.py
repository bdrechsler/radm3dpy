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
        ppdisk  = __import__('radmc3dPy.model_ppdisk', fromlist=['']) 
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

Generic warped protoplanetary disk model

"""
def getModelDesc():
    """
    A one line description of the model
    """

    return 'A generic warped disk model with analyitcal description for the warp'

def getDefaultParams():
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

    defpar = ppdisk.getDefaultParams()
    defpar.append(['warp_model', "'generic_warp'", ' Name of the warp model'])
    defpar.append(['warp_rin', '2.0*au', ' Inner radius of the warp'])
    defpar.append(['warp_powex', '-4.0', ' Power exponent of the inclination angle as a function of cylindrical radius'])
    defpar.append(['warp_iamp', '10./180.*pi', ' Maximum inclination angle of the warp with respect of the unperturbed disk midplane'])
    defpar.append(['warp_omprec0', '0.0', ' Precession angular velocity at the inner warp radius'])
    defpar.append(['warp_powexprec', '0.0', ' Power exponent of the precession angular velocity as a function of cylindrical radius'])
    defpar.append(['warp_t', '0.', ' Time at which the twist/precession term should be calculated'])

    return defpar




def getWarpAngles(grid=None, ppar=None, outer=False):
    """
    Function to calculate the vertical elevation of a disk warp (due to e.g. the presence of a companion
        in an inclinded orbit with respect to the unpertured disks' midplane). No precession is taken
        into account.

    The warp is described in a simple analytical formula:
        z0 = z0 * (R/Rin_warp)**(-p) * cos(phi - omega*t) 
        omega = omega0 * (R/Rin_warp)**(-q)

        using the inclination angle instead of the cylindrical z coorindate:
        i = arccos( cos(i_0) * (R/Rin_warp)**(p+1) * cos(phi - omega*t)


    INPUT:
    ------
        grid  : a radmc3dGrid class containing the spatial grid
        ppar  : a dictionary containing the parameters of the RADMC3D model
        outer : if set to True then the warp is calculated in the 'outer disk' (i.e. the outer
                parts of the disk will be the most distorted) oterwhise the warp is assumed to
                be in the inner disk edge.
    OUTPUT:
    -------
        This function returns a dictionary with the following keys:
        'z0'   : 2D Numpy array, vertical elevation of the disks midplane [cm]
        'incl' : 1D Numpy array, inclination angle of a disk 'ring' at a given radius [radian]
        'omt'  : 1D numpy array, twisting angle of the warp (omega_precession * t) [radian]
    """


#
# TODO: Check the gridstyle and if it's not spherical raise error/exception
# TODO: Check if all input parameters are set
# TODO: Check if we are in 3D (if nz>1) 
# 

    z0 = zeros([grid.nx, grid.nz], dtype=float64)
    if not outer:
        rcyl   = (grid.x / ppar['warp_rin'])
    else:
        rcyl = (grid.x / ppar['warp_rout'])

    omega = ppar['warp_omprec0'] * rcyl**(ppar['warp_powexprec'])
    incl = arctan(tan(ppar['warp_iamp']) * rcyl**(ppar['warp_powex']-1.)) 

    # Make sure that inside warp_rin no perturbation will be made
    incl[rcyl<1.] = 0.
    omega[rcyl<1.] = 0.

    for iz in range(grid.nz):
        z0[:,iz] = grid.x * tan(incl) * cos(grid.z[iz] - omega*ppar['warp_t'])

    fig = figure()
    ax  = fig.add_subplot(111, polar=True)
    c   = contourf(grid.z, grid.x/1.496e13, z0/1.496e13, 50)
    cb  = colorbar(c)
    xlabel('X [AU]')
    ylabel('Y [AU]')
    title('Z [AU')


    return {'z0':z0, 'incl':incl, 'omt':omega*ppar['warp_t']}




# ============================================================================================================================
#
# ============================================================================================================================
def getGasAbundance(grid=None, ppar=None, ispec=''):
    """
    Function to create the conversion factor from volume density to number density of molecule ispec.
    The number density of a molecule is rhogas * abun 
   
    INPUT:
    ------
        grid - An instance of the radmc3dGrid class containing the spatial and wavelength grid
        ppar - Dictionary containing all parameters of the model 
        ispec - The name of the gas species whose abundance should be calculated

    OUTPUT:
    -------
        returns the abundance as a Numpy array
    """

    # Read the dust density and temperature
    try: 
        data = readData(ddens=True, dtemp=True, binary=True)
    except:
        try: 
            data = readData(ddens=True, dtemp=True, binary=False)
        except:
            print 'WARNING!!'
            print 'No data could be read in binary or in formatted ascii format'
            print '  '
            return 0

    # Calculate continuum optical depth 
    data.getTau(axis='xy', wav=0.55)
    

    nspec = len(ppar['gasspec_mol_name'])
    ndust = data.dustemp.shape[3]

    if ppar['gasspec_mol_name'].__contains__(ispec):

        sid   = ppar['gasspec_mol_name'].index(ispec)
        # Check where the radial and vertical optical depth is below unity
        gasabun = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64)  
        
        for spec in range(nspec):
            gasabun[:,:,:] = ppar['gasspec_mol_abun'][sid]
           

        for iz in range(data.grid.nz):
            for iy in range(data.grid.ny):
                ii = (data.taux[:,iy,iz]<ppar['gasspec_mol_dissoc_taulim'][sid])
                gasabun[ii,iy,iz] = 1e-90

                ii = (data.dusttemp[:,iy,iz,ndust-1]<ppar['gasspec_mol_freezeout_temp'][sid])
                gasabun[ii,iy,iz] =  ppar['gasspec_mol_abun'][sid] * ppar['gasspec_mol_freezeout_dfact'][sid]
        
        #for iz in range(data.grid.nz):
            #for iy in range(data.grid.ny/2):

                #ii = (data.tauy[:,iy,iz]<ppar['gasspec_mol_dissoc_taulim'][sid])
                #gasabun[ii,iy,iz] = 1e-90
                #gasabun[ii,data.grid.ny-1-iy,iz] = 1e-90

    else:
        gasabun = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64) + 1e-10
        print 'WARNING !!!'
        print 'Molecule name "'+ispec+'" is not found in gasspec_mol_name'
        print 'A default 1e-10 abundance will be used'
        print ' ' 


    #gasabun = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64) 
    #gasabun[:,:,:] = ppar['gasspec_mol_abun'][0] / (2.4*mp)

    return gasabun
# -----------------------------------------------------------------------------------------------
# 
# -----------------------------------------------------------------------------------------------
def getGasDensity(grid=None, ppar=None):

#
# Apply perturbations 
#
 
    z0   = zeros([grid.nx, grid.nz, grid.ny], dtype=float64)
    dum  = getWarpAngles(grid=grid, ppar=ppar)
    for iy in range(grid.ny):
        z0[:,:,iy] = dum['z0']
#    z0   = dum['z0']

    rr, th = np.meshgrid(grid.x, grid.y)
    rcyl = rr * np.sin(th)
    zz   = rr * np.cos(th)

    for ix in range(grid.nx):
        for iy in range(grid.ny):
            if (rcyl[iy,ix]<ppar['warp_rin']):
                z0[ix,:,iy] = 0.0

    print dum['incl'].max()/pi*180.
    print z0.max() / au

    if ppar.has_key('warp_model'):
        if ppar['warp_model']=='generic_warp':
            rho = ppdisk.getGasDensity(z0=z0, grid=grid, ppar=ppar)
        else:
            rho = ppdisk.getGasDensity(grid=grid, ppar=ppar)
    else:
        rho = ppdisk.getDustDensity(grid=grid, ppar=ppar)



    return rho

# -----------------------------------------------------------------------------------------------
# 
# -----------------------------------------------------------------------------------------------
def getDustDensity(grid=None, ppar=None):

#
# Apply perturbations 
#
 
    z0   = zeros([grid.nx, grid.nz, grid.ny], dtype=float64)
    dum  = getWarpAngles(grid=grid, ppar=ppar)
    for iy in range(grid.ny):
        z0[:,:,iy] = dum['z0']
#    z0   = dum['z0']

    rr, th = np.meshgrid(grid.x, grid.y)
    rcyl = rr * np.sin(th)
    zz   = rr * np.cos(th)

    for ix in range(grid.nx):
        for iy in range(grid.ny):
            if (rcyl[iy,ix]<ppar['warp_rin']):
                z0[ix,:,iy] = 0.0

    print dum['incl'].max()/pi*180.
    print z0.max() / au

    if ppar.has_key('warp_model'):
        if ppar['warp_model']=='generic_warp':
            rho = ppdisk.getDustDensity(z0=z0, grid=grid, ppar=ppar)
        else:
            rho = ppdisk.getDustDensity(grid=grid, ppar=ppar)
    else:
        rho = ppdisk.get_density(grid=grid, ppar=ppar)



    return rho

# -----------------------------------------------------------------------------------------------
def getVTurb(grid=None, ppar=None):
    """
    Function to create the turbulent velocity field
    
    INPUT:
    ------
        grid - An instance of the radmc3dGrid class containing the spatial and wavelength grid
        ppar - Dictionary containing all parameters of the model 
    
    OUTPUT:
    -------
        returns the turbulent velocity in cm/s
    """

    vturb = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64) + ppar['gasspec_vturb']
    return vturb

# -----------------------------------------------------------------------------------------------
def getVelocity(grid=None, ppar=None):

    
    # Calculate the rotation angles of the precessing warp
    dum = getWarpAngles(grid=grid, ppar=ppar)
    rcyl = zeros([grid.nx, grid.ny, grid.nz], dtype=float64)
    for ix in range(grid.nx):
        for iz in range(grid.nz):
            rcyl[ix,:,iz] = abs(grid.x[ix]*sin(grid.y + dum['incl'][ix] * cos(grid.z[iz])))


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
        incl[ix,:,:]  = dum['incl'][ix]
        twist[ix,:,:] = dum['omt'][ix] 
        r[ix,:,:]     = grid.x[ix]
        
    for iz in range(grid.nz):
        phi[:,:,iz]   = grid.z[iz]
    for iy in range(grid.ny):
        theta[:,iy,:] = grid.y[iy]

    # Here comes the rotation of the vector components itself 
    vphi = vkep * cos(incl * sin(phi-twist))
    vtheta = vkep * (sin(incl * sin(phi-twist)))

    # Fill up the returning array and let' be DONE
    vel[:,:,:,1] = vtheta
    vel[:,:,:,2] = vphi


    return vel


