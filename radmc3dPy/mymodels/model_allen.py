try:
    import numpy as np
except:
    print 'ERROR'
    print ' Numpy cannot be imported '
    print ' To use the python module of RADMC-3D you need to install Numpy'


from radmc3dPy.natconst import *
try:
    from matplotlib.pylab import *
except:
    print ' WARNING'
    print ' matploblib.pylab cannot be imported ' 
    print ' To used the visualization functionality of the python module of RADMC-3D you need to install matplotlib'
    print ' Without matplotlib you can use the python module to set up a model but you will not be able to plot things or'
    print ' display images'

from radmc3dPy.crd_trans import vrot
import sys

"""
PYTHON module for RADMC3D 
(c) Attila Juhasz 2011/2012

Generic protoplanetary disk model

FUNCTIONS:
----------

    get_density()
    get_velocity()

"""
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

    defpar = {}

    defpar = ['r0', '10.0*au', 'Reference radius'], 
    ['n0', '1e-15', 'Volume density at the reference radius in the midplane']]

    return defpar

def get_density(rcyl=None, phi=None, z=None, grid=None, ppar=None):
    """
    Function to create the density distribution in a protoplanetary disk
    
    OUTPUT:
    -------
        returns the volume density in g/cm^3, whether the density is that of the gas
        or dust or both depends on what is specified in the surface density/mass
    """
    if (grid==None):
        print ' ***********************************************************'
        print 'This mode in model_ppdisk is still not yet finished'
        print 'Stopped'
        print ' ***********************************************************'
        sys.exit(0)
            
    else:

        rho = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64)

        for iy in range(ny):
            sint = sin(grid.y)
            for ix in range(nx):
                rho[ix,iy,:] = ppar['rho0'] * (1.8 * sint / (r[ir]/rc))**1.5


    
    # Split up the disk density distribution according to the given abundances
    if ppar.has_key('dustkappa_ext'):
        ngs  = len(ppar['dustkappa_ext'])
        if ppar.has_key('mfrac'):
            gsfact = ppar['mfrac'] / ppar['mfrac'].sum()
        else:
            ngs = 1
            
            
    else:
        ngs = ppar['ngs']
        #
        # WARNING!!!!!!
        # At the moment I assume that the multiple dust population differ from each other only in 
        # grain size but not in bulk density thus when I calculate the abundances / mass fractions 
        # they are independent of the grains bulk density since abundances/mass fractions are normalized
        # to the total mass. Thus I use 1g/cm^3 for all grain sizes.
        # TODO: Add the possibility to handle multiple dust species with different bulk densities and 
        # with multiple grain sizes.
        #
        gdens = zeros(ngs, dtype=float) + 1.0
        gs = ppar['gsmin'] * (ppar['gsmax']/ppar['gsmin']) ** (arange(ppar['ngs'], dtype=float64) / (float(ppar['ngs'])-1.))
        gmass = 4./3.*np.pi*gs**3. * gdens
        gsfact = gmass * gs**(ppar['gsdist_powex']+1)
        gsfact = gsfact / gsfact.sum()

    if (ngs>1):
       
        rho_old = array(rho)
        rho = np.zeros([grid.nx, grid.ny, grid.nz, ngs], dtype=np.float64)
        for igs in range(ngs):
            rho[:,:,:,igs] = rho_old[:,:,:] * gsfact[igs]

    return rho
# -----------------------------------------------------------------------------------------------
# 
# -----------------------------------------------------------------------------------------------
def get_velocity(rcyl=None, phi=None, z=None, z0=None, grid=None, ppar=None):
    """
    Function to create the velocity field in a protoplanetary disk
    """

   
    if (grid==None):
        nr   = rcyl.shape[1]
        nphi = phi.shpae[0]
        nz   = z.shape[0]

        vel = np.zeros([nr,nz,nphi, 3], dtype=np.float64)
        if (z0==None):
            vkep = np.sqrt(gg*ppar['mstar'][0]/rcyl)
            for iz in range(nz):
                for ip in range(nz):
                    vel[:,iz,ip, 2] = vkep
        else:
            rcyl_rot     = np.arange(nr, dtype=np.float64)
            for ir in range(nr):
                dum       = np.array(z0[ir,:])
                z0_max    = dum.max()
                rcyl_rot[ir]  = np.sqrt(rcyl[ir]**2. + z0_max**2.)
            
            vkep = np.sqrt(gg*ppar['mstar'][0]/rcyl_rot)
            for iz in range(nz):
                for ip in range(nz):
                    vel[:,iz,ip, 2] = vkep
    else:
        nr   = grid.nx
        nphi = grid.nz
        nz   = grid.ny
        rcyl = grid.x

        rr, th = np.meshgrid(grid.x, grid.y)
        rcyl_rot = rr * np.sin(th)
        
        vel = np.zeros([nr,nz,nphi,3], dtype=np.float64)
        vkep = np.sqrt(gg*ppar['mstar'][0]/rcyl)
        for iz in range(nz):
            for ip in range(nphi):
                vel[:,iz,ip,2] = vkep


    return vel
