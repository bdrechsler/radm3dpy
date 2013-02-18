try:
    import numpy as np
except:
    print 'ERROR'
    print ' Numpy cannot be imported '
    print ' To use the python module of RADMC-3D you need to install Numpy'


from radmc3d.natconst import *
try:
    from matplotlib.pylab import *
except:
    print ' WARNING'
    print ' matploblib.pylab cannot be imported ' 
    print ' To used the visualization functionality of the python module of RADMC-3D you need to install matplotlib'
    print ' Without matplotlib you can use the python module to set up a model but you will not be able to plot things or'
    print ' display images'

from radmc3d.crd_trans import vrot
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

    defpar = [['rin', '1.0*au', ' Inner radius of the disk'],
    ['rdisk', '200.0*au', ' Outer radius of the disk'],
    ['hrdisk', '0.1', ' Ratio of the pressure scale height over radius at hrpivot'],
    ['hrpivot', "200.0*au", ' Reference radius at which Hp/R is taken'],
    ['plh', '1./7.', ' Flaring index'],
    ['plsig1', '-1.0', ' Power exponent of the surface density distribution as a function of radius'],
    ['sig0', '0.0', ' Surface density at rdisk'],
    ['mdisk', '1e-4*ms', ' Mass of the disk (either sig0 or mdisk should be set to zero or commented out)'],
    ['bgdens', '1e-30', ' Background density (g/cm^3)'],
    ['gap_rin', '[0e0*au]', ' Inner radius of the gap'],
    ['gap_rout', '[0e0*au]', ' Outer radius of the gap'],
    ['gap_drfact', '[0e0]', ' Density reduction factor in the gap'],
    ['dusttogas', '0.01', ' Dust-to-gas mass ratio']]


    return defpar

def get_density(rcyl=None, phi=None, z=None, z0=None, hp=None, sigma=None, grid=None, ppar=None):
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

        nr   = rcyl.shape[1]
        nphi = phi.shape[0]
        nz   = z.shape[0]

        # Calculate the pressure scale height as a function of r, phi
        if (hp==None):
            hp  = np.zeros([nr,nphi], dtype=np.float64)
            dum = ppar['hrdisk'] * (rcyl[:,0]/ppar['rdisk'])**ppar['plh'] * rcyl[:,0]
            for ip in range(nphi):
                hp[:,ip] = dum
   
        # Calculate the surface density if it is not given (e.g from a hydrodynamic simulation)
        if (sigma==None):
            sigma  = np.zeros([nr,nphi], dtype=np.float64)
            dum    = ppar['sig0'] * (rcyl/ppar['rdisk'])**ppar['plsig1']
            for ip in range(nphi):
                ii = (rcyl[:,ip]>=ppar['rin']) & (rcyl[:,ip]<=ppar['rdisk'])
                sigma[ii,ip] = dum[ii,ip]
            
    else:
        rr, th = np.meshgrid(grid.x, grid.y)
        zz   = rr * np.cos(th)
        rcyl = rr * np.sin(th)

        # Calculate the pressure scale height as a function of r, phi
        if (hp==None):
            hp = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64)
            dum = ppar['hrdisk'] * (rcyl/ppar['rdisk'])**ppar['plh'] * rcyl
            dum = dum.swapaxes(0,1)
            for iz in range(grid.nz):
                hp[:,:,iz] = dum

        # Calculate the surface density if it is not given (e.g from a hydrodynamic simulation)
        if (sigma==None):
            sigma = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64)
            # Calculate sigma from sig0, rdisk and plsig1
            if ppar.has_key('sig0'):
                if ppar['sig0'] != 0.:
                    dum    = ppar['sig0'] * (rcyl/ppar['rdisk'])**ppar['plsig1']
                    dum = dum.swapaxes(0,1)
                    for iz in range(grid.nz):
                        sigma[:,:,iz] = dum
                else:
                    dum    = 1.0 * (rcyl/ppar['rdisk'])**ppar['plsig1']
                    dum = dum.swapaxes(0,1)
                    for iz in range(grid.nz):
                        sigma[:,:,iz] = dum

                for iy in range(grid.ny):
                    ii = (rcyl[iy,:]<ppar['rin'])|(rcyl[iy,:]>ppar['rdisk'])
                    sigma[ii,iy,:] = 0.0

    if (z0==None):
        z0 = np.zeros([grid.nx, grid.nz, grid.ny], dtype=np.float64)


    if ppar.has_key('do_incl_innerdisk'):
        if ppar['do_incl_innerdisk']:
            for ix in range(grid.nx):                
                if (grid.x[ix]<ppar['innerdisk_rout']):
                    for iy in range(grid.ny):
                        z0[ix,:,iy] = grid.x[ix] * np.sin(ppar['innerdisk_incl']) * np.cos(grid.z)


    rho = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64)
    for iz in range(grid.nz):
        for iy in range(grid.ny):
            rho[:,iy,iz] = sigma[:,iy,iz] / (hp[:,iy,iz] * np.sqrt(2.0*np.pi)) * \
                np.exp(-0.5 * ((zz[iy,:])-z0[:,iz,iy])*((zz[iy,:])-z0[:,iz,iy])/\
                (hp[:,iy,iz]*hp[:,iy,iz])) + ppar['bgdens']

    # Normalize the disk to mdisk if it is given instead of sig0
    if ppar.has_key('mdisk'):
        if ppar['mdisk']!=0.:
            # Calculate the volume of each grid cell
            vol  = grid.get_cell_volume()
            mass = (rho*vol).sum(0).sum(0).sum(0)
            rho  = rho * (ppar['mdisk']/mass)

    if (ppar['gap_rout']>ppar['rin']):
        for igap in range(len(ppar['gap_rout'])):
            for ix in range(grid.nx):
                if (grid.x[ix]>=ppar['gap_rin'][igap])&(grid.x[ix]<=ppar['gap_rout'][igap]):
                    rho[ix,:,:] = rho[ix,:,:] * ppar['gap_drfact'][igap]
#                for iy in range(grid.ny):
#                    if ((rcyl[iy,ix]>=ppar['gap_rin'][igap])&(rcyl[iy,ix]<=ppar['gap_rout'][igap])):
#                        rho[ix,iy,:] = rho[ix,iz,:] * ppar['gap_drfact'][igap] 


    
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
