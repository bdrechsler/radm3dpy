
try:
    import scipy.interpolate as itp
except:
    print 'ERROR'
    print ' scipy.interpolate cannot be imported '
    print ' To use the python module of RADMC-3D you need to install Scipy'

try:
    from fargo import fargo_data
except:
    print 'ERROR'
    print ' fargo module cannot be imported'

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
import sys, os
from radmc3dPy.analyze import readData

"""
PYTHON module for RADMC3D 
(c) Attila Juhasz 2011,2012,2013

Protoplanetary disk model with surface density taken from the 2D hydrodynamic code FARGO

"""

def getModelDesc():
    """
    Function to provide a brief description of the model
    """

    return "Protoplanetary disk model. The surface density is taken from the output "+\
            "of the hydrodynamic code FARGO. Vertical structure is calculated analytically."

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


    defpar = {}

    defpar = [['fargo_path', "'/disk2/juhasz/Projects/Asymmetric_disks/FARGO/ECCEntric_disks/out4'", ' '],
            ['fargo_frame', '625', 'Frame number in the FARGO simulation to be used'],
            ['fargo_linscale', '5.7', 'Scale factor to change linear scales in the model'],
            ['fargo_umass', 'ms', 'Unit mass'],
            ['fargo_ulength', 'fargo_linscale*au', 'Unit length'],
            ['fargo_uvelo', 'sqrt(gg*mstar[0]/(fargo_linscale*au))', 'Unit velocity'],
            ['gfargo', 'False', 'Set it to True if the output is coming from the GPU version of FARGO'],
            ['rin', '1.0*au', ' Inner radius of the disk'],
            ['rdisk', '200.0*au', ' Outer radius of the disk'],
            ['hrdisk', '0.1', ' Ratio of the pressure scale height over radius at hrpivot'],
            ['hrpivot', "200.0*au", ' Reference radius at which Hp/R is taken'],
            ['plh', '1./7.', ' Flaring index'],
            ['bgdens', '1e-30', ' Background density (g/cm^3)'],
            ['dusttogas', '0.01', ' Dust-to-gas mass ratio'], 
            ['gasspec_mol_dissoc_taulim', '[1.0]', 'Continuum optical depth limit below which all molecules dissociate'],
            ['gasspec_mol_freezeout_temp', '[19.0]', 'Freeze-out temperature of the molecules in Kelvin'],
            ['gasspec_mol_freezeout_dfact', '[1e-3]', 'Factor by which the molecular abundance should be decreased in the frezze-out zone']]

    return defpar

# -----------------------------------------------------------------------------------------------
def getDustDensity(grid=None, ppar=None):
    """
    Function to generate a 3D volume density distribtion from the surface density output frames of FARGO
    
    INPUT:
    ------
        grid : a radmc3dGgrid class containing the spatial grid onto which the FARGO output should
                be interpolated
        ppar : dictionary containing all parameters of the RADMC3D setup

    OUTPUT:
    -------
        Returns the volume density of the gas in g/cm^3
    """


# Get the gas surface density
    rhogas = getGasDensity(grid=grid, ppar=ppar)
    rho    = array(rhogas) * ppar['dusttogas']

# Split up the dust density if we have a grain size distribution
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
    else:
        rho_old = array(rho)
        rho = np.zeros([grid.nx, grid.ny, grid.nz, ngs], dtype=np.float64)
        rho[:,:,:,0] = rho_old

    return rho

# -----------------------------------------------------------------------------------------------
def getDustTemperature(grid=None, ppar=None):
    """
    Function to get the dust temperature
    The dust temperature is taken to be that of the gas (i.e. FARGO-ADSG output) and the disk is assumed to
    be vertically isothermal. 
    
    INPUT:
    ------
        grid : a radmc3dGgrid class containing the spatial grid onto which the FARGO output should
                be interpolated
        ppar : dictionary containing all parameters of the RADMC3D setup

    OUTPUT:
    -------
        Returns the dust temperature in Kelvin
    """

#
# Read the gas temperature from a FARGO frame
# WARNING!!!
# fdata['gastemp'] has a dimension [nphi, nr] 
# For other models I was using [nr, nphi] for sigma!!!!
#
    fdata = fargo_data(path=ppar['fargo_path'], frame=ppar['fargo_frame'], gfargo=ppar['gfargo'],\
                           gastemp=True, vel=False)

#
# Now map the FARGO surface density onto our (r, phi) grid
#
#
# TODO: 
# The RectBivariateSpline() function should be doublechecked as once it produced a weird wiggle
#  in the interpolated field similar to numerical instability in the hydrodinamic advenction schemes..
# 
    gastemp_sp = itp.RectBivariateSpline(fdata['phi'], fdata['r'], fdata['gastemp'])
    
    iirf     = ((grid.x>=fdata['r'].min()*ppar['fargo_ulength'])&(grid.x<=fdata['r'].max()*ppar['fargo_ulength']))  
    dr       = grid.x[iirf]/ppar['fargo_ulength']
    fgastemp   = gastemp_sp(grid.z, dr) 

#
# Add an outer disk 
#
    #if ppar['add_outer_disk']:
        #iiro =(grid.x>=fdata['r'].max()*ppar['fargo_ulength'])
        #od_r = grid.x[iiro]
        #od_sigma = np.zeros([grid.nz, od_r.shape[0]])
        #sig0 = fsigma[:,dr.shape[0]-1].mean()
        #print fsigma[:,dr.shape[0]-1]
        #dum = sig0 * (od_r / fdata['r'].max()/ppar['fargo_ulength'])**ppar['od_plsig1']
        #for iz in range(grid.nz):
            #od_sigma[iz,:] = dum
#
# Blow up the 2D surface density to a 3D volume density distribution
#
          
    dusttemp2D  = np.zeros([grid.nx,grid.nz], dtype=np.float64)
    for iz in range(grid.nz):
        dusttemp2D[iirf,iz] = fgastemp[iz,:]
    #if ppar['add_outer_disk']:
        #for iz in range(grid.nz):
            #sigma[iiro,iz] = od_sigma[iz,:]


    fig = figure()
    x = fig.add_subplot(111)
    plot(grid.x, dusttemp2D[:,0])
    x.set_xscale('log')
    x.set_yscale('log')

    show()
    

    dusttemp = np.zeros([grid.nx, grid.ny, grid.nz, ppar['ngs']], dtype=np.float64)
    z0  = np.zeros([grid.nx, grid.nz], dtype=np.float64)
    for iy in range(grid.ny):
        for iz in range(grid.nz):
            dusttemp[:,iy,iz,0] = dusttemp2D[:,iz]

    for igs in range(1,ppar['ngs']):
        dusttemp[:,:,:,igs] = dusttemp[:,:,:,0]

        
    
    return dusttemp
# -----------------------------------------------------------------------------------------------
def getGasDensity(grid=None, ppar=None):
    """
    Function to generate a 3D volume density distribtion from the surface density output frames of FARGO
    
    INPUT:
    ------
        grid : a radmc3dGgrid class containing the spatial grid onto which the FARGO output should
                be interpolated
        ppar : dictionary containing all parameters of the RADMC3D setup

    OUTPUT:
    -------
        Returns the volume density of the gas in g/cm^3
    """

#
# Read the surface density from a FARGO frame
# WARNING!!!
# fdata['sigma'] has a dimension [nphi, nr] 
# For other models I was using [nr, nphi] for sigma!!!!
#
    fdata = fargo_data(path=ppar['fargo_path'], frame=ppar['fargo_frame'], gfargo=ppar['gfargo'],\
                           sigma=True, vel=False)

#
# Now map the FARGO surface density onto our (r, phi) grid
#
#
# TODO: 
# The RectBivariateSpline() function should be doublechecked as once it produced a weird wiggle
#  in the interpolated field similar to numerical instability in the hydrodinamic advenction schemes..
# 
    sigma_sp = itp.RectBivariateSpline(fdata['phi'], fdata['r'], fdata['sigma'])
    
    iirf     = ((grid.x>=fdata['r'].min()*ppar['fargo_ulength'])&(grid.x<=fdata['r'].max()*ppar['fargo_ulength']))  
    dr       = grid.x[iirf]/ppar['fargo_ulength']
    fsigma   = sigma_sp(grid.z, dr) * ppar['fargo_umass'] / ppar['fargo_ulength']**2.

#
# Add an outer disk 
#
    #if ppar['add_outer_disk']:
        #iiro =(grid.x>=fdata['r'].max()*ppar['fargo_ulength'])
        #od_r = grid.x[iiro]
        #od_sigma = np.zeros([grid.nz, od_r.shape[0]])
        #sig0 = fsigma[:,dr.shape[0]-1].mean()
        #print fsigma[:,dr.shape[0]-1]
        #dum = sig0 * (od_r / fdata['r'].max()/ppar['fargo_ulength'])**ppar['od_plsig1']
        #for iz in range(grid.nz):
            #od_sigma[iz,:] = dum
#
# Blow up the 2D surface density to a 3D volume density distribution
#
    rr, th = np.meshgrid(grid.x, grid.y)
    zz   = rr * np.cos(th)
    rcyl = rr * np.sin(th)
    hp  = np.zeros([grid.nx, grid.nz], dtype=np.float64)
    dum = ppar['hrdisk'] * (grid.x/ppar['hrpivot'])**ppar['plh'] * grid.x
    for ip in range(grid.nz):
        hp[:,ip] = dum
        
    sigma  = np.zeros([grid.nx,grid.nz], dtype=np.float64)
    for iz in range(grid.nz):
        sigma[iirf,iz] = fsigma[iz,:]
    #if ppar['add_outer_disk']:
        #for iz in range(grid.nz):
            #sigma[iiro,iz] = od_sigma[iz,:]


    #fig = figure()
    #x = fig.add_subplot(111)
    #plot(grid.x, sigma[:,0])
    #x.set_xscale('log')
    #x.set_yscale('log')

    #show()
    

    rho = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64)
    z0  = np.zeros([grid.nx, grid.nz], dtype=np.float64)
    for iy in range(grid.ny):
        for iz in range(grid.nz):
            rho[:,iy,iz] = sigma[:,iz] / (hp[:,iz] * np.sqrt(2.0*np.pi)) * \
                np.exp(-0.5 * ((zz[iy,:])-z0[:,iz])*((zz[iy,:])-z0[:,iz]) /(hp[:,iz]*hp[:,iz]))


    #if (ppar['gap_rout']>ppar['rin']):
        #for igap in range(len(ppar['gap_rout'])):
            #for iy in range(grid.ny):
                #for ix in range(grid.nx):
                    #if ((rcyl[iy,ix]<ppar['gap_rout'][igap])&(rcyl[iy,ix]>ppar['gap_rin'][igap])):
                        #rho[ix,iy,:] = rho[ix,iy,:] * ppar['gap_drfact'][igap]
        
    
    return rho
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
            exit()
            return 0

    # Calculate continuum optical depth 
    data.getTau(axis='xy', wav=0.55)
    

    nspec = len(ppar['gasspec_mol_name'])
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

                ii = (data.dusttemp[:,iy,iz,0]<ppar['gasspec_mol_freezeout_temp'][sid])
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
    """
    Function to generate the gas velocity field input for RADMC3D from the output frames of FARGO
    
    INPUT:
    ------
        grid : a radmc3dGgrid class containing the spatial grid onto which the FARGO output should
                be interpolated
        ppar : dictionary containing all parameters of the RADMC3D setup

    OUTPUT:
    -------
        Returns the velocity field in spherical coordinates in cm/s
    """

#
# Read the surface density from a FARGO frame
# WARNING!!!
# fdata['sigma'] has a dimension [nphi, nr] 
# For other models I was using [nr, nphi] for sigma!!!!
#
    fdata = fargo_data(path=ppar['fargo_path'], frame=ppar['fargo_frame'], gfargo=ppar['gfargo'],\
                           sigma=False, vel=True)
    
    fdata['rotv'] = np.zeros([fdata['nphi'],fdata['nr']], dtype=np.float64)


#
# If the frame came from GFARGO
#
    if ppar['gfargo']:
        dum = sqrt(1.0/fdata['r'])
        for ip in range(fdata['nphi']):
            fdata['rotv'][ip,:] = dum
    else:
#
# If the frame came from FARGO
#
        calc_rotv = False
        try:
            rotvFile = open(ppar['fargo_path']+'/planet0.dat','r')
        except:
            print 'planet0.dat has not been found'
            print 'I will calculate the rotation velocity of the frames'
            calc_rotv = True

#
# Calculate the coordinate system rotation rate if no planet0.dat has been found
#
        if calc_rotv:
            dum = np.sqrt(gg*ppar['mstar'][0]/ppar['fargo_ulength']**3.) * \
                (fdata['r']*ppar['fargo_ulength'])/ ppar['fargo_uvelo']
            for ip in range(fdata['nphi']):
                fdata['rotv'][ip,:] = dum
        else:
#
# Read the coordinate system rotation rate from planet0.dat 
#
 
            dum = ' '
            while (dum != ''):
                dum = rotvFile.readline()
                dum = dum.split()
                if (int(dum[0]) == ppar['fargo_frame']):
                    omega = float(dum[8])
                    dum   = ''

            for ip in range(fdata['nphi']):
                fdata['rotv'][ip,:] = fdata['r'] * omega
            rotvFile.close()

#
# Now interpolate the velocity field onto our spatial grid
#

    vrad_sp = itp.RectBivariateSpline(fdata['phi'], fdata['r'], fdata['vrad'])
    vphi_sp = itp.RectBivariateSpline(fdata['phi'], fdata['r'], fdata['vphi'])
    rotv_sp = itp.RectBivariateSpline(fdata['phi'], fdata['r'], fdata['rotv'])

    iirf     = ((grid.x>=fdata['r'].min()*ppar['fargo_ulength'])&(grid.x<=fdata['r'].max()*ppar['fargo_ulength']))  
    dr       = grid.x[iirf]/ppar['fargo_ulength']
    dumvr    = vrad_sp(grid.z, dr).swapaxes(0,1) * ppar['fargo_uvelo']
    dumvp    = vphi_sp(grid.z, dr).swapaxes(0,1) * ppar['fargo_uvelo']
    dumvrot  = rotv_sp(grid.z, dr).swapaxes(0,1) * ppar['fargo_uvelo']

#
# Add an outer disk 
#
    #if ppar['add_outer_disk']:
        #iiro =(grid.x>=fdata['r'].max()*ppar['fargo_ulength'])
        #od_r = grid.x[iiro]
        #vkep = np.sqrt(gg*ppar['mstar']/od_r)
        #od_vphi = np.zeros([grid.nz, od_r.shape[0]])
        #for iz in range(grid.nz):
            #od_vphi[iz,:] = vkep

#
# Put everything together
#

    vel    = np.zeros([grid.nx, grid.ny, grid.nz, 3], dtype=np.float64)
    for iy in range(grid.ny):
        vel[iirf,iy,:,0] = dumvr
        vel[iirf,iy,:,2] = dumvp + dumvrot
        
    #if ppar['add_outer_disk']:
        #dum = od_vphi.swapaxes(0,1)
        #for iy in range(grid.ny):
            #vel[iiro,iy,:,1] = dum


    return vel


    

            

