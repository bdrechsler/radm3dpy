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


from radmc3dPy.crd_trans import vrot
import sys
"""
PYTHON module for RADMC3D 
(c) Attila Juhasz 2011,2012,2013

Circumstellar disk with spiral waves

"""

def getModelDesc():
    """
    Function to provide a brief description of the model
    """

    return "Circumstellar disk model with spiral waves (spiral wake equation is taken from Oglivie & Lubov 2002)"
           

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

    defpar.append(['dcomp', '10.0*au', 'Distance of the companion exciting the spiral waves from the central star'])
    defpar.append(['eps', '0.1', 'Inverse of the Mach-number'])
    defpar.append(['sp_sig', '0.5*au', ' Sets the width of the spiral in the radial direction (if omitted 4pi*eps^2*r/2.3548 will be used)'])
    defpar.append(['spa', '0.1', 'Maximum relative perturbation (perturbed variables are multiplied by 1+spiral_amp)'])
    defpar.append(['nfold', '1', 'Sets how many times the spiral is folded'])
    defpar.append(['spiral_location', '2', "Apply perturbation 0 - inside, 1 - outside, 2 - inside and outside of the companion's orbit"])
    defpar.append(['azim_shift', '0.', "Azimuthal angle of the perturber's location [in radian]"])

    return defpar



def initSpiralAmpOL02(rin=None, rout=None, dcomp=None, nr=None):
    """
    Function to initialize some variables for creating spiral arms in
      protoplanetary disks (Ogilvie & Lubow 2002, MNRAS, 330, 950)
      This initialization is important if the spiral arms are folded
      more than once. In this case it is not trivial to find the arm
      centers at a given azimuthal angle. This function intialize a
      grid over which linear interpolation is used to find the arm centers.
    
    SYNTAX
    -----
    INIT = initSpiralAmpOL02(rin=rin, rout=rout, dcomp=dcomp, nr=nr)
    
    INPUT
    -----
    
             rin              : Inner radius of the grid (<= Rin_disk)
             rout             : Outer radius of the grid (>= Rout_disk)
             dcomp            : Distance of the companion from the central star
             nr               : Number of grid points
    
    OUTPUT
    ------
    
             init             : A dictionary containing sompe pre-calculated functions
             r_in             : Cylindrical radius if the companion is inside of the disk
             fx_in            : Pre-calculated functions if the companion is inside of the disk
             r_out            : Cylindrical radius if the companion is outside of the disk
             fx_out           : Pre-calculated functions if the companion is outside of the disk
    
    """
    r_in   = rin * (dcomp/rin)**(np.arange(nr, dtype=np.float64)/(float(nr)-1.)) / dcomp
    fx_in  = r_in**(3.0/2.0) - 3.0/2.0 * np.log(r_in)

    r_out  = dcomp * (rout/dcomp)**(np.arange(nr, dtype=np.float64)/(float(nr)-1.)) / dcomp
    fx_out = r_out**(3.0/2.0) - 3.0/2.0 * np.log(r_out)

    res    = {'r_in':r_in*dcomp, 'r_out':r_out*dcomp, 'fx_in':fx_in, 'fx_out':fx_out}
    return res



def getSpiralAmpOL02(rc=None, phi=None, grid=None, init=None, eps=None, spa=None, sp_sig=None, nfold=None, azim_shift=None, \
                       location=None):
    """
    Function to calculate the spiral wake caused by a companion in a
    protoplanetary disk (Ogilvie & Lubow 2002, MNRAS, 330, 950)

    SYNTAX
    ------

    INIT = getSpiralAmpOL02(rin=rin, rout=rout, dcomp=dcomp, nr)

    INPUT
    -----
             CRD        : Two elements wector containing the cylindrical
                           coordinates [r,phi] at which the wake should
                           be calculated
             INIT       : Structure containing some gridded data to
                           calculate the spiral arm centers at a given
                           azimuthal angle inside and outside of
                           the companion's oribit
             EPS        : Inverse of the Mach number (M^-1) the sound
                          speed in the disk is assumed to be eps*sqrt(G*Mstar/r)
             SPA        : Amplitude of the spiral perturbation (shouldn't we set it to 1d0?)
             SP_SIG     : Spiral arms are assumed to be gaussian in the
                           radial direction. This parameter sets the standard
                           deviation of the gaussian. If not set the default value (dr/r=4*pi*eps) is used
             NFOLD      : How many times the spiral arm is folded
             AZIM_SHIFT : Azimuthal angular shift of the position of the
                           planet. If zero, the planet is assumed to lie
                           on the positive half of the x axis
             INNER      : Set this keyword if spiral arms inside of the
                           planetary orbit should be calculated
    
    OUTPUT
    ------

             spiral_amp : A Numpy array containing the amplitude of the spiral arms 
                          as a function of cylindrical radius and azimuthal angle

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
        # I have to reverse the phi array since the azimuthal velocity goes counterclockwise
        phi    = np.flipud(grid.z)
        nr     = grid.nx
        nz     = grid.ny
        nphi   = grid.nz
    else:
        nr     = rc.shape[0]
        nphi   = phi.shape[0]
        nz     = 1 
        zz      = array([0.0], dtype=float64)
        rcyl,z = np.meshgrid(rc, zz)

# Set the default value for the radial extent of the spiral  fwhm = 4*pi*eps^2*r
    if not sp_sig:
        sp_sig = 4.*np.pi*eps*eps*rcyl/2.35482004503
    else:
        sp_sig = np.zeros(rcyl.shape, dtype=np.float64) + sp_sig

#
# Calculate how many times the spiral is folded
#
    if (azim_shift==None): azim_shift=0.0
    infold = np.ceil(float(nfold) + azim_shift / (2.0*np.pi))
  

    spiral_amp_inner = np.zeros([nr,nphi,nz], dtype=np.float64)
    spiral_amp_outer = np.zeros([nr,nphi,nz], dtype=np.float64)
    spiral_amp = np.zeros([nr,nphi,nz], dtype=np.float64)
#
# Now calculate the spiral wake if the spiral is inwards the companion's orbit
#
    for i in range(int(infold)):
        phi3     = (2.*np.pi - np.array(phi)) + float(i) * 2.0 * np.pi - azim_shift
        c        = 1.0 + 3./2. * eps * phi3
        rc       = np.interp(c, np.flipud(init['fx_in']), np.flipud(init['r_in']))
        azim_amp = 1.0 - phi3 / (float(nfold) * 2.* np.pi)
#
# Prevent the amplitude of the spiral to go below zero
#   
        azim_amp = azim_amp.clip(0.0)
       

        azim_amp = azim_amp.T
        trcyl    = rcyl.T 


       

#
# Now put evertying together 
#
        for ip in range(nphi):
            spiral_amp_inner[:,ip,:] = spiral_amp_inner[:,ip,:] + \
                    spa * azim_amp[ip] * np.exp(-(trcyl - rc[ip])**2. / sp_sig.T**2.)   
# 
# Calculate the spiral wake if the spiral is outwards the companion's orbit
#
    phi2 = np.array(phi) - azim_shift
    for i in range(int(infold)):

        phi3 = (np.array(phi) - azim_shift) + float(i)*2.0*np.pi
        c      = 1.0 + 3./2. * eps * phi3 
        rc     = np.interp(c, init['fx_out'], init['r_out']) 
               
        azim_amp = (1.0 - phi3 / (float(nfold) * 2.0 * np.pi))
#
# Prevent the amplitude of the spiral to go below zero
#   
        azim_amp = azim_amp.clip(0.0)
            
        azim_amp = azim_amp.T
        trcyl    = rcyl.T
#           
# Now put evertying together 
#
        for ip in range(nphi):
            spiral_amp_outer[:,ip,:] = spiral_amp_outer[:,ip,:] + \
                    spa * azim_amp[ip] * np.exp(-(trcyl - rc[ip])**2. / sp_sig.T**2.)   

    if location==0:
        spiral_amp = spiral_amp_inner
    if location==1:
        spiral_amp = spiral_amp_outer
    if location==2:
        spiral_amp = spiral_amp_inner + spiral_amp_outer

    #fig = figure()
    #ax  = fig.add_subplot(111, polar=True)
    #pdat = np.squeeze(spiral_amp[:,:,nz/2.])
    #contour(grid.z, grid.x/au, pdat, 50)
    #show()

    #xx = raw_input('')
    
    return spiral_amp


# ***********************************************************************************************
# 
# ***********************************************************************************************
def getDustDensity(grid=None, ppar=None):

#
# Apply perturbations 
#

    rr, th = np.meshgrid(grid.x, grid.y)
    rcyl = rr * np.sin(th)    

    sp_amp  = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64)

    if ppar.has_key('do_spiral_ol02'):
        if ppar['do_spiral_ol02']:
            init = initSpiralAmpOL02(rin=ppar['rin'], rout=ppar['rdisk'], dcomp=ppar['dcomp'], \
                                            nr=ppar['nx'])
            if ppar.has_key('sp_sig'):
                sp_amp = getSpiralAmpOL02(grid=grid,  init=init, eps=ppar['eps'], \
                                         spa=ppar['spa'], sp_sig=ppar['sp_sig'], \
                                         nfold=ppar['nfold'], azim_shift=ppar['azim_shift'], \
                                         location=ppar['spiral_location'])
            else:
                sp_amp = getSpiralAmpOL02(grid=grid,  init=init, eps=ppar['eps'], \
                                         spa=ppar['spa'], \
                                         nfold=ppar['nfold'], azim_shift=ppar['azim_shift'], \
                                         location=ppar['spiral_location'])

            sp_amp = sp_amp.swapaxes(1,2)
# ------------------------------------------------------------------------------------------
# Do the perturbation in the pressure scale height only
# ------------------------------------------------------------------------------------------
            if (ppar['pert_var']==1):
                sigma = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64)
             
# Calculate the surface density
                if ppar.has_key('sig0'):
                    if ppar['sig0']!=0.:
                        sig0 = ppar['sig0']

                if ppar.has_key('mdisk'):
                    if (ppar['mdisk']!=0.0):
                        sig0 = 1.0
                        
                dum    = sig0 * (rcyl/ppar['rdisk'])**ppar['plsig1']
                dum = dum.swapaxes(0,1)
                for iz in range(grid.nz):
                    sigma[:,:,iz] = dum

                mass = 0.0
                for ix in range(grid.nx):
                    dum = sigma[ix,:,0] * (grid.xi[ix+1]**2. - grid.xi[ix]**2.) * np.pi * \
                        (grid.zi[1] - grid.zi[0]) / (2.0*np.pi)
                    mass = mass + dum.sum()
                    
                if ppar.has_key('mdisk'):
                    if (ppar['mdisk']!=0.0):
                        sigma = sigma * ppar['mdisk'] / mass
 
# Calculate the pressure scale height

                hp = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64)
                dum = ppar['hrdisk'] * (rcyl/ppar['rdisk'])**ppar['plh'] * rcyl
                dum = dum.swapaxes(0,1)
                for iz in range(grid.nz):
                    hp[:,:,iz] = dum

# Now do the perturbation

                hp = hp * (1.0 + sp_amp) 

# Calculate the volume density of the whole disk
                rho = ppdisk.getDustDensity(grid=grid, sigma=sigma, hp=hp, ppar=ppar)



# ------------------------------------------------------------------------------------------
# Do the perturbation in the surface density only
# ------------------------------------------------------------------------------------------


            if (ppar['pert_var']==2):

                sigma = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64)

# Calculate the surface density
                if ppar.has_key('sig0'):
                    if ppar['sig0']!=0.:
                        sig0 = ppar['sig0']

                if ppar.has_key('mdisk'):
                    if (ppar['mdisk']!=0.0):
                        sig0 = 1.0
                        
                dum    = sig0 * (rcyl/ppar['rdisk'])**ppar['plsig1']
                dum = dum.swapaxes(0,1)
                for iz in range(grid.nz):
                    sigma[:,:,iz] = dum

                mass = 0.0
                for ix in range(grid.nx):
                    dum = sigma[ix,:,0] * (grid.xi[ix+1]**2. - grid.xi[ix]**2.) * np.pi * \
                        (grid.zi[1] - grid.zi[0]) / (2.0*np.pi)
                    mass = mass + dum.sum()
                    
                if ppar.has_key('mdisk'):
                    if (ppar['mdisk']!=0.0):
                        sigma = sigma * ppar['mdisk'] / mass
                
# Now do the perturbation

                sigma = sigma * (1.0 + sp_amp)

# Calculate the volume density of the whole disk
                rho = ppdisk.getDustDensity(grid=grid, sigma=sigma,ppar=ppar)


# ------------------------------------------------------------------------------------------
# Do the perturbation both in the surface density and in the pressure scale height
# ------------------------------------------------------------------------------------------
            if (ppar['pert_var']==3):
                sigma = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64)

# Calculate the surface density
                if ppar.has_key('sig0'):
                    if ppar['sig0']!=0.:
                        sig0 = ppar['sig0']
                
                if ppar.has_key('mdisk'):
                    if (ppar['mdisk']!=0.0):
                        sig0 = 1.0
                        
                dum    = sig0 * (rcyl/ppar['rdisk'])**ppar['plsig1']
                dum = dum.swapaxes(0,1)
                for iz in range(grid.nz):
                    sigma[:,:,iz] = dum

                mass = 0.0
                for ix in range(grid.nx):
                    dum = sigma[ix,:,0] * (grid.xi[ix+1]**2. - grid.xi[ix]**2.) * np.pi * \
                        (grid.zi[1] - grid.zi[0]) / (2.0*np.pi)
                    mass = mass + dum.sum()
                    
                if ppar.has_key('mdisk'):
                    if (ppar['mdisk']!=0.0):
                        sigma = sigma * ppar['mdisk'] / mass
                
# Calculate the pressure scale height

                hp = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64)
                dum = ppar['hrdisk'] * (rcyl/ppar['rdisk'])**ppar['plh'] * rcyl
                dum = dum.swapaxes(0,1)
                for iz in range(grid.nz):
                    hp[:,:,iz] = dum

# Now do the perturbation

                sigma = sigma * (1.0 + sp_amp)
                hp = hp * (1.0 + sp_amp) 

# Calculate the volume density of the whole disk
                rho = ppdisk.getDustDensity(grid=grid, sigma=sigma, hp=hp, ppar=ppar)


    return rho
# ***********************************************************************************************
# 
# ***********************************************************************************************
def getGasDensity(grid=None, ppar=None):

#
# Apply perturbations 
#

    rr, th = np.meshgrid(grid.x, grid.y)
    rcyl = rr * np.sin(th)    

    sp_amp  = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64)

    if ppar.has_key('do_spiral_ol02'):
        if ppar['do_spiral_ol02']:
            init = initSpiralAmpOL02(rin=ppar['rin'], rout=ppar['rdisk'], dcomp=ppar['dcomp'], \
                                            nr=ppar['nx'])
            if ppar.has_key('sp_sig'):
                sp_amp = getSpiralAmpOL02(grid=grid,  init=init, eps=ppar['eps'], \
                                         spa=ppar['spa'], sp_sig=ppar['sp_sig'], \
                                         nfold=ppar['nfold'], azim_shift=ppar['azim_shift'], \
                                         location=ppar['spiral_location'])
            else:
                sp_amp = getSpiralAmpOL02(grid=grid,  init=init, eps=ppar['eps'], \
                                         spa=ppar['spa'], \
                                         nfold=ppar['nfold'], azim_shift=ppar['azim_shift'], \
                                         location=ppar['spiral_location'])

            sp_amp = sp_amp.swapaxes(1,2)
# ------------------------------------------------------------------------------------------
# Do the perturbation in the pressure scale height only
# ------------------------------------------------------------------------------------------
            if (ppar['pert_var']==1):
                sigma = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64)
             
# Calculate the surface density
                if ppar.has_key('sig0'):
                    if ppar['sig0']!=0.:
                        sig0 = ppar['sig0']

                if ppar.has_key('mdisk'):
                    if (ppar['mdisk']!=0.0):
                        sig0 = 1.0
                        
                dum    = sig0 * (rcyl/ppar['rdisk'])**ppar['plsig1']
                dum = dum.swapaxes(0,1)
                for iz in range(grid.nz):
                    sigma[:,:,iz] = dum

                mass = 0.0
                for ix in range(grid.nx):
                    dum = sigma[ix,:,0] * (grid.xi[ix+1]**2. - grid.xi[ix]**2.) * np.pi * \
                        (grid.zi[1] - grid.zi[0]) / (2.0*np.pi)
                    mass = mass + dum.sum()
                    
                if ppar.has_key('mdisk'):
                    if (ppar['mdisk']!=0.0):
                        sigma = sigma * ppar['mdisk'] / mass
 
# Calculate the pressure scale height

                hp = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64)
                dum = ppar['hrdisk'] * (rcyl/ppar['rdisk'])**ppar['plh'] * rcyl
                dum = dum.swapaxes(0,1)
                for iz in range(grid.nz):
                    hp[:,:,iz] = dum

# Now do the perturbation

                hp = hp * (1.0 + sp_amp) 

# Calculate the volume density of the whole disk
                rho = ppdisk.getGasDensity(grid=grid, sigma=sigma, hp=hp, ppar=ppar)



# ------------------------------------------------------------------------------------------
# Do the perturbation in the surface density only
# ------------------------------------------------------------------------------------------


            if (ppar['pert_var']==2):

                sigma = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64)

# Calculate the surface density
                if ppar.has_key('sig0'):
                    if ppar['sig0']!=0.:
                        sig0 = ppar['sig0']

                if ppar.has_key('mdisk'):
                    if (ppar['mdisk']!=0.0):
                        sig0 = 1.0
                        
                dum    = sig0 * (rcyl/ppar['rdisk'])**ppar['plsig1']
                dum = dum.swapaxes(0,1)
                for iz in range(grid.nz):
                    sigma[:,:,iz] = dum

                mass = 0.0
                for ix in range(grid.nx):
                    dum = sigma[ix,:,0] * (grid.xi[ix+1]**2. - grid.xi[ix]**2.) * np.pi * \
                        (grid.zi[1] - grid.zi[0]) / (2.0*np.pi)
                    mass = mass + dum.sum()
                    
                if ppar.has_key('mdisk'):
                    if (ppar['mdisk']!=0.0):
                        sigma = sigma * ppar['mdisk'] / mass
                
# Now do the perturbation

                sigma = sigma * (1.0 + sp_amp)

# Calculate the volume density of the whole disk
                rho = ppdisk.getGasDensity(grid=grid, sigma=sigma,ppar=ppar)


# ------------------------------------------------------------------------------------------
# Do the perturbation both in the surface density and in the pressure scale height
# ------------------------------------------------------------------------------------------
            if (ppar['pert_var']==3):
                sigma = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64)

# Calculate the surface density
                if ppar.has_key('sig0'):
                    if ppar['sig0']!=0.:
                        sig0 = ppar['sig0']
                
                if ppar.has_key('mdisk'):
                    if (ppar['mdisk']!=0.0):
                        sig0 = 1.0
                        
                dum    = sig0 * (rcyl/ppar['rdisk'])**ppar['plsig1']
                dum = dum.swapaxes(0,1)
                for iz in range(grid.nz):
                    sigma[:,:,iz] = dum

                mass = 0.0
                for ix in range(grid.nx):
                    dum = sigma[ix,:,0] * (grid.xi[ix+1]**2. - grid.xi[ix]**2.) * np.pi * \
                        (grid.zi[1] - grid.zi[0]) / (2.0*np.pi)
                    mass = mass + dum.sum()
                    
                if ppar.has_key('mdisk'):
                    if (ppar['mdisk']!=0.0):
                        sigma = sigma * ppar['mdisk'] / mass
                
# Calculate the pressure scale height

                hp = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64)
                dum = ppar['hrdisk'] * (rcyl/ppar['rdisk'])**ppar['plh'] * rcyl
                dum = dum.swapaxes(0,1)
                for iz in range(grid.nz):
                    hp[:,:,iz] = dum

# Now do the perturbation

                sigma = sigma * (1.0 + sp_amp)
                hp = hp * (1.0 + sp_amp) 

# Calculate the volume density of the whole disk
                rho = ppdisk.getGasDensity(grid=grid, sigma=sigma, hp=hp, ppar=ppar)


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
    """
    Function to calculate the velocity in a disk with a spiral wave in it.
        The velocity field is a simple Keplerian field with an additional
        perturbation on top of that.
        vphi = spira_vamp * exp(-r**2 / spiral
    """

# First calculate the cylindrical radii for each grid point and the corresponding circular Kepler velocity
    rcyl = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64)
    gx, gy = meshgrid(grid.x, grid.y)
    gx, gy = gx.T, gy.T
    rcyl2d = gx * sin(gy)
    for iz in range(grid.nz):
        rcyl[:,:,iz] = rcyl2d

    # OK, so we have to do a not really physical trick. Close to the rotation axis (pi/2-theta~ ~+-pi/2) 
    #  The cylindrical radius is very small, creating insanely high velocities. Thus I clip the cylindrical
    #  radius to be the inner disk radius. This is of course not physical, but the results should be fine as
    #  there shouldn't be much emission coming from high above the midplane, close to the pole. That would
    #  be rather a jet which is not part of this model description. 
    rcyl = rcyl.clip(ppar['rin'], rcyl.max())
    vkep = np.sqrt(gg * ppar['mstar'][0] / rcyl)
    

# Now calculate the sound speed in the disk midplane

    # This formula is merely for debugging purposes and it should be replaced by the real 
    #  calculated dust temperature

    # Calculate the planet's velocity
    v_planet = sqrt(gg*ppar['mstar'][0] / ppar['dcomp'])

    # Calculate the spiral wake
    init = initSpiralAmpOL02(rin=ppar['rin'], rout=ppar['rdisk'], dcomp=ppar['dcomp'], \
                                    nr=ppar['nx'])

    if ppar.has_key('sp_sig'):
        sp_amp = getSpiralAmpOL02(grid=grid,  init=init, eps=ppar['eps'], \
                                 spa=1.0, sp_sig=ppar['sp_sig'], \
                                 nfold=ppar['nfold'], azim_shift=ppar['azim_shift'], \
                                 location=ppar['spiral_location'])
    else:
        sp_amp = getSpiralAmpOL02(grid=grid,  init=init, eps=ppar['eps'], \
                                 spa=1.0, \
                                 nfold=ppar['nfold'], azim_shift=ppar['azim_shift'], \
                                 location=ppar['spiral_location'])

    sp_amp = sp_amp.swapaxes(1,2)

    ny = sp_amp.shape[1]


# Now calculate the velocity field so that in the spiral the azimuthal velocity component
# is v_planet and there is a smooth transition to a Keplerian velocity field outside 
# of the spiral 
   
    vel = np.zeros([grid.nx, grid.ny, grid.nz, 3], dtype=np.float64)
    vel[:,:,:,2] = vkep + (v_planet-vkep) * sp_amp
  
    fig = figure()
    ax  = fig.add_subplot(111, polar=True)
    pdat = np.squeeze(vel[:,ny/2,:,2])
    c = contourf(grid.z, grid.x, pdat, 50)
    cbar = colorbar(c)

    fig2 = figure()
    ax = fig.add_subplot(111)
    ir = 20
    plot(grid.z, vel[ir,ny/2,:,2]/1e5)
    title('R = '+str(grid.x[ir]/au))
    show()

    xx = raw_input('')
    return vel

