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

Warped circumstellar disk model of Terquem & Bertout 1993

"""

def getModelDesc():
    """
    A one line description of the model
    """

    return 'A warped disk model of Terquem & Bertout 1993'

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

    defpar.append(['warp_model', "'tb93'", ' Name of the warp model'])
    defpar.append(['warp_dcomp', '300.0*au', ' Distance of the companion causing the warp (warp_dcomp>rdisk (!!))'])
    defpar.append(['warp_icomp', '10.0', " Inclination angle of the companion's orbit with respect to the plane of the disk [degrees]"])
    defpar.append(['warp_mcomp', '1.0*mju', ' Mass of the companion'])
    defpar.append(['warp_t', '1e5*year', ' Time at which the warp should be calculated'])

    return defpar







def getWarpZ0TB93(rcyl=None, phi=None, grid=None, mstar=None, dcomp=None, rdisk=None, \
        icomp=None, mcomp=None, t=None):
    """
    Function to create a twist-free warp (Terquem & Bertout 1993, A&A, 274, 291)
    This model was developed to describe the warp perturbation caused
    by a companion outside of the disk (R_DISK<D_COMP)
    USAGE

    Z0 = getWarpZ0TB93(CRD=CRD, MSTAR=MSTAR, DCOMP=DCOMP, RDISK=RDISK, $
                             ICOMP=ICOMP, MCOMP=MCOMP, T=T)

    OUTPUT

          Z0     : Scalar, height of the midplane above the unperturbed value (z0)

    KEYWORDS 

          rcyl   : A numpy array containing cylindrical r coordinate (cell centers!)
          phi    : A numpy array containing azimuthal angles (cell centers!)
          mstar  : Mass of the central star [g] 
          rdisk  : Outer radius of the disk
          dcomp  : Distance of the companion from the central star [cm]
          icomp  : Inclination angle of the obit of the secondary
                   with respect to the plane of the disk [radian]
          mcomp  : Mass of the secondary [g]
          t      : Time at which the twisted warp should be calculated 
                   If t=0 no precession will be present only the steady term 
    """

    
    if (grid!=None):
        rr, th = np.meshgrid(grid.x, grid.y)
        rcyl = rr * np.sin(th)
        z    = rr * np.cos(th)
        phi  = grid.z
        nx   = grid.nx
        ny   = grid.nz
        nz   = grid.ny
    else:
        nx  = rcyl.shape[0]
        ny  = phi.shape[0]
        nz  = grid.ny
    
#
# Reduced time
#

    tau       = np.sqrt(np.pi*gg*mstar/rdisk**3.0) * t

#
# Steady-state solution
#
    g  = 3.0/8.0 * (mcomp/mstar) * np.sin(2.0*icomp) * rdisk**3./dcomp**3.

#
# Calculate the cylindrical coordinates
#
    z0 = np.zeros([nx, ny, nz], dtype=np.float64)
    for ix in range(nx):
        for iz in range(nz):
            hrr  = g * rcyl[iz,ix]**4.0/rdisk**4. 

#
# Precession term
#
            sr = np.sqrt(np.pi)/2.0 * (rcyl[iz,ix]/rdisk)**(1.5) * 3.0/(2.0*np.pi) * \
                mcomp / mstar * (rdisk/dcomp)**3.0 * np.sin(2.0*icomp)
            
            z0[ix,:,iz] = sr * tau * rcyl[iz,ix] * np.sin(phi) + hrr * rcyl[iz,ix] * np.cos(phi)

    return z0
# -----------------------------------------------------------------------------------------------
# 
# -----------------------------------------------------------------------------------------------
def getWarpOmegaPrec(rcyl=None, mstar=None, dcomp=None, rdisk=None, \
        icomp=None, mcomp=None, t=None):

    sr = np.sqrt(np.pi)/2.0 * (rcyl/rdisk)**(1.5) * 3.0/(2.0*np.pi) * \
        mcomp / mstar * (rdisk/dcomp)**3.0 * np.sin(2.0*icomp)

    return sr*sqrt(np.pi)

# -----------------------------------------------------------------------------------------------
# 
# -----------------------------------------------------------------------------------------------
def getDustDensity(grid=None, ppar=None):

#
# Apply perturbations 
#
    
    z0      = np.zeros([grid.nx, grid.nz, grid.ny], dtype=np.float64)

    if ppar.has_key('warp_model'):
        if ppar['warp_model'].lower()=='tb93':
            z0 =getWarpZ0TB93(grid=grid, mstar=ppar['mstar'][0], \
                                     dcomp=ppar['warp_dcomp'], icomp=ppar['warp_icomp'], \
                                     mcomp=ppar['warp_mcomp'], t=ppar['warp_t'], rdisk=ppar['rdisk'])

            
    rr, th = np.meshgrid(grid.x, grid.y)
    rcyl = rr * np.sin(th)
    zz   = rr * np.cos(th)

    for ix in range(grid.nx):
        for iy in range(grid.ny):
            if (rcyl[iy,ix]>ppar['warp_dcomp']):
                z0[ix,:,iy] = 0.0

    rho = ppdisk.getDustDensity(z0=z0, grid=grid, ppar=ppar)
    return rho
