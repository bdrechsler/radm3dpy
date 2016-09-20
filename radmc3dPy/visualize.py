"""This module contains classes and functions to read and write input/output data for RADMC-3D and
to do some simple analysis/diagnostics of the model.

"""
try:
    import numpy as np
except:
    print 'ERROR'
    print ' Numpy cannot be imported '
    print ' To use the python module of RADMC-3D you need to install Numpy'

try:
    import matplotlib.pylab as plb
except:
    print ' WARNING'
    print ' matploblib.pylab cannot be imported ' 
    print ' To used the visualization functionality of the python module of RADMC-3D you need to install matplotlib'
    print ' Without matplotlib you can use the python module to set up a model but you will not be able to plot things or'
    print ' display images'

import subprocess as sp
import sys
import os 
import copy 
import radmc3dPy.natconst as nc
import radmc3dPy.crd_trans  as crd_trans 
from staratm import StellarAtm

# --------------------------------------------------------------------------------------------------
def plotSpectrum(a,ev=False,kev=False,hz=False,micron=False,jy=False,lsun=False,
                lnu=False,nulnu=False,fnu=False,nufnu=False,dpc=1.e0,
                oplot=False,xlg=False,ylg=False,obs=False,
                mol=None,ilin=None):
   """Plot the spectrum / SED 

   Parameters
   ----------
   a               : ndarray
                    A 2D array of size [Nfreq,2] returned by readSpectrum(). 
                    [:,0] - wavelength in micrometer, or for line data the velocity in km/s
                    [:,1] - flux density in erg/s/cm/cm/Hz
   ev              : bool
                    True --> frequency in electronvolt (default=Hz)

   kev             : bool 
                    True --> frequency in kiloelectronvolt (default=Hz)

   micron          : bool
                    True --> wavelength in micron (default=Hz)

   jy              : bool
                    True --> Flux in Jansky

   lnu             : bool
                    True --> L_nu (default L_nu)

   nulnu           : bool
                    True --> nu*L_nu (default F_nu)

   lsun            : bool
                    True --> nu*L_nu in units of solar luminosity

   dpc             : bool
                    Distance of observer in units of parsec (Default: 1 pc)

   oplot           : bool
                    True --> Plot without refreshing subplot

   xlg             : bool
                    True --> logarithmic x-axis

   ylg             : bool
                    True --> logarithmic y-axis

   obs             : bool
                    True --> Treat the spectrum as an observation
                              (i.e. do not scale with dpc^(-2))

   mol             : bool
                    (optional) Molecule data (see radmc3dMolecule class)
                     This is required if you want to plot a line spectrum
                     with on the x-axis the radial velocity in km/s

   ilin            : bool
                    (if set) the index of the line (of mol; starting,
                     as in RADMC-3D, with the index 1) which shall act
                     as the 0 km/s wavelength reference. If ilin is set
                     the x axis will be in km/s (overriding other settings)

   """
   #
   # Basic
   #
   lam    = a[:,0]
   fluxnu = a[:,1]
   #
   # Calculate frequency in Hz
   #
   cc    = 2.9979245800000e10     # Light speed             [cm/s]
   freq  = 1e4*nc.cc/lam
   #
   # Default: frequency in Hz
   #
   xcoord = freq
   xtitle = '$\lambda [\mu\mathrm{m}]$'
   #
   # If ev: electronvolt
   #
   if ev:
       xcoord = 4.13568842841e-15 * freq
       xtitle = '$h\\nu [\mathrm{eV}]$'
   #
   # If kev: kiloelectronvolt
   #
   if kev:
       xcoord = 4.13568842841e-18 * freq
       xtitle = '$h\\nu [\mathrm{KeV}]$'
   #
   # If micron
   #
   if micron:
       xcoord = lam
       xtitle = '$h\\nu [\mathrm{KeV}]$'
   #
   # Plot nuFnu or Fnu (same with Lnu)? And what about Fnu vs Lnu?
   #
   # Default:
   sed=True
   ylum=False
   # The flags:
   if jy:
       sed=False
   if fnu:
       sed=False
       ylum=False
   if lnu:
       sed=False
       ylum=True
   if nulnu:
       sed=True
       ylum=True
   if fnu:
       sed=False
       ylum=False
   if nufnu:
       sed=True
       ylum=False
   if jy:
       ylum=False
   if lsun:
       ylum=True
       sed=True
   #
   # If ilin is set, then override the above and use instead the line
   # as a reference and use km/s as x-axis
   #
   if ilin is not None:
       if mol is None:
           print "Error in plotSpectrum(): if you specify ilin, you must give a molecule with mol=..."
           return
       else:
           freq0  = mol.freq[ilin-1]
           xcoord = 2.9979245800000e+10*(freq0-freq)/freq0/1.e5
           xtitle = '$\Delta v [\mathrm{km/h}]$'
   #
   # Which plot to make? Lum or flux?
   #
   if not ylum:
       #
       # Plot spectrum as flux at a certain distance
       #
       if not obs:
           distfact = 1.0 / (dpc**2)
       else:
           distfact = 1.0
       #
       # Set the vertical axis name
       #
       if not jy:
           if not sed:
               lumfact=1.0
               ytitle='$F_{\\nu}\; [\mathrm{erg}\,\mathrm{cm}^{-2}\,\mathrm{Hz}^{-1}\, \mathrm{s}^{-1}]$'
           else:
               lumfact=1.0*freq
               ytitle='$\\nu F_{\\nu}\; [\mathrm{erg}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}]$'
       else:
           if not sed:
               lumfact=1e+23
               ytitle='$F_{\\nu} [Jy]$'
           else:
               lumfact=1e+23*freq
               ytitle='$\\nu F_{\\nu} [JyHz]$'
   else:
       #
       # Plot spectrum as luminosity
       #
       if not obs:
           distfact = 1.1965280793e38   # = 4*pi*(1 parsec)^2 = 1.19d38 cm^2
       else:
           distfact = dpc**2 * 1.1965280793e38

       if not sed:
           lumfact=1.e0
           ytitle='L_{\\nu}\; [\mathrm{erg}\,\mathrm{Hz}^{-1}\, \mathrm{s}^{-1}]'
       else:
           if not lsun:
               lumfact = 1.0*freq
               ytitle  = '\\nu L_{\\nu}\; [\mathrm{erg}\, \mathrm{s}^{-1}]'
           else:
               lumfact = freq * 2.5956986e-34
               ytitle  = '\\nu L_{\\nu}\; [L_{\odot}]'

   #
   # The data on the y axis
   #
   ycoord = distfact*lumfact*fluxnu
   #
   # If not oplot, then reset the subplot and set the axes
   #
   if not oplot:
       plb.cla()
       if xlg:
           plb.xscale('log')
       if ylg:
           plb.yscale('log')
       plb.xlabel(xtitle)
       plb.ylabel(ytitle)
   #
   # Now plot
   #
   plb.plot(xcoord,ycoord)


# --------------------------------------------------------------------------------------------------
def plotStruct2D(d=None, var='ddens', scale='lin', iplane=0, icrd3=0., crd3=None, ispec=0, xlim=(), ylim=(), 
        minval=None, maxval=None, contours=False, clev=None, clmin=None, clmax=None, ncl=None, cllog=False, 
        clcol='k', cmap=None, overlay=False, alpha=1.0, ax=None, projection='spherical', au=True, lattitude=True, 
        deg=False):

    """ 
    Creates a contour plot of the crossection of the disk at various planes.

    Parameters:
    -----------
        
        d        : radmc3dData
                   An instance of the radmc3dData class, containing the variable to be plotted

        var      : {'ddens', 'dtemp', 'gdens', 'ndens', 'gtemp', 'vx', 'vy', 'vz', 'vturb', 'taux', 'tauy'}
                   The variable to be plotted

        scale    : {'lin', 'log'}
                   Linear or logarithmic scale to be used in the contours

    Options:
    --------

        contours : bool
                   If True contour lines are plotted, if False a colorscale plot will be created

        clev     : ndarray  
                   A numpy ndarray containing the levels to be displayed with contour lines. If clev is set
                   then clmin, clmax and ncl are omitted

        clmin    : float
                   Min. contour level (for setting auto-contours between clmin and clmax at ncl values)

        clmax    : float
                   Max. contour level (for setting auto-contours between clmin and clmax at ncl values)
        
        ncl      : float
                   Number of contour levels (for setting auto-contours between clmin and clmax at ncl values)

        cllog    : bool
                   If clmin, clmax and ncl are used to generate the contour levels, then if cllog is True
                   the contours will be log-scaled

        clcol    : str
                   Color-code for the contour lines for single color contours

        cmap     : matplotib color map
                   A color map used that could be used both for the color-scale image and also for the contour
                   lines if clcol is not set

        ax       : matplotlib.Figure.Axis
                   A matplotlib axis in which the plot should be made

        projection : {'cartesian', 'polar'}
                   Coordinate system of the plot

        au       : bool
                   If set to true, linear coordinates will have the unit of AU

        lattitude : bool
                   If the coordinate sytem used in RADMC-3D is spherical, then the 2nd coordiante is the co-lattitude.
                   If lattitude is set to True then the 2nd coordinate in the RADMC-3D grid will be transformet to true
                   lattitude (i.e. pi/2.-colattitude). If set to false the original co-lattitude will be used. 

        deg      : bool
                   If set to True angular coordinates will be displayed in degrees instead of radians. 
                
        overlay  : bool
                   If true the plot will be overlayed over an existing plot on the given axis

        alpha    : float
                   Transparency of the plot (0. - transparent, 1.0 - opaque)
             
    """

    #
    # Check the input consistency
    #
    if d==None:
        print 'ERROR'
        print ' No data to be plotted'

        return 


    #
    # Check what to plot
    #
    if var.strip().lower()=='ddens':
        if type(d.rhodust)==int:  
            print 'ERROR'
            print 'Dust density is not present in the passed radmc3dData instance'
            return
        else:
            data = d.rhodust
            d_label = r'$\rho_{\rm dust}$ [g/cm$^3$]'

    elif var.strip().lower()=='dtemp':
        if type(d.dusttemp)==int:  
            print 'ERROR'
            print 'Dust temperature is not present in the passed radmc3dData instance'
            return
        else:
            data = d.dusttemp
            d_label = r'$T_{\rm dust}$ [K]'
    
    elif var.strip().lower()=='gdens':
        if type(d.rhogas)==int:  
            print 'ERROR'
            print 'Gas density is not present in the passed radmc3dData instance'
            return
        else:
            data = d.rhogas
            d_label = r'$\rho_{\rm gas}$ [g/cm$^3$]'
    
    elif var.strip().lower()=='ndens':
        if type(d.ndens_mol)==int:  
            print 'ERROR'
            print 'Gas number density is not present in the passed radmc3dData instance'
            return
        else:
            data = d.ndens_mol
            d_label = r'$n_{\rm gas}$ [molecule/cm$^3$]'
    
    elif var.strip().lower()=='gtemp':
        if type(d.gastemp)==int:  
            print 'ERROR'
            print 'Gas temperture is not present in the passed radmc3dData instance'
            return
        else:
            data = d.gastemp
            d_label = r'$T_{\rm gas}$ [K]'
    
    elif var.strip().lower()=='vx':
        if type(d.gvel)==int:  
            print 'ERROR'
            print 'Gas velocity is not present in the passed radmc3dData instance'
            return
        else:
            data = d.gvel[:,:,:,0]
            d_label = r'$v_{\rm x}$ [cm/s]'
    
    elif var.strip().lower()=='vy':
        if type(d.gvel)==int:  
            print 'ERROR'
            print 'Gas velocity is not present in the passed radmc3dData instance'
            return
        else:
            data = d.gvel[:,:,:,1]
            d_label = r'$v_{\rm y}$ [cm/s]'
    elif var.strip().lower()=='vz':
        if type(d.gvel)==int:  
            print 'ERROR'
            print 'Gas velocity is not present in the passed radmc3dData instance'
            return
        else:
            data = d.gvel[:,:,:,2]
            d_label = r'$v_{\rm z}$ [cm/s]'
    
    elif var.strip().lower()=='vturb':
        if type(d.vturb)==int:  
            print 'ERROR'
            print 'Microturbulent velocity is not present in the passed radmc3dData instance'
            return
        else:
            data = d.vturb
            d_label = r'$v_{\rm turb}$ [cm/s]'
    
    if var.strip().lower()=='taux':
        if type(d.taux)==int:  
            print 'ERROR'
            print 'Optical depth is not present in the passed radmc3dData instance'
            return
        else:
            data = d.taux
            d_label = r'$\tau_{\rm r}$'
    
    if var.strip().lower()=='tauy':
        if type(d.tauy)==int:  
            print 'ERROR'
            print 'Optical depth is not present in the passed radmc3dData instance'
            return
        else:
            data = d.tauy
            d_label = r'$\tau_{\rm \theta}$'

    if scale == 'log':
        data = np.log10(data)

    
    if contours == True:
        if clev == None:
            if clmin == None:
                clmin = data.min()
            if clmax == None:
                clmax = data.max()
            if ncl == None:
                ncl = 12


    if d.grid.crd_sys == 'sph':
        # r-theta plane
        if iplane == 0:
            # Get the third coordinate index
            if ( (icrd3 == None) & (crd3 == None) ):
                print 'ERROR'
                print 'Either icrd3 or crd3 should be specificied'
                return 
            else:
                if crd3 != None:
                    icrd3 = abs(d.grid.z - crd3).argmin()

            # Select the data dimensions
            if len(data.shape)==4:
                if ispec>=0: 
                    data = data[:,:,icrd3,ispec]
                else:
                    data = data[:,:,icrd3,:]
                    data.sum(2)
            else:
                data = data[:,:,icrd3]



            if projection == 'spherical':
                xx,yy = np.meshgrid(d.grid.x, d.grid.y)
                x_label = 'r'
                y_label = r'$\theta$'
                if au == True:
                    xx /= nc.au
                    x_label += ' [AU]'
                else:
                    x_label += ' [cm]'
                
                if lattitude == True:
                    yy = np.pi/2.-yy
                    y_label = r'$\pi/2 - \theta$'

                if deg == True:
                    yy *= 180./np.pi
                    y_label += ' [deg]'
                else:
                    y_label += ' [rad]'

            elif projection == 'cartesian':
                rr,tt = np.meshgrid(d.grid.x, d.grid.y)
                xx    = rr*np.sin(tt)
                yy    = rr*np.cos(tt)
                x_label = 'x'
                y_label = 'y'

                if au == True:
                    xx /= nc.au
                    yy /= nc.au
                    x_label += ' [AU]'
                    y_label += ' [AU]'
                else:
                    x_label += ' [cm]'
                    y_label += ' [cm]'

            else:
                print 'ERROR'
                print 'Unkonwn projection ', projection
                print 'Projection should be either "spherical", or "cartesian"'
                return
        
        # r-phi plane
        elif iplane == 1:
            # Get the third coordinate index
            if ( (icrd3 == None) & (crd3 == None) ):
                print 'ERROR'
                print 'Either icrd3 or crd3 should be specificied'
                return 
            else:
                if crd3 != None:
                    if lattitude == True:
                        icrd3 = abs( (np.pi/2.-d.grid.y) - crd3).argmin()
                    else:
                        icrd3 = abs(d.grid.y - crd3).argmin()

            # Select the data dimensions
            if len(data.shape)==4:
                if ispec>=0: 
                    data = data[:,icrd3,:,ispec]
                else:
                    data = data[:,icrd3,:,:]
                    data.sum(2)
            else:
                data = data[:,icrd3,:]
           
            if projection == 'spherical':
                xx,yy = np.meshgrid(d.grid.x, d.grid.z)
                x_label = 'r'
                y_label = r'$\phi$'
                if au == True:
                    xx /= nc.au
                    x_label += ' [AU]'
                else:
                    x_label += ' [cm]'

                if deg == True:
                    xx *= 180./np.pi
                    yy *= 180./np.pi
                    y_label += ' [deg]'
                else:
                    y_label += ' [rad]'


            elif projection == 'cartesian':

                rr,pp = np.meshgrid(d.grid.x, d.grid.z)
                xx    = rr*np.sin(tt)
                yy    = rr*np.cos(tt)
                x_label = 'x'
                y_label = 'y'
               
                if au == True:
                    xx /= nc.au
                    yy /= nc.au
                    x_label += ' [AU]'
                    y_label += ' [AU]'
                else:
                    x_label += ' [cm]'
                    y_label += ' [cm]'

            else:
                print 'ERROR'
                print 'Unkonwn projection ', projection
                print 'Projection should be either "spherical", or "cartesian"'
                return



        # phi-theta plane
        elif iplane == 2:
           
            # Get the third coordinate index
            if ( (icrd3 == None) & (crd3 == None) ):
                print 'ERROR'
                print 'Either icrd3 or crd3 should be specificied'
                return 
            else:
                if crd3 != None:
                    icrd3 = abs(d.grid.x - crd3).argmin()

            # Select the data dimensions
            if len(data.shape)==4:
                if ispec>=0: 
                    data = data[icrd3,:,:,ispec]
                else:
                    data = data[icrd3,:,:,:]
                    data.sum(2)
            else:
                data = data[icrd3,:,:]

            if projection == 'spherical':
                xx,yy = np.meshgrid(d.grid.z, d.grid.y)
                data  = data.T
                x_label = r'$\phi$'
                y_label = r'$\theta$'

                if lattitude == True:
                    yy = np.pi/2. - yy
                    y_label = r'$\pi/2 - \theta$'

                if deg == True:
                    xx *= 180./np.pi
                    yy *= 180./np.pi
                    y_label += ' [deg]'
                    x_label += ' [deg]'
                else:
                    y_label += ' [rad]'
                    x_label += ' [rad]'

            elif projection == 'cartesian':

                tt,pp = np.meshgrid(d.grid.y, d.grid.z)
                rr = d.grid.x[icrd3] 
                xx = pp*rr
                yy = (pi/2.-tt)*rr
                x_label = 'x'
                y_label = 'y'
                
                if au == True:
                    xx /= nc.au
                    yy /= nc.au
                    x_label += ' [AU]'
                    y_label += ' [AU]'
                else:
                    x_label += ' [cm]'
                    y_label += ' [cm]'

            else:
                print 'ERROR'
                print 'Unkonwn projection ', projection
                print 'Projection should be either "spherical", or "cartesian"'
                return



    # Make the plot
    if ax == None:
        ax = plb.gca()

    # Clear the axis if overlay==True
    if overlay == False:
        ax.cla()


    if cmap == None:
        cmap = plb.cm.jet
       
    if contours == True:
        # Generate the contour levels
        if clev == None:
            if cllog == True:
                clev = clmin * (clmax/clmin)**(np.arange(ncl, dtype=float)/float(ncl-1))
            else:
                clev = clmin + (clmax-clmin)*(np.arange(ncl, dtype=float)/float(ncl-1))
        if clcol == 'none':
            if cmap != None:
                ax.contour(xx, yy, data, clev, cmap, alpha=alpha)
            else:
                print 'ERROR'
                print 'If clcol=="none" cmap should be specified'
                return
        else:
            ax.contour(xx, yy, data.T, clev, colors=clcol, alpha=alpha)
        
        cbar = None
    else:
        if minval == None:
            minval = data.min()
        if maxval == None:
            maxval = data.max()

        plb.axes(ax)
        plb.pcolormesh(xx,yy,data.T,cmap=cmap,alpha=alpha, vmin=minval, vmax=maxval)
        cbar = plb.colorbar(ax=ax)
        cbar.set_label(d_label)
        #cbar = None
    
    if len(xlim)==2:
        plb.xlim(xlim[0], xlim[1])
    else:
        plb.xlim(xx.min(), xx.max())
    if len(ylim)==2:
        plb.ylim(ylim[0], ylim[1])
    else:
        plb.xlim(xx.min(), xx.max())
   
    plb.xlabel(x_label)
    plb.ylabel(y_label)
    
    return {'ax':ax, 'cb':cbar}
     

