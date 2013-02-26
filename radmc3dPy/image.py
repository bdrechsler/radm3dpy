"""
PYTHON module for RADMC3D 
(c) Attila Juhasz 2011,2012,2013

This sub-module contains classes/functions to create and read images with radmc3d and to calculate
interferometric visibilities and write fits files
For help on the syntax or functionality of each function see the help of the individual functions

CLASSES:
--------
    radmc3dImage - RADMC3D image class
    radmc3dVisibility - Class of interferometric visibilities

FUNCTIONS:
----------

    get_psf() - Calculates a Gaussian PSF/beam
    get_visibility() - Calculates interferometric visiblities
    makeimage() - Runs RADMC3D to calculate images/channel maps
    plotimage() - Plots the image
    readimage() - Reads RADMC3D image(s)

"""
try:
    from matplotlib.pylab import *
except:
    print ' WARNING'
    print ' matploblib.pylab cannot be imported ' 
    print ' To used the visualization functionality of the python module of RADMC-3D you need to install matplotlib'
    print ' Without matplotlib you can use the python module to set up a model but you will not be able to plot things or'
    print ' display images'
try:
    from numpy import *
except:
    print 'ERROR'
    print ' Numpy cannot be imported '
    print ' To use the python module of RADMC-3D you need to install Numpy'

try:
    from scipy.interpolate import RectBivariateSpline as rbs
    from scipy.interpolate import BivariateSpline as bvspline
except:
    print 'WARNING'
    print ' Scipy.interpolate.BivariateSpline or Scipy.interpolate.RectBivariateSpline cannot be imported '

try:
    import pyfits as pf
except:
    print 'WARNING'
    print ' PyFits cannot be imported'
    print ' The PyFits module is needed to write RADMC-3D images to FITS format'
    print ' Without PyFits no fits file can be written'

from copy import deepcopy
import subprocess as sp
import sys, os


# **************************************************************************************************
class radmc3dVisibility():
    """
    Visiblity class to be returned by the get_vis function 
    """
    def __init__(self):
        self.fftImage = 0
        self.vis = 0
        self.amp = 0
        self.phase = 0
        self.u = 0
        self.v = 0
        self.bl = 0
        self.pa = 0
        self.blu = 0
        self.blv = 0 
        self.mu = 0
        self.mu_err = 0
        self.mv = 0 
        self.mv_err = 0 
        self.wav = 0

# **************************************************************************************************
class radmc3dImage():
    """
    RADMC3D image class
    """
    def __init__(self):
        self.image = 0
        self.imageJyppix = 0
        self.x = 0 
        self.y = 0
        self.nx = 0
        self.ny = 0
        self.sizepix_x = 0
        self.sizepix_y = 0
        self.nfreq = 0 
        self.freq = 0
        self.nwav = 0 
        self.wav = 0


# --------------------------------------------------------------------------------------------------
    def write_casafits(self, fname='', dpc=1., coord='03h10m05s -10d05m30s', bandwidthmhz=2000.0):
        """
        Function to write out a RADMC3D image data in fits format that is compatible with CASA
  
        INPUT:
        ------
         fname   : file name of the radmc3d output image (if omitted 'image.fits' is used)
         coord   : image center coordinates
         bandwidthmhz : if the image is in the continuum set the bandwidth of the image
        """
# --------------------------------------------------------------------------------------------------
        if fname=='':
            fname = 'image.fits'
        pc = 3.0857200e+18

        # Decode the image center cooridnates

        # Check first whether the format is OK
        dum = coord

        ra = []
        delim = ['h', 'm', 's']
        for i in delim:
            ind = dum.find(i)
            if ind<=0:
                print 'ERROR'
                print 'coord keyword has a wrong format'
                print 'coord="0h10m05s -10d05m30s"'
                print ra
                print dum
                return
            ra.append(float(dum[:ind]))
            dum = dum[ind+1:]

        dec = []
        delim = ['d', 'm', 's']
        for i in delim:
            ind = dum.find(i)
            if ind<=0:
                print 'ERROR'
                print 'coord keyword has a wrong format'
                print 'coord="0h10m05s -10d05m30s"'
                return
            dec.append(float(dum[:ind]))
            dum = dum[ind+1:]



        target_ra = (ra[0] + ra[1]/60. + ra[2]/3600.) * 15.
        target_dec = (dec[0] + dec[1]/60. + dec[2]/3600.) 

        # Conversion from erg/s/cm/cm/ster to Jy/pixel
        conv = self.sizepix_x * self.sizepix_y / (dpc * pc)**2. * 1e23

        # Create the data to be written
        data = zeros([1, self.nfreq, self.ny, self.nx], dtype=float)
        for inu in range(self.nfreq):
            data[0,inu,:,:] = self.image[:,:] * conv

        hdu     = pf.PrimaryHDU(data)
        hdulist = pf.HDUList([hdu])
        
        hdulist[0].header.update('CRPIX1', (self.nx+1.)/2., ' ')
        hdulist[0].header.update('CDELT1', self.sizepix_x/1.496e13/dpc/3600., '')
        hdulist[0].header.update('CRVAL1', self.sizepix_x/1.496e13/dpc*0.5+target_ra, '')
        hdulist[0].header.update('CUNIT1', '     DEG', '')
        hdulist[0].header.update('CTYPE1', 'RA---SIN', '')
       
        hdulist[0].header.update('CRPIX2', (self.ny+1.)/2., '')
        hdulist[0].header.update('CDELT2', self.sizepix_y/1.496e13/dpc/3600., '')
        hdulist[0].header.update('CRVAL2', self.sizepix_y/1.496e13/dpc*0.5+target_dec, '')
        hdulist[0].header.update('CUNIT2', '     DEG', '')
        hdulist[0].header.update('CTYPE2', 'DEC--SIN', '')
     
        hdulist[0].header.update('CRPIX3', 1., '')
        hdulist[0].header.update('CDELT3', 1., '')
        hdulist[0].header.update('CRVAL3', 1., '')
        hdulist[0].header.update('CUNIT3', '        ','')
        hdulist[0].header.update('CTYPE3', 'STOKES  ','')

        if self.nwav==1:
            hdulist[0].header.update('CRPIX4', 1.0, '')
            hdulist[0].header.update('CDELT4', bandwidthmhz*1e6, '')
            hdulist[0].header.update('CRVAL4', self.freq[0], '')
            hdulist[0].header.update('CUNIT4', '      HZ', '')
            hdulist[0].header.update('CTYPE4', 'FREQ-LSR', '')
        else:
            hdulist[0].header.update('CRPIX4', 1.0, '')
            hdulist[0].header.update('CDELT4', (self.freq[1]-self.freq[0]), '')
            hdulist[0].header.update('CRVAL4', self.freq[0], '')
            hdulist[0].header.update('CTYPE4', '      HZ', '')
            hdulist[0].header.update('CTYPE4', 'FREQ-LSR', '')
        
        hdulist[0].header.update('BUNIT', 'JY/PIXEL', '')
        hdulist[0].header.update('BTYPE', 'INTENSITY', '')
        hdulist[0].header.update('BZERO', 0.0, '')
        hdulist[0].header.update('BSCALE', 1.0, '')
        
        hdulist[0].header.update('EPOCH', 2000.0, '')
        hdulist[0].header.update('LONPOLE', 180.0, '')

        if os.path.exists(fname):
            print fname+' already exists'
            dum = raw_input('Do you want to overwrite it (yes/no)?')
            if (dum.strip()[0]=='y')|(dum.strip()[0]=='Y'):
                os.remove(fname)
                hdu.writeto(fname)
            else:
                print 'No image has been written'
        else:
            hdu.writeto(fname)
# --------------------------------------------------------------------------------------------------
    def writefits(self, fname='', dpc=1.):
        """
        Function to write out a RADMC3D image data in fits format  
  
        INPUT:
        ------
         fname   : file name of the radmc3d output image (if omitted 'image.fits' is used)
 
        """
# --------------------------------------------------------------------------------------------------
        if fname=='':
            fname = 'image.fits'
        pc = 3.0857200e+18

        # Conversion from erg/s/cm/cm/ster to Jy/pixel
        conv = self.sizepix_x * self.sizepix_y / (dpc * pc)**2. * 1e23
       
        hdu     = pf.PrimaryHDU(self.image*conv)
        hdulist = pf.HDUList([hdu])
        
        hdulist[0].header.update('CRPIX1', '%f'%((self.nx+1.)/2.), ' ')
        hdulist[0].header.update('CDELT1', '%f'%(self.sizepix_x/1.496e13/dpc), '')
        hdulist[0].header.update('CRVAL1', '%f'%(self.sizepix_x/1.496e13/dpc*0.5), '')
        hdulist[0].header.update('CUNIT1', '  ARCSEC', '')
       
        hdulist[0].header.update('CRPIX2', '%f'%((self.ny+1.)/2.), '')
        hdulist[0].header.update('CDELT2', '%f'%(self.sizepix_y/1.496e13/dpc), '')
        hdulist[0].header.update('CRVAL2', '%f'%(self.sizepix_y/1.496e13/dpc*0.5), '')
        hdulist[0].header.update('CUNIT2', '  ARCSEC', '')
      
        if self.nwav==1:
            hdulist[0].header.update('CRPIX3', '%f'%1, '')
            hdulist[0].header.update('CDELT3', '%f'%0, '')
            hdulist[0].header.update('CRVAL3', '%f'%self.wav, '')
            hdulist[0].header.update('CUNIT3', '  MICRON', '')
        else:
            hdulist[0].header.update('CRPIX3', '%f'%1, '')
            hdulist[0].header.update('CDELT3', '%f'%(self.wav[1]-self.wav[0]), '')
            hdulist[0].header.update('CRVAL3', '%f'%self.wav[0], '')
            hdulist[0].header.update('CUNIT3', '  MICRON', '')
        
        hdulist[0].header.update('BUNIT', 'JY/PIXEL', '')
        hdulist[0].header.update('BTYPE', 'INTENSITY', '')

        if os.path.exists(fname):
            print fname+' already exists'
            dum = raw_input('Do you want to overwrite it (yes/no)?')
            if (dum.strip()[0]=='y')|(dum.strip()[0]=='Y'):
                os.remove(fname)
                hdu.writeto(fname)
            else:
                print 'No image has been written'
        else:
            hdu.writeto(fname)
# --------------------------------------------------------------------------------------------------
    def plot_momentmap(self, moment=0, nu0=0, wav0=0, zmapnorm=False, dpc=1., au=False, arcsec=False, cmap=None, vclip=None):
        """
        Function to plot moment maps

        INPUT:
        ------
            moment : moment of the channel maps to be calculated 
            nu0    : rest frequency of the line in Hz
            wav0   : rest wavelength of the line in micron
            dpc    : distance of the source in pc

        OUTPUT:
        -------
            map : Numpy array with the same dimension as the individual channel maps
        """

        # I/O error handling
        if nu0==0:
            if wav0==0:
                print 'ERROR'
                print 'Neither rest frequency (nu0) nor rest wavelength (wav0) of the line is specified'
                return
            else:
                nu0 = 2.99792458e10/wav0*1e4
        

        if len(self.image.shape)!=3:
            print 'ERROR'
            print ' Channel map calculation requires a three dimensional array with Nx * Ny * Nnu dimensions'
            print ' The current image array contains '+str(len(self.image.shape))+' dimensions'
            return

        mmap = self.momentmap(moment=moment, nu0=nu0, wav0=wav0)

        if moment>0:
            if zmapnorm:
                mmap0 = self.momentmap(moment=0, nu0=nu0, wav0=wav0)
                mmap = mmap / mmap0


# Select the coordinates of the data
        if au:
            x = self.x/1.496e13
            y = self.y/1.496e13
            xlab = 'X [AU]'
            ylab = 'Y [AU]'
        elif arcsec:
            x = self.x/1.496e13/dpc
            y = self.y/1.496e13/dpc
            xlab = 'RA offset ["]'
            ylab = 'DEC offset ["]'
        else:
            x = self.x
            y = self.y
            xlab = 'X [cm]'
            ylab = 'Y [cm]'

        ext = (x[0], x[self.nx-1], y[0], y[self.ny-1])

        if moment==0:
            mmap = mmap / (dpc*dpc) 
            cb_label = 'I'+r'$_\nu$'+' [erg/s/cm/cm/Hz/ster*km/s]'
        if moment==1:
            mmap = mmap / (dpc*dpc) 
            if not zmapnorm:
                cb_label = 'I'+r'$_\nu$'+' [erg/s/cm/cm/Hz/ster'+r'$\rm (km/s)^2$'+']'
            else:
                cb_label = 'v [km/s]'
        if moment>1:
            mmap = mmap / (dpc*dpc) 
            powex = str(moment)
            if not zmapnorm:
                cb_label = 'I'+r'$_\nu$'+' [erg/s/cm/cm/Hz/ster'+r'$(km/s)^'+powex+'$'+']'
            else:
                cb_label = r'v$^'+powex+'$ [(km/s)$^'+powex+'$]'

        close()

        if vclip!=None:
            if len(vclip)!=2:
                print 'ERROR'
                print 'vclip should be a two element list with (clipmin, clipmax)'
                return
            else:
                mmap = mmap.clip(vclip[0], vclip[1])
        implot = imshow(mmap, extent=ext, cmap=cmap)
        cbar = colorbar(implot)
        cbar.set_label(cb_label)
        xlabel(xlab)
        ylabel(ylab)
# --------------------------------------------------------------------------------------------------
    def momentmap(self, moment=0, nu0=0, wav0=0):
        """
        Function to calculate moment maps

        INPUT:
        ------
            moment : moment of the channel maps to be calculated 
            nu0    : rest frequency of the line in Hz
            wav0   : rest wavelength of the line in micron

        OUTPUT:
        -------
            map : Numpy array with the same dimension as the individual channel maps
        """

        # I/O error handling
        if nu0==0:
            if wav0==0:
                print 'ERROR'
                print 'Neither rest frequency (nu0) nor rest wavelength (wav0) of the line is specified'
                return
            else:
                nu0 = 2.99792458e10/wav0*1e4
        

        if len(self.image.shape)!=3:
            print 'ERROR'
            print ' Channel map calculation requires a three dimensional array with Nx * Ny * Nnu dimensions'
            print ' The current image array contains '+str(len(self.image.shape))+' dimensions'
            return
        
        # First calculate the velocity field
        v = 2.99792458e10*(nu0-self.freq)/nu0/1e5


        vmap = zeros([self.nx, self.ny, self.nfreq], dtype=float64)
        for ifreq in range(self.nfreq):
            vmap[:,:,ifreq] = v[ifreq]

        # Now calculate the moment map
        y = self.image * (vmap**moment)


        dum = (vmap[:,:,1:] - vmap[:,:,:-1]) * (y[:,:,1:] + y[:,:,:-1]) * 0.5

        return dum.sum(2)

# --------------------------------------------------------------------------------------------------
    def readimage(self, fname=None):
        """
        Function to read an image calculated by RADMC3D 
     
        INPUT:
        ------
         fname   : file name of the radmc3d output image (if omitted 'image.out' is used)
 
        """
# --------------------------------------------------------------------------------------------------
        pc   = 3.08572e18

# Look for the image file

        if (fname==None): 
            fname = 'image.out'

        try:
            rfile = open(fname, 'r')
        except:
            print 'ERROR!'
            print 'No '+fname+' file has been found!'
            return -1
    
        dum = ''
    
# Format number
        dum = rfile.readline()

# Nr of pixels
        dum = rfile.readline()
        dum = dum.split()
        self.nx  = int(dum[0])
        self.ny  = int(dum[1])
# Nr of frequencies
        self.nfreq = int(rfile.readline())
        self.nwav  = self.nfreq 
# Pixel sizes
        dum = rfile.readline()
        dum = dum.split()
        self.sizepix_x = float(dum[0])
        self.sizepix_y = float(dum[1])
# Wavelength of the image
        self.wav = []
        for iwav in range(self.nwav):
            self.wav.append(float(rfile.readline()))
        self.wav = array(self.wav)
        self.freq = 2.99792458e10 / self.wav * 1e4
    
        
        self.image = zeros([self.nx,self.ny,self.nwav], dtype=float64)
        for iwav in range(self.nwav):
# Blank line
            dum = rfile.readline()
            for ix in range(self.nx):
                for iy in range(self.ny):
                    self.image[ix,iy,iwav] = float(rfile.readline())
       
        self.image = squeeze(self.image)
        rfile.close()
    
# Conversion from erg/s/cm/cm/Hz/ster to Jy/pixel
    
        conv  = self.sizepix_x * self.sizepix_y / pc**2. * 1e23 
        self.imageJyppix = self.image * conv

        self.x = ((arange(self.nx, dtype=float64) + 0.5) - self.nx/2) * self.sizepix_x
        self.y = ((arange(self.ny, dtype=float64) + 0.5) - self.ny/2) * self.sizepix_y

# --------------------------------------------------------------------------------------------------
    def get_vis(self, bl=None, pa=None, dpc=None, int_type='nearest'):
        """
        Function to calculate visibilities 
   
        SYNTAX:
        -------
            vis = get_visibility(image=image, wav=[10.0], bl=[15.0], pa=[40.], dpc=103.)
 
        INPUT:
        ------
            image    : a radmc3dImage class returned by readimage
            bl       : list of projected baselines in meter
            pa       : list of position angles (start from north counterclockwise)
            dpc      : distance of the source in parsec
            int_type : interpolation type 'spline' or 'nearest'

            NOTE!!!! bl and pa should have the same number of elements!

        OUTPUT:
        -------
            result          : a radmc3dVisibility class 
            result.fftImage : Fourier-transform of the image (shifted, but not normalized!)
            result.u        : spatial frequency calculated from the image axes
            result.v        : spatial frequency calculated from the image axes
            result.bl       : projected baseline in meter
            result.pa       : position angle of the projected baseline
            result.blu      : spatial frequency corresponding to the projected baselines
            result.blv      : spatial frequency corresponding to the projected baselines
            result.mu       : spatial frequency at which the fourier transform was taken (nearest neighbour to result.blu)
            result.mv       : spatial frequency at which the fourier transform was taken (nearest neighbour to result.blv)
            result.mu_err   : relative error of in the u coordinate ((result.mu-result.blu)/result.blu)
            result.mv_err   : relative error of in the u coordinate ((result.mv-result.blv)/result.blv)
            result.wav      : wavelength of the image/visibility
        """
# --------------------------------------------------------------------------------------------------

#
# Safety check
#
        if (len(bl) != len(pa)):
            print 'ERROR'
            print ' Number of baselines and number of position angles differ!'
            return -1   

#
# Natural constants
#
        au = 1.496e13

#
# Get the dimensions
#
        nbl  = len(bl)
        npa  = len(pa)
           

# Create the radmc3dVisibility that will be returned at the end

        res = radmc3dVisibility()

        res.bl     = bl
        res.pa     = pa
        res.wav    = self.wav
        res.cvis   = zeros([self.nwav, nbl], dtype=complex128)
        res.vis    = zeros([self.nwav, nbl], dtype=float64)
        res.amp    = zeros([self.nwav, nbl], dtype=float64)
        res.phase  = zeros([self.nwav, nbl], dtype=float64)
        res.mu     = zeros(nbl, dtype=float64)
        res.mv     = zeros(nbl, dtype=float64)
        res.mu_err = zeros([self.nwav, nbl], dtype=float64)
        res.mv_err = zeros([self.nwav, nbl], dtype=float64)
        res.blu    = zeros(nbl, dtype=float64)
        res.blv    = zeros(nbl, dtype=float64)
#
# Now do two big loops to get the visibilites for each baseline-PA pair and each wavelength point
#
        for ilam in range(self.nwav):
            image  = rot90(rot90(self.image))
# Calculate the ra, dec axes
            r_ra  = (arange(self.nx, dtype=float64) - (self.nx*0.5) + 0.5) * (-self.sizepix_x / au / dpc / 3600.)
            r_dec = (arange(self.ny, dtype=float64) - (self.ny*0.5) + 0.5) * (self.sizepix_y / au / dpc / 3600.)
# Calculate the Fourier transform of the image

            f_image = fft.fft2(image)
#            f_image = fft.rfft2(image)
#
# Shift the image along both axes and normalize it to u,v=(0,0)
#  NOTE: I make the following assumption fft(image)[0,0] = max(fft(image)) 
#
            res.fftImage = fft.fftshift(f_image)
            nsf_image = res.fftImage / self.image.sum()#res.fftImage.max() 
            nsf_image_real = real(nsf_image)
            nsf_image_imag = imag(nsf_image)
#
# Calculate the u-v points 
#

            res.u = (arange(self.nx, dtype=float64) - self.nx/2. + 1.) / \
                (-(self.sizepix_x / au / dpc / 3600.) /180.*pi * self.nx)
        
            res.v = (arange(self.ny, dtype=float64) - self.ny/2. + 1.) / \
                ((self.sizepix_x / au / dpc / 3600.) /180.*pi * self.ny)
        

            if int_type.strip().lower()=='spline':

                #sp_real = bvspline(res.v, res.u[::-1], nsf_image_real[:,::-1])#,sx=0,sy=0)
                #sp_imag = bvspline(res.v, res.u[::-1], nsf_image_imag[:,::-1])#,sx=0,sy=0)
                
                sp_real = rbs(res.v, res.u[::-1], nsf_image_real[:,::-1],kx=2,ky=2,s=0)
                sp_imag = rbs(res.v, res.u[::-1], nsf_image_imag[:,::-1],kx=2,ky=2,s=0)
                
                for ibl in range(nbl):
#
# Calculate the u,v, coordinates from the basline length and position angle
#

                    res.blu[ibl] = bl[ibl] * 1e6 / self.wav[ilam] * cos(pa[ibl]/180.*pi + pi/2.)
                    res.blv[ibl] = bl[ibl] * 1e6 / self.wav[ilam] * sin(pa[ibl]/180.*pi + pi/2.)

#
# Do spline interpolation
#

                    dum1 = sp_real(res.blv[ibl], res.blu[ibl])
                    dum2 = sp_imag(res.blv[ibl], res.blu[ibl])
                    
                    res.cvis[ilam,ibl]   = complex(dum1, dum2)
                    res.vis[ilam,ibl]    = abs(res.cvis[ilam,ibl])
                    res.amp[ilam,ibl]    = real(sqrt(res.cvis[ilam,ibl] * conj(res.cvis[ilam,ibl])))
                    res.phase[ilam,ibl]  = arccos(real(res.cvis[ilam,ibl]) / res.amp[ilam,ibl])
                    if imag(res.cvis[ilam,ibl])<0.0 : 
                        res.phase[ilam,ibl] = 2.0*pi -res.phase[ilam,ibl]
                    res.mu[ibl]          = res.blu[ibl]
                    res.mv[ibl]          = res.blv[ibl]
                    res.mu_err[ilam, ibl] = (res.mu[ibl]-res.blu[ibl])/res.blu[ibl]
                    res.mv_err[ilam, ibl] = (res.mu[ibl]-res.blv[ibl])/res.blv[ibl]


            elif int_type.strip().lower()=='nearest':
                
                for ibl in range(nbl):
#
# Calculate the u,v, coordinates from the basline length and position angle
#

                    res.blu[ibl] = bl[ibl] * 1e6 / self.wav[ilam] * cos(pa[ibl]/180.*pi + pi/2.)
                    res.blv[ibl] = bl[ibl] * 1e6 / self.wav[ilam] * sin(pa[ibl]/180.*pi + pi/2.)

#
# Take the nearest neighbour of the points (m_u, m_v) in the fourier transform 
#

                    diff = abs(res.u - res.blu[ibl])
                    ii   = (diff==diff.min())
                    
                
                
                    diff = abs(res.v - res.blv[ibl])
                    jj   = (diff==diff.min())

                    print (res.u[1]-res.u[0])/1e6*self.wav[ilam], res.blu[ibl]/1e6*self.wav[ilam]

                    res.cvis[ilam,ibl]   = nsf_image[ii,jj][0]
                    res.vis[ilam,ibl]    = abs(nsf_image[ii,jj])
                    res.amp[ilam,ibl]    = real(sqrt(nsf_image[ii,jj] * conj(nsf_image[ii,jj])))
                    res.phase[ilam,ibl]  = arccos(real(res.cvis[ilam,ibl]) / res.amp[ilam,ibl])
                    if imag(res.cvis[ilam,ibl])<0.0 : 
                        res.phase[ilam,ibl] = 2.0*pi -res.phase[ilam,ibl]
                    res.mu[ibl]          = res.u[ii][0]
                    res.mv[ibl]          = res.v[jj][0]
                    res.mu_err[ilam, ibl] = (res.mu[ibl]-res.blu[ibl])/res.blu[ibl]
                    res.mv_err[ilam, ibl] = (res.mv[ibl]-res.blv[ibl])/res.blv[ibl]


        um = res.u/1e6*self.wav[0]
        print um.min(), um.max()

        clf()
        imshow(abs(imag(nsf_image)), extent=(um.min(), um.max(), um.min(), um.max()), \
                interpolation='nearest')

        dum = raw_input()

        plot(res.blu/1e6*self.wav[0], res.blv/1e6*self.wav[0], 'ko')
        #dum = raw_input()
        #exit()

        return res
# --------------------------------------------------------------------------------------------------
    def imconv(self, fwhm=None, pa=None, dpc=None):
        """
        Function to convolve a radmc3d image with a two dimensional Gaussian psf 
    
        INPUT:
        ------
              fwhm    : A list of two numbers; the FWHM of the two dimensional psf along the two principal axes
                            The unit is assumed to be arcsec if dpc keyword is set, otherwise the unit is pixel
              pa      : Position angle of the psf ellipse (counts from North counterclockwise)
              dpc     : Distance of the source in pc, if omitted the unit of FWHM is assumed to be pixel
    
        OUTPUT:
        -------
              result  : same  
              'cimage': The convolved image with the psf (unit is erg/s/cm/cm/Hz/ster)
              'image' : The original unconvolved image (unit is erg/s/cm/cm/Hz/ster)
              'psf'   : Two dimensional psf
              'x'     : first coordinate axis of the psf/image
              'y'     : second coordinate axis of the psf/image
        """
# --------------------------------------------------------------------------------------------------
# Natural constants    
        au = 1.496e13
    
        imag = rot90(rot90(self.image))
        nx = self.nx
        ny = self.ny
        dx = self.sizepix_x / au / dpc
        dy = self.sizepix_y/au/dpc
        nfreq = self.nfreq
    
    
# Calculate the Gaussian psf
        dum   = get_psf(nx=self.nx, ny=self.ny, fwhm=fwhm, pa=pa, pscale=[dx,dy])
        psf   = dum['psf']
        f_psf = fft.fft2(psf)
        cimage = zeros([self.nx,self.ny,self.nfreq], dtype=float64)
        for ifreq in range(nfreq):
            f_imag  = fft.fft2(imag)
            f_cimag = f_psf * f_imag
            cimage[:,:,ifreq] = abs(fft.ifftshift(fft.ifft2(f_cimag)))
            
        cimage = squeeze(cimage)
  
# Return the convolved image (copy the image class and replace the image attribute to the convolved image)

        res = deepcopy(self)
        res.image = cimage
        res.psf   = psf
        res.fwhm  = fwhm
        res.pa    = pa
        res.dpc   = dpc


        return res

# --------------------------------------------------------------------------------------------------
def get_psf(nx=None, ny=None, fwhm=None, pa=None, pscale=None):
    """
    Function to generate a two dimensional Gaussian PSF
    
    INPUT:
    ------
          nx      : image size in the first dimension
          ny      : image size in the second dimension
          fwhm    : full width at half maximum of the psf in each dimension [fwhm_x, fwhm_y]
          pa      : position angle of the gaussian if the gaussian is not symmetric
          pscale  : pixelscale of the image, if set fwhm should be in the same unit, if not set unit of fwhm is pixels

    OUTPUT:
    -------
          result  : dictionary containing the following keys
          'psf'   : two dimensional numpy array containing the normalized psf
          'x'     : first coordinate axis of the psf
          'y'     : seonc coordinate axis of the psf
          
    """
# --------------------------------------------------------------------------------------------------

# Create the two axes

    if (pscale!=None):
        dx,dy = pscale[0], pscale[1]
    else:
        dx,dy = 1., 1.

    x = ((arange(nx, dtype=float64) + 0.5) - nx/2) * dx
    y = ((arange(ny, dtype=float64) + 0.5) - ny/2) * dy

# Calculate the standard deviation of the Gaussians
    sigmax = fwhm[0] / (2.0 * sqrt(2.0 * log(2.)))
    sigmay = fwhm[1] / (2.0 * sqrt(2.0 * log(2.)))
    norm   = 1./(2. * pi * sigmax * sigmay)


# Pre-compute sin and cos angles

    sin_pa = sin(pa/180.*pi - pi/2.)
    cos_pa = cos(pa/180.*pi - pi/2.)

# Define the psf
    psf = zeros([nx,ny], dtype=float64)
    for ix in range(nx):
        for iy in range(ny):
            xx = cos_pa * x[ix] - sin_pa * y[iy]
            yy = sin_pa * x[ix] + cos_pa * y[iy]

            psf[ix,iy] = exp(-0.5*xx*xx/sigmax/sigmax - 0.5*yy*yy/sigmay/sigmay)

    psf = psf / psf.sum()
    res = {'psf':psf, 'x':x, 'y':y}

    return res

# --------------------------------------------------------------------------------------------------
def readimage(fname=None):
    """
    Function to read an image calculated by RADMC3D 
     
    INPUT:
    ------
        fname   : file name of the radmc3d output image (if omitted 'image.out' is used)
 
    """

    dum = radmc3dImage()
    dum.readimage(fname=fname)
    return dum

# ***************************************************************************************************************
def plotimage(image=None, arcsec=False, au=False, log=False, dpc=None, maxlog=None, saturate=None, bunit=None, \
                  ifreq=None, cmap=None, cmask_rad=None, interpolation='none'):
    """
    Function to plot a radmc3d image
    
    SYNTAX:
    -------
          result = plotimage(image='image.out', arcsec=True, au=False, log=True, dpc=140, maxlog=-6., 
                             saturate=0.1, bunit='Jy')

    INPUT:
    ------
          image    : A radmc3dImage class returned by readimage   
          arcsec   : If True image axis will have the unit arcsec (NOTE: dpc keyword should also be set!)
          au       : If True image axis will have the unit AU
          log      : If True image scale will be logarithmic, otherwise linear
          dpc      : Distance to the source in parsec (This keywords should be set if arcsec=True, or bunit!=None)
          maxlog   : Logarithm of the lowest pixel value to be plotted, lower pixel values will be clippde
          saturate : Highest pixel values to be plotted in terms of the peak value, higher pixel values will be clipped
          bunit    : Unit of the image, can be 'Jy' = [Jy/pixel], 'abs'-[erg/s/cm/cm/Hz/ster], 'None'-I_nu/max(I_nu)
          ifreq    : If the image file/array consists of multiple frequencies/wavelengths ifreq denotes the index
                     of the frequency/wavelength in the image array to be plotted
          cmask_rad : Simulates coronographyic mask : sets the image values to zero within this radius of the image center
                      The unit is the same as the image axis (au, arcsec, cm)
                      NOTE: this works only on the plot, the image array is not changed (for that used the cmask() function)
    """
# ***************************************************************************************************************

# Natural constants
    pc   = 3.08572e18

# Check whether or not we need to mask the image

    if cmask_rad!=None:
        dum_image = cmask(image, rad=cmask_rad, au=au, arcsec=arcsec, dpc=dpc) 
    else:
        dum_image = image

    if (image.nfreq>1):
        if (ifreq==None):
            ifreq = 0
        data = squeeze(dum_image.image[:,:,ifreq])
        data = rot90(rot90(data))
    else:
        data = rot90(rot90(dum_image.image))

    norm  = data.max()
    if (bunit==None):
        data = data/norm

    clipnorm = data.max()
# Check if the data should be plotted on a log scale
    if log:
        clipmin = log10(data[data>0.].min())
        data = log10(data.clip(1e-90))
        
# Clipping the data
        if (maxlog!=None):
            clipmin = -maxlog + log10(clipnorm)
    else:
        clipmin  = data.min()

    if (saturate!=None):
        if (saturate>1.): 
            saturate = 1.0
        if log:
            clipmax = log10(saturate) + log10(clipnorm)
        else:
            clipmax = clipnorm * saturate
    else:
        clipmax = clipnorm

    data = data.clip(clipmin, clipmax)

# Select the unit of the data

    if (bunit==None):
        if log:
            cb_label = 'log(F'+r'$_\nu$'+'/max(F'+r'$_\nu$'+'))'
        else:
            cb_label = 'F'+r'$_\nu$'+'/max(F'+r'$_\nu$'+')'
    elif (bunit=='abs'):
        if log:
            data    = data - log10(dpc*dpc) 
            cb_label = 'log(F'+r'$_\nu$'+' [erg/s/cm/cm/Hz])'
        else:
            data    = data / (dpc*dpc) 
            cb_label = 'F'+r'$_\nu$'+' [erg/s/cm/cm/Hz/]'

    elif (bunit=='Jy'):
        if dpc==None:
            print 'ERROR'
            print ' If Jy/pixel is selected for the image unit the dpc keyword should also be set'
            return
        else:
            if log:
                data    = data + log10(image.sizepix_x * image.sizepix_y / (dpc*pc)**2. * 1e23) 
                cb_label = 'log(I'+r'$_\nu$'+ '[Jy/pixel])'
            else:
                data    = data * (image.sizepix_x * image.sizepix_y / (dpc*pc)**2. * 1e23) 
                cb_label = 'I'+r'$_\nu$'+' [Jy/pixel]'

# Set the color bar boundaries
    if log:
        cb_bound = (data.max(), data.min())
    else:
        cb_bound = (data.min(), data.max())

# Select the coordinates of the data
    if au:
        x = image.x/1.496e13
        y = image.y/1.496e13
        xlab = 'X [AU]'
        ylab = 'Y [AU]'
    elif arcsec:
        x = image.x/1.496e13/dpc
        y = image.y/1.496e13/dpc
        xlab = 'RA offset ["]'
        ylab = 'DEC offset ["]'
    else:
        x = image.x
        y = image.y
        xlab = 'X [cm]'
        ylab = 'Y [cm]'

    ext = (x[0], x[image.nx-1], y[0], y[image.ny-1])


# Now finally put everything together and plot the data
    delaxes()
    delaxes()

    if (cmap==None): 
        cmap = cm.gist_gray
#    implot = imshow(data, extent=ext, cmap=cm.gist_gray)
    implot = imshow(data, extent=ext, cmap=cmap)
    xlabel(xlab)
    ylabel(ylab)
    title(r'$\lambda$='+str(image.wav)+r'$\mu$m')
    cbar = colorbar(implot)
    cbar.set_label(cb_label)
    show()
# ***************************************************************************************************************

def makeimage(npix=None, incl=None, wav=None, sizeau=None, phi=None, posang=None, pointau=None, \
                  fluxcons=True, nostar=False, noscat=False):
    """
    Function to call RADMC3D to calculate a rectangular image
    
    SYNTAX:
    -------
           makeimage(npix=100, incl=60.0, wav=10.0, sizeau=300., phi=0., posang=15., 
                     pointau=[0., 0.,0.], fluxcons=True, nostar=False, noscat=False)
           
    INPUT:
    ------
           npix   : number of pixels on the rectangular images
           sizeau : diameter of the image in au
           incl   : inclination angle of the source
           dpc    : distance of the source in parsec
           phi    : azimuthal rotation angle of the source in the model space
           posang : position angle of the source in the image plane
           pointau: three elements list of the cartesian coordinates of the image center
    
    KEYWORDS:
    ---------
           fluxcons : this should not even be a keyword argument, it ensures flux conservation 
           (adaptive subpixeling) in the rectangular images
           nostar   : if True the calculated images will not contain stellar emission
           noscat   : if True, scattered emission will be neglected in the source function, however, 
                          extinction will contain scattering if kappa_scat is not zero.  
    """
# **************************************************************************************************
# 
# The basic keywords that should be set
#
    if (npix==None):
        print 'ERROR!'
        print ' npix keyword is not set'
        return -1
    if (incl==None):
        print 'ERROR!'
        print ' incl keyword is not set'
        return -1
    if (wav==None):
        print 'ERROR!'
        print ' wav/lambda keyword is not set'
        return -1
    if (sizeau==None):
        print 'ERROR!'
        print ' sizeau keyword is not set'
        return -1

    com = 'radmc3d image' 
    com = com + ' npix ' + str(npix)
    com = com + ' lambda ' + str(wav)
    com = com + ' incl ' + str(incl)
    com = com + ' sizeau ' + str(sizeau)
    
#
# Now add additional optional keywords/arguments
#
    if (phi!=None):
        com = com + ' phi ' + str(phi)

    if (posang!=None):
        com = com + ' posang ' + str(posang)

    if (pointau!=None):
        if (len(pointau)!=3):
            print 'ERROR'
            print ' pointau should be a list of 3 elements corresponding to the '
            print ' cartesian coordinates of the image center'
            return -1
        else:
            com = com + ' pointau ' + str(posang[0]) + ' ' + str(posang[1]) + ' ' + str(posang[2]) 
    else:
        com = com + ' pointau 0.0  0.0  0.0'

    if fluxcons:
        com = com + ' fluxcons'

#
# Write the RADMC3D control file
#
#    print 'Writing radmc3d.inp'
#    wfile = open('radmc3d.inp', 'w')
#    wfile.write('%s\n'%('nphot='+str('%d'%nphot)))
#    wfile.write('%s\n'%('nphot_scat='+str('%d'%nphot_scat)))
#    wfile.write('%s\n'%('scattering_mode_max='+str('%d'%scattering_mode_max)))
#    wfile.write('%s\n'%('tgas_eq_tdust='+str('%d'%tgas_eq_tdust)))
#    wfile.close()
#    if (nostar==False):
#        com = com + ' inclstar'
#
#    if noscat:
#        com = com + ' noscat'

#
# Now finally run radmc3d and calculate the image
#

    print com
    #dum = sp.Popen([com], stdout=sp.PIPE, shell=True).wait()
    dum = sp.Popen([com], shell=True).wait()

    return 0
# **************************************************************************************************

def get_visibility(image=None, bl=None, pa=None, dpc=None):
    """
    Function to calculate visibilities 

    SYNTAX:
    -------
        vis = get_visibility(image=image, wav=[10.0], bl=[15.0], pa=[40.], dpc=103.)
 
    INPUT:
    ------
        image  : a radmc3dImage class returned by readimage
        bl     : list of projected baselines in meter
        pa     : list of position angles (start from north counterclockwise)
        dpc    : distance of the source in parsec

        NOTE!!!! bl and pa should have the same number of elements!

    OUTPUT:
    -------
         result          : a radmc3dVisibility class 
         result.fftImage : Fourier-transform of the image (shifted, but not normalized!)
         result.u        : spatial frequency calculated from the image axes
         result.v        : spatial frequency calculated from the image axes
         result.bl       : projected baseline in meter
         result.pa       : position angle of the projected baseline
         result.blu      : spatial frequency corresponding to the projected baselines
         result.blv      : spatial frequency corresponding to the projected baselines
         result.mu       : spatial frequency at which the fourier transform was taken (nearest neighbour to result.blu)
         result.mv       : spatial frequency at which the fourier transform was taken (nearest neighbour to result.blv)
         result.mu_err   : relative error of in the u coordinate ((result.mu-result.blu)/result.blu)
         result.mv_err   : relative error of in the u coordinate ((result.mv-result.blv)/result.blv)
         result.wav      : wavelength of the image/visibility
    """
    res = image.get_vis(bl=bl, pa=pa, dpc=dpc)
    return res
           

def cmask(im=None, rad=0.0, au=False, arcsec=False, dpc=None):
    """
    Function to simulate a coronographic mask by
    setting the image values to zero within circle of a given radius around the
    image center

    INPUT:
    ------
        im     : a radmc3dImage class
        rad    : radius of the mask 
        au     : if true the radius is taken to have a unit of AU
        arcsec : if true the radius is taken to have a unit of arcsec (dpc
                  should also be set)
        dpc    : distance of the source (required if arcsec = True)

        NOTE: if arcsec=False and au=False rad is taken to have a unit of pixel

    OUTPUT:
    -------
        res    : a radmc3dImage class containing the masked image
    """

    if au:
        if arcsec:
            print 'ERROR'
            print ' Either au or arcsec should be set, but not both of them'
            return
       
        crad = rad*1.496e+13
    else:
        if arcsec:
            crad = rad * 1.496e13 * dpc
        else:
            crad = rad* im.sizepix_x

    res = deepcopy(im)
    if im.nfreq!=1:
        for ix in range(im.nx):
            r = sqrt(im.y**2 + im.x[ix]**2)
            ii = r<=crad
            res.image[ix,ii,:] = 0.0
    else:
        for ix in range(im.nx):
            r = sqrt(im.y**2 + im.x[ix]**2)
            ii = r<=crad
            res.image[ix,ii] = 0.0

    return res

