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
    print ' Scipy.interpolate.BivariateSpline or Scipy.interpolate.RectBivariateSpline cannot be imported'

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
class radmc3dImage():
    """
    RADMC3D image class

    ATTRIBUTES:
    -----------
        image       - The image as calculated by radmc3d (the values are intensities in erg/s/cm^2/Hz/ster)
        imageJyppix - The image with pixel units of Jy/pixel
        x           - x coordinate of the image [cm]
        y           - y coordinate of the image [cm]
        nx          - Number of pixels in the horizontal direction
        ny          - Number of pixels in the vertical direction
        sizepix_x   - Pixel size in the horizontal direction [cm]
        sizepix_y   - Pixel size in the vertical direction [cm]
        nfreq       - Number of frequencies in the image cube
        freq        - Frequency grid in the image cube
        nwav        - Number of wavelengths in the image cube (same as nfreq)
        wav         - Wavelength grid in the image cube

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
    def get_closure_phase(self, bl=None, pa=None, dpc=None):
        """
        Function to calculate clusure phases for a given model image for any arbitrary baseline triplet

        INPUT:
        ------
            bl       - a list or Numpy array containing the length of projected baselines in meter(!)
            pa       - a list or Numpy array containing the position angles of projected baselines in degree(!)
            dpc      - distance of the source in parsec

            NOTE, bl and pa should either be an array with dimension [N,3] or if they are lists each element of
                the list should be a list of length 3, since closure phases are calculated only for closed triangles

        OUTPUT:
        -------
            returns a dictionary with the following keys:
                
                bl     - projected baseline in meter
                pa     - position angle of the projected baseline in degree
                nbl    - number of baselines
                u      - spatial frequency along the x axis of the image
                v      - spatial frequency along the v axis of the image
                vis    - complex visibility at points (u,v)
                amp    - correlation amplitude 
                phase  - Fourier phase
                cp     - closure phase
                wav    - wavelength 
                nwav   - number of wavelengths

        """

        res = {}
        res['bl']    = array(bl, dtype=float)
        res['pa']    = array(pa, dtype=float)
        res['ntri']  = res['bl'].shape[0]
        res['nbl']   = 3
        res['nwav']  = self.nwav
        res['wav']   = self.wav

        res['ntri']  = res['bl'].shape[0]
        res['u']     = zeros([res['ntri'], 3, res['nwav']], dtype=float)
        res['v']     = zeros([res['ntri'], 3, res['nwav']], dtype=float)
        res['vis']   = zeros([res['ntri'], 3, res['nwav']], dtype=complex64)
        res['amp']   = zeros([res['ntri'], 3, res['nwav']], dtype=float)
        res['phase'] = zeros([res['ntri'], 3, res['nwav']], dtype=float)
        res['cp']    = zeros([res['ntri'], res['nwav']], dtype=float)
        
        
        l   = self.x / 1.496e13 / dpc / 3600. / 180.*pi 
        m   = self.y / 1.496e13 / dpc / 3600. / 180.*pi 
        dl  = l[1]-l[0]
        dm  = m[1]-m[0]
        
        for itri in range(res['ntri']):
            print 'Calculating baseline triangle # : ', itri
            
            dum = self.get_visibility(bl=res['bl'][itri,:], pa=res['pa'][itri,:], dpc=dpc)
            res['u'][itri,:,:] = dum['u']
            res['v'][itri,:,:] = dum['v']
            res['vis'][itri,:,:] = dum['vis']
            res['amp'][itri,:,:] = dum['amp']
            res['phase'][itri,:,:] = dum['phase']
            res['cp'][itri,:] = (dum['phase'].sum(0) / pi * 180.)%360.
            ii = res['cp'][itri,:]>180.
            if (res['cp'][itri,ii]).shape[0]>0:
                res['cp'][itri,ii] = res['cp'][itri,ii]-360.

        return res


# --------------------------------------------------------------------------------------------------
    def get_visibility(self, bl=None, pa=None, dpc=None):
        """
        Function to calculate visibilities for a given set of projected baselines and position angles
        with the Discrete Fourier Transform

        INPUT:
        ------
            bl       - a list or Numpy array containing the length of projected baselines in meter(!)
            pa       - a list or Numpy array containing the position angles of projected baselines in degree(!)
            dpc      - distance of the source in parsec

        OUTPUT:
        -------
            returns a dictionary with the following keys:
                
                bl     - projected baseline in meter
                pa     - position angle of the projected baseline in degree
                nbl    - number of baselines
                u      - spatial frequency along the x axis of the image
                v      - spatial frequency along the v axis of the image
                vis    - complex visibility at points (u,v)
                amp    - correlation amplitude 
                phase  - phase
                wav    - wavelength 
                nwav   - number of wavelengths

        """

        res   = {}
        res['bl']    = array(bl, dtype=float)
        res['pa']    = array(pa, dtype=float)
        res['nbl']   = res['bl'].shape[0]
        
        res['wav']   = array(self.wav)
        res['nwav']  = self.nwav

        res['u']     = zeros([res['nbl'], self.nwav], dtype=float)
        res['v']     = zeros([res['nbl'], self.nwav], dtype=float)

        res['vis']   = zeros([res['nbl'], self.nwav], dtype=complex64)
        res['amp']   = zeros([res['nbl'], self.nwav], dtype=float64)
        res['phase'] = zeros([res['nbl'], self.nwav], dtype=float64)


        l   = self.x / 1.496e13 / dpc / 3600. / 180.*pi 
        m   = self.y / 1.496e13 / dpc / 3600. / 180.*pi 
        dl  = l[1]-l[0]
        dm  = m[1]-m[0]

        for iwav in range(res['nwav']):
        
            # Calculate spatial frequencies 
            res['u'][:,iwav] = res['bl'] * cos(res['pa']) * 1e6 / self.wav[iwav]
            res['v'][:,iwav] = res['bl'] * sin(res['pa']) * 1e6 / self.wav[iwav]


            for ibl in range(res['nbl']):
                dum = complex(0.)
                imu = complex(0., 1.)

                dmv = zeros(self.image.shape[1])*0. + dm
                

                for il in range(len(l)):
                    phase = 2.*pi * (res['u'][ibl]*l[il] + res['v'][ibl]*m)
                    cterm = cos(phase)
                    sterm = -sin(phase)
                    dum   = dum + (self.image[il,:]*(cterm + imu*sterm)).sum()*dl*dm
               
                res['vis'][ibl, iwav]     = dum
                res['amp'][ibl, iwav]     = sqrt(abs(dum * conj(dum)))
                res['phase'][ibl, iwav]   = arccos(real(dum)/res['amp'][ibl])
                if imag(dum)<0.:
                    res['phase'][ibl, iwav] = 2.*pi - res['phase'][ibl, iwav]
                     
                print 'Calculating baseline # : ', ibl, ' wavelength # : ', iwav 

        return res
# --------------------------------------------------------------------------------------------------
    def writefits(self, fname='', dpc=1., coord='03h10m05s -10d05m30s', bandwidthmhz=2000.0, casa=False):
        """
        Function to write out a RADMC3D image data in fits format (CASA compatible)
  
        INPUT:
        ------
         fname        : File name of the radmc3d output image (if omitted 'image.fits' is used)
         coord        : Image center coordinates
         bandwidthmhz : Bandwidth of the image in MHz (equivalent of the CDELT keyword in the fits header)
         casa         : If set to True a Stokes axis will be added to the image cube to make it compatible
                        with the casa simulator
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

        if casa:
            data = zeros([self.nfreq, 1, self.ny, self.nx], dtype=float)
            if self.nfreq==1:
                data[0,0,:,:] = self.image[:,:] * conv

            else:
                for inu in range(self.nfreq):
                    data[inu,0,:,:] = self.image[:,:,inu] * conv
        else:
            data = zeros([self.nfreq, self.ny, self.nx], dtype=float)
            if self.nfreq==1:
                data[0,:,:] = self.image[:,:] * conv

            else:
                for inu in range(self.nfreq):
                    data[inu,:,:] = self.image[:,:,inu] * conv


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
    

        if casa:
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
        else:
            if self.nwav==1:
                hdulist[0].header.update('CRPIX3', 1.0, '')
                hdulist[0].header.update('CDELT3', bandwidthmhz*1e6, '')
                hdulist[0].header.update('CRVAL3', self.freq[0], '')
                hdulist[0].header.update('CUNIT3', '      HZ', '')
                hdulist[0].header.update('CTYPE3', 'FREQ-LSR', '')
            else:
                hdulist[0].header.update('CRPIX3', 1.0, '')
                hdulist[0].header.update('CDELT3', (self.freq[1]-self.freq[0]), '')
                hdulist[0].header.update('CRVAL3', self.freq[0], '')
                hdulist[0].header.update('CTYPE3', '      HZ', '')
                hdulist[0].header.update('CTYPE3', 'FREQ-LSR', '')


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
    def plot_momentmap(self, moment=0, nu0=0, wav0=0, dpc=1., au=False, arcsec=False, cmap=None, vclip=None):
        """
        Function to plot moment maps

        INPUT:
        ------
            moment : moment of the channel maps to be calculated 
            nu0    : rest frequency of the line in Hz
            wav0   : rest wavelength of the line in micron
            dpc    : distance of the source in pc
            au     : If true displays the image with AU as the spatial axis unit
            arcsec : If true displays the image with arcsec as the spatial axis unit (dpc should also be set!)
            cmap   : matplotlib colormap
            vclip  : two element list / Numpy array containin the lower and upper limits for the values in the moment
                      map to be displayed

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

        mmap = self.get_momentmap(moment=moment, nu0=nu0, wav0=wav0)

        if moment>0:
            mmap0 = self.get_momentmap(moment=0, nu0=nu0, wav0=wav0)
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
            cb_label = 'v [km/s]'
        if moment>1:
            mmap = mmap / (dpc*dpc) 
            powex = str(moment)
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
    def get_momentmap(self, moment=0, nu0=0, wav0=0):
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
        ifreq=None, cmask_rad=None, interpolation='nearest', cmap=cm.gist_gray, **kwargs):
                  #ifreq=None, cmap=None, cmask_rad=None, interpolation='nearest'):
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
          bunit    : Unit of the image, (None - Inu/max(Inu), 'inu' - Inu, fnu - Jy/pixel) 
          ifreq    : If the image file/array consists of multiple frequencies/wavelengths ifreq denotes the index
                     of the frequency/wavelength in the image array to be plotted
          cmask_rad : Simulates coronographyic mask : sets the image values to zero within this radius of the image center
                      The unit is the same as the image axis (au, arcsec, cm)
                      NOTE: this works only on the plot, the image array is not changed (for that used the cmask() function)

          cmap     : matplotlib color map
          interpolation: interpolation keyword for imshow (e.g. 'nearest', 'bilinear', 'bicubic')
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
        data = squeeze(dum_image.image[::-1,:,ifreq])
    else:
        data = dum_image.image[::-1,:] 

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
            cb_label = 'log(I'+r'$_\nu$'+'/max(I'+r'$_\nu$'+'))'
        else:
            cb_label = 'I'+r'$_\nu$'+'/max(I'+r'$_\nu$'+')'
    elif (bunit=='inu'):
        if log:
            cb_label = 'log(I'+r'$_\nu$'+' [erg/s/cm/cm/Hz/ster])'
        else:
            cb_label = 'I'+r'$_\nu$'+' [erg/s/cm/cm/Hz/ster]'
    elif (bunit=='fnu'):
        if dpc==None:
            print 'ERROR'
            print ' If Jy/pixel is selected for the image unit the dpc keyword should also be set'
            return
        else:
            if log:
                data    = data + log10(image.sizepix_x * image.sizepix_y / (dpc*pc)**2. * 1e23) 
                cb_label = 'log(F'+r'$_\nu$'+ '[Jy/pixel])'
            else:
                data    = data * (image.sizepix_x * image.sizepix_y / (dpc*pc)**2. * 1e23) 
                cb_label = 'F'+r'$_\nu$'+' [Jy/pixel]'

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

    #if (cmap==None): 
        #cmap = cm.gist_gray
#    implot = imshow(data, extent=ext, cmap=cm.gist_gray)
    implot = imshow(data, extent=ext, cmap=cmap, interpolation=interpolation, **kwargs)
    xlabel(xlab)
    ylabel(ylab)
    title(r'$\lambda$='+("%.5f"%image.wav[ifreq])+r'$\mu$m')
    cbar = colorbar(implot)
    cbar.set_label(cb_label)
    show()
# ***************************************************************************************************************

def makeimage(npix=None, incl=None, wav=None, sizeau=None, phi=None, posang=None, pointau=None, \
                  fluxcons=True, nostar=False, noscat=False, \
                  widthkms=None, linenlam=None, vkms=None, iline=None):
    """
    Function to call RADMC3D to calculate a rectangular image
    
    SYNTAX:
    -------
           makeimage(npix=100, incl=60.0, wav=10.0, sizeau=300., phi=0., posang=15., 
                     pointau=[0., 0.,0.], fluxcons=True, nostar=False, noscat=False)
           
    INPUT:
    ------
           npix    : number of pixels on the rectangular images
           sizeau  : diameter of the image in au
           incl    : inclination angle of the source
           dpc     : distance of the source in parsec
           phi     : azimuthal rotation angle of the source in the model space
           posang  : position angle of the source in the image plane
           pointau : three elements list of the cartesian coordinates of the image center
           widthkms: width of the frequency axis of the channel maps
           linenlam: number of wavelengths to calculate images at
           vkms    : a single velocity value at which a channel map should be calculated 
           iline   : line transition index
    
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

    if widthkms:
        if vkms:
            print ' ERROR '
            print ' Either widthkms or vkms keyword should be set but not both'
            return -1
        com = com + ' widthkms '+("%.5e"%widthkms)

    if vkms:
        if widthkms:
            print ' ERROR '
            print ' Either widthkms or vkms keyword should be set but not both'
            return -1
        com = com + ' vkms '+("%.5e"%vkms)

    if linenlam:
        com = com + ' linenlam '+("%d"%linenlam)

    if iline:
        com = com + ' iline '+("%d"%iline)


#
# Now finally run radmc3d and calculate the image
#

    #dum = sp.Popen([com], stdout=sp.PIPE, shell=True).wait()
    dum = sp.Popen([com], shell=True).wait()

    return 0
# **************************************************************************************************

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

