"""This module contains classes for handling dust opacities
"""
from __future__ import absolute_import
from __future__ import print_function
import traceback
import subprocess as sp
import os

try:
    import numpy as np
except ImportError:
    np = None
    print(' Numpy cannot be imported ')
    print(' To use the python module of RADMC-3D you need to install Numpy')
    print(traceback.format_exc())

try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None
    print('Warning')
    print('matplotlib.pyplot cannot be imported')
    print('Without matplotlib you can use the python module to set up a model but you will not be able to plot things')
    print('or display images')

from . import natconst as nc
from . import miescat
from . reggrid import *


class radmc3dDustOpac(object):
    """
    Class to handle dust opacities.


    Attributes
    ----------

    wav     : list
                Each element of the list contains an ndarray with the wavelength grid

    freq    : list
                Each element of the list contains an ndarray with the frequency grid

    nwav    : list
                Each element of the list contains an integer with the number of wavelengths

    kabs    : list
                Each element of the list contains an ndarray with the absorption coefficient per unit mass

    ksca    : list
                Each element of the list contains an ndarray with the scattering coefficient per unit mass

    phase_g : list
                Each element of the list contains an ndarray with the hase function

    ext     : list
                Each element of the list contains a string wht the file name extension of the duskappa_ext.Kappa file

    therm   : list
                Each element of the list contains a bool, if it is set to False the dust grains are quantum-heated
                (default: True)

    idust   : lisintt
                Each element of the list contains an integer with the index of the dust species in the dust density
                distribution array

    scatmat : list
                Each element is a boolean indicating whether the dust opacity table includes (True) the full scattering
                matrix or not (False)

    nang    : list
                Each element is a string, containing the number of scattering angles in the scattering matrix if its
                given

    scatang : list
                Each element is a numpy ndarray containing the scattering angles in the scattering matrix if its given

    z11     : list
                Each element is a numpy ndarray containing the (1,1) element of the scattering angles in the scattering
                matrix if its given

    z12     : list
                Each element is a numpy ndarray containing the (1,2) element of the scattering angles in the scattering
                matrix if its given

    z22     : list
                Each element is a numpy ndarray containing the (2,2) element of the scattering angles in the scattering
                matrix if its given

    z33     : list
                Each element is a numpy ndarray containing the (3,3) element of the scattering angles in the scattering
                matrix if its given

    z34     : list
                Each element is a numpy ndarray containing the (3,4) element of the scattering angles in the scattering
                matrix if its given

    z44     : list
                Each element is a numpy ndarray containing the (4,4) element of the scattering angles in the scattering
                matrix if its given

    """

    # --------------------------------------------------------------------------------------------------
    def __init__(self):

        self.wav = []
        self.freq = []
        self.nwav = []
        self.nfreq = []
        self.kabs = []
        self.ksca = []
        self.phase_g = []
        self.ext = []
        self.idust = []
        self.therm = []
        self.scatmat = []
        self.z11 = []
        self.z12 = []
        self.z22 = []
        self.z33 = []
        self.z34 = []
        self.z44 = []
        self.scatang = []
        self.nang = []

    # --------------------------------------------------------------------------------------------------
    def readOpac(self, ext=None, idust=None, scatmat=None, old=False):
        """Reads the dust opacity files.

        Parameters
        ----------

        ext  : list
                File name extension (file names should look like 'dustkappa_ext.inp')

        idust: list
                Indices of the dust species in the master opacity file (dustopac.inp') - starts at 0

        scatmat: list
                If specified, its elements should be booleans indicating whether the opacity file
                contains also the full scattering matrix (True) or only dust opacities (False)

        old   : bool, optional
                If set to True the file format of the previous, 2D version of radmc will be used
        """

        # Check the input keywords and if single strings are given convert them to lists
        # This assumes, though, that there is a single dust opacity file or dust species, though!!
        if ext is None:
            if idust is None:
                raise ValueError('Unknown ext and idust. '
                                 'File name extension must be given to be able to read the opacity from file.')
            else:
                if isinstance(idust, int):
                    idust = [idust]
        else:
            if isinstance(ext, str):
                ext = [ext]

            if (len(ext) == 1) & (ext[0] != ''):
                if idust is not None:
                    raise ValueError('Either idust or ext should be specified, but not both')

        if scatmat is None:
            # If the scatmat keyword is not given (i.e. if it is None) then assume that
            # it is False for all dust species
            scatmat = []
            if idust is None:
                for i in range(len(ext)):
                    scatmat.append(False)

            else:
                for i in range(len(idust)):
                    scatmat.append(False)
        else:
            if isinstance(scatmat, str):
                scatmat = [scatmat]

        # Find the file name extensions in the master opacity file if idust is specified instead of ext
        if idust:
            # Read the master dust opacity file to get the dust indices and dustkappa file name extensions
            mopac = self.readMasterOpac()

            ext = []
            for ispec in idust:
                if (ispec + 1) > len(mopac['ext']):
                    raise ValueError('No dust species found at index ' + ("%d" % ispec))
                else:
                    ext.append(mopac['ext'][ispec])

        # If only the extension is specified look for the master opacity file and find the index of this dust species
        #  or set the index to -1 if no such dust species is present in the master opacity file
        else:
            # # Read the master dust opacity file to get the dust indices and dustkappa file name extensions
            idust = [i for i in range(len(ext))]

        # Now read all dust opacities
        for i in range(len(ext)):
            if scatmat[i]:
                with open('dustkapscatmat_' + ext[i] + '.inp', 'r') as rfile:

                    print('Reading dustkapscatmat_' + ext[i] + '.inp ....')

                    self.ext.append(ext[i])

                    # Read the header/comment field
                    dum = rfile.readline()
                    while dum.strip()[0] == '#':
                        dum = rfile.readline()

                    # for j in range(6):
                    # dum = rfile.readline()

                    # Read the file format
                    iformat = int(dum)
                    # iformat = int(rfile.readline())
                    if iformat != 1:
                        rfile.close()
                        raise ValueError('Format number of the file dustkapscatmat_' + ext[i] + '.inp (iformat=' + (
                            "%d" % iformat) + ') is unkown')

                    # Read the number of wavelengths in the file
                    dum = int(rfile.readline())
                    self.nwav.append(dum)
                    self.nfreq.append(dum)
                    self.idust.append(idust[i])
                    idu = len(self.nwav) - 1

                    # Read the scattering angular grid
                    self.nang.append(int(rfile.readline()))
                    wav = np.zeros(self.nwav[idu], dtype=np.float64)
                    kabs = np.zeros(self.nwav[idu], dtype=np.float64)
                    ksca = np.zeros(self.nwav[idu], dtype=np.float64)
                    phase_g = np.zeros(self.nwav[idu], dtype=np.float64)
                    scatang = np.zeros(self.nang[idu], dtype=np.float64)
                    z11 = np.zeros([self.nwav[idu], self.nang[idu]], dtype=np.float64)
                    z12 = np.zeros([self.nwav[idu], self.nang[idu]], dtype=np.float64)
                    z22 = np.zeros([self.nwav[idu], self.nang[idu]], dtype=np.float64)
                    z33 = np.zeros([self.nwav[idu], self.nang[idu]], dtype=np.float64)
                    z34 = np.zeros([self.nwav[idu], self.nang[idu]], dtype=np.float64)
                    z44 = np.zeros([self.nwav[idu], self.nang[idu]], dtype=np.float64)

                    print('Reading the opacities..')
                    dum = rfile.readline()
                    for ilam in range(self.nwav[idu]):
                        dum = rfile.readline().split()
                        wav[ilam] = float(dum[0])
                        kabs[ilam] = float(dum[1])
                        ksca[ilam] = float(dum[2])
                        phase_g[ilam] = float(dum[3])

                    print('Reading the angular grid..')
                    dum = rfile.readline()
                    for iang in range(self.nang[idu]):
                        dum = rfile.readline()
                        scatang[iang] = float(dum)

                    print('Reading the scattering matrix..')
                    for ilam in range(self.nwav[idu]):
                        dum = rfile.readline()
                        for iang in range(self.nang[idu]):
                            dum = rfile.readline().split()
                            z11[ilam, iang] = float(dum[0])
                            z12[ilam, iang] = float(dum[1])
                            z22[ilam, iang] = float(dum[2])
                            z33[ilam, iang] = float(dum[3])
                            z34[ilam, iang] = float(dum[4])
                            z44[ilam, iang] = float(dum[5])

                    self.wav.append(wav)
                    self.freq.append(nc.cc / wav * 1e4)
                    self.kabs.append(kabs)
                    self.ksca.append(ksca)
                    self.phase_g.append(phase_g)
                    self.scatang.append(scatang)
                    self.z11.append(z11)
                    self.z12.append(z12)
                    self.z22.append(z22)
                    self.z33.append(z33)
                    self.z34.append(z34)
                    self.z44.append(z44)

            else:
                if not old:
                    fname = 'dustkappa_' + ext[i] + '.inp'
                    with open(fname, 'r') as rfile:

                        self.ext.append(ext[i])

                        # Read the file format
                        iformat = int(rfile.readline())
                        if (iformat < 1) | (iformat > 3):
                            rfile.close()
                            raise ValueError('Unknown file format in the dust opacity file ' + fname)

                        # Read the number of wavelengths in the file
                        dum = rfile.readline()
                        self.nwav.append(int(dum))
                        self.nfreq.append(int(dum))
                        self.idust.append(idust[i])
                        idu = len(self.nwav) - 1

                        # If only the absorption coefficients are specified
                        if iformat == 1:
                            wav = np.zeros(self.nwav[idu], dtype=np.float64)
                            kabs = np.zeros(self.nwav[idu], dtype=np.float64)
                            for ilam in range(self.nwav[idu]):
                                dum = rfile.readline().split()
                                wav[ilam] = float(dum[0])
                                kabs[ilam] = float(dum[1])
                            self.wav.append(wav)
                            self.freq.append(nc.cc / wav * 1e4)
                            self.kabs.append(kabs)
                            self.ksca.append([-1])
                            self.phase_g.append([-1])
                        # If the absorption and scattering coefficients are specified
                        elif iformat == 2:
                            wav = np.zeros(self.nwav[idu], dtype=np.float64)
                            kabs = np.zeros(self.nwav[idu], dtype=np.float64)
                            ksca = np.zeros(self.nwav[idu], dtype=np.float64)
                            for ilam in range(self.nwav[idu]):
                                dum = rfile.readline().split()
                                wav[ilam] = float(dum[0])
                                kabs[ilam] = float(dum[1])
                                ksca[ilam] = float(dum[2])
                            self.wav.append(wav)
                            self.freq.append(nc.cc / wav * 1e4)
                            self.kabs.append(kabs)
                            self.ksca.append(ksca)
                            self.phase_g.append([-1])

                        # If the absorption and scattering coefficients and also the scattering phase
                        # function are specified
                        elif iformat == 3:
                            wav = np.zeros(self.nwav[idu], dtype=np.float64)
                            kabs = np.zeros(self.nwav[idu], dtype=np.float64)
                            ksca = np.zeros(self.nwav[idu], dtype=np.float64)
                            phase_g = np.zeros(self.nwav[idu], dtype=np.float64)
                            for ilam in range(self.nwav[idu]):
                                dum = rfile.readline().split()
                                wav[ilam] = float(dum[0])
                                kabs[ilam] = float(dum[1])
                                ksca[ilam] = float(dum[2])
                                phase_g[ilam] = float(dum[3])

                            self.wav.append(wav)
                            self.freq.append(nc.cc / wav * 1e4)
                            self.kabs.append(kabs)
                            self.ksca.append(ksca)
                            self.phase_g.append(phase_g)

                else:
                    fname = 'dustopac_' + ext[i] + '.inp'
                    with open(fname, 'r') as rfile:
                        print('Reading ' + fname)
                        freq = np.fromfile('frequency.inp', count=-1, sep="\n", dtype=float)
                        nfreq = int(freq[0])
                        freq = freq[1:]

                        self.ext.append(ext[i])
                        dum = rfile.readline().split()
                        if int(dum[0]) != nfreq:
                            rfile.close()
                            raise ValueError(fname + ' contains a different number of frequencies than frequency.inp')

                        wav = nc.cc / freq * 1e4
                        kabs = np.zeros(nfreq, dtype=float)
                        ksca = np.zeros(nfreq, dtype=float)

                        dum = rfile.readline()
                        for ilam in range(nfreq):
                            kabs[ilam] = float(rfile.readline())
                        dum = rfile.readline()
                        for ilam in range(nfreq):
                            ksca[ilam] = float(rfile.readline())

                    self.wav.append(wav[::-1])
                    self.freq.append(freq[::-1])
                    self.kabs.append(kabs[::-1])
                    self.ksca.append(ksca[::-1])
                    self.phase_g.append([-1])

        return 0

    def makeOpac(self, ppar=None, wav=None, old=False, code='python',
                 theta=None, logawidth=None, wfact=3.0, na=20, chopforward=0., errtol=0.01,
                 verbose=False, extrapolate=False):
        """Createst the dust opacities using a Mie code distributed with RADMC-3D.

        Parameters
        ----------

        ppar        : dictionary
                      Parameters of the simulations

        wav         : ndarray, optional
                      Wavelength grid on which the mass absorption coefficients should be calculated

        code        : {'python', 'fortran'}
                      Version of the mie scattering code BHMIE to be used. 'fortran' - use the original fortran77
                      code of Bruce Drain (should be downloaded separately, compiled and its path added to the PATH
                      environment variable), 'python' a python version of BHMIE by Kees Dullemond (radmc3dPy.miescat).

        theta       : ndarray, optional
                      Angular grid (a numpy array) between 0 and 180
                      which are the scattering angle sampling points at
                      which the scattering phase function is computed.

        logawidth   : float, optional
                     If set, the size agrain will instead be a
                     sample of sizes around agrain. This helps to smooth out
                     the strong wiggles in the phase function and opacity
                     of spheres at an exact size. Since in Nature it rarely
                     happens that grains all have exactly the same size, this
                     is quite natural. The value of logawidth sets the width
                     of the Gauss in ln(agrain), so for logawidth<<1 this
                     give a real width of logawidth*agraincm.

        wfact       : float
                      Grid width of na sampling points in units
                      of logawidth. The Gauss distribution of grain sizes is
                      cut off at agrain * exp(wfact*logawidth) and
                      agrain * exp(-wfact*logawidth). Default = 3


        na          : int
                      Number of size sampling points (if logawidth set, default=20)

        chopforward : float
                      If >0 this gives the angle (in degrees from forward)
                      within which the scattering phase function should be
                      kept constant, essentially removing the strongly peaked
                      forward scattering. This is useful for large grains
                      (large ratio 2*pi*agraincm/lamcm) where the forward
                      scattering peak is extremely strong, yet extremely
                      narrow. If we are not interested in very forward-peaked
                      scattering (e.g. only relevant when modeling e.g. the
                      halo around the moon on a cold winter night), this will
                      remove this component and allow a lower angular grid
                      resolution for the theta grid.


        errtol      : float
                      Tolerance of the relative difference between kscat
                      and the integral over the zscat Z11 element over angle.
                      If this tolerance is exceeded, a warning is given.

        verbose     : bool
                      If set to True, the code will give some feedback so
                      that one knows what it is doing if it becomes slow.

        extrapolate : bool
                      If set to True, then if the wavelength grid lamcm goes
                      out of the range of the wavelength grid of the
                      optical constants file, then it will make a suitable
                      extrapolation: keeping the optical constants constant
                      for lamcm < minimum, and extrapolating log-log for
                      lamcm > maximum.


        old         : bool, optional
                      If set to True the file format of the previous, 2D version of radmc will be used
        """

        #
        # Create the wavelength grid if it is not specified
        #
        if wav is None:
            grid = radmc3dGrid()
            grid.makeWavelengthGrid(ppar=ppar)
            wav = grid.wav

            #
            # Do we need to mix the opacities?
            #
        if ppar is None:
            raise ValueError('Unknown ppar. The parameter dictionary is required to get the lnk file names.')

        if isinstance(ppar['lnk_fname'], str):
            ppar['lnk_fname'] = [ppar['lnk_fname']]

        if len(ppar['lnk_fname']) > 1:
            ext = []
            for idust in range(len(ppar['lnk_fname'])):

                # makedust needs the lnk file to be sorted in wavelength so create a dummy file
                # which contains the sorted optical constants
                with open(ppar['lnk_fname'][idust], 'r') as rfile:
                    w = []
                    n = []
                    k = []
                    dum = rfile.readline()
                    while len(dum) > 0:
                        dum = dum.split()
                        w.append(dum[0])
                        n.append(dum[1])
                        k.append(dum[2])
                        dum = rfile.readline()

                w = np.array(w, dtype=float)
                n = np.array(n, dtype=float)
                k = np.array(k, dtype=float)

                if float(w[0]) > float(w[w.shape[0] - 1]):
                    w = w[::-1]
                    n = n[::-1]
                    k = k[::-1]

                # Write out the dummy file containing the sorted optical constants
                with open('opt_const.dat', 'w') as wfile:
                    for iwav in range(w.shape[0]):
                        wfile.write("%s %s %s \n" % (w[iwav], n[iwav], k[iwav]))

                if code.lower().strip() == 'fortran':
                    # Run makedust
                    self.runMakedust(freq=nc.cc / wav * 1e4, gmin=ppar['gsmin'], gmax=ppar['gsmax'], ngs=ppar['ngs'],
                                     lnk_fname='opt_const.dat', gdens=ppar['gdens'][idust])

                    # Change the name of makedust's output
                    for igs in range(ppar['ngs']):
                        dum = sp.Popen('mv dustkappa_' + str(igs + 1) + '.inp dustkappa_idust_' + str(idust + 1)
                                       + '_igsize_' + str(igs + 1) + '.inp', shell=True).wait()
                        ext.append('idust_' + str(idust + 1) + '_igsize_' + str(igs + 1))

                elif code.lower().strip() == 'python':

                    if 'nscatang' in ppar:
                        nang = ppar['nscatang']
                    else:
                        nang = 180
                    theta = 180. * np.arange(nang, dtype=np.float) / np.float(nang - 1)

                    if 'logawidth' in ppar:
                        logawidth = ppar['logawidth']
                    else:
                        logawidth = None

                    if 'wfact' in ppar:
                        wfact = ppar['wfact']
                    else:
                        wfact = 3.0

                    if 'chopforward' in ppar:
                        if ppar['chopforward'] > 0.:
                            chopforward = ppar['chopforward']
                        else:
                            chopforward = None
                    else:
                        chopforward = 0.0

                    if 'errtol' in ppar:
                        errtol = ppar['errtol']
                    else:
                        errtol = 0.01

                    if 'miescat_verbose' in ppar:
                        verbose = ppar['miescat_verbose']
                    else:
                        verbose = False

                    if 'extrapolate' in ppar:
                        extrapolate = ppar['extrapolate']
                    else:
                        extrapolate = False

                    # Get the grain sizes in micrometer
                    gsize = ppar['gsmin'] + (ppar['gsmax'] / ppar['gsmin'])**(
                        np.arange(ppar['ngs'], dtype=np.float64) / (float() - 1.))

                    for igs in range(ppar['ngs']):
                        o = miescat.compute_opac_mie(fname='opt_const.dat', matdens=ppar['gdens'][idust],
                                                     agraincm=gsize[igs] * 1e-4, lamcm=wav * 1e-4, theta=theta,
                                                     logawidth=logawidth, wfact=wfact, na=na, chopforward=chopforward,
                                                     errtol=errtol, verbose=verbose, extrapolate=extrapolate)

                        if ppar['scattering_mode_max'] <= 2:
                            miescat.write_radmc3d_kappa_file(package=o, name='idust_1_igsize_' + str(igs + 1))
                        else:
                            miescat.write_radmc3d_scatmat_file(package=o, name='idust_1_igsize_' + str(igs + 1))

                os.remove('opt_const.dat')

            # Mix the opacity of different dust species for a given grain size if mixing is requested
            if 'mixabun' in ppar:
                if len(ppar['mixabun']) == len(ppar['lnk_fname']):
                    ext = []
                    for igs in range(ppar['ngs']):
                        mixnames = ['dustkappa_igsize_' + str(igs + 1) + '.inp']
                        mixspecs = [['dustkappa_idust_' + str(idust + 1) + '_igsize_' + str(igs + 1) + '.inp'
                                     for idust in range(len(ppar['lnk_fname']))]]
                        self.mixOpac(mixnames=mixnames, mixspecs=mixspecs, mixabun=[ppar['mixabun']])
                        ext.append('igsize_' + str(igs + 1))
                else:
                    raise ValueError('ppar["mixabun"] or ppar["lnk_fname"] has the wrong shape. They both should have '
                                     + 'the same number of elements, but the number of elements are different.')

            therm = [True for i in range(len(ext))]
            self.writeMasterOpac(ext=ext, therm=therm, scattering_mode_max=ppar['scattering_mode_max'], old=old)

            if old:
                self.makeopacRadmc2D(ext=ext)

        else:
            # makedust needs the lnk file to be sorted in wavelength so create a dummy file
            # which contains the sorted optical constants
            with open(ppar['lnk_fname'][0], 'r') as rfile:
                w = []
                n = []
                k = []
                dum = rfile.readline()
                while len(dum) > 0:
                    dum = dum.split()
                    w.append(dum[0])
                    n.append(dum[1])
                    k.append(dum[2])
                    dum = rfile.readline()

            w = np.array(w, dtype=float)
            n = np.array(n, dtype=float)
            k = np.array(k, dtype=float)

            if float(w[0]) > float(w[w.shape[0] - 1]):
                w = w[::-1]
                n = n[::-1]
                k = k[::-1]

            # Write out the dummy file containing the sorted optical constants
            with open('opt_const.dat', 'w') as wfile:
                for iwav in range(w.shape[0]):
                    wfile.write("%s %s %s \n" % (w[iwav], n[iwav], k[iwav]))

            if code.lower().strip() == 'fortran':
                # Run makedust
                self.runMakedust(freq=nc.cc / wav * 1e4, gmin=ppar['gsmin'], gmax=ppar['gsmax'], ngs=ppar['ngs'],
                                 lnk_fname='opt_const.dat', gdens=ppar['gdens'][0])

                # Change the name of makedust's output
                ext = []
                therm = []
                for igs in range(ppar['ngs']):
                    dum = sp.Popen('mv dustkappa_' + str(igs + 1) + '.inp dustkappa_idust_1_igsize_' + str(igs + 1)
                                   + '.inp', shell=True).wait()
                    ext.append('idust_1_igsize_' + str(igs + 1))
                    therm.append(True)

            elif code.lower().strip() == 'python':

                if 'nscatang' in ppar:
                    nang = ppar['nscatang']
                else:
                    nang = 180
                theta = 180. * np.arange(nang, dtype=np.float) / np.float(nang - 1)

                if 'logawidth' in ppar:
                    logawidth = ppar['logawidth']
                else:
                    logawidth = None

                if 'wfact' in ppar:
                    wfact = ppar['wfact']
                else:
                    wfact = 3.0

                if 'chopforward' in ppar:
                    if ppar['chopforward'] > 0.:
                        chopforward = ppar['chopforward']
                    else:
                        chopforward = None
                else:
                    chopforward = 0.0

                if 'errtol' in ppar:
                    errtol = ppar['errtol']
                else:
                    errtol = 0.01

                if 'miescat_verbose' in ppar:
                    verbose = ppar['miescat_verbose']
                else:
                    verbose = False

                if 'extrapolate' in ppar:
                    extrapolate = ppar['extrapolate']
                else:
                    extrapolate = False

                o = miescat.compute_opac_mie(fname='opt_const.dat', matdens=ppar['gdens'][0],
                                             agraincm=ppar['gsmin'] * 1e-4, lamcm=wav * 1e-4, theta=theta,
                                             logawidth=logawidth, wfact=wfact, na=na, chopforward=chopforward,
                                             errtol=errtol, verbose=verbose, extrapolate=extrapolate)

                if ppar['scattering_mode_max'] <= 2:
                    miescat.write_radmc3d_kappa_file(package=o, name='idust_1_igsize_1')
                else:
                    miescat.write_radmc3d_scatmat_file(package=o, name='idust_1_igsize_1')

                therm = [True]
                ext = ['idust_1_igsize_1']

            else:
                msg = 'Unknown mie scattering code version ' + code
                raise ValueError(msg)

            self.writeMasterOpac(ext=ext, therm=therm, scattering_mode_max=ppar['scattering_mode_max'], old=old)
            if old:
                self.makeopacRadmc2D(ext=ext)

        # Clean up and remove dust.inp and frequency.inp
        os.remove('dust.inp')
        if not old:
            os.remove('frequency.inp')

    @staticmethod
    def mixOpac(ppar=None, mixnames=None, mixspecs=None, mixabun=None, writefile=True):
        """Mixes dust opacities.


        Parameters
        -----------
        ppar      : dictionary, optional
                    All parameters of the actual model setup.

        mixnames  : list, optional
                    Names of the files into which the mixed dust opacities will be written
                    (not needed if writefile=False)

        mixspecs  : list, optional
                    Names of the files from which the dust opacities are read (not needed if readfile=False)

        mixabun   : list, optional
                    Abundances of different dust species

        writefile : bool
                    If False the mixed opacities will not be written out to files given in mixnames.

        NOTE, either ppar or  mixname, mixspecs, and mixabun should be set.

        """

        if writefile:
            if mixnames is None:
                if ppar is None:
                    raise ValueError('Neither ppar nor mixnames are set in mixOpac')
                else:
                    mixnames = ppar['mixnames']

        if mixspecs is None:
            if ppar is None:
                raise ValueError(' Neither ppar nor mixspecs are set in mixOpac ')
            else:
                mixspecs = ppar['mixspecs']

        if mixabun is None:
            if ppar is None:
                raise ValueError(' Neither ppar nor mixabun are set in mixOpac ')
            else:
                mixabun = ppar['mixabun']

        for i in range(len(mixnames)):
            #
            # Read the dust opacities to be mixed for composite dust species #1
            #
            ocabs = []
            ocsca = []
            ogsym = []
            oform = 0
            for j in range(len(mixspecs[i])):
                with open(mixspecs[i][j], 'r') as rfile:
                    form = int(rfile.readline())
                    nwav = int(rfile.readline())
                    dw = np.zeros(nwav, dtype=float)
                    dcabs = np.zeros(nwav, dtype=float)
                    dcsca = np.zeros(nwav, dtype=float)
                    gsym = np.zeros(nwav, dtype=float)
                    if form == 1:
                        if (oform == 0) | (oform == 1):
                            oform = 1
                        else:
                            print(' ')
                            print('WARNING')
                            print(' You are trying to mix opacity tables with different formats. Some of the tables \n'
                                  + ' contain scattering coefficients while (format>=2) while others do not '
                                  + ' (format=1).\n'
                                  + ' If you wish to continue mixing will only be done for the absorption and the \n'
                                  + 'output opacity table will have a format number of 1.')
                            dum = input('Do you wish to continue (1-yes, 0-no) ?')
                            if dum.strip() != '1':
                                return

                        for iwav in range(nwav):
                            dum = rfile.readline().split()
                            dw[iwav], dcabs[iwav] = float(dum[0]), float(dum[1])
                    if form == 2:
                        if (oform == 0) | (oform == 2):
                            oform = 2
                        else:
                            print(' ')
                            print('WARNING')
                            print(' You are trying to mix opacity tables with different formats. Some of the tables \n'
                                  + ' contain scattering coefficients while (format>=2) while other do not '
                                  + '(format=1). \n'
                                  + ' If you wish to continue mixing will only be done for the absorption and the \n'
                                  + 'output opacity table will have a format number of 1.')

                            dum = input('Do you wish to continue (1-yes, 0-no) ?')
                            if dum.strip() != '1':
                                return
                        for iwav in range(nwav):
                            dum = rfile.readline().split()
                            dw[iwav], dcabs[iwav], dcsca[iwav] = float(dum[0]), float(dum[1]), float(dum[2])
                    if form == 3:
                        if (oform == 0) | (oform == 3):
                            oform = 3
                        else:
                            print(' ')
                            print('WARNING')
                            print(' You are trying to mix opacity tables with different formats. Some of the tables \n'
                                  + ' contain scattering coefficients while (format>=2) while other do not '
                                  + '(format=1) \n'
                                  + ' If you wish to continue mixing will only be done for the absorption and the '
                                  + 'output opacity table will have a format number of 1.')
                            dum = input('Do you wish to continue (1-yes, 0-no) ?')
                            if dum.strip() != '1':
                                return
                        for iwav in range(nwav):
                            dum = rfile.readline().split()
                            dw[iwav], dcabs[iwav], dcsca[iwav], gsym[iwav] = float(dum[0]), float(dum[1]), float(
                                dum[2]), float(dum[3])
                    if form > 3:
                        raise ValueError(' Unsupported dust opacity table format (format number: ' + ("%d" % form) + ')'
                                         + ' Currently only format number 1 and 2 are supported')

                    if dw[1] < dw[0]:
                        print(' Dust opacity table seems to be sorted in frequency instead of wavelength')
                        print(' Reversing the arrays')
                        dw = dw[::-1]
                        dcabs = dcabs[::-1]
                        dcsca = dcsca[::-1]

                if j == 0:
                    ocabs = np.array(dcabs) * mixabun[i][j]
                    ocsca = np.array(dcsca) * mixabun[i][j]
                    ogsym = np.array(gsym) * mixabun[i][j]
                    nwav0 = dw.shape[0]
                    owav = np.array(dw)
                else:
                    #
                    # Interpolate dust opacities to the wavelength grid of the first dust species
                    #
                    ii = ((owav >= dw[0]) & (owav <= dw[nwav - 1]))
                    il = (owav < dw[0])
                    ih = (owav > dw[nwav - 1])
                    dum = np.zeros(nwav0, dtype=float)
                    dum[ii] = 10. ** np.interp(np.log10(owav[ii]), np.log10(dw), np.log10(dcabs))

                    # Edwtrapolate the absorption coefficients using linear fit in log-log space
                    # (i.e. fitting a polinomial) for short wavelengths
                    # der = np.log10(dcabs[1] / dcabs[0]) / np.log10(dw[1] / dw[0])
                    dum[il] = 10. ** (np.log10(dcabs[0]) + np.log10(dw[0] / owav[il]))

                    # Edwtrapolate the absorption coefficients using linear fit in log-log space
                    # (i.e. fitting a polinomial) for long wavelengths
                    # der = np.log10(dcabs[nwav - 1] / dcabs[nwav - 2]) / np.log10(dw[nwav - 1] / dw[nwav - 2])
                    dum[ih] = 10. ** (np.log10(dcabs[nwav - 1]) + np.log10(owav[il] / dw[nwav - 1]))

                    ocabs = ocabs + np.array(dum) * mixabun[i][j]

                    if oform == 2:
                        # Do the inter-/extrapolation of for the scattering coefficients
                        dum = np.zeros(nwav0, dtype=float)
                        dum[ii] = 10. ** np.interp(np.log10(owav[ii]), np.log10(dw), np.log10(dcsca))

                        # der = np.log10(dcsca[1] / dcsca[0]) / np.log10(dw[1] / dw[0])
                        dum[il] = 10. ** (np.log10(dcsca[0]) + np.log10(dw[0] / owav[il]))

                        # der = np.log10(dcsca[nwav - 1] / dcsca[nwav - 2]) / np.log10(dw[nwav - 1] / dw[nwav - 2])
                        dum[ih] = 10. ** (np.log10(dcsca[nwav - 1]) + np.log10(owav[il] / dw[nwav - 1]))

                        ocsca = ocsca + np.array(dum) * mixabun[i][j]

                    if oform == 3:
                        # Do the inter-/extrapolation of for the scattering phase function
                        dum = np.zeros(nwav0, dtype=float)
                        dum[ii] = 10. ** np.interp(np.log10(owav[ii]), np.log10(dw), np.log10(gsym))

                        # der = np.log10(gsym[1] / gsym[0]) / np.log10(dw[1] / dw[0])
                        dum[il] = 10. ** (np.log10(gsym[0]) + np.log10(dw[0] / owav[il]))

                        # der = np.log10(gsym[nwav - 1] / gsym[nwav - 2]) / np.log10(dw[nwav - 1] / dw[nwav - 2])
                        dum[ih] = 10. ** (np.log10(gsym[nwav - 1]) + np.log10(owav[il] / dw[nwav - 1]))

                        ogsym = ogsym + np.array(dum) * mixabun[i][j]

            #
            # Write out the mixed dust opacities
            #
            with open(mixnames[i], 'w') as wfile:
                wfile.write("%d\n" % oform)
                wfile.write("%d\n" % owav.shape[0])
                if oform == 1:
                    for iwav in range(owav.shape[0]):
                        wfile.write("%.9e %.9e\n" % (owav[iwav], ocabs[iwav]))
                if oform == 2:
                    for iwav in range(owav.shape[0]):
                        wfile.write("%.9e %.9e %.9e\n" % (owav[iwav], ocabs[iwav], ocsca[iwav]))
                if oform == 3:
                    for iwav in range(owav.shape[0]):
                        wfile.write("%.9e %.9e %.9e %.9e\n" % (owav[iwav], ocabs[iwav], ocsca[iwav], ogsym[iwav]))

        return

    @staticmethod
    def readMasterOpac():
        """Reads the master opacity file 'dustopac.inp'.
        It reads the dustkappa filename extensions (dustkappa_ext.inp) corresponding to dust species indices

        Returns
        -------

        Returns a dictionary with the following keys:

            *ext   : list of dustkappa file name extensions

            *therm : a list of integers specifying whether the dust grain is thermal or quantum heated
            (0 - thermal, 1 - quantum heated)
        """

        with open('dustopac.inp', 'r') as rfile:

            # file format
            dum = rfile.readline()
            # nr of dust species
            ndust = int(rfile.readline().split()[0])
            # Comment line
            dum = rfile.readline()

            ext = []
            therm = []
            scatmat = []
            for idust in range(ndust):
                # Check if we have dust opacities also for the full scattering matrix
                dum = rfile.readline().split()
                if int(dum[0]) == 1:
                    scatmat.append(False)
                elif int(dum[0]) == 10:
                    scatmat.append(True)

                # Check if the dust grain is thermal or quantum heated
                dum = int(rfile.readline().split()[0])
                if dum == 0:
                    therm.append(True)
                else:
                    therm.append(False)
                # Dustkappa filename extension
                dum = rfile.readline().split()[0]
                ext.append(dum)
                # Comment line
                dum = rfile.readline()

        return {'ext': ext, 'therm': therm, 'scatmat': scatmat}

    @staticmethod
    def writeMasterOpac(ext=None, therm=None, scattering_mode_max=1, old=False):
        """Writes the master opacity file 'dustopac.inp'.

        Parameters
        ----------

        ext                 : list
                              List of dustkappa file name extensions

        therm               : list
                              List of integers specifying whether the dust grain is thermal or quantum heated
                              (0-thermal, 1-quantum)

        scattering_mode_max : int
                              Scattering mode code in radmc3d : 0 - no scattering, 1 - isotropic scattering,
                              2 - anisotropic scattering with Henyei-Greenstein phase function, 5 - anisotropic
                              scattering using the full scattering matrix and stokes vectors.

        old                 : bool, optional
                              If set to True the file format of the previous, 2D version of radmc will be used
        """

        print('Writing dustopac.inp')

        if not ext:
            raise ValueError('Unknown ext. '
                             + ' No file name extension is specified. Without it dustopac.inp cannot be written')
        else:
            if isinstance(ext, str):
                ext = [ext]

        if therm is None:
            # If therm is not specified it is assumed that all grains are thermal, no quantum heating
            therm = [True for i in range(len(ext))]
        else:
            if isinstance(therm, int):
                therm = [therm]
            if len(ext) != len(therm):
                raise ValueError(' The number of dust species in ext and in therm are different')

        with open('dustopac.inp', 'w') as wfile:

            # File format
            wfile.write('%-15s %s\n' % ('2', 'Format number of this file'))
            # Number of dust species
            wfile.write('%-15s %s\n' % (str(len(ext)), 'Nr of dust species'))
            # Separator
            wfile.write('%s\n' % '============================================================================')

            if not old:
                for idust in range(len(ext)):
                    # Dust opacity will be read from a file
                    if scattering_mode_max < 5:
                        wfile.write('%-15s %s\n' % ('1', 'Way in which this dust species is read'))
                    else:
                        wfile.write('%-15s %s\n' % ('10', 'Way in which this dust species is read'))

                    # Check if the dust grain is thermal or quantum heated
                    if therm:
                        if therm[idust]:
                            wfile.write('%-15s %s\n' % ('0', '0=Thermal grain, 1=Quantum heated'))
                    else:
                        wfile.write('%-15s %s\n' % ('1', '0=Thermal grain, 1=Quantum heated'))

                    # Dustkappa filename extension
                    wfile.write('%s %s %s\n' % (ext[idust], '    ', 'Extension of name of dustkappa_***.inp file'))
                    # Separator
                    wfile.write('%s\n' % '----------------------------------------------------------------------------')
            else:
                for idust in range(len(ext)):
                    # Dust opacity will be read from a file
                    wfile.write('%-15s %s\n' % ('-1', 'Way in which this dust species is read (-1=File)'))

                    # Check if the dust grain is thermal or quantum heated
                    wfile.write('%-15s %s\n' % ('0', '0=Thermal grain, 1=Quantum heated'))
                    # Dustkappa filename extension
                    wfile.write('%d %s %s\n' % ((idust + 1), '    ', 'Extension of name of dustopac_***.inp file'))
                    # Separator
                    wfile.write('%s\n' % '----------------------------------------------------------------------------')

    def makeopacRadmc2D(self, ext=None):
        """
        Creates dust opacities (dustopac_*.inp files) for the previous 2D version of radmc
        It takes the input dust opacity files and interpolates them onto the used frequency grid

        Parameters
        ----------

            ext : list
                  List of dustkappa file name extensions, i.e. the input file name has to be named
                  as dustkappa_ext[i].inp

        """

        if ext is None:
            raise ValueError('Unknown ext. Dust opacity file name extensions are mandatory.')

        else:
            if isinstance(ext, str):
                ext = [ext]

        self.readOpac(ext=ext, old=False)
        #
        # Read the frequency.inp file
        #
        freq = np.fromfile('frequency.inp', count=-1, sep="\n", dtype=float)
        nfreq = int(freq[0])
        freq = freq[1:]
        freq = freq[::-1]
        wav = nc.cc / freq * 1e4
        #
        # Check if the frequency grid is ordered in frequency or in wavelength
        #
        worder = False
        if freq[-1] < freq[0]:
            worder = True

        for i in range(len(ext)):
            kabs = np.zeros(nfreq, dtype=float)
            ksca = np.zeros(nfreq, dtype=float)
            ish = (wav < self.wav[i][0])
            ilo = (wav > self.wav[i][-1])
            ii = ((wav >= self.wav[i][0]) & (wav <= self.wav[i][-1]))

            #
            # Do logarithmic interpolation for the overlapping wavelenght domain
            #
            kabs[ii] = 10. ** np.interp(np.log10(wav[ii]), np.log10(self.wav[i]), np.log10(self.kabs[i]))
            if len(self.ksca[i]) > 1:
                ksca[ii] = 10. ** np.interp(np.log10(wav[ii]), np.log10(self.wav[i]), np.log10(self.ksca[i]))

            #
            # Do the long wavelength part
            #
            if True in ilo:
                x1 = np.log10(self.wav[i][-1])
                x0 = np.log10(self.wav[i][-2])

                y1 = np.log10(self.kabs[i][-1])
                y0 = np.log10(self.kabs[i][-2])
                der = (y1 - y0) / (x1 - x0)
                kabs[ilo] = 10. ** (y1 + der * (np.log10(wav[ilo]) - x1))

                y1 = np.log10(self.ksca[i][-1])
                y0 = np.log10(self.ksca[i][-2])
                der = (y1 - y0) / (x1 - x0)
                ksca[ilo] = 10. ** (y1 + der * (np.log10(wav[ilo]) - x1))

            #
            # Do the shorter wavelength
            #
            if True in ish:
                kabs[ish] = self.kabs[0][0]
                ksca[ish] = self.ksca[0][0]

            #
            # Now write the results to file
            #
            fname = 'dustopac_' + ("%d" % (i + 1)) + '.inp'
            with open(fname, 'w') as wfile:
                print('Writing ' + fname)
                wfile.write("%d 1\n" % nfreq)
                wfile.write(" \n")
                #
                # Reverse the order of kabs,ksca as they are ordered in frequency in radmc
                #
                if worder:
                    x = kabs[::-1]
                else:
                    x = kabs
                for ilam in range(nfreq):
                    wfile.write("%.7e\n" % x[ilam])

                wfile.write(" \n")
                if worder:
                    x = ksca[::-1]
                else:
                    x = ksca
                for ilam in range(nfreq):
                    wfile.write("%.7e\n" % x[ilam])

    @staticmethod
    def runMakedust(freq=None, gmin=None, gmax=None, ngs=None, lnk_fname=None, gdens=None):
        """Interface function to the F77 code makedust to calculate mass absorption coefficients.

        Parameters
        ----------
        freq       : ndarray
                    Contains the frequency grid on which the opacities should be calculated

        gmin       : float
                    Minimum grain size

        gmax       : float
                    Maximum grain size

        ngs        : int
                    Number of grain sizes

        gdens      : float
                    Density of the dust grain in g/cm^3

        lnk_fname  : str
                    Name of the file in which the optical constants are stored

        Returns
        -------

        Returns an ndarray with [nfreq,ngs] dimensions containing the resulting opacities
        """

        #
        # Calculate the grain sizes
        #
        if ngs > 1:
            gsize = gmin * (gmax / gmin) ** (np.arange(ngs, dtype=np.float64) / (float(ngs) - 1.))
        else:
            gsize = [gmin]

        #
        # Write the frequency.inp file
        #
        with open('frequency.inp', 'w') as wfile:
            wfile.write("%d\n" % freq.shape[0])
            wfile.write("  \n")
            for i in range(freq.shape[0]):
                wfile.write("%.10e\n" % freq[i])

        #
        # Write the dust.inp file (makedust main control file)
        #
        with open('dust.inp', 'w') as wfile:
            for igs in range(ngs):
                wfile.write("%s\n" % lnk_fname)
                wfile.write("%s\n" % "MIE")
                wfile.write("%d %f %f %f %d %f %f %f\n" %
                            (1, 0.0, np.log10(gsize[igs]), np.log10(gsize[igs]), 1., -3.5, gdens, -2.0))

        #
        # Run the Mie-code
        #
        dum = sp.Popen('makedust', shell=True).wait()
