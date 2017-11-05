"""This module contains classes to handle molecular data
"""
from __future__ import absolute_import
from __future__ import print_function
import traceback

try:
    import numpy as np
except ImportError:
    np = None
    print(' Numpy cannot be imported ')
    print(' To use the python module of RADMC-3D you need to install Numpy')
    print(traceback.format_exc())


from . import natconst as nc


class radmc3dMolecule(object):
    """
    RADMC-3D molecule class
    Based on the Leiden LAMDA database, but is in principle generic

    NOTE: For now only the levels and lines are included, not the
          collision rates.

    Attributes
    ----------
    name            : str
                     The name as listed in the molecule file

    molweight       : float
                     Molecular weight in units of proton mass

    nlev            : int
                     Nr of energy levels

    nlin            : int
                     Nr of lines

    energycminv     : float
                     Energy[ilev] of level ilev in 1/cm

    energy          : float
                     Energy[ilev] of level ilev in erg

    wgt             : float
                     Statistical weight[ilev] of level ilev

    jrot            : float
                     Quantum rotational J[ilev] of level ilev

    iup             : int
                     ilev of upper level of line ilin (starting with 0)

    ilow            : int
                     ilev of lower level of line ilin (starting with 0)

    aud             : float
                     Einstein A up low of line ilin in 1/second

    freq            : float
                     Frequency of line ilin in Hz

    lam             : float
                     Wavelength of line ilin in micron

    """

    def __init__(self):
        self.name = ""
        self.molweight = 0.0
        self.nlev = 0
        self.nlin = 0
        self.energycminv = 0.0
        self.energy = 0.0
        self.wgt = 0.0
        self.jrot = 0.0
        self.iup = 0
        self.ilow = 0
        self.aud = 0.0
        self.freq = 0.0
        self.lam = 0.0

    # --------------------------------------------------------------------------------------------------
    def read(self, mol=None, fname=None):
        """Read the molecule_<mol>.inp file

        The file format is the format of the Leiden LAMDA molecular database

        Parameters
        ----------
        mol             : str
                         molecule name (e.g. 'co') if the file name is in the form of 'molecule_<mol>.inp'

        fname           : str
                         full file name
        """

        if fname is None:
            if mol is None:
                raise ValueError('Unknown fname and mol. Either fname or mol must be set.')
            else:
                fname = 'molecule_' + mol + '.inp'

        with open(fname, 'r') as rfile:
            print('Reading ' + fname + '...')
            # with open(fname,'r') as f:
            dum = rfile.readline()
            dum = rfile.readline().split()
            self.name = dum[0]
            dum = rfile.readline()
            self.molweight = float(rfile.readline())
            dum = rfile.readline()
            self.nlev = int(rfile.readline())
            dum = rfile.readline()
            self.energycminv = np.zeros(self.nlev, dtype=np.float64)
            self.energy = np.zeros(self.nlev, dtype=np.float64)
            self.wgt = np.zeros(self.nlev, dtype=np.float64)
            self.jrot = np.zeros(self.nlev, dtype=np.float64)
            for i in range(self.nlev):
                dum = rfile.readline().split()
                self.energycminv[i] = float(dum[1])
                self.energy[i] = float(dum[1]) * nc.hh * nc.cc
                self.wgt[i] = float(dum[2])
                self.jrot[i] = float(dum[3])
            dum = rfile.readline()
            self.nlin = int(rfile.readline())
            dum = rfile.readline()
            self.iup = np.zeros(self.nlin, dtype=np.int)
            self.ilow = np.zeros(self.nlin, dtype=np.int)
            self.aud = np.zeros(self.nlin, dtype=np.float64)
            self.freq = np.zeros(self.nlin, dtype=np.float64)
            self.lam = np.zeros(self.nlin, dtype=np.float64)
            for i in range(self.nlin):
                dum = rfile.readline().split()
                self.iup[i] = int(dum[1])  # Use as index: [iup-1]
                self.ilow[i] = int(dum[2])  # Use as index: [ilow-1]
                self.aud[i] = float(dum[3])
                self.freq[i] = float(dum[4]) * 1e9
                self.lam[i] = nc.cc / self.freq[i]

        return True
