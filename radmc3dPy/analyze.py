"""This module contains classes and functions to read and write input/output data for RADMC-3D and
to do some simple analysis/diagnostics of the model.
"""
from __future__ import print_function
import traceback
import subprocess as sp
from multiprocessing import Pool
from functools import partial
import inspect
import os
import importlib

try:
    import numpy as np
except ImportError:
    np = None
    print(' Numpy cannot be imported ')
    print(' To use the python module of RADMC-3D you need to install Numpy')
    print(traceback.format_exc())

try:
    import matplotlib.pylab as plb
except ImportError:
    plb = None
    print('Warning')
    print('matplotlib.pylab cannot be imported')
    print('Without matplotlib you can use the python module to set up a model but you will not be able to plot things')
    print('or display images')

from matplotlib.colors import LogNorm
import matplotlib.patches as patches
import matplotlib.lines as ml

from . natconst import *
from . import natconst as nc
from . import crd_trans
from . import staratm


class radmc3dOctree(object):
    """
    Octree-like object with switchable resolution in each dimension

    Attributes
    ----------

    xi                : ndarray
                        Base grid cell interface grid in the first dimension

    yi                : ndarray
                        Base grid cell interface grid in the second dimension

    zi                : ndarray
                        Base grid cell interface grid in the third dimension

    xc                : ndarray
                        Base grid cell center grid in the first dimension

    yc                : ndarray
                        Base grid cell center grid in the second dimension

    zc                : ndarray
                        Base grid cell center grid in the third dimension

    x                 : ndarray
                        Tree cell center array in the first dimension

    y                 : ndarray
                        Tree cell center array in the second dimension

    z                 : ndarray
                        Tree cell center array in the third dimension

    dx                : ndarray
                        Tree cell halfwidth array in the first dimension

    dy                : ndarray
                        Tree cell halfwidth array in the second dimension

    dz                : ndarray
                        Tree cell halfwidth array in the third dimension

    leafID            : ndarray
                        Leaf index array,  mapping between a full tree and an array containing only the leaves

    isLeaf            : ndarray
                        Boolean array for the cell type (True - leaf, False - branch)

    level             : ndarray
                        Level array (base grid is level 0)

    parentID          : ndarray
                        Array containing the index of the parent cell (currently unused, only needed if we go up 
                        in the tree)

    childID           : list
                        List of children indices. Each list element is an ndarray with nChild elements containing the 
                        child indices

    act_dim           : list
                        A three element array to indicate which dimension is active, i.e. which dimensions are the 
                        cells resolved (0 - inactive, 1 - active)

    nCell             : int
                        Nr of cells (both branch and leaf) in the tree

    nxRoot            : int
                        Nr of cells in the base grid in the first dimension

    nyRoot            : int
                        Nr of cells in the base grid in the second dimension

    nzRoot            : int
                        Nr of cells in the base grid in the third dimension

    nLeaf             : int
                        Nr of leaf cells (i.e. true, unresolved grid cells)

    nBranch           : int
                        Nr of branches (i.e. resolved cells)

    nChild            : int
                        Nr of children (i.e. 8, 4, or 2 for 3, 2, 1 active dimensions, respectively)

    levelMax          : int
                        Highest actual level in the tree (Base grid has a level of 0 and the level increases)

    levelMaxLimit     : int
                        Highest allowed level in the tree (only used in tree building)

    crd_sys           : {'car', 'sph'}
                        Coordinate system type cartesian or spherical


    """

    def __init__(self):

        #
        # Arrays for the base grid (redundant information but improves speed and the memory overhead is not that large)
        #
        self.xi = np.zeros(0, dtype=np.float64)
        self.yi = np.zeros(0, dtype=np.float64)
        self.zi = np.zeros(0, dtype=np.float64)
        self.xc = np.zeros(0, dtype=np.float64)
        self.yc = np.zeros(0, dtype=np.float64)
        self.zc = np.zeros(0, dtype=np.float64)
        #
        # Grid cell center arrays
        #
        self.x = np.zeros(0, dtype=np.float64)
        self.y = np.zeros(0, dtype=np.float64)
        self.z = np.zeros(0, dtype=np.float64)
        #
        # Grid cell half(!)widths (this is in principle also redundant information as the only relevant information is
        # the cell width in the base grid, which is a single float for each dimension. The cell width at any level can
        # be simply calculated as cellwidth_base * 2^level. I still need to check how the performance would be affected
        # if I'd calculate the cell width on the fly instead of storing it in an array
        #
        self.dx = np.zeros(0, dtype=np.float64)
        self.dy = np.zeros(0, dtype=np.float64)
        self.dz = np.zeros(0, dtype=np.float64)
        #
        # Leaf index array - mapping between a full array and an array containing only leaves
        #
        self.leafID = np.zeros(0, dtype=np.int)
        #
        # Boolean array for the cell type (True - leaf, False - branch)
        #
        self.isLeaf = np.zeros(0, dtype=bool)
        #
        # Level array (base grid is level 0)
        #
        self.level = np.zeros(0, dtype=np.int)
        #
        # Array containing the index of the parent cell (currently unused, only needed if we go up in the tree)
        #
        self.parentID = np.zeros(0, dtype=np.int)
        #
        # List of children indices
        #
        self.childID = []
        self.cID = np.zeros(0, dtype=int)
        #
        # Number of cells in the whole tree and for the base grid
        #
        self.nCell = 0
        self.nxRoot = 0
        self.nyRoot = 0
        self.nzRoot = 0
        #
        # Highest actual level in the tree (Base grid has a level of 0 and the level increases)
        #
        self.levelMax = 0
        #
        # Highest allowed level in the tree (only used in tree building)
        #
        self.levelMaxLimit = 0
        #
        # Nr of branches (i.e. resolved cells)
        #
        self.nBranch = 0
        #
        # Nr of leaf cells (i.e. true, unresolved grid cells)
        #
        self.nLeaf = 0
        #
        # Coordinate sytem and active (resolvable) dimensions
        #
        self.crd_sys = 'car'
        self.act_dim = [1, 1, 1]
        self.grid_style = 1
        self.nChild = 0
        #
        # Stuff used for tree building
        #
        self.model = None
        #
        # Variable meant to be used internally only
        #
        self.cellIDCur = -1

        self.nwav = 0
        self.nfreq = 0
        self.wav = np.zeros(0, dtype=np.int)
        self.freq = np.zeros(0, dtype=np.int)

        self.counter = -1

    def getCellVolume(self, fullTree=False):
        """
        Calculates the grid cell volume

        Parameters
        ----------

        fullTree    : bool, optional
                      If True the cell volumes of the full tree (including both branches and leaves) will be 
                      calculated, while if set to False (default) the volume of only the leaf cells will be calculated

        Returns
        -------
        An linear ndarray containing the cell volumes
        """
        if self.crd_sys == 'car':
            vol = (2. * self.dx) * (2. * self.dy) * (2. * self.dz)
        else:
            vol = 1. / 3. * (self.x + self.dx)**3 - (self.x - self.dx)**3 \
                                                      * (np.cos(self.y - self.dy) - np.cos(self.y + self.dy)) \
                                                      * self.dz

        if not fullTree:
            ii = (self.leafID >= 0)
            dummy_vol = np.array(vol)
            vol = np.zeros(self.nLeaf, dtype=np.float64)
            vol[self.leafID[ii]] = dummy_vol[ii]

        return vol

    def _getContainerLeafIDRec(self, crd=(), cellID=-1):
        """
        Recursive function to find the tree index of a leaf that contains a given coordinate

        Parameters
        ----------
        crd         : tuple
                      List/tuple/ndarray containing the tree dimensional coordinates of the point 

        cellID      : int
                      Cell index
        """

        xmin = self.x[cellID] - self.dx[cellID]
        xmax = self.x[cellID] + self.dx[cellID]
        ymin = self.y[cellID] - self.dy[cellID]
        ymax = self.y[cellID] + self.dy[cellID]
        zmin = self.z[cellID] - self.dz[cellID]
        zmax = self.z[cellID] + self.dz[cellID]

        if self.isLeaf[cellID]:
            if (((crd[0] >= xmin) & (crd[0] < xmax)) &
                    ((crd[1] >= ymin) & (crd[1] < ymax)) &
                    ((crd[2] >= zmin) & (crd[2] < zmax))):
                return cellID
            else:
                return None
        else:
            dum = None
            for i in range(self.nChild):
                dum = self._getContainerLeafIDRec(crd, self.childID[cellID][i])
                if dum is not None:
                    break
            return dum

    def getContainerLeafID(self, crd=()):
        """
        Finds the tree index of a leaf that contains a given coordinate

        Parameters
        ----------
        crd         : tuple
                      List/tuple/ndarray containing the tree dimensional coordinates of the point 
        """

        leafID = -1

        if (crd[0] < self.xi[0]) | (crd[0] > self.xi[-1]):
            return leafID
        if (crd[1] < self.yi[0]) | (crd[1] > self.yi[-1]):
            return leafID
        if (crd[2] < self.zi[0]) | (crd[2] > self.zi[-1]):
            return leafID

        ix = np.searchsorted(self.xi, crd[0])
        iy = np.searchsorted(self.yi, crd[1])
        iz = np.searchsorted(self.zi, crd[2])

        if self.xi[ix] != crd[0]:
            ix -= 1
        if self.yi[iy] != crd[1]:
            iy -= 1
        if self.zi[iz] != crd[2]:
            iz -= 1

        if crd[0] == self.xi[-1]:
            ix = self.nxRoot - 1
        if crd[1] == self.yi[-1]:
            iy = self.nyRoot - 1
        if crd[2] == self.zi[-1]:
            iz = self.nzRoot - 1

        ind = iz * self.nyRoot * self.nxRoot + iy * self.nxRoot + ix
        leafID = self._getContainerLeafIDRec(crd, ind)

        return leafID

    # --------------------------------------------------------------------------------------------------
    def makeWavelengthGrid(self, wbound=None, nw=None, ppar=None):
        """Creates the wavelength/frequency grid.

        Parameters
        ----------

        wbound : list 
                 Contains the wavelength boundaries of the wavelength grid (should contain at least two elements)

        nw     : list 
                 Contains len(wbound)-1 elements containing the number of wavelengths between the bounds
                 set by wbound

        ppar   : dictionary, optional
                 Contains all input parameters with the parameter names as keys  
        """

        if ppar:
            if wbound is None:
                wbound = ppar['wbound']
            if nw is None:
                nw = ppar['nw']

        if wbound is None:
            raise ValueError('Unknown wbound. Wavelength grid bondaries are not specified.')
        if nw is None:
            raise ValueError('Unknown nw. Number of wavelength grid points are not specified.')

        self.nwav = nw[0]
        self.wav = wbound[0] * (wbound[1] / wbound[0])**(np.arange(nw[0], dtype=np.float64) / nw[0])

        for ipart in range(1, len(nw) - 1):
            dum = wbound[ipart] * (wbound[ipart + 1]
                                   / wbound[ipart])**(np.arange(nw[ipart], dtype=np.float64) / nw[ipart])
            self.wav = np.append(self.wav, dum)

        ipart = len(nw) - 1
        dum = wbound[ipart] * (wbound[ipart + 1]
                               / wbound[ipart])**(np.arange(nw[ipart], dtype=np.float64) / (nw[ipart] - 1.))
        self.wav = np.append(self.wav, dum)
        self.nwav = self.wav.shape[0]
        self.freq = nc.cc / self.wav * 1e4
        self.nfreq = self.nwav

    def readGrid(self):
        """
        Reads the spatial and wavelength grids from files
        """
        self.readWavelengthGrid()
        self.readSpatialGrid()

    # --------------------------------------------------------------------------------------------------
    def readWavelengthGrid(self, fname='wavelength_micron.inp'):
        """
        Function to read the wavelength/frequency grid

        Parameters
        ----------

        fname       : str, optional
                      Name of the file to read the wavelength grid from (if not specified wavelenth_micron.inp will 
                      be used)

        """
        #
        # Read the frequency grid
        #

        with open(fname, 'r') as rfile:
            self.nwav = int(rfile.readline())
            self.nfreq = self.nwav
            self.wav = np.zeros(self.nwav, dtype=np.float64)

            for i in range(self.nwav):
                self.wav[i] = float(rfile.readline())

            self.freq = nc.cc / self.wav * 1e4

    def writeWavelengthGrid(self, fname='', old=False):
        """Wriites the wavelength grid to a file (e.g. wavelength_micron.inp).

        Parameters
        ----------

        fname  : str, optional
                 File name into which the wavelength grid should be written. If omitted 'wavelength_micron.inp' will 
                 be used

        old    : bool, optional
                 If set to True the file format of the previous, 2D version of radmc will be used
        """

        if not old:
            if fname == '':
                fname = 'wavelength_micron.inp'

            print('Writing ' + fname)

            with open(fname, 'w') as wfile:
                wfile.write('%d\n' % self.nwav)
                for ilam in range(self.nwav):
                    wfile.write('%.9e\n' % self.wav[ilam])
        else:
            if fname == '':
                fname = 'frequency.inp'
            with open(fname, 'w') as wfile:
                print('Writing ' + fname)
                wfile.write("%d\n" % self.nfreq)
                wfile.write(" \n")
                #
                # Reverse the order of the frequency grid as it is ordered in frequency in radmc
                #
                freq = self.freq[::-1]
                for i in range(self.nfreq):
                    wfile.write("%.7e\n" % freq[i])

    def putNode(self, crd=(), cellsize=(), level=None, parentID=-1, cellID=None):
        """
        Function to put the data of a single node into the tree. This funcion assumes that all the arrays
        have already been allocated for the tree so input cell indices must refer to already existing array elements.

        Parameters
        ----------
        crd      : tuple
                   Cell center coordinates of the node

        cellsize : tuple
                   Full size of the cell in each dimension

        level    : int
                   Level of the cell in the tree

        parentID : int
                   Tree index of the parent cell

        cellID   : int
                   Tree index of the cell to be added
        """

        #
        # Add the cell centre and cell half width to the arrays
        #
        self.x[cellID] = crd[0]
        self.y[cellID] = crd[1]
        self.z[cellID] = crd[2]

        self.dx[cellID] = cellsize[0] * 0.5
        self.dy[cellID] = cellsize[1] * 0.5
        self.dz[cellID] = cellsize[2] * 0.5

        self.isLeaf[cellID] = True
        self.level[cellID] = level
        self.parentID[cellID] = parentID
        self.childID.append(np.zeros(self.nChild, dtype=np.int))

        return

    def resolveNodes(self, rsIDs=None):
        """
        Resolve multiple nodes simultaneously and add the children of the resolved node to the tree arrays extending
        the tree array

        Parameters
        ----------
        rsIDs       : list  
                      List/tuple/array of indices of the resolvable cell in the tree array
        """

        if isinstance(rsIDs, np.ndarray):
            ncell = rsIDs.shape[0]
        else:
            ncell = len(rsIDs)

        self.nLeaf -= ncell
        self.nBranch += ncell
        self.nLeaf += ncell * self.nChild

        x = np.zeros(ncell * self.nChild, dtype=np.float64)
        y = np.zeros(ncell * self.nChild, dtype=np.float64)
        z = np.zeros(ncell * self.nChild, dtype=np.float64)
        dx = np.zeros(ncell * self.nChild, dtype=np.float64)
        dy = np.zeros(ncell * self.nChild, dtype=np.float64)
        dz = np.zeros(ncell * self.nChild, dtype=np.float64)
        isLeaf = np.zeros(ncell * self.nChild, dtype=bool)
        level = np.zeros(ncell * self.nChild, dtype=np.int)
        parentID = np.zeros(ncell * self.nChild, dtype=np.int)
        ind = np.arange(ncell, dtype=np.int) * self.nChild
        nx = self.x.shape[0]

        xc_offset = None
        yc_offset = None
        zc_offset = None
        #
        # Generate the cell center offsets for a proper octree
        #
        if self.nChild == 8:
            xc_offset = np.array([-0.5, 0.5, -0.5, 0.5, -0.5, 0.5, -0.5, 0.5], dtype=np.float64) * self.dx[rsIDs][0]
            yc_offset = np.array([-0.5, -0.5, 0.5, 0.5, -0.5, -0.5, 0.5, 0.5], dtype=np.float64) * self.dy[rsIDs][0]
            zc_offset = np.array([-0.5, -0.5, -0.5, -0.5, 0.5, 0.5, 0.5, 0.5], dtype=np.float64) * self.dz[rsIDs][0]

        #
        # If we are resolving only two dimensions
        #
        elif self.nChild == 4:

            if self.act_dim[0] == 0:
                xc_offset = np.array([0., 0., 0., 0.], dtype=np.float64)
                yc_offset = np.array([-0.5, 0.5, -0.5, 0.5], dtype=np.float64) * self.dy[rsIDs][0]
                zc_offset = np.array([-0.5, -0.5, 0.5, 0.5], dtype=np.float64) * self.dz[rsIDs][0]
            elif self.act_dim[1] == 0:
                xc_offset = np.array([-0.5, 0.5, -0.5, 0.5], dtype=np.float64) * self.dx[rsIDs][0]
                yc_offset = np.array([0., 0., 0., 0.], dtype=np.float64)
                zc_offset = np.array([-0.5, -0.5, 0.5, 0.5], dtype=np.float64) * self.dz[rsIDs][0]
            elif self.act_dim[2] == 0:
                xc_offset = np.array([-0.5, 0.5, -0.5, 0.5], dtype=np.float64) * self.dx[rsIDs][0]
                yc_offset = np.array([-0.5, -0.5, 0.5, 0.5], dtype=np.float64) * self.dy[rsIDs][0]
                zc_offset = np.array([0., 0., 0., 0.], dtype=np.float64)

        #
        # Generate the cell center offsets if only two dimensions are resolved
        #
        elif self.nChild == 2:

            if self.act_dim[0] == 1:
                xc_offset = np.array([-0.5, 0.5], dtype=np.float64) * self.dx[rsIDs][0]
                yc_offset = np.array([0., 0.], dtype=np.float64)
                zc_offset = np.array([0., 0.], dtype=np.float64)

            if self.act_dim[1] == 1:
                xc_offset = np.array([0., 0.], dtype=np.float64)
                yc_offset = np.array([-0.5, 0.5], dtype=np.float64) * self.dy[rsIDs][0]
                zc_offset = np.array([0., 0.], dtype=np.float64)

            if self.act_dim[2] == 1:
                xc_offset = np.array([0., 0.], dtype=np.float64)
                yc_offset = np.array([0., 0.], dtype=np.float64)
                zc_offset = np.array([-0.5, 0.5], dtype=np.float64) * self.dz[rsIDs][0]

        #
        # Generate the cell center offsets if only one dimension is resolved
        #
        for i in range(self.nChild):
            x[ind + i] = self.x[rsIDs] + xc_offset[i]
            y[ind + i] = self.y[rsIDs] + yc_offset[i]
            z[ind + i] = self.z[rsIDs] + zc_offset[i]
            if self.act_dim[0] == 1:
                dx[ind + i] = np.zeros(ncell, dtype=np.float64) + self.dx[rsIDs][0] * 0.5
            else:
                dx[ind + i] = np.zeros(ncell, dtype=np.float64) + self.dx[rsIDs][0]

            if self.act_dim[1] == 1:
                dy[ind + i] = np.zeros(ncell, dtype=np.float64) + self.dy[rsIDs][0] * 0.5
            else:
                dy[ind + i] = np.zeros(ncell, dtype=np.float64) + self.dy[rsIDs][0]

            if self.act_dim[2] == 1:
                dz[ind + i] = np.zeros(ncell, dtype=np.float64) + self.dz[rsIDs][0] * 0.5
            else:
                dz[ind + i] = np.zeros(ncell, dtype=np.float64) + self.dz[rsIDs][0]

            isLeaf[ind + i] = np.ones(ncell, dtype=bool)
            level[ind + i] = np.zeros(ncell, dtype=np.int) + self.level[rsIDs][0] + 1
            parentID[ind + i] = rsIDs

        childID = []
        cid = np.arange(self.nChild, dtype=np.int)
        for i in range(ncell * self.nChild):
            childID.append(cid)

        #
        # Add the child cells to the tree
        #
        self.x = np.append(self.x, x)
        self.y = np.append(self.y, y)
        self.z = np.append(self.z, z)
        self.dx = np.append(self.dx, dx)
        self.dy = np.append(self.dy, dy)
        self.dz = np.append(self.dz, dz)
        self.isLeaf = np.append(self.isLeaf, isLeaf)
        self.level = np.append(self.level, level)
        self.parentID = np.append(self.parentID, parentID)
        self.childID.extend(childID)

        #
        # Now add the child indices to the parents
        #
        for i in range(ncell):
            ind = rsIDs[i]
            self.childID[ind] = cid + nx + i * self.nChild
            self.isLeaf[ind] = False
            #
        # Update the array length variables
        #
        self.nCell = self.x.shape[0]
        self.cID = np.arange(self.nCell, dtype=np.int)

        return

    def makeSpatialGrid(self, ppar=None, levelMaxLimit=None, dfunc=None, model='', **kwargs):
        """
        Function to create an octree-like AMR grid

        Parameters
        ----------

        ppar             : dictionary
                            Dictionary containing all input parameters of the model (from the problem_params.inp file)

        model            : str
                            Name of the model to be used in the tree building

        dfunc            : function
                            A user defined function that decides whether an AMR grid cell should be refined

        levelMaxLimit    : int, optional
                            Highest allowable level of the tree. This keyword is optional. If not specified at input as
                            a separate keyword, levelMaxLimit should be present in the problem_params.inp file.

        """

        #
        # Set the model
        #
        self.setModel(model)

        if dfunc is None:
            fnamelist = [f[0] for f in inspect.getmembers(self.model) if inspect.isfunction(f[1])]
            if 'decisionFunction' in fnamelist:
                dfunc = self.model.decisionFunction

            else:
                raise NameError('decisionFunction is not defined. \n '
                                + 'It is required to decide when a node / cell should be resolved. \n '
                                + 'Decision function should be given either as a dfunc keyword argument in \n'
                                + 'setup function call or should be implemented in the model as decisionFunction()')

        if levelMaxLimit is None:
            if 'levelMaxLimit' not in ppar.keys():
                raise KeyError('Unknown levelMaxLimit. \n It is required for AMR-style grid generation. levelMaxLimit '
                               + 'should be defined either in problem_params.inp or '
                               + 'in the call of the makeSpatialGrid/problemSetupDust/problemSetupGas functions.')
            else:
                self.levelMaxLimit = ppar['levelMaxLimit']
        else:
            self.levelMaxLimit = levelMaxLimit
        self.crd_sys = ppar['crd_sys']
        self.act_dim = [1, 1, 1]
        if ppar['nx'] == 0:
            self.act_dim[0] = 0
        if ppar['ny'] == 0:
            self.act_dim[1] = 0
        if ppar['nz'] == 0:
            self.act_dim[2] = 0

            #
        # Generate the base grid
        #
        self.xi = ppar['xbound'][0] + np.arange(ppar['nx'][0] + 1, dtype=float) / float(ppar['nx'][0]) \
                  * (ppar['xbound'][1] - ppar['xbound'][0])
        self.yi = ppar['ybound'][0] + np.arange(ppar['ny'][0] + 1, dtype=float) / float(ppar['ny'][0]) \
                  * (ppar['ybound'][1] - ppar['ybound'][0])
        self.zi = ppar['zbound'][0] + np.arange(ppar['nz'][0] + 1, dtype=float) / float(ppar['nz'][0]) \
                  * (ppar['zbound'][1] - ppar['zbound'][0])
        self.xc = 0.5 * (self.xi[1:] + self.xi[:-1])
        self.yc = 0.5 * (self.yi[1:] + self.yi[:-1])
        self.zc = 0.5 * (self.zi[1:] + self.zi[:-1])

        cellsize_x = self.xi[1] - self.xi[0]
        cellsize_y = self.yi[1] - self.yi[0]
        cellsize_z = self.zi[1] - self.zi[0]
        self.nxRoot = ppar['nx'][0]
        self.nyRoot = ppar['ny'][0]
        self.nzRoot = ppar['nz'][0]

        self.levelMax = 0

        ind = 0
        for iz in range(self.nzRoot):
            for iy in range(self.nyRoot):
                for ix in range(self.nxRoot):
                    #
                    # Add nodes to the tree
                    #
                    self.x = np.append(self.x, self.xc[ix])
                    self.y = np.append(self.y, self.yc[iy])
                    self.z = np.append(self.z, self.zc[iz])

                    self.dx = np.append(self.dx, [cellsize_x * 0.5])
                    self.dy = np.append(self.dy, [cellsize_y * 0.5])
                    self.dz = np.append(self.dz, [cellsize_z * 0.5])

                    self.isLeaf = np.append(self.isLeaf, True)
                    self.level = np.append(self.level, 0)
                    self.parentID = np.append(self.parentID, -1)
                    self.childID.append(np.zeros(self.nChild, dtype=np.int))

                    self.nLeaf += 1
                    ind += 1
        self.nCell = self.x.shape[0]

        #
        # Now build the tree
        #
        if 1 in self.act_dim:
            self.nChild = 2**(np.array(self.act_dim, dtype=int).sum())

        if self.nChild > 0:
            print('Adaptive Mesh Refinement (AMR) is active')
        txt = 'Active dimensions : '
        for i in range(3):
            if self.act_dim[i] == 1:
                txt += ("%d " % i)
        print(txt)

        #
        # Now go level by level and check which cells are to be resolved and resolve what's necessary
        #
        for ilev in range(self.levelMaxLimit):
            print('Resolving level ' + ("%d" % ilev))
            #
            # Select the cells at the current level
            #
            cID = np.arange(self.x.shape[0], dtype=int)
            ii = (self.level == ilev)

            if True in ii:
                #
                # Check which cells to resolve
                #
                resolve = dfunc(self.x[ii], self.y[ii], self.z[ii], self.dx[ii][0], self.dy[ii][0], self.dz[ii][0],
                                model=self.model, ppar=ppar, **kwargs)

                jj = resolve[self.level[ii] < self.levelMaxLimit]
                #
                # If there are some to resolve do so
                #
                if True in jj:
                    ncell2resolve = cID[ii][jj].shape[0]
                    print('Cells to resolve at this level : ', ncell2resolve)
                    self.resolveNodes(rsIDs=cID[ii][jj])
                    self.levelMax += 1
                else:
                    print('No cells to resolve at this level')
            else:
                print('No cells to resolve at this level')

        self.childID = np.array(self.childID)
        #
        # Print out some statistics
        #
        print('Tree building done')
        print('Maximum tree depth : ', self.levelMax)
        print('Nr of branches     : ', self.nBranch)
        print('Nr of leaves       : ', self.nLeaf)
        # ncells_fullgrid = self.nChild**self.levelMax * self.nxRoot * self.nyRoot * self.nzRoot
        # cell_fraction = float(self.nLeaf + self.nBranch) / ncells_fullgrid
        # print 'Using '+("%.3f"%(cell_fraction*100))+'% memory of a regular grid at max resolution'

        self.generateLeafID()
        return

    def setModel(self, model=''):
        """
        Sets the model to be used for tree building

        Parameters
        ----------
        model       : str
                      Name of the model
        """

        try:
            self.model = importlib.import_module(model)
        except ImportError:
            try:
                # self.model = __import__('radmc3dPy.models.' + model, fromlist=[''])
                self.model = importlib.import_module('radmc3dPy.models.' + model)
            except ImportError:
                print(model + '.py could not be imported. \n '
                      + 'The model files should either be in the current working directory \n '
                      + 'or in the radmc3d python module directory')
                print(traceback.format_exc())

    def _selfCheckCounterRec(self, cellID=None):
        """
        Recursive function for consistency check of the tree
        """

        if self.isLeaf[cellID]:
            self.counter[0] += 1
            if self.childID[cellID].max() > self.counter[2]:
                self.counter[2] = self.childID[cellID].max()
        else:
            self.counter[1] += 1
            for i in range(self.nChild):
                self._selfCheckCounterRec(cellID=self.childID[cellID][i])

        return

    def selfCheck(self):
        """
        Performs a self-check of the tree allocation and report it to the screen
        """

        self.counter = np.zeros([3], dtype=np.int)
        nRoot = self.nxRoot * self.nyRoot * self.nzRoot
        for i in range(nRoot):
            self._selfCheckCounterRec(cellID=i)

        print('Tree consistency check')
        print('Tree depth      : ' + ("%d" % self.levelMax))
        print('Nr of leaves    : ' + ("%d" % self.counter[0]) + " should be " + ("%d" % self.nLeaf))
        print('Nr of branches  : ' + ("%d" % self.counter[1]) + " should be " + ("%d" % self.nBranch))
        print('Nr of cells     : ' + ("%d" % self.nCell) + " should be " + ("%d" % (self.nBranch + self.nLeaf)))
        print('Leaf array      : ' + ("%d" % self.isLeaf.shape[0]))
        print('Level array     : ' + ("%d" % self.level.shape[0]))
        print('ParentID array  : ' + ("%d" % self.parentID.shape[0]))
        print('ChildID list    : ' + ("%d" % len(self.childID)))
        print('Max childID     : ' + ("%d" % self.counter[2]))
        print('x array         : ' + ("%d" % self.x.shape[0]))
        print('y array         : ' + ("%d" % self.y.shape[0]))
        print('z array         : ' + ("%d" % self.z.shape[0]))
        print('dx array        : ' + ("%d" % self.dx.shape[0]))
        print('dy array        : ' + ("%d" % self.dy.shape[0]))
        print('dz array        : ' + ("%d" % self.dz.shape[0]))

        return

    def _generateLeafIDRec(self, cellID=None):
        """
        Recursive function to generate the leaf indices 
        """

        if self.isLeaf[cellID]:
            self.cellIDCur += 1
            self.leafID[cellID] = self.cellIDCur
        else:
            for i in range(self.nChild):
                self._generateLeafIDRec(self.childID[cellID][i])

    def generateLeafID(self):
        """
        Function to generate the cell index mapping from arrays containing the full tree and those containing 
        only the leaves
        """

        print('Generating leaf indices')
        self.leafID = np.zeros(self.nCell, dtype=np.int) - 1
        self.cellIDCur = -1
        nRoot = self.nxRoot * self.nyRoot * self.nzRoot
        for i in range(nRoot):
            self._generateLeafIDRec(i)

        print('Done')

    def convArrLeaf2Tree(self, var=None):
        """
        Converts a leaf array to full tree size. 

        Parameters
        ----------
        var     : ndarray
                  A one or two dimensional ndarray with the first dimension is the size of the full tree

        Returns
        -------
        A one or two dimensional ndarray with size of of the full tree in the first dimension
        """

        ii = (self.leafID >= 0)
        ndim = len(var.shape)
        treeVar = None

        if ndim == 1:
            treeVar = np.zeros(self.nLeaf + self.nBranch, dtype=var.dtype)
            treeVar[ii] = var[self.leafID[ii]]
        elif ndim == 2:
            treeVar = np.zeros([self.nLeaf + self.nBranch, var.shape[1]], dtype=var.dtype)
            for i in range(var.shape[1]):
                treeVar[ii, i] = var[self.leafID[ii], i]
        else:
            assert (ndim >= 2), 'Incorrect shape for input array. Input array has too many dimensions. ' \
                                + 'Octree AMR only supports one or two dimensional arrays with the first dimension ' \
                                + 'being the cell/ spatial dimension'

            assert (ndim < 1), 'Incorrect shape for input array. Input array has too few dimensions (i.e. < 1). '

        return treeVar

    def convArrTree2Leaf(self, var=None):
        """
        Converts a data array to leaf size. The input is a scalar or vector variable defined at all nodes and the 
        returned variable will only represent values at leaf nodes thereby reduced in length compared to the input.

        Parameters
        ----------
        var     : ndarray
                  A one or two dimensional ndarray with the first dimension is the size of the full tree

        Returns
        -------
        A one or two dimensional ndarray with size of nLeaf in the first dimension
        """

        ii = (self.leafID >= 0)
        ndim = len(var.shape)
        leafVar = None

        if ndim == 1:
            leafVar = np.zeros(self.nLeaf, dtype=var.dtype)
            leafVar[self.leafID[ii]] = var[ii]
        elif ndim == 2:
            leafVar = np.zeros([self.nLeaf, var.shape[1]], dtype=var.dtype)
            for i in range(var.shape[1]):
                leafVar[self.leafID[ii], i] = var[ii, i]
        else:
            assert (ndim >= 2), 'Incorrect shape for input array. Input array has too many dimensions. ' \
                                'Octree AMR only supports one or two dimensional arrays with the first dimension being'\
                                'the cell/ spatial dimension'

            assert (ndim < 1), 'Incorrect shape for input array. Input array has too few dimensions (i.e. < 1). '

        return leafVar

    def writeSpatialGrid(self, fname=''):
        """
        Writes the wavelength grid to a file (e.g. amr_grid.inp).

        Parameters
        ----------

        fname : str, optional
                File name into which the spatial grid should be written. If omitted 'amr_grid.inp' will be used. 

        """

        if fname == '':
            fname = 'amr_grid.inp'

        print('Writing ' + fname)
        with open(fname, 'w') as wfile:
            wfile.write('%d\n' % 1)  # Format number

            wfile.write('\n')

            wfile.write('%d\n' % 1)  # AMR self.style (0=regular, 1 - Octree, 10 - Layered)
            if self.crd_sys == 'car':
                wfile.write('%d\n' % 0)  # Coordinate system (0-99 cartesian, 100-199 spherical, 200-299 cylindrical)
            if self.crd_sys == 'sph':
                wfile.write('%d\n' % 100)  # Coordinate system (0-99 cartesian, 100-199 spherical, 200-299 cylindrical)
            if self.crd_sys == 'cyl':
                wfile.write('%d\n' % 200)  # Coordinate system (0-99 cartesian, 100-199 spherical, 200-299 cylindrical)
            wfile.write('%d\n' % 0)  # Gridinfo

            wfile.write('\n')

            # Active dimensions
            wfile.write('%d %d %d \n' % (self.act_dim[0], self.act_dim[1], self.act_dim[2]))
            # Grid size (x,y,z or r,phi,theta, or r,phi,z)
            wfile.write(
                '%d %d %d \n' % (self.nxRoot, self.nyRoot, self.nzRoot))

            wfile.write('\n')

            # Max. refinement level, Nr of leaves and Nr of branches
            wfile.write('%d %d %d \n' % (self.levelMax, self.nLeaf, self.nBranch + self.nLeaf))
            wfile.write('\n')
            for i in range(self.nxRoot + 1):
                wfile.write('%.9e\n' % self.xi[i])
            wfile.write('\n')
            for i in range(self.nyRoot + 1):
                wfile.write('%.9e\n' % self.yi[i])
            wfile.write('\n')
            for i in range(self.nzRoot + 1):
                wfile.write('%.9e\n' % self.zi[i])
            wfile.write('\n')
            wfile.write('\n')

            nRoot = self.nxRoot * self.nyRoot * self.nzRoot
            for i in range(nRoot):
                self._writeOcTreeNodeTypeRec(cellID=i, wfile=wfile)

        return

    def _writeOcTreeNodeTypeRec(self, cellID=None, wfile=None):
        """
        Recursive function to write the node type to file

        Parameters
        ----------

        cellID      : int
                      Tree index of the cell to be written

        wfile       : file
                      File object to write to
        """
        if self.isLeaf[cellID]:
            wfile.write("0\n")
        else:
            wfile.write("1\n")
            for i in range(self.childID[cellID].shape[0]):
                self._writeOcTreeNodeTypeRec(cellID=self.childID[cellID][i], wfile=wfile)

        return

    def readSpatialGrid(self, fname=''):
        """
        Reads the spatial grid from amr_grid.inp

        Parameters
        ----------

        fname : str, optional
                File name from which the spatial grid should be read. If omitted 'amr_grid.inp' will be used. 

        """

        if fname == '':
            fname = 'amr_grid.inp'

        with open(fname, 'r') as rfile:
            print('Reading ' + fname)

            form = float(rfile.readline())
            dum = rfile.readline()
            grid_style = float(rfile.readline())
            if grid_style != 1:
                raise ValueError('Unsupported AMR style in ' + fname + '. Current only Octree AMR is supported.')

            crd_system = int(rfile.readline())
            if crd_system < 100:
                self.crd_sys = 'car'
                print("Reading cartesian grid")
            elif (crd_system >= 100) & (crd_system < 200):
                self.crd_sys = 'sph'
                print("Reading spherical grid")
            elif (crd_system >= 200) & (crd_system < 300):
                self.crd_sys = 'cyl'
                print("Reading cylindrical grid")
            else:
                raise ValueError('Unsupported coordinate system in ' + fname)

            grid_info = float(rfile.readline())

            dum = rfile.readline()

            dum = rfile.readline().split()
            self.act_dim = [int(dum[i]) for i in range(len(dum))]
            print("Active dimensions : ", self.act_dim[0], self.act_dim[1], self.act_dim[2])
            if 1 in self.act_dim:
                self.nChild = 2**(np.array(self.act_dim, dtype=int).sum())
            dum = rfile.readline().split()
            self.nxRoot, self.nyRoot, self.nzRoot = int(dum[0]), int(dum[1]), int(dum[2])
            print("Base grid size : ", self.nxRoot, self.nyRoot, self.nzRoot)
            dum = rfile.readline()

            dum = rfile.readline().split()
            levelMax, nLeaf, nCell = int(dum[0]), int(dum[1]), int(dum[2])
            nBranch = nCell - nLeaf

            print("Tree depth : ", levelMax)
            print("Nr of leaves : ", nLeaf)
            print("Nr of cells : ", nCell)
            dum = rfile.readline()

            self.levelMax = 0
            self.nLeaf = 0
            self.nBranch = 0

            self.x = np.zeros(nCell, dtype=np.float64)
            self.y = np.zeros(nCell, dtype=np.float64)
            self.z = np.zeros(nCell, dtype=np.float64)
            self.dx = np.zeros(nCell, dtype=np.float64)
            self.dy = np.zeros(nCell, dtype=np.float64)
            self.dz = np.zeros(nCell, dtype=np.float64)

            self.isLeaf = np.ones(nCell, dtype=bool)
            self.level = np.zeros(nCell, dtype=np.int)
            self.parentID = np.zeros(nCell, dtype=np.int)
            self.childID = []

            #
            # First of all read the base grid
            #
            self.xi = np.zeros(self.nxRoot + 1, dtype=np.float64)
            self.yi = np.zeros(self.nyRoot + 1, dtype=np.float64)
            self.zi = np.zeros(self.nzRoot + 1, dtype=np.float64)

            for i in range(self.nxRoot + 1):
                self.xi[i] = float(rfile.readline())
            dum = rfile.readline()
            for i in range(self.nyRoot + 1):
                self.yi[i] = float(rfile.readline())
            dum = rfile.readline()
            for i in range(self.nzRoot + 1):
                self.zi[i] = float(rfile.readline())

            dum = rfile.readline()
            dum = rfile.readline()

            self.xc = (self.xi[0:self.nxRoot] + self.xi[1:self.nxRoot + 1]) * 0.5
            self.yc = (self.yi[0:self.nyRoot] + self.yi[1:self.nyRoot + 1]) * 0.5
            self.zc = (self.zi[0:self.nzRoot] + self.zi[1:self.nzRoot + 1]) * 0.5

            dx = self.xi[1] - self.xi[0]
            dy = self.yi[1] - self.yi[0]
            dz = self.zi[1] - self.zi[0]

            self.cellIDCur = 0
            for iz in range(self.nzRoot):
                for iy in range(self.nyRoot):
                    for ix in range(self.nxRoot):
                        self.putNode(crd=(self.xc[ix], self.yc[iy], self.zc[iz]), cellsize=(dx, dy, dz), level=0,
                                     parentID=-1, cellID=self.cellIDCur)
                        self.cellIDCur += 1
                        self.nLeaf += 1

            nRoot = self.nxRoot * self.nyRoot * self.nzRoot

            for i in range(nRoot):
                self._readGridNodeTypeOcTreeRec(cellID=i, rfile=rfile)
        #
        # Now read which of the cells should be resolved
        #
        self.nCell = self.nLeaf + self.nBranch
        self.generateLeafID()

        self.selfCheck()
        return

    def _readGridNodeTypeOcTreeRec(self, cellID=None, rfile=None):
        """
        Recursive function to write the node type to file

        Parameters
        ----------

        cellID      : int
                      Tree index of the cell to be read

        rfile       : file
                      File object to read from

        """
        dum = int(rfile.readline())
        if dum == 1:

            self.nLeaf -= 1
            self.nBranch += 1
            self.nLeaf += self.nChild
            if self.levelMax < self.level[cellID] + 1:
                self.levelMax = self.level[cellID] + 1
                print('Tree depth : ', self.levelMax)

            #
            # Generate the cell center offsets for a proper octree
            #
            if self.nChild == 8:
                xc_offset = np.array([-0.5, 0.5, -0.5, 0.5, -0.5, 0.5, -0.5, 0.5], dtype=np.float64) * self.dx[cellID]
                yc_offset = np.array([-0.5, -0.5, 0.5, 0.5, -0.5, -0.5, 0.5, 0.5], dtype=np.float64) * self.dy[cellID]
                zc_offset = np.array([-0.5, -0.5, -0.5, -0.5, 0.5, 0.5, 0.5, 0.5], dtype=np.float64) * self.dz[cellID]
            #
            # If we are resolving only two dimensions
            #
            elif self.nChild == 4:

                if self.act_dim[0] == 0:
                    xc_offset = np.array([0., 0., 0., 0.], dtype=np.float64)
                    yc_offset = np.array([-0.5, 0.5, -0.5, 0.5], dtype=np.float64) * self.dy[cellID]
                    zc_offset = np.array([-0.5, -0.5, 0.5, 0.5], dtype=np.float64) * self.dz[cellID]
                elif self.act_dim[1] == 0:
                    xc_offset = np.array([-0.5, 0.5, -0.5, 0.5], dtype=np.float64) * self.dx[cellID]
                    yc_offset = np.array([0., 0., 0., 0.], dtype=np.float64)
                    zc_offset = np.array([-0.5, -0.5, 0.5, 0.5], dtype=np.float64) * self.dz[cellID]
                elif self.act_dim[2] == 0:
                    xc_offset = np.array([-0.5, 0.5, -0.5, 0.5], dtype=np.float64) * self.dx[cellID]
                    yc_offset = np.array([-0.5, -0.5, 0.5, 0.5], dtype=np.float64) * self.dy[cellID]
                    zc_offset = np.array([0., 0., 0., 0.], dtype=np.float64)

            #
            # Generate the cell center offsets if only two dimensions are resolved
            #
            elif self.nChild == 2:

                if self.act_dim[0] == 1:
                    xc_offset = np.array([-0.5, 0.5], dtype=np.float64) * self.dx[cellID]
                    yc_offset = np.array([0., 0.], dtype=np.float64)
                    zc_offset = np.array([0., 0.], dtype=np.float64)

                if self.act_dim[1] == 1:
                    xc_offset = np.array([0., 0.], dtype=np.float64)
                    yc_offset = np.array([-0.5, 0.5], dtype=np.float64) * self.dy[cellID]
                    zc_offset = np.array([0., 0.], dtype=np.float64)

                if self.act_dim[2] == 1:
                    xc_offset = np.array([0., 0.], dtype=np.float64)
                    yc_offset = np.array([0., 0.], dtype=np.float64)
                    zc_offset = np.array([-0.5, 0.5], dtype=np.float64) * self.dz[cellID]

            else:
                raise ValueError('Wrong number of child/leaf cells. The number of leaf cells should be 2,4,8 in '
                                 + '1,2,3 dimensions, respectively. Insted ' + ("%d" % self.nChild) + ' have been '
                                 + 'claimed.')

            self.isLeaf[cellID] = False
            self.childID[cellID] = np.zeros(self.nChild, dtype=np.int)
            for i in range(self.nChild):
                self.childID[cellID][i] = self.cellIDCur
                dx = self.dx[cellID] * (2.0 - self.act_dim[0])
                dy = self.dy[cellID] * (2.0 - self.act_dim[1])
                dz = self.dz[cellID] * (2.0 - self.act_dim[2])

                self.putNode(
                    crd=(self.x[cellID] + xc_offset[i], self.y[cellID] + yc_offset[i], self.z[cellID] + zc_offset[i]),
                    cellsize=(dx, dy, dz), level=self.level[cellID] + 1,
                    parentID=cellID, cellID=self.cellIDCur)
                self.cellIDCur += 1
                self._readGridNodeTypeOcTreeRec(cellID=self.childID[cellID][i], rfile=rfile)

        return


class radmc3dData(object):
    """RADMC-3D data class.
        Reading and writing dust density/temperature, gas density/temperature/velocity,
        generating a legacy vtk file for visualization.


    Attributes
    ----------

    grid      : radmc3dGrid, radmc3dOctree
                Instance of the radmc3dGrid class, contains the spatial and frequency grids

    rhodust   : ndarray
                Dust density in g/cm^3 

    dusttemp  : ndarray
                Dust temperature in K 

    rhogas    : ndarray
                Gas density in g/cm^3

    ndens_mol : ndarray
                Number density of the molecule [molecule/cm^3]

    ndens_cp  : ndarray
                Number density of the collisional partner [molecule/cm^3]

    gasvel    : ndarray
                Gas velocity in cm/s 

    gastemp   : ndarray
                Gas temperature in K

    vturb     : ndarray
                Mictroturbulence in cm/s

    taux      : ndarray
                Optical depth along the x (cartesian) / r (cylindrical) / r (spherical) dimension

    tauy      : ndarray
                Optical depth along the y (cartesian) / theta (cylindrical) / theta (spherical) dimension

    tauz      : ndarray
                Optical depth along the z (cartesian) / z (cylindrical) / phi (spherical) dimension

    sigmadust : ndarray
                Dust surface density in g/cm^2

    sigmagas  : ndarray
                Gas surface density in molecule/cm^2 (or g/cm^2 depending on the dimension of rhogas)
    """

    def __init__(self, grid=None):

        if grid:
            self.grid = grid
        else:
            self.grid = None

        self.rhodust = np.zeros(0, dtype=np.float64)
        self.dusttemp = np.zeros(0, dtype=np.float64)
        self.rhogas = np.zeros(0, dtype=np.float64)
        self.ndens_mol = np.zeros(0, dtype=np.float64)
        self.ndens_cp = np.zeros(0, dtype=np.float64)
        self.gasvel = np.zeros(0, dtype=np.float64)
        self.gastemp = np.zeros(0, dtype=np.float64)
        self.vturb = np.zeros(0, dtype=np.float64)
        self.taux = np.zeros(0, dtype=np.float64)
        self.tauy = np.zeros(0, dtype=np.float64)
        self.tauz = np.zeros(0, dtype=np.float64)
        self.sigmadust = np.zeros(0, dtype=np.float64)
        self.sigmagas = np.zeros(0, dtype=np.float64)

    def _scalarfieldWriter(self, data=None, fname='', binary=True, octree=False):
        """Writes a scalar field to a file.

        Parameters
        ----------

        data   : ndarray
                Scalar variable to be written

        fname  : str
                Name of the file containing a scalar variable

        binary : bool
                If True the file will be in binary format, if False the file format is formatted ASCII text

        octree : bool
                If True data will be written for an octree model (data should be a 1D numpy array)

        """

        with open(fname, 'w') as wfile:
            if binary:
                if octree:
                    if len(data.shape) == 1:
                        hdr = np.array([1, 8, data.shape[0], 1], dtype=np.int)
                    elif len(data.shape) == 2:
                        hdr = np.array([1, 8, data.shape[0], data.shape[1]], dtype=np.int)
                    else:
                        raise ValueError('Incorrect shape. Data stored in an Octree-type grid should have 1 or 2 '
                                         + 'dimensions, with the second dimension being the dust species. '
                                         + ' The data to be written has a dimension of ' + ("%d" % len(data.shape))
                                         + '\n No data has been written')
                    hdr.tofile(wfile)
                    data.flatten(order='f').tofile(wfile)
                else:
                    if len(data.shape) == 3:
                        hdr = np.array([1, 8, self.grid.nx * self.grid.ny * self.grid.nz], dtype=int)
                        hdr.tofile(wfile)
                    elif len(data.shape) == 4:
                        hdr = np.array([1, 8, self.grid.nx * self.grid.ny * self.grid.nz, data.shape[3]], dtype=int)
                        hdr.tofile(wfile)
                    # Now we need to flatten the dust density array since the Ndarray.tofile function writes the
                    # array always in C-order while we need Fortran-order to be written
                    if len(data.shape) == 4:
                        data = np.swapaxes(data, 0, 3)
                        data = np.swapaxes(data, 1, 2)
                        data.tofile(wfile)
                    elif len(data.shape) == 3:
                        data = np.swapaxes(data, 0, 2)
                        data.tofile(wfile)
                    else:
                        raise ValueError('Incorrect shape. Data stored in a regular grid should have 3 or 4 dimensions'
                                         + ' with the fourth dimension being the dust species. '
                                         + ' The data to be written has a dimension of ' + ("%d" % len(data.shape))
                                         + '\n No data has been written')
            else:
                if octree:
                    if len(data.shape) == 1:
                        hdr = np.array([1, data.shape[0], 1], dtype=np.int)
                    elif len(data.shape) == 2:
                        hdr = np.array([1, data.shape[0], data.shape[1]], dtype=np.int)
                    else:
                        raise ValueError('Incorrect shape. Data stored in an Octree-type grid should have 1 or 2 '
                                         + 'dimensions, with the second dimension being the dust species. '
                                         + ' The data to be written has a dimension of ' + ("%d" % len(data.shape))
                                         + '\n No data has been written')

                    hdr.tofile(wfile, sep=" ", format="%d\n")
                    data.flatten(order='f').tofile(wfile, sep=" ", format="%.9e\n")
                else:
                    if len(data.shape) == 3:
                        hdr = np.array([1, self.grid.nx * self.grid.ny * self.grid.nz], dtype=int)
                        hdr.tofile(wfile, sep=" ", format="%d\n")
                        # Now we need to flatten the dust density array since the Ndarray.tofile function writes the
                        # array always in C-order while we need Fortran-order to be written
                        data = np.swapaxes(data, 0, 2)
                        data.tofile(wfile, sep=" ", format="%.9e\n")

                    elif len(data.shape) == 4:
                        hdr = np.array([1, self.grid.nx * self.grid.ny * self.grid.nz, data.shape[3]], dtype=int)
                        hdr.tofile(wfile, sep=" ", format="%d\n")
                        # Now we need to flatten the dust density array since the Ndarray.tofile function writes the
                        # array always in C-order while we need Fortran-order to be written
                        data = np.swapaxes(data, 0, 3)
                        data = np.swapaxes(data, 1, 2)
                        data.tofile(wfile, sep=" ", format="%.9e\n")
                    else:
                        raise ValueError('Incorrect shape. Data stored in a regular grid should have 3 or 4 dimensions'
                                         + ' with the fourth dimension being the dust species. '
                                         + ' The data to be written has a dimension of ' + ("%d" % len(data.shape))
                                         + '\n No data has been written')

    def _scalarfieldReader(self, fname='', binary=True, octree=False):
        """Reads a scalar field from file.

        Parameters
        ----------

        fname  : str 
                Name of the file containing a scalar variable

        binary : bool
                If True the file is in binary format, if False the file format is formatted ASCII text

        octree : bool
                If True data will be read from an octree model and will be stored in a 1D numpy array 

        Returns
        -------

        Returns a numpy Ndarray with the scalar field
        """

        if binary:
            if octree:
                hdr = np.fromfile(fname, count=4, dtype=int)
                if hdr[2] != self.grid.nLeaf:
                    print(hdr[1], self.grid.nLeaf)
                    raise ValueError('Number of cells in ' + fname + ' is different from that in amr_grid.inp'
                                     + ' nr cells in ' + fname + ' : ' + ("%d" % hdr[2]) + '\n '
                                     + ' nr of cells in amr_grid.inp : '
                                     + ("%d" % self.grid.nLeaf))

                if hdr[1] == 8:
                    data = np.fromfile(fname, count=-1, dtype=np.float64)
                elif hdr[1] == 4:
                    data = np.fromfile(fname, count=-1, dtype=float)
                else:
                    raise TypeError('Unknown datatype/precision in ' + fname + '. RADMC-3D binary files store 4 byte '
                                    + 'floats or 8 byte doubles. The precision in the file header is '
                                    + ("%d" % hdr[1]))

                if data.shape[0] == (hdr[2] + 3):
                    data = np.reshape(data[3:], [self.grid.nLeaf, 1], order='f')
                elif data.shape[0] == (hdr[2] * hdr[3] + 4):
                    data = np.reshape(data[4:], [self.grid.nLeaf, hdr[3]], order='f')

            else:
                # hdr[0] = format number
                # hdr[1] = data precision (4=single, 8=double)
                # hdr[2] = nr of cells
                # hdr[3] = nr of dust species
                hdr = np.fromfile(fname, count=4, dtype=int)
                if hdr[2] != (self.grid.nx * self.grid.ny * self.grid.nz):
                    raise ValueError('Number of grid cells in ' + fname + ' is different from that in amr_grid.inp '
                                     + ' nr cells in ' + fname + ' : ' + ("%d" % hdr[2]) + '\n '
                                     + ' nr of cells in amr_grid.inp : '
                                     + ("%d" % (self.grid.nx * self.grid.ny * self.grid.nz)))

                if hdr[1] == 8:
                    data = np.fromfile(fname, count=-1, dtype=np.float64)
                elif hdr[1] == 4:
                    data = np.fromfile(fname, count=-1, dtype=float)
                else:
                    raise TypeError('Unknown datatype/precision in ' + fname + '. RADMC-3D binary files store 4 byte '
                                    + 'floats or 8 byte doubles. The precision in the file header is '
                                    + ("%d" % hdr[1]))

                if data.shape[0] == (hdr[2] + 3):
                    data = np.reshape(data[3:], [1, self.grid.nz, self.grid.ny, self.grid.nx])
                elif data.shape[0] == (hdr[2] * hdr[3] + 4):
                    data = np.reshape(data[4:], [hdr[3], self.grid.nz, self.grid.ny, self.grid.nx])

                # data = reshape(data, [hdr[3],self.grid.nz,self.grid.ny,self.grid.nx])
                # We need to change the axis orders as Numpy always writes binaries in C-order while RADMC-3D
                # uses Fortran-order
                data = np.swapaxes(data, 0, 3)
                data = np.swapaxes(data, 1, 2)

        else:
            if octree:
                hdr = np.fromfile(fname, count=3, sep='\n', dtype=int)

                if hdr[1] != self.grid.nLeaf:
                    raise ValueError('Number of cells in ' + fname + ' is different from that in amr_grid.inp'
                                     + ' nr cells in ' + fname + ' : ' + ("%d" % hdr[1]) + '\n '
                                     + ' nr of cells in amr_grid.inp : '
                                     + ("%d" % self.grid.nLeaf))

                data = np.fromfile(fname, count=-1, sep='\n', dtype=np.float64)[3:]
                data = data.reshape([hdr[1], hdr[2]], order='f')

            else:
                hdr = np.fromfile(fname, count=3, sep='\n', dtype=int)

                if (self.grid.nx * self.grid.ny * self.grid.nz) != hdr[1]:
                    raise ValueError('Number of grid cells in ' + fname + ' is different from that in amr_grid.inp '
                                     + ' nr cells in ' + fname + ' : ' + ("%d" % hdr[2]) + '\n '
                                     + ' nr of cells in amr_grid.inp : '
                                     + ("%d" % (self.grid.nx * self.grid.ny * self.grid.nz)))
                else:

                    data = np.fromfile(fname, count=-1, sep="\n", dtype=np.float64)
                    if data.shape[0] == hdr[1] + 2:
                        data = np.reshape(data[2:], [1, self.grid.nz, self.grid.ny, self.grid.nx])
                    elif data.shape[0] == hdr[1] * hdr[2] + 3:
                        data = np.reshape(data[3:], [hdr[2], self.grid.nz, self.grid.ny, self.grid.nx])
                    # We need to change the axis orders as Numpy always reads  in C-order while RADMC-3D
                    # uses Fortran-order
                    data = np.swapaxes(data, 0, 3)
                    data = np.swapaxes(data, 1, 2)

        return data

    def getTauOneDust(self, idust=0, axis='', kappa=0.):
        """Calculates the optical depth of a single dust species along any given combination of the axes.

        Parameters
        ----------

        idust : int
                Index of the dust species whose optical depth should be calculated

        axis  : str
                Name of the axis/axes along which the optical depth should be calculated 
                (e.g. 'x' for the first dimension or 'xyz' for all three dimensions)

        kappa : float
                Mass extinction coefficients of the dust species at the desired wavelength

        Returns
        -------

        Returns a dictionary with the following keys

            taux  : ndarray
                    optical depth along the first dimension
            tauy  : ndarray
                    optical depth along the second dimension

            (tauz is not yet implemented)
        """

        # Check along which axis should the optical depth be calculated
        do_taux = False
        do_tauy = False

        if 'x' in axis:
            do_taux = True
        if 'y' in axis:
            do_tauy = True

        # Calculate the optical depth along the x-axis (r in spherical coordinates)
        if do_taux:
            taux = np.zeros([self.grid.nx, self.grid.ny, self.grid.nz], dtype=np.float64)
            diff_x = self.grid.xi[1:] - self.grid.xi[:-1]
            taux[0, :, :] = self.rhodust[0, :, :, idust] * kappa * diff_x[0]
            for ix in range(1, self.grid.nx):
                taux[ix, :, :] = taux[ix - 1, :, :] + self.rhodust[ix, :, :, idust] * kappa * diff_x[ix]
        else:
            taux = [-1.]

        # Calculate the optical depth along the theta in spherical coordinates
        # Warning the formulation below is valid only in spherical coordinate sytem

        dum_x = np.zeros([self.grid.nx, self.grid.nz], dtype=np.float64)
        for iz in range(self.grid.nz):
            dum_x[:, iz] = self.grid.x

        if do_tauy:
            tauy = np.zeros([self.grid.nx, self.grid.ny, self.grid.nz], dtype=np.float64)
            diff_y = self.grid.yi[1:] - self.grid.yi[:-1]
            tauy[:, 0, :] = self.rhodust[:, 0, :, idust] * kappa * diff_y[0] * dum_x
            for iy in range(1, self.grid.ny):
                tauy[:, iy, :] = tauy[:, iy - 1, :] + self.rhodust[:, iy, :, idust] * kappa * diff_y[iy] * dum_x
        else:
            tauy = [-1.]

        return {'taux': taux, 'tauy': tauy}

    def getTau(self, idust=None, axis='xy', wav=0., kappa=None, old=False):
        """Calculates the optical depth along any given combination of the axes.

        Parameters
        ----------

        idust : list
                List of dust component indices whose optical depth should be calculated
                If multiple indices are set the total optical depth is calculated summing 
                over all dust species in idust

        axis  : str
                Name of the axis/axes along which the optical depth should be calculated 
                (e.g. 'x' for the first dimension or 'xyz' for all three dimensions)

        wav   : float
                Wavelength at which the optical depth should be calculated

        kappa : list, tuple
                If set it should be a list of mass extinction coefficients at the desired wavelength
                The number of elements in the list should be equal to that in the idust keyword

        old   : bool, optional
                If set to True the file format of the previous, 2D version of radmc will be used
        """

        # Check if the input idust indices can be found in rhoudust
        if len(self.rhodust.shape) == 3:
            ndust = 1
        else:
            ndust = self.rhodust.shape[3]

        if idust is None:
            idust = np.arange(ndust)

        if max(idust) > ndust:
            raise ValueError(' There are less number of dust species than some of the indices in idust')

        scatmat = None
        # If the kappa keyword is set it should be used during the optical depth calculation
        if kappa is not None:
            # Safety check
            if len(kappa) != len(idust):
                raise ValueError(' The number of kappa values should be identical to the number of specified '
                                 + 'dust species ')
        else:
            if not old:
                # Read the master opacity file to get the dustkappa file name extensions
                dum = radmc3dDustOpac()
                mo = dum.readMasterOpac()
                dummy_ext = mo['ext']
                scatmat = mo['scatmat']
                if len(dummy_ext) <= max(idust):
                    raise ValueError('There are less dust species specified in dustopac.inp than '
                                     + 'some of the specified idust indices')
                else:
                    ext = [dummy_ext[i] for i in idust]

        if 'x' in axis:
            self.taux = np.zeros([self.grid.nx, self.grid.ny, self.grid.nz], dtype=np.float64)
        if 'y' in axis:
            self.tauy = np.zeros([self.grid.nx, self.grid.ny, self.grid.nz], dtype=np.float64)

        for i in idust:
            if kappa is None:
                if old:
                    opac = readOpac(ext=[("%d" % (i + 1))], old=True)
                else:
                    opac = readOpac(ext=ext[i], scatmat=scatmat)

                if not opac.ext:
                    return -1
                else:
                    kabs = 10. ** np.interp(np.log10(np.array(wav)), np.log10(opac.wav[0]), np.log10(opac.kabs[0]))
                if opac.ksca[0][0] > 0:
                    ksca = 10. ** np.interp(np.log10(np.array(wav)), np.log10(opac.wav[0]), np.log10(opac.ksca[0]))
                else:
                    ksca = np.array(kabs) * 0.

                print('Opacity at ' + ("%.2f" % wav) + 'um : ', kabs + ksca)
                dum = self.getTauOneDust(i, axis=axis, kappa=kabs + ksca)
            else:
                dum = self.getTauOneDust(i, axis=axis, kappa=kappa[i])

            if 'x' in axis:
                self.taux = self.taux + dum['taux']
            if 'y' in axis:
                self.tauy = self.tauy + dum['tauy']

    def getGasMass(self, mweight=2.3, rhogas=False):
        """Calculates the gas mass in radmc3dData.ndens_mol or radmc3dData.rhogas

        Parameters
        ----------
            mweight   : float
                        Molecular weight [atomic mass unit / molecule, i.e. same unit as mean molecular weight]

            rhogas    : bool, optional
                        If True the gas mass will be calculated from radmc3dData.rhogas, while if set to False
                        the gas mass will be calculated from radmc3dData.ndens_mol. The mweight parameter is only
                        required for the latter. 

        Returns
        -------
            A single float being the gas mass in gramm
        """

        vol = self.grid.getCellVolume()
        if not rhogas:
            #
            # I'm not sure if this is the right way of doing it but right now I don't have a better idea
            #
            if isinstance(self.grid, radmc3dOctree):
                return (vol * self.ndens_mol[:, 0] * mweight * nc.mp).sum()
            else:
                return (vol * self.ndens_mol[:, :, :, 0] * mweight * nc.mp).sum()

        else:
            if isinstance(self.grid, radmc3dOctree):
                return (vol * self.rhogas[:, 0]).sum()
            else:
                return (vol * self.rhogas[:, :, :, 0]).sum()


    def getDustMass(self, idust=-1):
        """Calculates the dust mass in radmc3dData.rhodust

        Parameters
        ----------
            idust   : int
                      Dust index whose dust should be calculated. If it is set to -1 (default) the total
                      dust mass is calculated summing up all dust species

        Returns
        -------
            A single float being the dust mass in gramm
        """

        vol = self.grid.getCellVolume()
        if idust > 0:
            #
            # I'm not sure if this is the right way of doing it but right now I don't have a better idea
            #
            if isinstance(self.grid, radmc3dOctree):
                dmass = (vol * self.rhodust[:, idust]).sum()
            else:
                dmass = (vol * self.rhodust[:, :, :, idust]).sum()
        else:
            dmass = 0.
            if isinstance(self.grid, radmc3dOctree):
                for i in range(self.rhodust.shape[1]):
                    dmass += (vol * self.rhodust[:, i]).sum()
            else:
                for i in range(self.rhodust.shape[3]):
                    dmass += (vol * self.rhodust[:, :, :, i]).sum()

        return dmass

    def readDustDens(self, fname='', binary=True, old=False, octree=False):
        """Reads the dust density.

        Parameters
        ----------

        fname   : str, optional
                  Name of the file that contains the dust density. If omitted 'dust_density.inp' is used
                  (or if binary=True the 'dust_density.binp' is used).

        binary  : bool, optional 
                  If true the data will be read in binary format, otherwise the file format is ascii

        old     : bool, optional
                  If set to True the file format of the previous, 2D version of radmc will be used

        octree  : bool, optional
                  If the data is defined on an octree-like AMR
        """

        if self.grid is None:
            if octree:
                self.grid = radmc3dOctree()
            else:
                self.grid = radmc3dGrid()
                self.grid.readGrid(old=old)

        print('Reading dust density')
        #
        # Read radmc3d output
        #
        if not old:

            if binary:
                if fname == '':
                    fname = 'dust_density.binp'
            else:
                if fname == '':
                    fname = 'dust_density.inp'

            self.rhodust = self._scalarfieldReader(fname=fname, binary=binary, octree=octree)
        #
        # Read the output of the previous 2d version of the code
        #
        else:
            with open('dustdens.inp', 'r') as rfile:

                dum = rfile.readline().split()
                ndust = int(dum[0])
                nr = int(dum[1])
                nt = int(dum[2])

                rfile.readline()
                self.rhodust = np.zeros([nr, nt * 2, 1, ndust], dtype=float)

                for idust in range(ndust):
                    for ix in range(nr):
                        for iy in range(nt):
                            self.rhodust[ix, iy, 0, idust] = float(rfile.readline())
                            self.rhodust[ix, self.grid.ny - 1 - iy, 0, idust] = self.rhodust[ix, iy, 0, idust]

        return True

    def readDustTemp(self, fname='', binary=True, octree=False, old=False):
        """Reads the dust temperature.

        Parameters
        ----------

        fname   : str, optional
                  Name of the file that contains the dust temperature. 

        binary  : bool, optional 
                  If true the data will be read in binary format, otherwise the file format is ascii

        octree  : bool, optional
                  If the data is defined on an octree-like AMR

        old     : bool, optional
                  If True dust temperature will be written in the old RADMC format
        """

        if self.grid is None:
            if octree:
                self.grid = radmc3dOctree()
            else:
                self.grid = radmc3dGrid()
                self.grid.readGrid(old=old)

        print('Reading dust temperature')

        if not old:
            if binary:
                if fname == '':
                    fname = 'dust_temperature.bdat'
            else:
                if fname == '':
                    fname = 'dust_temperature.dat'

            self.dusttemp = self._scalarfieldReader(fname=fname, binary=binary, octree=octree)
        else:
            with open('dusttemp_final.dat', 'r') as rfile:
                dum = rfile.readline().split()
                ndust = int(dum[0])
                nr = int(dum[1])
                nt = int(dum[2])

                rfile.readline()
                self.dusttemp = np.zeros([nr, nt * 2, 1, ndust], dtype=float)

                for idust in range(ndust):
                    rfile.readline()
                    for ix in range(nr):
                        for iy in range(nt):
                            self.dusttemp[ix, iy, 0, idust] = float(rfile.readline())
                            self.dusttemp[ix, self.grid.ny - 1 - iy, 0, idust] = self.dusttemp[ix, iy, 0, idust]

        return True

    def readGasVel(self, fname='', binary=True, octree=False):
        """Reads the gas velocity.  

        Parameters
        -----------

        fname : str, optional
                Name of the file that contains the gas velocity
                If omitted 'gas_velocity.inp' (if binary=True 'gas_velocity.binp')is used.

        binary : bool
                If true the data will be read in binary format, otherwise the file format is ascii

        octree  : bool, optional
                  If the data is defined on an octree-like AMR
        """

        if self.grid is None:
            if octree:
                self.grid = radmc3dOctree()
            else:
                self.grid = radmc3dGrid()
                self.grid.readGrid(old=False)

        if binary:
            if fname == '':
                fname = 'gas_velocity.binp'

            print('Reading gas velocity')
            if os.path.isfile(fname):
                # If we have an octree grid
                if isinstance(self.grid, radmc3dOctree):
                    hdr = np.fromfile(fname, count=3, dtype=int)
                    if hdr[2] != self.grid.nLeaf:
                        raise ValueError('Number of cells in ' + fname + ' is different from that in amr_grid.inp'
                                         + ' nr cells in ' + fname + ' : ' + ("%d" % hdr[2]) + '\n '
                                         + ' nr of cells in amr_grid.inp : '
                                         + ("%d" % self.grid.nLeaf))

                    if hdr[1] == 8:
                        self.gasvel = np.fromfile(fname, count=-1, dtype=np.float64)
                    elif hdr[1] == 4:
                        self.gasvel = np.fromfile(fname, count=-1, dtype=float)
                    else:
                        raise TypeError(
                            'Unknown datatype/precision in ' + fname + '. RADMC-3D binary files store 4 byte '
                            + 'floats or 8 byte doubles. The precision in the file header is '
                            + ("%d" % hdr[1]))
                    self.gasvel = np.reshape(self.gasvel[3:], [self.nLeaf, 3], order='f')

                else:
                    hdr = np.fromfile(fname, count=3, dtype=int)
                    if hdr[2] != self.grid.nx * self.grid.ny * self.grid.nz:
                        raise ValueError('Number of grid cells in ' + fname + ' is different from that in amr_grid.inp '
                                         + ' nr cells in ' + fname + ' : ' + ("%d" % hdr[2]) + '\n '
                                         + ' nr of cells in amr_grid.inp : '
                                         + ("%d" % (self.grid.nx * self.grid.ny * self.grid.nz)))

                    if hdr[1] == 8:
                        self.gasvel = np.fromfile(fname, count=-1, dtype=np.float64)
                    elif hdr[1] == 4:
                        self.gasvel = np.fromfile(fname, count=-1, dtype=float)
                    else:
                        raise TypeError(
                            'Unknown datatype/precision in ' + fname + '. RADMC-3D binary files store 4 byte '
                            + 'floats or 8 byte doubles. The precision in the file header is '
                            + ("%d" % hdr[1]))
                    self.gasvel = np.reshape(self.gasvel[3:], [self.grid.nz, self.grid.ny, self.grid.nx, 3])
                    self.gasvel = np.swapaxes(self.gasvel, 0, 2)
            else:
                raise FileNotFoundError(fname + 'was not found')

        else:

            if fname == '':
                fname = 'gas_velocity.inp'

            with open(fname, 'r') as rfile:

                print('Reading gas velocity')
                dum = rfile.readline()
                dum = int(rfile.readline())

                # If we have an octree grid
                if octree:
                    if self.grid.nLeaf != dum:
                        raise ValueError('Number of cells in ' + fname + ' is different from that in amr_grid.inp'
                                         + ' nr cells in ' + fname + ' : ' + ("%d" % dum) + '\n '
                                         + ' nr of cells in amr_grid.inp : '
                                         + ("%d" % self.grid.nLeaf))

                    self.gasvel = np.zeros([dum, 3], dtype=np.float64)
                    for i in range(dum):
                        val = rfile.readline().split()
                        self.gasvel[i, 0] = float(dum[0])
                        self.gasvel[i, 1] = float(dum[1])
                        self.gasvel[i, 2] = float(dum[2])
                else:

                    if (self.grid.nx * self.grid.ny * self.grid.nz) != dum:
                        raise ValueError('Number of grid cells in ' + fname + ' is different from that in amr_grid.inp '
                                         + ' nr cells in ' + fname + ' : ' + ("%d" % dum) + '\n '
                                         + ' nr of cells in amr_grid.inp : '
                                         + ("%d" % (self.grid.nx * self.grid.ny * self.grid.nz)))

                    self.gasvel = np.zeros([self.grid.nx, self.grid.ny, self.grid.nz, 3], dtype=np.float64)

                    for k in range(self.grid.nz):
                        for j in range(self.grid.ny):
                            for i in range(self.grid.nx):
                                dum = rfile.readline().split()
                                self.gasvel[i, j, k, 0] = float(dum[0])
                                self.gasvel[i, j, k, 1] = float(dum[1])
                                self.gasvel[i, j, k, 2] = float(dum[2])

        return True

    def readVTurb(self, fname='', binary=True, octree=False):
        """Reads the turbulent velocity field. 

        Parameters
        ----------

        fname   : str, optional 
                  Name of the file that contains the turbulent velocity field
                  If omitted 'microturbulence.inp' (if binary=True 'microturbulence.binp') is used.

        binary  : bool 
                  If true the data will be read in binary format, otherwise the file format is ascii

        octree  : bool, optional
                  If the data is defined on an octree-like AMR
        """
        if self.grid is None:
            if octree:
                self.grid = radmc3dOctree()
            else:
                self.grid = radmc3dGrid()
                self.grid.readGrid(old=False)

        print('Reading microturbulence')

        if binary:
            if fname == '':
                fname = 'microturbulence.binp'
        else:
            if fname == '':
                fname = 'microturbulence.inp'

        self.vturb = self._scalarfieldReader(fname=fname, binary=binary, octree=octree)
        if octree:
            self.vturb = np.squeeze(self.vturb)

        return True

    def readGasDens(self, ispec='', binary=True, octree=False):
        """Reads the gas density.

        Parameters
        ----------

        ispec   : str 
                  File name extension of the 'numberdens_ispec.inp' (or if binary=True 'numberdens_ispec.binp') file.

        binary  : bool 
                  If true the data will be read in binary format, otherwise the file format is ascii

        octree  : bool, optional
                  If the data is defined on an octree-like AMR
        """

        if self.grid is None:
            if octree:
                self.grid = radmc3dOctree()
            else:
                self.grid = radmc3dGrid()
                self.grid.readGrid(old=False)

        if binary:
            fname = 'numberdens_' + ispec + '.binp'
        else:
            fname = 'numberdens_' + ispec + '.inp'

        print('Reading gas density (' + fname + ')')
        self.ndens_mol = self._scalarfieldReader(fname=fname, binary=binary, octree=octree)
        if octree:
            self.ndens_mol = np.squeeze(self.ndens_mol)

        return True

    def readGasTemp(self, fname='', binary=True, octree=False):
        """Reads the gas temperature.

        Parameters
        ----------

        fname   : str,optional
                  Name of the file that contains the gas temperature. If omitted 'gas_temperature.inp' 
                  (or if binary=True 'gas_tempearture.binp') is used.

        binary  : bool
                  If true the data will be read in binary format, otherwise the file format is ascii

        octree  : bool, optional
                  If the data is defined on an octree-like AMR
        """

        if self.grid is None:
            if octree:
                self.grid = radmc3dOctree()
            else:
                self.grid = radmc3dGrid()
                self.grid.readGrid(old=False)

        print('Reading gas temperature')

        if binary:
            if fname == '':
                fname = 'gas_temperature.binp'
        else:
            if fname == '':
                fname = 'gas_temperature.inp'

        self.gastemp = self._scalarfieldReader(fname=fname, binary=binary, octree=octree)
        if octree:
            self.gastemp = np.squeeze(self.gastemp)

        return True

    def writeDustDens(self, fname='', binary=True, old=False, octree=False):
        """Writes the dust density.

        Parameters
        ----------

        fname   : str, optional
                  Name of the file into which the dust density should be written. If omitted 'dust_density.inp' is used.

        binary  : bool
                  If true the data will be written in binary format, otherwise the file format is ascii

        old     : bool, optional
                  If set to True the file format of the previous, 2D version of radmc will be used

        octree  : bool, optional
                  If the data is defined on an octree-like AMR
        """

        #
        # Write dust density for radmc3d
        #
        if not old:
            if fname == '':
                if binary:
                    fname = 'dust_density.binp'
                else:
                    fname = 'dust_density.inp'

            print('Writing ' + fname)

            if octree:
                self._scalarfieldWriter(data=self.grid.convArrTree2Leaf(self.rhodust), fname=fname, binary=binary,
                                        octree=True)
            else:
                self._scalarfieldWriter(data=self.rhodust, fname=fname, binary=binary, octree=False)

        #
        # Write dust density for the previous 2D version of the code
        #
        else:
            if self.rhodust.shape[2] > 1:
                raise ValueError('Wrong dimensions for dust density. RADMC (the predecessor code to RADMC-3D is '
                                 + 'strictly 2D, yet the dust density to be written in the RADMC format is 3D.')

            if fname == '':
                fname = 'dustdens.inp'

            with open(fname, 'w') as wfile:

                ntheta = int(self.grid.ny / 2)
                wfile.write("%d %d %d 1\n" % (self.rhodust.shape[3], self.grid.nx, ntheta))
                wfile.write(" \n")
                for idust in range(self.rhodust.shape[3]):
                    for ix in range(self.grid.nx):
                        for iy in range(ntheta):
                            wfile.write("%.7e\n" % self.rhodust[ix, iy, 0, idust])

    def writeDustTemp(self, fname='', binary=True, octree=False):
        """Writes the dust density.

        Parameters
        ----------

        fname   : str, optional
                Name of the file into which the dust density should be written. If omitted 'dust_density.inp' is used.

        binary  : bool
                  If true the data will be written in binary format, otherwise the file format is ascii

        octree  : bool, optional
                  If the data is defined on an octree-like AMR
        """
        if fname == '':
            if binary:
                fname = 'dust_temperature.bdat'
            else:
                fname = 'dust_temperature.dat'

        print('Writing ' + fname)
        if octree:
            self._scalarfieldWriter(data=self.grid.convArrTree2Leaf(self.dusttemp), fname=fname, binary=binary,
                                    octree=True)
        else:
            self._scalarfieldWriter(data=self.dusttemp, fname=fname, binary=binary, octree=False)

    def writeGasDens(self, fname='', ispec='', binary=True, octree=False):
        """Writes the gas density.

        Parameters
        ----------

        fname   : str, optional
                  Name of the file into which the data will be written. If omitted "numberdens_xxx.inp" and
                  "numberdens_xxx.binp" will be used for ascii and binary format, respectively 
                  (xxx is the name of the molecule).

        ispec   : str
                  File name extension of the 'numberdens_ispec.inp' (if binary=True 'numberdens_ispec.binp') 
                  file into which the gas density should be written

        binary  : bool
                  If true the data will be written in binary format, otherwise the file format is ascii

        octree  : bool, optional
                  If the data is defined on an octree-like AMR
        """
        if ispec == '':
            raise ValueError('Unknown ispec. Without knowing the name of the gas species the gas density '
                             + 'cannot be written, as the filename "numberdens_ispec.inp" cannot be generated.')
        else:
            if fname == '':
                if binary:
                    fname = 'numberdens_' + ispec + '.binp'
                else:
                    fname = 'numberdens_' + ispec + '.inp'

            print('Writing ' + fname)
            if octree:
                self._scalarfieldWriter(data=self.grid.convArrTree2Leaf(self.ndens_mol), fname=fname, binary=binary,
                                        octree=True)
            else:
                self._scalarfieldWriter(data=self.ndens_mol, fname=fname, binary=binary, octree=False)

    def writeGasTemp(self, fname='', binary=True, octree=False):
        """Writes the gas temperature.

        Parameters
        ----------

        fname : str, optional 
                Name of the file into which the gas temperature should be written. If omitted 
                'gas_temperature.inp' (if binary=True 'gas_tempearture.binp') is used.

        binary : bool
                If true the data will be written in binary format, otherwise the file format is ascii

        octree  : bool, optional
                  If the data is defined on an octree-like AMR
        """
        if fname == '':
            if binary:
                fname = 'gas_temperature.binp'
            else:
                fname = 'gas_temperature.inp'

        print('Writing ' + fname)

        if octree:
            self._scalarfieldWriter(data=self.grid.convArrTree2Leaf(self.gastemp), fname=fname, binary=binary,
                                    octree=True)
        else:
            self._scalarfieldWriter(data=self.gastemp, fname=fname, binary=binary, octree=False)

    def writeGasVel(self, fname='', binary=True, octree=False):
        """Writes the gas velocity.

        Parameters
        ----------

        fname  : str, optional
                Name of the file into which the gas temperature should be written. 
                If omitted 'gas_velocity.inp' (if binary=True 'gas_velocity.binp') is used.

        binary : bool
                If true the data will be written in binary format, otherwise the file format is ascii

        octree  : bool, optional
                  If the data is defined on an octree-like AMR
        """

        if binary:
            if fname == '':
                fname = 'gas_velocity.binp'

            print('Writing ' + fname)
            with open(fname, 'w') as wfile:

                if octree:
                    #
                    # Check if the gas velocity contains the full tree or only the leaf nodes
                    #
                    if self.gasvel.shape[0] == self.grid.nLeaf:
                        hdr = np.array([1, 8, self.gasvel.shape[0]], dtype=int)
                        hdr.tofile(wfile)
                        self.gasvel.flatten(order='f').tofile(wfile)
                    else:
                        hdr = np.array([1, 8, self.grid.nLeaf], dtype=int)
                        hdr.tofile(wfile)
                        dummy = self.grid.convArrTree2Leaf(self.gasvel)
                        dummy.flatten(order='f').tofile(wfile)
                        # self.gasvel.flatten(order='f').tofile(wfile)

                else:
                    hdr = np.array([1, 8, self.grid.nx * self.grid.ny * self.grid.nz], dtype=int)
                    hdr.tofile(wfile)
                    # Now we need to change the axis orders since the Ndarray.tofile function writes the
                    # array always in C-order while we need Fortran-order to be written
                    self.gasvel = np.swapaxes(self.gasvel, 0, 2)
                    self.gasvel.tofile(wfile)

                    # Switch back to the original axis order
                    self.gasvel = np.swapaxes(self.gasvel, 0, 2)
        else:
            if fname == '':
                fname = 'gas_velocity.inp'

            print('Writing ' + fname)

            # If we have an octree grid
            if octree:
                #
                # Check if the gas velocity contains the full tree or only the leaf nodes
                #
                if self.gasvel.shape[0] == self.grid.nLeaf:
                    hdr = '1\n'
                    hdr += ("%d\n" % self.gasvel.shape[0])
                    try:
                        np.savetxt(fname, self.gasvel, fmt="%.9e %.9e %.9e", header=hdr, comments='')
                    except Exception as e:
                        print(e)
                else:
                    hdr = '1\n'
                    hdr += ("%d\n" % self.grid.nLeaf)
                    dummy = self.grid.convArrTree2Leaf(self.gasvel)
                    try:
                        np.savetxt(fname, dummy, fmt="%.9e %.9e %.9e", header=hdr, comments='')
                    except Exception as e:
                        print(e)

            else:
                # hdr = "1\n"
                # hdr += ('%d'%(self.grid.nx*self.grid.ny*self.grid.nz))
                # np.savetxt(fname, self.gasvel, fmt="%.9 %.9 %.9", header=hdr, comments='')

                with open(fname, 'w') as wfile:

                    wfile.write('%d\n' % 1)
                    wfile.write('%d\n' % (self.grid.nx * self.grid.ny * self.grid.nz))

                    for iz in range(self.grid.nz):
                        for iy in range(self.grid.ny):
                            for ix in range(self.grid.nx):
                                wfile.write("%9e %9e %9e\n" % (self.gasvel[ix, iy, iz, 0], self.gasvel[ix, iy, iz, 1],
                                                               self.gasvel[ix, iy, iz, 2]))

    def writeVTurb(self, fname='', binary=True, octree=False):
        """Writes the microturbulence file.

        Parameters
        ----------

        fname   : str, optional
                  Name of the file into which the turubulent velocity field should be written. 
                  If omitted 'microturbulence.inp' (if binary=True 'microturbuulence.binp') is used.

        binary  : bool
                  If true the data will be written in binary format, otherwise the file format is ascii

        octree  : bool, optional
                  If the data is defined on an octree-like AMR
        """

        if fname == '':
            if binary:
                fname = 'microturbulence.binp'
            else:
                fname = 'microturbulence.inp'

        print('Writing ' + fname)

        if octree:
            self._scalarfieldWriter(data=self.grid.convArrTree2Leaf(self.vturb), fname=fname, binary=binary,
                                    octree=True)
        else:
            self._scalarfieldWriter(data=self.vturb, fname=fname, binary=binary, octree=False)

    def writeVTK(self, vtk_fname='', ddens=False, dtemp=False, idust=None, gdens=False, gvel=False, gtemp=False):
        """Writes physical variables to a legacy vtk file.

        Parameters
        ----------

        vtk_fname : str
                    Name of the file to be written, if not specified 'radmc3d_data.vtk' will be used

        ddens     : bool
                    If set to True the dust density will be written to the vtk file

        dtemp     : bool
                    If set to True the dust temperature will be written to the vtk file

        idust     : list
                    List of indices that specifies which dust component should be written 
                    if not set then the first dust species (zero index) will be used

        gdens     : bool
                    If set to True the gas density will be written to the vtk file

        gtemp     : bool
                    If set to True the gas temperature will be written to the vtk file

        gvel      : bool
                    If set to True the gas velocity will be written to the vtk file
        """

        if isinstance(idust, int):
            idust = [idust]

        if (ddens is True) | (dtemp is True):
            if idust is None:
                msg = 'Unknown dust species. Dust density or temperature is to be written but  ' \
                      + 'it has not been specified for which dust species.'
                raise ValueError(msg)


        if vtk_fname == '':
            vtk_fname = 'radmc3d_data.vtk'
        else:
            vtk_fname = str(vtk_fname)

        #
        # Get the grid
        #

        x = self.grid.xi
        # For the theta axis I leave out the poles
        #  The current cell type is hexahedron and the cells near the pole are
        #    rather tetrahedra than hexahedra and this is not yet implemented
        y = np.array(self.grid.yi[1:self.grid.nyi - 1])
        z = self.grid.zi
        nxi = x.shape[0]
        nyi = y.shape[0]
        nzi = z.shape[0]

        #
        # Gas velocity field (Should be corner centered)
        # TODO
        #  The lines below should be double checked and re-implemented as
        #  the re-mapping of the cell centered velocity field to the cell corners
        #  is physically not correct and very messy...
        #
        if gvel:
            vgas = np.zeros([nxi, nyi, nzi, 3], dtype=np.float64)
            vgas[0:nxi - 1, 0:nyi, 0:nzi - 1, :] = self.gasvel[:, 1:nyi + 1, :, :]
            vgas[nxi - 1, :, :, :] = vgas[nxi - 2, :, :, :]
            vgas[:, nyi - 1, :, :] = vgas[:, nyi - 2, :, :]
            vgas[:, :, nzi - 1, :] = vgas[:, :, nzi - 2, :]

        #
        # Header
        #
        with open(vtk_fname, 'w') as wfile:
            wfile.write('%s\n' % '# vtk DataFile Version 3.0')
            wfile.write('%s\n' % 'RADMC-3D Data')
            wfile.write('%s\n' % 'ASCII')
            wfile.write('%s\n' % 'DATASET UNSTRUCTURED_GRID')

            #
            # Write out the coordinates of the cell corners
            #
            wfile.write('%s\n' % ('POINTS ' + str(nxi * nyi * nzi).strip() + ' double'))
            print('Writing POINTS: ')
            for ix in range(nxi):
                print(ix, nxi)
                for iy in range(nyi):
                    for iz in range(nzi):
                        crd = crd_trans.ctransSph2Cart([x[ix], z[iz], y[iy]])
                        wfile.write('%.9e %9e %9e\n' % (crd[0], crd[1], crd[2]))

                        # ---------------------------------------------------------------------------------------------
                        # Write out the indices of the cell interface mesh that define a
                        # hexahedron (VTK cell type #12)
                        #
                        # The indexing of a hexahedron is as follows
                        #
                        #                  7________6
                        #                 /|      / |
                        #                / |     /  |
                        #               4_------5   |             z ^   ^ y
                        #               |  3____|___2               |  /
                        #               | /     |  /                | /
                        #               |/      | /                 |/
                        #               0-------1                   0-----> x
                        #
                        # ---------------------------------------------------------------------------------------------

            wfile.write('%s %d %d\n' % ('CELLS ', ((nxi - 1) * (nyi - 1) * (nzi - 1)),
                                        ((nxi - 1) * (nyi - 1) * (nzi - 1)) * 9))

            for ix in range(nxi - 1):
                print('Writing CELL COORDINATES: ', ix, self.grid.nxi - 2)
                for iy in range(nyi - 1):
                    for iz in range(nzi - 1):
                        id1 = nzi * nyi * ix + nzi * iy + iz
                        id2 = nzi * nyi * ix + nzi * (iy + 1) + iz
                        id4 = nzi * nyi * ix + nzi * iy + ((iz + 1) % (nzi - 1))
                        id3 = nzi * nyi * ix + nzi * (iy + 1) + ((iz + 1) % (nzi - 1))
                        id5 = nzi * nyi * (ix + 1) + nzi * iy + iz
                        id6 = nzi * nyi * (ix + 1) + nzi * (iy + 1) + iz
                        id7 = nzi * nyi * (ix + 1) + nzi * (iy + 1) + ((iz + 1) % (nzi - 1))
                        id8 = nzi * nyi * (ix + 1) + nzi * iy + ((iz + 1) % (nzi - 1))

                        line = np.array([8, id1, id2, id3, id4, id5, id6, id7, id8])
                        line.tofile(wfile, sep=' ', format='%d')
                        wfile.write('\n')
                        #
                        # Now write out the type of each cell (#12)
                        #
            wfile.write('%s %d\n' % ('CELL_TYPES', ((nxi - 1) * (nyi - 1) * (nzi - 1))))

            for ix in range(nxi - 1):
                for iy in range(nyi - 1):
                    for iz in range(nzi - 1):
                        wfile.write('%d\n' % 12)
                        #
                        # Now write out the corner centered velocities
                        #

            if gvel:
                wfile.write('%s %d\n' % ('POINT_DATA', (nxi * nyi * nzi)))
                wfile.write('%s\n' % 'VECTORS gas_velocity double')
                for ix in range(nxi):
                    print('Writing velocity : ', ix, nxi - 1)
                    for iy in range(nyi):
                        for iz in range(nzi):
                            vsph = np.array([vgas[ix, iy, iz, 0], vgas[ix, iy, iz, 2], vgas[ix, iy, iz, 1]])
                            vxyz = crd_trans.vtransSph2Cart([x[ix], z[iz], y[iy]], vsph)

                            wfile.write('%.9e %.9e %.9e\n' % (vxyz[0], vxyz[1], vxyz[2]))

            #
            # Write out the cell centered scalars
            #
            wfile.write('%s %d\n' % ('CELL_DATA', ((nxi - 1) * (nyi - 1) * (nzi - 1))))

            if ddens:
                for ids in idust:
                    wfile.write('%s\n' % ('SCALARS dust_density_' + str(int(ids)) + ' double'))
                    wfile.write('%s\n' % 'LOOKUP_TABLE default')

                    for ix in range(nxi - 1):
                        print('Writing dust density : ', ix, nxi - 2)
                        for iy in range(nyi - 1):
                            for iz in range(nzi - 1):
                                wfile.write('%.9e\n' % self.rhodust[ix, iy, iz, ids])

            if dtemp:
                for ids in idust:
                    wfile.write('%s\n' % ('SCALARS dust_temperature_' + str(int(ids)) + ' double'))
                    wfile.write('%s\n' % 'LOOKUP_TABLE default')

                    for ix in range(nxi - 1):
                        print('Writing dust temperature : ', ix, nxi - 2)
                        for iy in range(nyi - 1):
                            for iz in range(nzi - 1):
                                wfile.write('%.9e\n' % self.dusttemp[ix, iy, iz, ids])

            if gdens:
                wfile.write('%s\n' % 'SCALARS gas_numberdensity double')
                wfile.write('%s\n' % 'LOOKUP_TABLE default')

                for ix in range(nxi - 1):
                    print('Writing gas density : ', ix, nxi - 2)
                    for iy in range(nyi - 1):
                        for iz in range(nzi - 1):
                            wfile.write('%.9e\n' % self.ndens_mol[ix, iy, iz])

            if gtemp:
                wfile.write('%s\n' % 'SCALARS gas_temperature double')
                wfile.write('%s\n' % 'LOOKUP_TABLE default')

                for ix in range(nxi - 1):
                    print('Writing dust temperature : ', ix, nxi - 2)
                    for iy in range(nyi - 1):
                        for iz in range(nzi - 1):
                            wfile.write('%.9e\n' % self.gastemp[ix, iy, iz])

    def getSigmaDust(self, idust=-1):
        """Calculates the dust surface density.

        Parameters
        ----------

        idust : int, optional
                Index of the dust species for which the surface density should be calculated 
                if omitted the calculated surface density will be the sum over all dust species
        """

        # Calculate the volume of each grid cell
        vol = self.grid.getCellVolume()
        # Dustmass in each grid cell
        if len(self.rhodust) > 3:
            if idust >= 0:
                mass = vol * self.rhodust[:, :, :, idust]
            else:
                mass = vol * self.rhodust.sum(3)
        else:
            mass = vol * self.rhodust

        # Calculate the surface of each grid facet in the midplane
        surf = np.zeros([self.grid.nx, self.grid.nz], dtype=np.float64)
        diff_r2 = (self.grid.xi[1:] ** 2 - self.grid.xi[:-1] ** 2) * 0.5
        diff_phi = self.grid.zi[1:] - self.grid.zi[:-1]
        for ix in range(self.grid.nx):
            surf[ix, :] = diff_r2[ix] * diff_phi

        # Now get the surface density
        dum = np.squeeze(mass.sum(1))
        self.sigmadust = dum / np.squeeze(surf)

    def getSigmaGas(self):
        """Calculates the gas surface density.
        This method uses radmc3dData.rhogas to calculate the surface density, thus the 
        unit of surface density depends on the unit of radmc3dData.rhogas (g/cm^2 or molecule/cm^2)
        """

        # Calculate the volume of each grid cell
        vol = self.grid.getCellVolume()
        # Total number of molecules / gas mass in each grid cell
        mass = vol * self.rhogas
        # Calculate the surface are of each grid facet in the midplane
        surf = np.zeros([self.grid.nx, self.grid.nz], dtype=np.float64)
        diff_r2 = (self.grid.xi[1:] ** 2 - self.grid.xi[:-1] ** 2) * 0.5
        diff_phi = self.grid.zi[1:] - self.grid.zi[:-1]
        for ix in range(self.grid.nx):
            surf[ix, :] = diff_r2[ix] * diff_phi

        # Now get the surface density
        dum = np.squeeze(mass.sum(1))
        self.sigmagas = dum / np.squeeze(surf)


class radmc3dGrid(object):
    """ Class for spatial and frequency grid used by RADMC-3D.

    Attributes
    ----------

    act_dim    : ndarray 
                A three element vector the i-th element is 1 if the i-th dimension is active, 
                otherwize the i-th element is zero

    crd_sys    : {'sph', 'cyl', 'car'} 
                coordinate system of the spatial grid

    nx         : int
                Number of grid points in the x (cartesian) / r (cylindrical) / r (spherical) dimension

    ny         : int
                Number of grid points in the y (cartesian) / theta (cylindrical) / theta (spherical) dimension

    nz         : int
                Number of grid points in the z (cartesian) / z (cylindrical) / phi (spherical) dimension

    nxi        : int
                Number of cell interfaces in the x (cartesian) / r (cylindrical) / r (spherical) dimension

    nyi        : int
                Number of cell interfaces in the y (cartesian) / theta (cylindrical) / theta (spherical) dimension

    nzi        : int
                Number of cell interfaces in the z (cartesian) / z (cylindrical) / phi (spherical) dimension

    nwav       : int
                Number of wavelengths in the wavelength grid

    nfreq      : int
                Number of frequencies in the grid (equal to nwav)

    x          : ndarray
                Cell centered x (cartesian) / r (cylindrical) / r (spherical)  grid points

    y          : ndarray
                Cell centered y (cartesian) / theta (cylindrical) / theta (spherical)  grid points

    z          : ndarray
                Cell centered z (cartesian) / z (cylindrical) / phi (spherical)  grid points

    xi         : ndarray
                Cell interfaces in the x (cartesian) / r (cylindrical) / r (spherical)  dimension

    yi         : ndarray
                Cell interfaces in the y (cartesian) / theta (cylindrical) / theta (spherical)  dimension

    zi         : ndarray
                Cell interfaces in the z (cartesian) / z (cylindrical) / phi (spherical)  dimension

    wav        : ndarray
                Wavelengh  grid

    freq       : ndarray
                Frequency  grid


    """

    # --------------------------------------------------------------------------------------------------

    def __init__(self):

        self.crd_sys = 'sph'
        self.act_dim = [1, 1, 1]
        self.grid_style = 0
        self.nx = 0
        self.ny = 0
        self.nz = 0
        self.nxi = 0
        self.nyi = 0
        self.nzi = 0
        self.x = np.zeros(0, dtype=np.float64)
        self.y = np.zeros(0, dtype=np.float64)
        self.z = np.zeros(0, dtype=np.float64)
        self.xi = np.zeros(0, dtype=np.float64)
        self.yi = np.zeros(0, dtype=np.float64)
        self.zi = np.zeros(0, dtype=np.float64)

        self.nwav = 0
        self.nfreq = 0
        self.wav = np.zeros(0, dtype=np.float64)
        self.freq = np.zeros(0, dtype=np.float64)

    # --------------------------------------------------------------------------------------------------
    def makeWavelengthGrid(self, wbound=None, nw=None, ppar=None):
        """Creates the wavelength/frequency grid.

        Parameters
        ----------

        wbound : list 
                 Contains the wavelength boundaries of the wavelength grid (should contain at least two elements)

        nw     : list 
                 Contains len(wbound)-1 elements containing the number of wavelengths between the bounds
                 set by wbound

        ppar   : dictionary, optional
                 Contains all input parameters with the parameter names as keys  
        """
        if ppar:
            if wbound is None:
                wbound = ppar['wbound']
            if nw is None:
                nw = ppar['nw']
        else:
            if wbound is None:
                raise ValueError('Unknown wbound. Without the grid boundaries the wavelength grid cannot be generated')
            if nw is None:
                raise ValueError('Unknown nw. Without the number of grid points the wavelength grid cannot be '
                                 + 'generated')

        self.nwav = nw[0]
        self.wav = wbound[0] * (wbound[1] / wbound[0]) ** (np.arange(nw[0], dtype=np.float64) / nw[0])

        for ipart in range(1, len(nw) - 1):
            dum = wbound[ipart] * (wbound[ipart + 1]
                                   / wbound[ipart])**(np.arange(nw[ipart], dtype=np.float64) / nw[ipart])
            self.wav = np.append(self.wav, dum)

        ipart = len(nw) - 1
        dum = wbound[ipart] * (wbound[ipart + 1]
                               / wbound[ipart])**(np.arange(nw[ipart], dtype=np.float64) / (nw[ipart] - 1.))
        self.wav = np.append(self.wav, dum)
        self.nwav = self.wav.shape[0]
        self.freq = nc.cc / self.wav * 1e4
        self.nfreq = self.nwav

    # --------------------------------------------------------------------------------------------------
    def writeWavelengthGrid(self, fname='', old=False):
        """Wriites the wavelength grid to a file (e.g. wavelength_micron.inp).

        Parameters
        ----------

        fname  : str, optional
                 File name into which the wavelength grid should be written. If omitted 'wavelength_micron.inp' 
                 will be used

        old    : bool, optional
                 If set to True the file format of the previous, 2D version of radmc will be used
        """

        if not old:
            if fname == '':
                fname = 'wavelength_micron.inp'

            print('Writing ' + fname)
            with open(fname, 'w') as wfile:
                wfile.write('%d\n' % self.nwav)
                for ilam in range(self.nwav):
                    wfile.write('%.9e\n' % self.wav[ilam])
        else:
            if fname == '':
                fname = 'frequency.inp'

            with open(fname, 'w') as wfile:

                print('Writing ' + fname)
                wfile.write("%d\n" % self.nfreq)
                wfile.write(" \n")
                #
                # Reverse the order of the frequency grid as it is ordered in frequency in radmc
                #
                freq = self.freq[::-1]
                for i in range(self.nfreq):
                    wfile.write("%.7e\n" % freq[i])

            # --------------------------------------------------------------------------------------------------

    def makeSpatialGrid(self, crd_sys=None, xbound=None, ybound=None, zbound=None, nxi=None, nyi=None, nzi=None,
                        ppar=None):
        """Calculates the spatial grid.

        Parameters
        ----------

        crd_sys : {'sph','car'}  
                    Coordinate system of the spatial grid

        xbound  : list 
                    (with at least two elements) of boundaries for the grid along the first dimension

        ybound  : list 
                    (with at least two elements) of boundaries for the grid along the second dimension

        zbound  : list 
                    (with at least two elements) of boundaries for the grid along the third dimension

        nxi     : int
                    Number of grid points along the first dimension. List with len(xbound)-1 elements with 
                    nxi[i] being the number of grid points between xbound[i] and xbound[i+1]

        nyi     : int
                    Same as nxi but for the second dimension

        nzi     : int
                    Same as nxi but for the third dimension

        ppar    : Dictionary containing all input parameters of the model (from the problem_params.inp file)
                   if ppar is set all keyword arguments that are not set will be taken from this dictionary
        """

        self.act_dim = [1, 1, 1]
        if ppar:
            if not crd_sys:
                crd_sys = ppar['crd_sys']
            self.crd_sys = crd_sys

            if not xbound:
                if 'xbound' in ppar:
                    xbound = ppar['xbound']
                else:
                    print(' No boundary for the first dimension is given, first dimension is deactivated.')
                    self.act_dim[0] = 0
            if not nxi:
                if 'nx' in ppar:
                    if not isinstance(ppar['nx'], list):
                        ppar['nx'] = [ppar['nx']]
                    nxi = [i + 1 for i in ppar['nx']]
                    if ppar['nx'][0] == 0:
                        self.act_dim[0] = 0
                else:
                    self.act_dim[0] = 0

            if not ybound:
                if 'ybound' in ppar:
                    ybound = ppar['ybound']
                else:
                    print(' No boundary for the second dimension is given, second dimension is deactivated.')
                    self.act_dim[1] = 0
            if not nyi:
                if 'ny' in ppar:
                    if not isinstance(ppar['ny'], list):
                        nyi = [ppar['ny'] + 1]
                        ppar['ny'] = [ppar['ny']]
                    else:
                        ppar['ny'] = ppar['ny']
                        nyi = [i + 1 for i in ppar['ny']]

                    if ppar['ny'][0] == 0:
                        self.act_dim[1] = 0
                else:
                    self.act_dim[1] = 0

            if not zbound:
                if 'zbound' in ppar:
                    zbound = ppar['zbound']
                else:
                    print(' No boundary for the third dimension is given, third dimension is deactivated.')
                    self.act_dim[2] = 0
            if not nzi:
                if 'nz' in ppar:
                    if not isinstance(ppar['nz'], list):
                        ppar['nz'] = [ppar['nz']]
                    if ppar['nz'][0] >  0:
                        nzi = [i + 1 for i in ppar['nz']]
                    if ppar['nz'][0] == 0:
                        self.act_dim[2] = 0
                else:
                    self.act_dim[2] = 0
                    nzi = [0]
        #
        # Type checking
        #
        if not isinstance(nxi, list):
            nxi = [nxi]
        if not isinstance(nyi, list):
            nyi = [nyi]
        if not isinstance(nzi, list):
            nzi = [nzi]

        if crd_sys == 'car':
            #
            # First check whether the grid boundaries are specified
            #
            if xbound is None:
                raise ValueError('Unknown xbound. Boundaries for the cartesian x-axis are not specified. '
                                 + 'Without the boundaries the grid cannot be generated')

            if ybound is None:
                raise ValueError('Unknown ybound. Boundaries for the cartesian y-axis are not specified. '
                                 + 'Without the boundaries the grid cannot be generated')
            if zbound is None:
                raise ValueError('Unknown zbound. Boundaries for the cartesian z-axis are not specified. '
                                 + 'Without the boundaries the grid cannot be generated')

            if nxi is None:
                raise ValueError('Unknown nxi. Number of grid points for the cartesian x-axis is not specified. '
                                 + 'The grid cannot be generated')

            if nyi is None:
                raise ValueError('Unknown nyi. Number of grid points for the cartesian y-axis is not specified. '
                                 + 'The grid cannot be generated')

            if nzi is None:
                raise ValueError('Unknown nzi. Number of grid points for the cartesian z-axis is not specified. '
                                 + 'The grid cannot be generated')

            #
            # Create the x-axis
            #
            if len(nxi) > 1:
                self.nxi = sum(nxi)
                self.nx = self.nxi - 1
                self.xi = xbound[0] + (xbound[1] - xbound[0]) * (np.arange(nxi[0], dtype=np.float64) / float(nxi[0]))
                for ipart in range(1, len(nxi) - 1):
                    dum = xbound[ipart] + (xbound[ipart + 1] - xbound[ipart]) \
                                          * (np.arange(nxi[ipart], dtype=np.float64) / float(nxi[ipart]))
                    self.xi = np.append(self.xi, dum)

                ipart = len(nxi) - 1
                dum = xbound[ipart] + (xbound[ipart + 1] - xbound[ipart]) * (
                    np.arange(nxi[ipart], dtype=np.float64) / float(nxi[ipart] - 1))
                self.xi = np.append(self.xi, dum)
                self.x = 0.5 * (self.xi[0:self.nx] + self.xi[1:self.nx + 1])
            else:
                if self.act_dim[0] == 1:
                    self.nxi = nxi[0]
                    self.xi = xbound[0] + (xbound[1] - xbound[0]) * (
                        np.arange(self.nxi, dtype=np.float64) / float(self.nxi - 1.))
                    self.nx = self.nxi - 1
                    self.x = 0.5 * (self.xi[0:self.nx] + self.xi[1:self.nx + 1])
                else:
                    self.x = [0.]
                    self.xi = [0., 0., ]
                    self.nx = 1
                    self.nxi = 2

                    #
                    # Create the y-ayis
                    #
            if len(nyi) > 1:
                self.nyi = sum(nyi)
                self.ny = self.nyi - 1
                self.yi = ybound[0] + (ybound[1] - ybound[0]) * (np.arange(nyi[0], dtype=np.float64) / float(nyi[0]))
                for ipart in range(1, len(nyi) - 1):
                    dum = ybound[ipart] + (ybound[ipart + 1] - ybound[ipart]) * (
                        np.arange(nyi[ipart], dtype=np.float64) / float(nyi[ipart]))
                    self.yi = np.append(self.yi, dum)

                ipart = len(nyi) - 1
                dum = ybound[ipart] + (ybound[ipart + 1] - ybound[ipart]) * (
                    np.arange(nyi[ipart], dtype=np.float64) / float(nyi[ipart] - 1))
                self.yi = np.append(self.yi, dum)
                self.y = 0.5 * (self.yi[0:self.ny] + self.yi[1:self.ny + 1])
            else:
                if self.act_dim[0] == 1:
                    self.nyi = nyi[0]
                    self.yi = ybound[0] + (ybound[1] - ybound[0]) * (
                        np.arange(self.nyi, dtype=np.float64) / float(self.nyi - 1.))
                    self.ny = self.nyi - 1
                    self.y = 0.5 * (self.yi[0:self.ny] + self.yi[1:self.ny + 1])
                else:
                    self.y = [0.]
                    self.yi = [0., 0., ]
                    self.ny = 1
                    self.nyi = 2

                    #
                    # Create the z-azis
                    #
            if len(nzi) > 1:
                self.nzi = sum(nzi)
                self.nz = self.nzi - 1
                self.zi = zbound[0] + (zbound[1] - zbound[0]) * (np.arange(nzi[0], dtype=np.float64) / float(nzi[0]))
                for ipart in range(1, len(nzi) - 1):
                    dum = zbound[ipart] + (zbound[ipart + 1] - zbound[ipart]) \
                                          * (np.arange(nzi[ipart], dtype=np.float64) / float(nzi[ipart]))
                    self.zi = np.append(self.zi, dum)

                ipart = len(nzi) - 1
                dum = zbound[ipart] + (zbound[ipart + 1] - zbound[ipart]) * (
                    np.arange(nzi[ipart], dtype=np.float64) / float(nzi[ipart] - 1))
                self.zi = np.append(self.zi, dum)
                self.z = 0.5 * (self.zi[0:self.nz] + self.zi[1:self.nz + 1])
            else:
                if self.act_dim[0] == 1:
                    self.nzi = nzi[0]
                    self.zi = zbound[0] + (zbound[1] - zbound[0]) * (
                        np.arange(self.nzi, dtype=np.float64) / float(self.nzi - 1.))
                    self.nz = self.nzi - 1
                    self.z = 0.5 * (self.zi[0:self.nz] + self.zi[1:self.nz + 1])
                else:
                    self.z = [0.]
                    self.zi = [0., 0.0]
                    self.nz = 1
                    self.nzi = 2

        if crd_sys == 'sph':
            #
            # r->x, theta->y, phi-z
            #
            if xbound is None:
                raise ValueError('Unknown xbound. Boundaries for the spherical radial grid are not specified. '
                                 + 'Without the boundaries the grid cannot be generated')

            if ybound is None:
                print('Unknown ybound. Setting the spherical co-lattitude grid to [0, pi]')
                ybound = [0.0, np.pi]

            if zbound is None:
                print('Unknown zbound. Setting the spherical azimuth angle grid to [0, 2*pi]')
                zbound = [0.0, 2.0*np.pi]

            if nxi is None:
                raise ValueError('Unknown nxi. Number of grid points for the spherical radial grid is not specified. '
                                 + 'The grid cannot be generated')

            if nyi is None:
                raise ValueError('Unknown nyi. Number of grid points for the spherical co-lattitude grid is not '
                                 + ' specified. The grid cannot be generated')

            if nzi is None:
                raise ValueError('Unknown nzi. Number of grid points for the spherical azimuthal angle grid is not '
                                 + 'specified. The grid cannot be generated')

            #
            # Create the x axis
            #
            if len(nxi) > 1:
                self.nxi = sum(nxi)
                self.nx = self.nxi - 1
                self.xi = xbound[0] * (xbound[1] / xbound[0])**(np.arange(nxi[0], dtype=np.float64) / float(nxi[0]))
                for ipart in range(1, len(nxi) - 1):
                    dum = xbound[ipart] * (xbound[ipart + 1] / xbound[ipart])**(
                        np.arange(nxi[ipart], dtype=np.float64) / float(nxi[ipart]))
                    self.xi = np.append(self.xi, dum)

                ipart = len(nxi) - 1
                dum = xbound[ipart] * (xbound[ipart + 1] / xbound[ipart])**(
                    np.arange(nxi[ipart], dtype=np.float64) / float(nxi[ipart] - 1))
                self.xi = np.append(self.xi, dum)
                self.x = np.sqrt(self.xi[0:self.nx] * self.xi[1:self.nx + 1])
            else:
                if self.act_dim[0] == 1:
                    self.nxi = nxi[0]
                    self.xi = xbound[0] * (xbound[1] / xbound[0])**(
                        np.arange(self.nxi, dtype=np.float64) / float(self.nxi - 1.))
                    self.nx = self.nxi - 1
                    self.x = np.sqrt(self.xi[0:self.nx] * self.xi[1:self.nx + 1])
                else:
                    self.x = [0.]
                    self.xi = [0., 0., ]
                    self.nx = 1
                    self.nxi = 2


            # Refinement of the inner edge of the grid
            # This has to be done properly
            if 'xres_nlev' in ppar:
                if ppar['xres_nlev'] > 0:
                    ri_ext = np.array([self.xi[0], self.xi[ppar['xres_nspan']]])
                    for i in range(ppar['xres_nlev']):
                        dum_ri = ri_ext[0] + (ri_ext[1] - ri_ext[0]) * np.arange(ppar['xres_nstep'] + 1,
                                                                                 dtype=np.float64) / float(
                            ppar['xres_nstep'])
                        # print ri_ext[0:2]/au
                        # print dum_ri/au
                        ri_ext_old = np.array(ri_ext)
                        ri_ext = np.array(dum_ri)
                        ri_ext = np.append(ri_ext, ri_ext_old[2:])

                    r_ext = (ri_ext[1:] + ri_ext[:-1]) * 0.5

                    self.xi = np.append(ri_ext, self.xi[ppar['xres_nspan'] + 1:])
                    self.x = np.append(r_ext, self.x[ppar['xres_nspan']:])
                    self.nx = self.x.shape[0]
                    self.nxi = self.xi.shape[0]

                    #
                    # Create the y axis
                    #
            if len(nyi) > 1:

                # Check if we go to the full [0,pi] interval or only use the upper half-plane [0, pi/2]

                if ybound[len(ybound) - 1] != np.pi / 2.:
                    self.nyi = sum(nyi) + 1
                    self.ny = self.nyi - 1
                    self.yi = ybound[0] + (ybound[1] - ybound[0]) * (
                        np.arange(nyi[0], dtype=np.float64) / float(nyi[0]))

                    for ipart in range(1, len(nyi) - 1):
                        # Now make sure that pi/2 will be a cell interface
                        #
                        # BUGFIX! 16-05-2012
                        # The grid was not symmetric to pi/2 when the grid contained multiple sections (i.e. len(nyi)>1)
                        # This is now fixed
                        if ybound[ipart] < np.pi / 2.:
                            dum = ybound[ipart] + (ybound[ipart + 1] - ybound[ipart]) * (
                                np.arange(nyi[ipart], dtype=np.float64) / float(nyi[ipart]))
                        else:
                            if ybound[ipart] == np.pi / 2.:
                                dum = ybound[ipart] + (ybound[ipart + 1] - ybound[ipart]) * (
                                    (np.arange(nyi[ipart] + 1, dtype=np.float64)) / (float(nyi[ipart])))
                            else:
                                dum = ybound[ipart] + (ybound[ipart + 1] - ybound[ipart]) * (
                                    (np.arange(nyi[ipart], dtype=np.float64) + 1.) / float(nyi[ipart]))

                        self.yi = np.append(self.yi, dum)

                    ipart = len(nyi) - 1
                    if len(nyi) == 2:
                        dum = ybound[ipart] + (ybound[ipart + 1] - ybound[ipart]) * (
                            (np.arange(nyi[ipart] + 1, dtype=np.float64)) / (float(nyi[ipart])))
                    else:
                        dum = ybound[ipart] + (ybound[ipart + 1] - ybound[ipart]) * (
                            (np.arange(nyi[ipart], dtype=np.float64) + 1.) / float(nyi[ipart]))

                else:
                    self.nyi = sum(nyi) + 1
                    self.ny = self.nyi - 1
                    self.yi = ybound[0] + (ybound[1] - ybound[0]) * (
                            np.arange(nyi[0], dtype=np.float64) / float(nyi[0]))
                    for ipart in range(1, len(nyi) - 1):
                        # Now make sure that pi/2 will be a cell interface
                        #
                        # BUGFIX! 16-05-2012
                        # The grid was not symmetric to pi/2 when the grid contained multiple sections (i.e. len(nyi)>1)
                        # This is now fixed
                        if ybound[ipart] < np.pi / 2.:
                            dum = ybound[ipart] + (ybound[ipart + 1] - ybound[ipart]) * (
                                np.arange(nyi[ipart], dtype=np.float64) / float(nyi[ipart]))
                        else:
                            dum = ybound[ipart] + (ybound[ipart + 1] - ybound[ipart]) * (
                                (np.arange(nyi[ipart] + 1, dtype=np.float64)) / (float(nyi[ipart])))

                        self.yi = np.append(self.yi, dum)

                    ipart = len(nyi) - 1

                    if len(nyi) == 2:
                        dum = ybound[ipart] + (ybound[ipart + 1] - ybound[ipart]) * (
                            (np.arange(nyi[ipart] + 1, dtype=np.float64)) / (float(nyi[ipart])))
                    else:
                        dum = ybound[ipart] + (ybound[ipart + 1] - ybound[ipart]) * (
                            (np.arange(nyi[ipart], dtype=np.float64) + 1.) / float(nyi[ipart]))

                self.yi = np.append(self.yi, dum)
                self.y = 0.5 * (self.yi[0:self.ny] + self.yi[1:self.ny + 1])

            else:
                if self.act_dim[1] == 1:
                    self.nyi = nyi[0]
                    self.yi = ybound[0] + (ybound[1] - ybound[0]) * (
                        np.arange(self.nyi, dtype=np.float64) / float(self.nyi - 1.))
                    self.ny = self.nyi - 1
                    self.y = 0.5 * (self.yi[0:self.ny] + self.yi[1:self.ny + 1])
                else:
                    self.y = [0.]
                    self.yi = [0., 0., ]
                    self.ny = 1
                    self.nyi = 2
                    #
                    # Create the z axis

            if len(nzi) > 1:
                self.nzi = sum(nzi)
                self.nz = self.nzi - 1

                self.zi = zbound[0] + (zbound[1] - zbound[0]) * (np.arange(nzi[0], dtype=np.float64) / float(nzi[0]))
                for ipart in range(1, len(nzi) - 1):
                    dum = zbound[ipart] + (zbound[ipart + 1] - zbound[ipart]) * (
                        np.arange(nzi[ipart], dtype=np.float64) / float(nzi[ipart]))
                    self.zi = np.append(self.zi, dum)
                ipart = len(nzi) - 1
                dum = zbound[ipart] + (zbound[ipart + 1] - zbound[ipart]) * (
                    np.arange(nzi[ipart], dtype=np.float64) / float(nzi[ipart] - 1))
                self.zi = np.append(self.zi, dum)
                self.z = 0.5 * (self.zi[0:self.nz] + self.zi[1:self.nz + 1])
            else:
                if self.act_dim[2] == 1:
                    self.nzi = nzi[0]
                    self.zi = zbound[0] + (zbound[1] - zbound[0]) * (
                        np.arange(self.nzi, dtype=np.float64) / float(self.nzi - 1))
                    self.nz = self.nzi - 1
                    self.z = 0.5 * (self.zi[0:self.nz] + self.zi[1:self.nz + 1])
                else:
                    self.z = np.array([0.])
                    self.zi = np.array([0., np.pi * 2.])
                    self.nz = 1
                    self.nzi = 2

    def writeSpatialGrid(self, fname='', old=False):
        """Writes the wavelength grid to a file (e.g. amr_grid.inp).

        Parameters
        ----------

        fname : str, optional
                File name into which the spatial grid should be written. If omitted 'amr_grid.inp' will be used. 

        old   : bool, optional
                If set to True the file format of the previous, 2D version of radmc will be used
        """

        #
        # Write the spatial grid for radmc3d
        #
        if not old:
            if fname == '':
                fname = 'amr_grid.inp'

            print('Writing ' + fname)
            with open(fname, 'w') as wfile:
                # Format number
                wfile.write('%d\n' % 1)
                # AMR self.style (0=regular self. NO AMR)
                wfile.write('%d\n' % 0)
                # Coordinate system (0-99 cartesian, 100-199 spherical, 200-299 cylindrical)
                if self.crd_sys == 'car':
                    wfile.write('%d\n' % 0)
                # Coordinate system (0-99 cartesian, 100-199 spherical, 200-299 cylindrical)
                if self.crd_sys == 'sph':
                    wfile.write('%d\n' % 100)
                # Coordinate system (0-99 cartesian, 100-199 spherical, 200-299 cylindrical)
                if self.crd_sys == 'cyl':
                    wfile.write('%d\n' % 200)
                # Gridinfo
                wfile.write('%d\n' % 0)

                # Active dimensions
                wfile.write('%d %d %d \n' % (self.act_dim[0], self.act_dim[1], self.act_dim[2]))
                # Grid size (x,y,z or r,phi,theta, or r,phi,z)
                wfile.write('%d %d %d \n' % (self.nx, self.ny, self.nz))
                for i in range(self.nxi):
                    wfile.write('%.9e\n' % self.xi[i])
                for i in range(self.nyi):
                    wfile.write('%.9e\n' % self.yi[i])
                for i in range(self.nzi):
                    wfile.write('%.9e\n' % self.zi[i])
            wfile.close()
        #
        # Write the spatial grid for radmc
        #
        else:

            fname = 'radius.inp'
            with open(fname, 'w') as wfile:

                print('Writing ' + fname)
                x = np.sqrt(self.xi[1:] * self.xi[:-1])
                wfile.write("%d\n" % self.nx)
                wfile.write(" \n")
                for i in range(self.nx):
                    wfile.write("%.7e\n" % x[i])

            fname = 'theta.inp'
            with open(fname, 'w') as wfile:
                print('Writing ' + fname)
                wfile.write("%d 1\n" % (self.ny / 2))
                wfile.write(" \n")
                for i in range(int(self.ny / 2)):
                    wfile.write("%.7e\n" % self.y[i])

    def readWavelengthGrid(self, fname='', old=False):
        """Reads the wavelength grid 

        Parameters
        ----------

        fname : str, optional
                File name from which the spatial grid should be read. If omitted 'wavelength_micron.inp' will be used. 

        old   : bool, optional
                If set to True the file format of the previous, 2D version of radmc will be used
        """
        #
        # Read the radmc3d format
        #
        if not old:
            if fname == '':
                fname = 'wavelength_micron.inp'
            #
            # Read the frequency grid
            #
            print('Reading '+fname)
            with open(fname, 'r') as rfile:

                self.nwav = int(rfile.readline())
                self.nfreq = self.nwav
                self.wav = np.zeros(self.nwav, dtype=np.float64)

                for i in range(self.nwav):
                    self.wav[i] = float(rfile.readline())

                self.freq = nc.cc / self.wav * 1e4

        #
        # Read the old radmc format
        #
        else:
            if fname == '':
                fname = 'frequency.inp'

            with open(fname, 'r') as rfile:

                self.nfreq = int(rfile.readline())
                self.nwav = self.nfreq
                dum = rfile.readline()
                self.freq = np.zeros(self.nfreq, dtype=float)
                self.wav = np.zeros(self.nfreq, dtype=float)
                for i in range(self.nfreq):
                    self.freq[i] = float(rfile.readline())
                    self.wav[i] = nc.cc / self.freq[i] * 1e4

        return

    def readSpatialGrid(self, fname='', old=False):
        """Reads the spatial grid 

        Parameters
        ----------

        fname : str, optional
                File name from which the spatial grid should be read. If omitted 'amr_grid.inp' will be used. 

        old   : bool, optional
                If set to True the file format of the previous, 2D version of radmc will be used
        """
        #
        # Read the radmc3d format
        #
        if not old:
            if fname == '':
                fname = 'amr_grid.inp'
                #
                # Read the spatial grid

            with open(fname, 'r') as rfile:
                form = float(rfile.readline())
                grid_style = float(rfile.readline())
                crd_system = int(rfile.readline())
                if crd_system < 100:
                    self.crd_sys = 'car'
                elif (crd_system >= 100) & (crd_system < 200):
                    self.crd_sys = 'sph'
                elif (crd_system >= 200) & (crd_system < 300):
                    self.crd_sys = 'cyl'
                else:
                    rfile.close()
                    raise ValueError('Unsupported coordinate system identification in the ' + fname + ' file.')

                grid_info = float(rfile.readline())
                dum = rfile.readline().split()
                self.act_dim = [int(dum[i]) for i in range(len(dum))]
                dum = rfile.readline().split()
                self.nx, self.ny, self.nz = int(dum[0]), int(dum[1]), int(dum[2])
                self.nxi, self.nyi, self.nzi = self.nx + 1, self.ny + 1, self.nz + 1

                self.xi = np.zeros(self.nx + 1, dtype=np.float64)
                self.yi = np.zeros(self.ny + 1, dtype=np.float64)
                self.zi = np.zeros(self.nz + 1, dtype=np.float64)

                for i in range(self.nxi):
                    self.xi[i] = float(rfile.readline())
                for i in range(self.nyi):
                    self.yi[i] = float(rfile.readline())
                for i in range(self.nzi):
                    self.zi[i] = float(rfile.readline())

                if self.crd_sys == 'car':
                    self.x = (self.xi[0:self.nx] + self.xi[1:self.nx + 1]) * 0.5
                    self.y = (self.yi[0:self.ny] + self.yi[1:self.ny + 1]) * 0.5
                    self.z = (self.zi[0:self.nz] + self.zi[1:self.nz + 1]) * 0.5
                else:
                    self.x = np.sqrt(self.xi[0:self.nx] * self.xi[1:self.nx + 1])
                    self.y = (self.yi[0:self.ny] + self.yi[1:self.ny + 1]) * 0.5
                    self.z = (self.zi[0:self.nz] + self.zi[1:self.nz + 1]) * 0.5

        #
        # Read the old radmc format
        #
        else:
            self.crd_sys = 'sph'
            self.act_dim = [1, 1, 0]

            #
            # Read the radial grid
            #
            with open('radius.inp', 'r') as rfile:

                self.nx = int(rfile.readline())
                self.nxi = self.nx + 1
                dum = rfile.readline()
                self.x = np.zeros(self.nx, dtype=float)
                self.xi = np.zeros(self.nxi, dtype=float)
                for i in range(self.nx):
                    self.x[i] = float(rfile.readline())
                self.xi[1:-1] = 0.5 * (self.x[1:] + self.x[:-1])
                self.xi[0] = self.x[0] - (self.xi[1] - self.x[0])
                self.xi[-1] = self.x[-1] + (self.x[-1] - self.xi[-2])

            #
            # Read the poloidal angular grid
            #
            with open('theta.inp', 'r') as rfile:

                ny = int(rfile.readline().split()[0])
                self.ny = ny * 2
                self.nyi = self.ny + 1
                dum = rfile.readline()

                self.y = np.zeros(self.ny, dtype=float)
                self.yi = np.zeros(self.nyi, dtype=float)
                self.yi[0] = 0.
                self.yi[-1] = np.pi
                self.yi[ny] = np.pi * 0.5

                for i in range(int(self.ny / 2)):
                    self.y[i] = float(rfile.readline())
                    self.y[self.ny - 1 - i] = np.pi - self.y[i]

                self.yi[1:-1] = 0.5 * (self.y[1:] + self.y[:-1])
                self.yi[ny] = np.pi * 0.5

            #
            # Create the azimuthal grid
            #

            self.nz = 1
            self.zi = np.array([0., 2. * np.pi], dtype=float)

        return

    def readGrid(self, old=False):
        """Reads the spatial (amr_grid.inp) and frequency grid (wavelength_micron.inp).

        Parameters
        ----------

        old   : bool, optional
                If set to True the file format of the previous, 2D version of radmc will be used
        """

        self.readSpatialGrid(old=old)
        self.readWavelengthGrid(old=old)

        return

    def getCellVolume(self):
        """Calculates the volume of grid cells.

        """
        if self.crd_sys == 'sph':

            if self.act_dim[0] == 0:
                raise ValueError('The first (r) dimension of a shserical grid is switched off')
            elif self.act_dim[1] == 0:
                if self.act_dim[2] == 0:
                    vol = np.zeros([self.nx, self.ny, self.nz], dtype=np.float64)
                    diff_r3 = self.xi[1:]**3 - self.xi[:-1]**3
                    diff_cost = 2.0
                    diff_phi = 2. * np.pi
                    for ix in range(self.nx):
                        vol[ix, 0, 0] = 1. / 3. * diff_r3[ix] * diff_cost * diff_phi

                else:
                    vol = np.zeros([self.nx, self.ny, self.nz], dtype=np.float64)
                    diff_r3 = self.xi[1:]**3 - self.xi[:-1]**3
                    diff_cost = 2.0
                    diff_phi = self.zi[1:] - self.zi[:-1]
                    for ix in range(self.nx):
                        for iz in range(self.nz):
                            vol[ix, 0, iz] = 1. / 3. * diff_r3[ix] * diff_cost * diff_phi[iz]

            elif self.act_dim[2] == 0:
                vol = np.zeros([self.nx, self.ny, self.nz], dtype=np.float64)
                diff_r3 = self.xi[1:]**3 - self.xi[:-1]**3
                diff_cost = np.cos(self.yi[:-1]) - np.cos(self.yi[1:])
                diff_phi = 2. * np.pi
                for ix in range(self.nx):
                    for iy in range(self.ny):
                        vol[ix, iy, :] = 1. / 3. * diff_r3[ix] * diff_cost[iy] * diff_phi

            else:
                vol = np.zeros([self.nx, self.ny, self.nz], dtype=np.float64)
                diff_r3 = self.xi[1:]**3 - self.xi[:-1]**3
                diff_cost = np.cos(self.yi[:-1]) - np.cos(self.yi[1:])
                diff_phi = self.zi[1:] - self.zi[:-1]
                for ix in range(self.nx):
                    for iy in range(self.ny):
                        vol[ix, iy, :] = 1. / 3. * diff_r3[ix] * diff_cost[iy] * diff_phi
        else:
            raise ValueError('Coordinate system ' + self.crd_sys + ' is not yet supported.')

        return vol


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
                            raise ValueError('Unknown file format in the dust opacity file '+fname)

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
                        print('Reading '+fname)
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

    # --------------------------------------------------------------------------------------------------------------------
    def makeOpac(self, ppar=None, wav=None, old=False):
        """Createst the dust opacities using a Mie code distributed with RADMC-3D.

        Parameters
        ----------

        ppar  : dictionary
                Parameters of the simulations

        wav   : ndarray, optional
                Wavelength grid on which the mass absorption coefficients should be calculated

        old   : bool, optional
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
            # if old:
            # print 'ERROR'
            # print 'Making old (RADMC) style opacities is not finished for more than one dust species'
            # exit()

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

                # Run makedust
                self.runMakedust(freq=nc.cc / wav * 1e4, gmin=ppar['gsmin'], gmax=ppar['gsmax'], ngs=ppar['ngs'],
                                 lnk_fname='opt_const.dat', gdens=ppar['gdens'][idust])

                # Change the name of makedust's output
                for igs in range(ppar['ngs']):
                    dum = sp.Popen('mv dustkappa_' + str(igs + 1) + '.inp dustkappa_idust_' + str(idust + 1)
                                   + '_igsize_' + str(igs + 1) + '.inp', shell=True).wait()
                    ext.append('idust_' + str(idust + 1) + '_igsize_' + str(igs + 1))

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

            # # Change the name of makedust's output
            #            dum = sp.Popen('mv dustkappa_1.inp dustkappa_idust_1_igsize_1.inp', shell=True).wait()
            #            os.remove('opt_const.dat')

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
                    dum[ii] = 10.**np.interp(np.log10(owav[ii]), np.log10(dw), np.log10(dcabs))

                    # Edwtrapolate the absorption coefficients using linear fit in log-log space
                    # (i.e. fitting a polinomial) for short wavelengths
                    # der = np.log10(dcabs[1] / dcabs[0]) / np.log10(dw[1] / dw[0])
                    dum[il] = 10.**(np.log10(dcabs[0]) + np.log10(dw[0] / owav[il]))

                    # Edwtrapolate the absorption coefficients using linear fit in log-log space
                    # (i.e. fitting a polinomial) for long wavelengths
                    # der = np.log10(dcabs[nwav - 1] / dcabs[nwav - 2]) / np.log10(dw[nwav - 1] / dw[nwav - 2])
                    dum[ih] = 10.**(np.log10(dcabs[nwav - 1]) + np.log10(owav[il] / dw[nwav - 1]))

                    ocabs = ocabs + np.array(dum) * mixabun[i][j]

                    if oform == 2:
                        # Do the inter-/extrapolation of for the scattering coefficients
                        dum = np.zeros(nwav0, dtype=float)
                        dum[ii] = 10.**np.interp(np.log10(owav[ii]), np.log10(dw), np.log10(dcsca))

                        # der = np.log10(dcsca[1] / dcsca[0]) / np.log10(dw[1] / dw[0])
                        dum[il] = 10.**(np.log10(dcsca[0]) + np.log10(dw[0] / owav[il]))

                        # der = np.log10(dcsca[nwav - 1] / dcsca[nwav - 2]) / np.log10(dw[nwav - 1] / dw[nwav - 2])
                        dum[ih] = 10.**(np.log10(dcsca[nwav - 1]) + np.log10(owav[il] / dw[nwav - 1]))

                        ocsca = ocsca + np.array(dum) * mixabun[i][j]

                    if oform == 3:
                        # Do the inter-/extrapolation of for the scattering phase function
                        dum = np.zeros(nwav0, dtype=float)
                        dum[ii] = 10.**np.interp(np.log10(owav[ii]), np.log10(dw), np.log10(gsym))

                        # der = np.log10(gsym[1] / gsym[0]) / np.log10(dw[1] / dw[0])
                        dum[il] = 10.**(np.log10(gsym[0]) + np.log10(dw[0] / owav[il]))

                        # der = np.log10(gsym[nwav - 1] / gsym[nwav - 2]) / np.log10(dw[nwav - 1] / dw[nwav - 2])
                        dum[ih] = 10.**(np.log10(gsym[nwav - 1]) + np.log10(owav[il] / dw[nwav - 1]))

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
            kabs[ii] = 10.**np.interp(np.log10(wav[ii]), np.log10(self.wav[i]), np.log10(self.kabs[i]))
            if len(self.ksca[i]) > 1:
                ksca[ii] = 10.**np.interp(np.log10(wav[ii]), np.log10(self.wav[i]), np.log10(self.ksca[i]))

            #
            # Do the long wavelength part
            #
            if True in ilo:
                x1 = np.log10(self.wav[i][-1])
                x0 = np.log10(self.wav[i][-2])

                y1 = np.log10(self.kabs[i][-1])
                y0 = np.log10(self.kabs[i][-2])
                der = (y1 - y0) / (x1 - x0)
                kabs[ilo] = 10.**(y1 + der * (np.log10(wav[ilo]) - x1))

                y1 = np.log10(self.ksca[i][-1])
                y0 = np.log10(self.ksca[i][-2])
                der = (y1 - y0) / (x1 - x0)
                ksca[ilo] = 10.**(y1 + der * (np.log10(wav[ilo]) - x1))

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
            gsize = gmin * (gmax / gmin)**(np.arange(ngs, dtype=np.float64) / (float(ngs) - 1.))
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


class radmc3dPar(object):
    """Parameter class for a RADMC-3D model.


    Attributes
    ----------

    ppar   : dictionary
            Contains parameter values with parameter names as keys 

    pdesc  : dictionary
            Contains parameter description (comments in the parameter file) with parameter names as keys

    pblock : dictionary
            Contains the block names in the parameter file and parameter names as values 

    pvalstr: dictionary
            Contains parameter values as strings with parameter names as keys

    """

    def __init__(self):

        self.ppar = {}
        self.pdesc = {}
        self.pblock = {}
        self.pvalstr = {}

    def readPar(self, fname=''):
        """Reads a parameter file.
        The parameters in the files should follow python syntax


        Parameters
        ----------

        fname  : str, optional
                File name to be read (if omitted problem_params.inp is used)

        Returns
        -------
            Returns a dictionary with the parameter names as keys

        """

        if fname == '':
            fname = 'problem_params.inp'

        with open(fname, 'r') as rfile:

            cchar = '#'
            lbchar = ""

            # Add pi to the local variables
            pi = np.pi
            # --------------------------------------------------------------------------------------------------------
            # First read every line that is not commented (i.e. does not begin with a comment character)
            # --------------------------------------------------------------------------------------------------------
            dumlist = []
            dumline = '-'

            dumline = rfile.readline()
            while dumline != '':
                # First check if the line is commented out, in which case ignore it
                comment = False
                if dumline[0] == cchar:
                    if dumline.find('Block') < 0:
                        comment = True

                # OK, the line is not commented out, now check if it contains a '=' sign (if not ignore the line)
                if not comment:
                    # Check if we have an empty line in which case also ignore it
                    if dumline.strip() != '':
                        dumlist.append(dumline)

                # Read the next line
                dumline = rfile.readline()

        # -------------------------------------------------------------------------------------------------------------
        # After every line in the file was read try to decode the lines to
        #  [variable name] = [variable value] # [comment]
        # also try to catch if an expression has been broken into multiple lines
        # -------------------------------------------------------------------------------------------------------------

        varlist = []
        iline = 0
        while iline < len(dumlist):
            # First check if the line contains an '=' sign if not we have a problem
            #  expression broken into multiple lines are should already be dealt with
            ind = dumlist[iline].find('=')
            if ind <= 0:
                if dumlist[iline].find('Block') <= 0:
                    raise SyntaxError(' Invalid expression in line ' + ("%d" % iline))
                else:
                    if dumlist[iline].find(':') <= 0:
                        raise SyntaxError('Invalid block identified. \n'
                                          + ' The syntax of the block name field is : \n'
                                          + ' # Block : BlockName ')
                    else:
                        blockname = dumlist[iline].split(':')[1].strip()

            else:
                # The line contains a '=' sign and a variable name, so let's check if the
                #  value expression is broken into multiple lines
                vlist = dumlist[iline].split('=')
                lbind = vlist[1].find('\\')
                cind = vlist[1].find('#')

                # The line is full not broken
                if lbind == -1:
                    # Check if there is a comment field
                    if cind >= 0:
                        vlist = [vlist[0], vlist[1][:cind], vlist[1][cind + 1:], blockname]
                    else:
                        vlist = [vlist[0], vlist[1][:cind], ' ', blockname]

                    varlist.append(vlist)
                # The value expression is broken into multiple lines; take all lines and join the pieces
                else:
                    # Check if there is any comment in the line
                    inBrokenLine = False
                    if cind >= 0:
                        # Part of the line is commented, now check if the line break is before or after the
                        # comment character
                        if lbind > cind:
                            # The line break is in the comment field so there is no real line break
                            vlist = [vlist[0], vlist[1][:cind], vlist[1][cind + 1:], blockname]
                        else:
                            # The line break is before the comment character
                            inBrokenLine = True
                            expr = vlist[1][:lbind]
                            com = vlist[1][cind + 1:]
                    else:
                        inBrokenLine = True
                        expr = vlist[1][:lbind]
                        com = ' '

                    if inBrokenLine:
                        # Now gather all other pieces of this line
                        iline2 = 0
                        while inBrokenLine:
                            iline2 = iline2 + 1
                            dummy = dumlist[iline + iline2]
                            # Search for comments
                            cind2 = dummy.find('#')
                            # Search for another line break
                            lbind2 = dummy.find('\\')

                            # TODO:
                            # At the moment I neglect the possiblity that the second line in a broken long line begins
                            # with a linebreak or commented out

                            # There is comment
                            if cind2 > 0:

                                # There is line break
                                if lbind2 > 0:
                                    # The line break is commented out
                                    if lbind2 > cind:
                                        expr = expr + dummy[:cind2].strip()
                                        com = com + dummy[cind2 + 1:]
                                        inBrokenLine = False
                                    else:
                                        # The line break is not commented out
                                        expr = expr + dummy[:lbind2].strip()
                                        com = com + dummy[cind2 + 1:]
                                else:
                                    # There is no line break
                                    expr = expr + dummy[:cind2].strip()
                                    com = com + dummy[cind2 + 1:]
                                    inBrokenLine = False

                            # There is no comment
                            else:
                                # There is a line break
                                if lbind2 > 0:
                                    expr = expr + dummy[:lbind2].strip()
                                    com = com + dummy[cind2 + 1:]

                                # There is no line break
                                else:
                                    expr = expr + dummy[:cind2].strip()
                                    com = com + ' '
                                    inBrokenLine = False
                        iline = iline + iline2
                        vlist = [vlist[0], expr, com, blockname]
                        varlist.append(vlist)

            iline = iline + 1
        # -------------------------------------------------------------------------------------------------------
        # Now evaluate the expressions in the value field and make the final dictionary
        # -------------------------------------------------------------------------------------------------------
        self.ppar = {}
        glob = globals()
        glob['pi'] = np.pi
        loc = locals()
        loc['pi'] = np.pi
        for i in range(len(varlist)):
            try:
                val = eval(varlist[i][1], glob)
                glob[varlist[i][0].strip()] = val
            except:
                try:
                    val = eval(varlist[i][1], loc)
                    loc[varlist[i][0].strip()] = val
                except:
                    raise SyntaxError('Unknown expression "' + varlist[i][1] + '"')
            self.ppar[varlist[i][0].strip()] = val
            self.pvalstr[varlist[i][0].strip()] = varlist[i][1].strip()
            self.pdesc[varlist[i][0].strip()] = varlist[i][2].strip()
            self.pblock[varlist[i][0].strip()] = varlist[i][3].strip()
        return

    # --------------------------------------------------------------------------------------------------
    def setPar(self, parlist=None):
        """Sets a parameter in the radmc3DPar class.
        If the paramter is already defined its value will be modified

        Parameters
        ----------

        parlist : list
                  If the parameter is already defined parlist should be a two element
                  list 1) parameter name, 2) parameter expression/value as a string

                  If the parameter is not yet defined parlist should be a four element
                  list 1) parameter name, 2) parameter expression/value as a string
                  3) Parameter description (= comment field in the parameter file)
        """

        if parlist is None:
            raise ValueError('Unknown parlist. \n No parameters to set.')
        else:
            parname = parlist[0].strip()

        #
        # Check whether or not the parameter is already defined
        #
        new_par = False
        if len(parlist) == 2:
            if parname not in self.ppar:
                raise ValueError(' parlist has too few elements.\n'
                                 + ' The argument of radmc3dPar.setPar() should be a four element list if a new \n'
                                 + ' parameter is defined 1) parameter name, 2) parameter expression/value as a '
                                 + ' string\n'
                                 + ' 3) Parameter description (= comment field in the parameter file) \n'
                                 + ' 4) The name of the block in which the parameter must be placed in the \n'
                                 + ' problem_params.inp file')
        else:
            new_par = True

        # Add pi to the local variables
        pi = np.pi
        #
        # Add the parameter to the dictionaries /change its value
        #
        glob = globals()
        glob['pi'] = np.pi
        loc = locals()
        loc['pi'] = np.pi

        try:
            self.ppar[parname] = eval(parlist[1].strip(), glob)
            glob[parname] = self.ppar[parname]
        except:
            try:
                self.ppar[parname] = eval(parlist[1].strip(), loc)
                loc[parname] = self.ppar[parname]
            except:
                raise SyntaxError('Unknown expression ' + parlist[1].strip())

        self.pvalstr[parname] = parlist[1].strip()

        if new_par:
            if parname not in self.pdesc:
                self.pdesc[parname] = parlist[2].strip()
            if len(parlist) == 4:
                if parname not in self.pblock:
                    self.pblock[parname] = parlist[3].strip()

    def loadDefaults(self, model='', ppar=None, reset=True):
        """Sets all parameters to default values.

        Parameters
        ----------
        model : str
                Model name whose paraemters should also be loaded

        ppar  : dictionary
                Contains parameter values as string and parameter names as keys
                Default values will be re-set to the values in this dictionary

        reset : bool
                If True the all class attributes will be re-initialized before
                the default values would be loaded. I.e. it will remove all entries
                from the dictionary that does not conain default values either in this
                function or in the optional ppar keyword argument
        """

        if reset:
            self.ppar = {}
            self.pvalstr = {}
            self.pdesc = {}
            self.pblock = {}

        #
        # Radiation sources
        #
        self.setPar(['incl_disc_stellarsrc', 'True', '# Switches on (True) or off (False) discrete stellar sources)',
                     'Radiation sources'])
        self.setPar(['mstar', '[1.0*ms]', '# Mass of the star(s)', 'Radiation sources'])
        self.setPar(['rstar', '[2.0*rs]', '# Radius of the star(s)', 'Radiation sources'])
        self.setPar(['tstar', '[4000.0]', '# Effective temperature of the star(s) [K]', 'Radiation sources'])
        self.setPar(
            ['pstar', '[0.0, 0.0, 0.0]', '# Position of the star(s) (cartesian coordinates)', 'Radiation sources'])
        self.setPar(['staremis_type', '["blackbody"]', '# Stellar emission type ("blackbody", "kurucz", "nextgen")',
                     'Radiation sources'])
        self.setPar(
            ['incl_cont_stellarsrc', 'False', '# Switches on (True) or off (False) continuous stellar sources )',
             'Radiation sources'])
        #
        # Grid parameters
        #
        self.setPar(['crd_sys', "'sph'", '  Coordinate system used (car/sph)', 'Grid parameters'])
        self.setPar(
            ['grid_style', '0', '  0 - Regular grid, 1 - Octree AMR, 10 - Layered/nested grid (not yet supported)',
             'Grid parameters'])
        self.setPar(
            ['nx', '50', '  Number of grid points in the first dimension (to switch off this dimension set it to 0)',
             'Grid parameters'])
        self.setPar(
            ['ny', '30', '  Number of grid points in the second dimension (to switch off this dimension set it to 0)',
             'Grid parameters'])
        self.setPar(
            ['nz', '36', '  Number of grid points in the third dimension (to switch off this dimension set it to 0)',
             'Grid parameters'])
        self.setPar(['xbound', '[1.0*au, 100.*au]', '  Boundaries for the x grid', 'Grid parameters'])
        self.setPar(['ybound', '[0.0, pi]', '  Boundaries for the y grid', 'Grid parameters'])
        self.setPar(['zbound', '[0.0, 2.0*pi]', '  Boundraries for the z grid', 'Grid parameters'])
        self.setPar(['xres_nlev', '3', 'Number of refinement levels (spherical coordinates only', 'Grid parameters'])
        self.setPar(['xres_nspan', '3', 'Number of the original grid cells to refine (spherical coordinates only)',
                     'Grid parameters'])
        self.setPar(
            ['xres_nstep', '3', 'Number of grid cells to create in a refinement level (spherical coordinates only)',
             'Grid parameters'])
        self.setPar(['wbound', '[0.1, 7.0, 25., 1e4]', '  Boundraries for the wavelength grid', 'Grid parameters'])
        self.setPar(['nw', '[19, 50, 30]', '  Number of points in the wavelength grid', 'Grid parameters'])
        self.setPar(['levelMaxLimit', '5', '  Highest refinement level in octree AMR', 'Grid parameters'])

        #
        # Dust opacity
        #
        self.setPar(['lnk_fname',
                     "['/disk2/juhasz/Data/JPDOC/astrosil/astrosil_WD2001_new.lnk', "
                     "'/disk2/juhasz/Data/JPDOC/carbon/A/cel600.lnk']",
                     ' ', 'Dust opacity'])
        self.setPar(['gdens', '[3.6, 1.8]', ' Bulk density of the materials in g/cm^3', 'Dust opacity'])
        self.setPar(['gsmin', '0.1', ' Minimum grain size', 'Dust opacity'])
        self.setPar(['gsmax', '10.0', ' Maximum grain size', 'Dust opacity'])
        self.setPar(['ngs', '1', ' Number of grain sizes', 'Dust opacity'])
        self.setPar(['gsdist_powex', '-3.5', ' Grain size distribution power exponent', 'Dust opacity'])
        self.setPar(['mixabun', '[0.75, 0.25]', ' Mass fractions of the dust componetns to be mixed', 'Dust opacity'])
        self.setPar(['dustkappa_ext', "['silicate']", ' ', 'Dust opacity'])

        #
        # Gas line RT
        #
        self.setPar(
            ['gasspec_mol_name', "['co']", '  Name of the gas species - the extension of the molecule_EXT.inp file',
             'Gas line RT'])
        self.setPar(['gasspec_mol_abun', '[1e-4]', '  Abundance of the molecule', 'Gas line RT'])
        self.setPar(['gasspec_mol_dbase_type', "['leiden']", '  leiden or linelist', 'Gas line RT'])
        self.setPar(
            ['gasspec_colpart_name', "['h2']", '  Name of the gas species - the extension of the molecule_EXT.inp file',
             'Gas line RT'])
        self.setPar(['gasspec_colpart_abun', '[1e0]', '  Abundance of the molecule', 'Gas line RT'])
        self.setPar(['gasspec_vturb', '0.1e5', '  Microturbulence', 'Gas line RT'])
        #
        # Code parameters
        #
        self.setPar(['nphot', 'int(1e5)', '  Nr of photons for the thermal Monte Carlo', 'Code parameters'])
        self.setPar(['nphot_scat', 'int(3e5)', '  Nr of photons for the scattering Monte Carlo (for images)',
                     'Code parameters'])
        self.setPar(['nphot_spec', 'int(1e5)', '  Nr of photons for the scattering Monte Carlo (for spectra)',
                     'Code parameters'])
        self.setPar(
            ['scattering_mode_max', '1', '  0 - no scattering, 1 - isotropic scattering, 2 - anizotropic scattering',
             'Code parameters'])
        self.setPar(['lines_mode', '-1', '  Line raytracing mode', 'Code parameters'])
        self.setPar(['istar_sphere', '0',
                     '  1 - take into account the finite size of the star, 0 - take the star to be point-like',
                     'Code parameters'])
        self.setPar(['itempdecoup', '1', '  Enable for different dust components to have different temperatures',
                     'Code parameters'])
        self.setPar(['tgas_eq_tdust', '1', '  Take the dust temperature to identical to the gas temperature',
                     'Code parameters'])
        self.setPar(
            ['rto_style', '1', '  Format of outpuf files (1-ascii, 2-unformatted f77, 3-binary', 'Code parameters'])
        self.setPar(
            ['modified_random_walk', '0', '  Switched on (1) and off (0) modified random walk', 'Code parameters'])
        #
        # Model parameters
        #
        if model != '':
            try:
                mdl = importlib.import_module(model)
            except:
                try:
                    mdl = importlib.import_module('radmc3dPy.models.' + model)
                except:
                    print(model + '.py could not be imported. \n '
                          + 'The model files should either be in the current working directory \n '
                          + 'or in the radmc3d python module directory')
                    print(traceback.format_exc())

            modpar = mdl.getDefaultParams()
            for i in range(len(modpar)):
                dum = modpar[i]
                if len(dum) == 3:
                    dum.append('Model ' + model)
                self.setPar(dum)

    def printPar(self):
        """Prints the parameters of the current model.

        """

        #
        # First get the unique block names
        #

        blocknames = ['Radiation sources', 'Grid parameters', 'Dust opacity', 'Gas line RT', 'Code parameters']
        for key in self.pblock.keys():
            dum = self.pblock[key]
            if not blocknames.__contains__(dum):
                blocknames.append(dum)

        #
        # Get the parameter block names and distionary keys
        #
        par_keys = list(self.pblock.keys())
        par_block = list(self.pblock.values())

        #
        # Print the parameters by blocks
        #
        for iblock in blocknames:
            print(
                '%s' % '# ------------------------------------------------------------------------------------------'
                       + '-------------------------------')
            txt = '# Block: ' + iblock
            print('%s' % txt)
            print(
                '%s' % '# ------------------------------------------------------------------------------------------'
                       + '-------------------------------')

            keys = []
            for i in range(len(par_block)):
                if par_block[i] == iblock:
                    keys.append(par_keys[i])

            keys.sort()
            for key in keys:
                print(key.ljust(25) + ' = ' + self.pvalstr[key].strip() + '  # ' + self.pdesc[key].strip())
                # --------------------------------------------------------------------------------------------------

    def writeParfile(self, fname=''):
        """Writes a parameter file.

        Parameters
        ----------

        fname  : str, optional
                File name to be read (if omitted problem_params.inp is used)

        """

        if fname == '':
            fname = 'problem_params.inp'

        print('Writing ' + fname)

        #
        # First get the uniq block names
        #

        blocknames = ['Radiation sources', 'Grid parameters', 'Dust opacity', 'Gas line RT', 'Code parameters']
        for key in self.pblock:
            dum = self.pblock[key]
            if not blocknames.__contains__(dum):
                blocknames.append(dum)

        with open(fname, 'w') as wfile:
            #
            # Write header
            #

            wfile.write(
                '%s\n' % ('#########################################################################################'
                         + '##################################'))
            wfile.write('%s\n' % '# RADMC-3D PARAMETER SETUP')
            wfile.write('%s\n' % '# Created by the python module of RADMC-3D')
            wfile.write(
                '%s\n' % ('#########################################################################################'
                         + '##################################'))

            #
            # Get the parameter block names and distionary keys
            #
            par_keys = list(self.pblock.keys())
            par_block = list(self.pblock.values())

            #
            # Write the parameterfile
            #
            for iblock in blocknames:
                wfile.write(
                    '%s\n' % ('# -----------------------------------------------------------------------------------'
                             + '--------------------------------------'))
                txt = '# Block: ' + iblock
                wfile.write('%s\n' % txt)
                wfile.write(
                    '%s\n' % ('# -----------------------------------------------------------------------------------'
                             + '--------------------------------------'))

                keys = []
                for i in range(len(par_block)):
                    if par_block[i] == iblock:
                        keys.append(par_keys[i])

                keys.sort()
                for key in keys:
                    wfile.write(key.ljust(25) + ' = ' + self.pvalstr[key].strip() + '  # '
                                + self.pdesc[key].strip() + '\n')


class radmc3dRadSources(object):
    """Class of the radiation sources.
    Currently discrete stars and continuous starlike source, the latter only in spherical coordinates.


    Attributes
    ----------

    wav          : ndarray
                    Wavelength for the stellar spectrum

    freq         : ndarray
                    Frequency for the stellar spectrum

    nwav         : int
                    Number of wavelenghts in the stellar spectrum

    nfreq        : int
                    Number of frequencies in the stellar spectrum

    mstar        : list
                    List of stellar masses

    tstar        : list
                    List of stellar effective temperatures

    rstar        : list
                    List of stellar radii

    lstar        : list 
                    List of stellar luminosities 

    nstar        : int
                    Number of stars

    pstar        : list
                    Each element of the list contains a three element list, the cartesian coordinates of the stars

    fnustar      : ndarray
                    Stellar spectrum (flux@1pc in erg/s/cm/cm/Hz)

    csdens       : ndarray
                    Stellar density for continuous starlike source

    csntemplate  : int
                    Number of stellar templates

    cstemp       : ndarray
                    Stellar template

    cstemptype   : int
                    Stellar template type 1 - Blackbody given by the effective temperature 
                    2 - Frequency dependent spectrum

    cststar      : ndarray
                    Stellar effective temperature

    csmstar      : ndarray
                    Stellar mass

    csrstar      : ndarray
                    Stellar radius

    tacc         : ndarray
                    Effective temperature of a viscous accretion disk as a function of cylindrical radius

    accrate      : float
                    Accretion rate of the viscous accretion disk [g/s]

    fnuaccdisk   : ndarray
                    Spatially integrated frequency-dependent flux density of the accretion disk @ 1pc distance

    tspot        : float
                    Temperature of the hot spot / boundary layer on the stellar surface

    starsurffrac : float
                    Fraction of the stellar surface covered by the hot spot / boundary layer

    fnustpot     : ndarray
                    Frequency-dependent flux density of the hot spot / boundary layer @ 1pc distance

    """

    def __init__(self, ppar=None, grid=None):

        # Spatial and frequency grid
        self.grid = grid

        # Discrete stars

        self.mstar = []
        self.tstar = []
        self.rstar = []
        self.lstar = []
        self.nstar = 0
        self.pstar = []
        self.fnustar = []
        self.staremis_type = []

        # Continuous starlike source
        self.csdens = []
        self.csntemplate = 0
        self.cstemp = []
        self.cstemptype = 1
        self.cststar = []
        self.csmstar = []
        self.csrstar = []

        # Viscous accretion disk
        self.incl_accretion = False
        self.tacc = []
        self.accrate = 0.
        self.fnuaccdisk = []

        # Hot spot / boundary layer - to model accretion in YSOs
        self.tspot = 0.
        self.starsurffrac = 0.
        self.fnuspot = []

        if ppar:
            if isinstance(ppar['mstar'], list):
                self.mstar = ppar['mstar']
            else:
                self.mstar = [ppar['mstar']]

            if isinstance(ppar['tstar'], list):
                self.tstar = ppar['tstar']
            else:
                self.tstar = [ppar['tstar']]

            if isinstance(ppar['rstar'], list):
                self.rstar = ppar['rstar']
            else:
                self.rstar = [ppar['rstar']]

            for istar in range(self.nstar):
                self.lstar.append(4. * np.pi * self.rstar[istar]**2. * nc.ss * self.tstar[istar]**4.)
            self.pstar = ppar['pstar']

            if 'incl_cont_stellarsrc' in ppar:
                self.incl_accretion = ppar['incl_cont_stellarsrc']
            else:
                self.incl_accretion = False

            if 'accrate' in ppar:
                self.accrate = ppar['accrate']
            else:
                self.accrate = 0.

    def findPeakStarspec(self):

        """Calculates the peak wavelength of the stellar spectrum.

        Returns
        -------

        The peak wavelength of the stellar spectrum in nu*Fnu for all 
            stars as a list
        """

        pwav = []

        for istar in range(self.nstar):
            ii = (self.fnustar[:, istar] * self.grid.freq).argmax()
            pwav.append(self.grid.wav[ii])

        return pwav

    def readStarsinp(self, fname=''):
        """Reads the data of discrete stellar sources from the stars.inp file.

        Parameters
        ----------

        fname : str, optional
                File name of the file that should be read (if omitted stars.inp will be used)
        """

        if fname == '':
            fname = 'stars.inp'

        with open(fname, 'r') as rfile:

            dum = rfile.readline()
            iformat = int(dum)
            if iformat != 2:
                raise ValueError(' Unknown file format. Format number ', iformat, ' is unknown.')

            dum = rfile.readline().split()
            self.nstar = int(dum[0])
            # self.grid  = radmc3dGrid()
            self.grid = readGrid(sgrid=False)
            self.grid.nwav = int(dum[1])
            self.grid.nfreq = self.grid.nwav
            self.rstar = []
            self.mstar = []
            self.tstar = []
            for istar in range(self.nstar):
                dum = rfile.readline().split()
                self.rstar.append(float(dum[0]))
                self.mstar.append(float(dum[1]))
                self.pstar.append([float(dum[2]), float(dum[3]), float(dum[4])])

            dum = rfile.readline()
            wav = []
            for ilam in range(self.grid.nwav):
                dum = rfile.readline()
                wav.append(float(dum))

            self.grid.wav = np.array(wav, dtype=float)
            self.grid.freq = nc.cc / self.grid.wav * 1e4

            # Empty line
            dum = rfile.readline()
            self.fnustar = np.zeros([self.grid.nfreq, self.nstar], dtype=float)
            for istar in range(self.nstar):
                # First data line containing the stellar spectrum or the effective temperature (if negative)
                dum = rfile.readline()
                # If we have only the stellar temperature
                if float(dum) < 0:
                    self.tstar.append(-float(dum))
                    #
                    # Now calculates the stellar spectrum
                    #
                    self.fnustar[:, istar] = 2. * nc.hh * self.grid.freq**3. / nc.cc**2 / (
                                         np.exp(nc.hh * self.grid.freq / nc.kk / self.tstar[istar]) - 1.0) \
                                         * np.pi * self.rstar[istar]**2. / nc.pc**2.

                # If we have the full frequency-dependent spectrum
                else:
                    self.tstar.append(0.0)
                    self.fnustar[0, istar] = float(dum)
                    for ifreq in range(1, self.grid.nfreq):
                        self.fnustar[ifreq, istar] = float(rfile.readline())

    def writeStarsinp(self, ppar=None, wav=None, freq=None, old=False):
        """Writes the input file for discrete stellar sources (stars.inp).

        Parameters
        ----------

        ppar  : dictionary
                Dictionary containing all parameters of the model (only mandatory if accretion is switched on)

        wav   : ndarray, optional
                Wavelength grid for the stellar spectrum

        freq  : ndarray, optional
                Frequency grid for the stellar spectrum (either freq or wav should be set)

        old   : bool, optional
                If set to True the file format of the previous, 2D version of radmc will be used
        """
        if not old:
            if freq is not None:
                if wav is not None:
                    raise ValueError('Either wav or freq should be set but not both.')
                else:
                    self.grid.wav = nc.cc / np.array(freq) * 1e4
                    self.grid.freq = np.array(freq)
                    self.grid.nwav = self.grid.wav.shape[0]
                    self.grid.nfreq = self.grid.nwav

            elif wav is not None:
                if freq is not None:
                    raise ValueError('Either wav or freq should be set but not both.')
                else:
                    self.grid.wav = np.array(wav)
                    self.grid.freq = nc.cc / self.grid.wav * 1e4
                    self.grid.nwav = self.grid.wav.shape[0]
                    self.grid.nfreq = self.grid.nwav


            self.nstar = len(self.rstar)
            # self.pstar = ppar['pstar']

            print('Writing stars.inp')
            with open('stars.inp', 'w') as wfile:
                wfile.write('%d\n' % 2)
                wfile.write('%d %d\n' % (self.nstar, self.grid.nwav))

                if self.nstar > 1:
                    for istar in range(self.nstar):
                        wfile.write('%.9e %.9e %.9e %.9e %.9e\n' % (self.rstar[istar], self.mstar[istar],
                                                                    self.pstar[istar][0], self.pstar[istar][1],
                                                                    self.pstar[istar][2]))
                else:
                    wfile.write('%.9e %.9e %.9e %.9e %.9e\n' % (self.rstar[0], self.mstar[0],
                                                                self.pstar[0], self.pstar[1], self.pstar[2]))

                wfile.write('%s\n' % ' ')
                for ilam in range(self.grid.nwav):
                    wfile.write('%.9e\n' % self.grid.wav[ilam])
                wfile.write('%s\n' % ' ')

                # If accretion is active write the sum of the spot and the star emission
                # NOTE, for now the spot emission is only added to the first star in the list
                if self.incl_accretion:
                    # Get the stellar spectrum
                    self.getStarSpectrum(tstar=self.tstar, rstar=self.rstar)
                    # Get the spot spectrum
                    self.getSpotSpectrum(ppar=ppar)

                    # Write out the spectrum of the first star with the additional spot luminosity
                    for ilam in range(self.grid.nwav):
                        wfile.write('%.9e\n' % (self.fnustar[ilam, 0] + self.fnuspot[ilam]))

                    # Now write the spectrum of all other discrete stars without the spot emission
                    for istar in range(1, self.nstar):
                        for ilam in range(self.grid.nwav):
                            wfile.write('%.9e\n' % (self.fnustar[ilam, istar]))

                else:
                    self.getStarSpectrum(ppar=ppar)
                    for istar in range(self.nstar):
                        if self.staremis_type[istar].strip().lower() == "blackbody":
                            wfile.write('%.9e\n' % (-self.tstar[istar]))
                        else:
                            for ilam in range(self.grid.nwav):
                                wfile.write('%.9e\n' % (self.fnustar[ilam, istar]))
        else:
            if freq is not None:
                self.grid.wav = nc.cc / np.array(freq) * 1e4
                self.grid.freq = np.array(freq)
                self.grid.nwav = self.grid.wav.shape[0]
                self.grid.nfreq = self.grid.nwav

            if wav is not None:
                self.grid.wav = np.array(wav)
                self.grid.freq = nc.cc / self.grid.wav * 1e4
                self.grid.nwav = self.grid.wav.shape[0]
                self.grid.nfreq = self.grid.nwav

            else:
                raise ValueError('Either wav or freq should be set but not both.')

            self.nstar = len(self.rstar)
            self.getStarSpectrum(ppar=ppar)

            if self.grid.freq[-1] < self.grid.freq[0]:
                self.grid.freq = self.grid.freq[::-1]
                self.grid.wav = self.grid.wav[::-1]
                for istar in range(self.nstar):
                    self.fnustar[:, istar] = self.fnustar[::-1, istar]

            print('Writing starinfo.inp')
            fname = 'starinfo.inp'
            with open(fname, 'w') as wfile:
                wfile.write("1\n")
                wfile.write("%.7e\n" % ppar['rstar'][0])
                wfile.write("%.7e\n" % ppar['mstar'][0])
                wfile.write("%.7e\n" % ppar['tstar'][0])

            print('Writing starspectrum.inp')
            fname = 'starspectrum.inp'
            with open(fname, 'w') as wfile:
                wfile.write("%d\n" % self.grid.nfreq)
                wfile.write(" \n")
                for i in range(self.grid.nfreq):
                    wfile.write("%.7e %.7e\n" % (self.grid.freq[i], self.fnustar[i, 0]))

    def getStarSpectrum(self, tstar=None, rstar=None, lstar=None, mstar=None, ppar=None, grid=None):
        """Calculates a blackbody stellar spectrum.

        Parameters
        ----------

        tstar : list
                Effective temperature of the stars in [K]

        rstar : list
                Radius of the stars in [cm]

        lstar : list
                Bolometric luminosity of the star [erg/s] (either rstar or lstar should be given)

        mstar : list
                Stellar mass in [g] (only required if an atmosphere model is used to calculate logg)

        ppar  : dictionary
                Dictionary containing all input parameters

        grid  : radmc3dGrid, optional
                An instance of a radmc3dGrid class containing the spatial and wavelength grid
        """
        #
        # Check the input which parameters are set and which should be calculated
        #

        if grid is not None:
            self.grid = grid

        if ppar is not None:
            if tstar is None:
                tstar = ppar['tstar']
            if rstar is None:
                rstar = ppar['rstar']
            if mstar is None:
                self.mstar = ppar['mstar']

        if mstar:
            self.mstar = mstar

        if tstar:
            if not isinstance(tstar, list):
                tstar = [tstar]
            dum1 = len(tstar)
            if lstar and rstar:
                raise ValueError(' Only two of the input variables tstar, rstar, lstar should be set not all three')
            elif lstar:
                if len(lstar) != dum1:
                    raise ValueError('lstar and tstar have different number of elements')
                else:
                    self.tstar = np.array(tstar)
                    self.lstar = np.array(lstar)
                    self.nstar = self.lstar.shape[0]
                    self.rstar = np.sqrt(self.lstar / (4. * np.pi * nc.ss * self.tstar**4.))
        else:
            if lstar and rstar:
                if len(lstar) != len(rstar):
                    raise ValueError('lstar and rstar have different number of elements')
                else:
                    self.lstar = np.array(lstar)
                    self.rstar = np.array(rstar)
                    self.nstar = self.rstar.shape[0]
                    self.tstar = (self.lstar / (4. * np.pi * nc.ss * self.rstar**2.))**0.25

        #
        # If we take blackbody spectrum
        #

        if ppar:
            if 'staremis_type' in ppar:
                self.staremis_type = ppar['staremis_type']
            else:
                self.staremis_type = []
                for i in range(self.nstar):
                    self.staremis_type.append("blackbody")
        else:
            self.staremis_type = []
            for i in range(self.nstar):
                self.staremis_type.append("blackbody")

        self.fnustar = np.zeros([self.grid.nwav, len(self.tstar)], dtype=np.float64)
        for istar in range(self.nstar):
            if self.staremis_type[istar].strip().lower() == "blackbody":
                self.fnustar[:, istar] = 2. * nc.hh * self.grid.freq**3. / nc.cc**2 / (
                                         np.exp(nc.hh * self.grid.freq / nc.kk / self.tstar[istar]) - 1.0) \
                                         * np.pi * self.rstar[istar]**2. / nc.pc**2.
            elif self.staremis_type[istar].strip().lower() == "kurucz":
                dum = staratm.getAtmModel(teff=self.tstar[istar], mstar=self.mstar[istar], rstar=self.rstar[istar],
                                          iwav=self.grid.wav, model="kurucz")
                self.fnustar[:, istar] = dum['lnu'] / (4. * np.pi * nc.pc**2)

            elif self.staremis_type[istar].strip().lower() == "nextgen":
                dum = staratm.getAtmModel(teff=self.tstar[istar], mstar=self.mstar[istar], rstar=self.rstar[istar],
                                          iwav=self.grid.wav, model="nextgen")
                self.fnustar[:, istar] = dum['lnu'] / (4. * np.pi * nc.pc**2)

            elif self.staremis_type[istar].strip().lower() == "ames-dusty":
                dum = staratm.getAtmModel(teff=self.tstar[istar], mstar=self.mstar[istar], rstar=self.rstar[istar],
                                          iwav=self.grid.wav, model="ames-dusty")
                self.fnustar[:, istar] = dum['lnu'] / (4. * np.pi * nc.pc**2)

            else:
                raise ValueError('Unknown stellar atmosphere model : ', self.staremis_type[istar],
                                 ' Only "kurucz" or "nextgen" are supported.')

    def getAccdiskTemperature(self, ppar=None, grid=None):
        """Calculates the effective temperature of a viscous accretion disk.

        Parameters
        ----------

        ppar : dictionary
               Dictionary containing all input parameters keys should include
               * mstar   : stellar mass
               * rstar   : stellar radius
               * accrate : accretion rate

               NOTE, that for the calculation of the effective disk temperature only the first
               star is used if more than one values are given in mstar and rstar. 

        grid : radmc3dGrid, optional 
               An instance of a radmc3dGrid class containing the spatial and wavelength grid

        """

        if grid is not None:
            self.grid = grid

        self.tacc = ((3.0 * nc.gg * ppar['mstar'][0] * ppar['accrate']) / (8.0 * np.pi * self.grid.x**3 * nc.ss)
                     * (1.0 - (ppar['rstar'][0] / self.grid.x)**0.5))**0.25

    # --------------------------------------------------------------------------------------------------
    def getSpotSpectrum(self, ppar=None, grid=None):
        """Calculates the spectrum of a hot spot / boundary layer on the stellar surface

        Parameters
        ----------

        ppar : dictionary
               Dictionary containing all input parameters keys should include
               * mstar   : stellar mass
               * rstar   : stellar radius
               * accrate : accretion rate

               NOTE, that for the calculation of the effective disk temperature only the first
               star is used if more than one values are given in mstar and rstar. 

        grid : radmc3dGrid, optional 
               An instance of a radmc3dGrid class containing the spatial and wavelength grid
        """
        #
        # Check the input which parameters are set and which should be calculated
        #

        if grid is not None:
            self.grid = grid

        # Calculate the total spot luminosity assuming boundary layer accretion, i.e.
        # that half of the total accretion luminosity is released from the boundary layer

        if 'accrate' in ppar:
            tot_acclum = 0.5 * nc.gg * ppar['mstar'][0] * ppar['accrate'] / ppar['rstar'][0]
            spotsurf = 4. * np.pi * ppar['rstar'][0]**2 * ppar['starsurffrac']
            self.starsurffrac = ppar['starsurffrac']
            if spotsurf == 0.:
                self.tspot = 0.

                # Now calculate the spot spectrum (i.e. the flux density @ 1pc distance)
                self.fnuspot = np.zeros(self.grid.nfreq, dtype=float)
            else:
                self.tspot = (tot_acclum / spotsurf / nc.ss)**0.25

                # Now calculate the spot spectrum (i.e. the flux density @ 1pc distance)
                self.fnuspot = np.pi * ppar['rstar'][0]**2 * ppar['starsurffrac'] / nc.pc**2 \
                               * 2. * nc.hh * self.grid.freq**3 / nc.cc**2 \
                               / (np.exp(nc.hh * self.grid.freq / nc.kk / self.tspot) - 1.0)
        else:
            self.fnuspot = np.zeros(self.grid.nfreq, dtype=float)

    def getTotalLuminosities(self, readInput=True):
        """Calcultes the frequency integrated luminosities of all radiation sources.


        Parameters
        ----------

        readInput : bool, optional
                    If true the input files of the radiation sources are read and the the total luminosities
                    are calculated from them. If readInput is set to False, the luminosities are calculated
                    by semi-analytic spectra.

        Returns
        -------

        Returns a dictionary with the following keys
            * lnu_star    : Luminosity of the discrete stars
            * lnu_accdisk : Luminosity of the accretion disk
            * lnu_spot    : Luminosity of the hot spot / boundary layer on the stellar surface
        """

        res = {'lnu_star': np.zeros(self.nstar, dtype=float),
               'lnu_accdisk': 0.,
               'lnu_spot': 0.}

        # If readIpnut
        if readInput:
            self.readStarsinp()

            # Note the negative sign in dnu is there because the frequency array is ordered in wavelength not in
            # frequency
            dnu = -(self.grid.freq[1:] - self.grid.freq[:-1])

            for istar in range(self.nstar):
                res['lnu_star'][istar] = 0.5 * (
                    (self.fnustar[1:, istar] + self.fnustar[:-1, istar]) * dnu).sum() * 4. * np.pi * nc.pc**2

            # Calculate the luminosity in the continuous stellar sources
            csDensFound = False
            csTempFound = False
            if os.path.exists('stellarsrc_density.inp'):
                self.readStellarsrcDensity(fname='stellarsrc_density.inp', binary=False)
                csDensFound = True
            if os.path.exists('stellarsrc_density.binp'):
                self.readStellarsrcDensity(fname='stellarsrc_density.binp', binary=True)
                csDensFound = True
            if os.path.exists('stellarsrc_templates.inp'):
                self.readStellarsrcTemplates()
                csTempFound = True
            if csDensFound & csTempFound:
                vol = self.grid.getCellVolume()

                dnu = abs(self.grid.freq[1:] - self.grid.freq[:-1])
                lum = 0.
                for itemp in range(self.csntemplate):
                    for ix in range(self.grid.nx):
                        for iy in range(self.grid.ny):
                            for iz in range(self.grid.nz):
                                if self.csdens[ix, iy, iz, itemp] > 0.:
                                    expterm = (nc.hh * self.grid.freq / nc.kk / (-self.cststar[ix])).clip(-600, 600)
                                    bb = 2. * nc.hh * self.grid.freq**3 / nc.cc**2 / (np.exp(expterm) - 1.)
                                    dum = bb * np.pi * vol[ix, iy, iz] * 4. * np.pi * self.csdens[ix, iy, iz, itemp]
                                    lum = lum + ((dum[1:] + dum[:-1]) * 0.5 * dnu).sum()

                res['lnu_accdisk'] = lum

            else:
                res['lnu_spot'] = 0.
                res['lnu_accdisk'] = 0.
        else:

            for istar in range(self.nstar):
                res['lnu_star'][istar] = 4. * np.pi * self.rstar[istar]**2 * nc.ss * self.tstar[istar]**4

            if self.accrate == 0.:
                print('Viscsous accretion is switched off')
                res['lnu_spot'] = 0.
                res['lnu_accdisk'] = 0.
            if not self.incl_accretion:
                print('Viscsous accretion is switched off')
                res['lnu_spot'] = 0.
                res['lnu_accdisk'] = 0.
            else:
                # Calculate the spot luminosity
                res['lnu_spot'] = 0.5 * nc.gg * self.mstar[0] * self.accrate / self.rstar[0]
                # Calculate the accretion disk luminosity
                res['lnu_accdisk'] = 0.5 * nc.gg * self.mstar[0] * self.accrate / self.rstar[0]

        return res

    def getAccdiskSpectra(self, ppar=None, grid=None, incl=0.):
        """Calculates the emergent spectra of an optically thick accretion disk at face-on orientation (incl=0deg).

        Parameters
        ----------

        ppar : dictionary
               Dictionary containing all input parameters keys should include
               * mstar   : stellar mass
               * rstar   : stellar radius
               * accrate : accretion rate

               NOTE, that for the calculation of the effective disk temperature only the first
               star is used if more than one values are given in mstar and rstar. 

        incl : float, optional
               Inclination angle in degrees at which the spectrum be calculated  (default - 0deg)

        grid : radmc3dGrid, optional 
               An instance of a radmc3dGrid class containing the spatial and wavelength grid

        """

        if 'accrate' in ppar:
            if ppar['accrate'] > 0.:
                self.getAccdiskTemperature(ppar=ppar, grid=grid)
                fnuaccdisk = np.zeros([self.grid.nx, self.grid.nwav], dtype=float)
                for ix in range(self.grid.nx):
                    dum = nc.hh * self.grid.freq / nc.kk / self.tacc[ix]
                    dum = dum.clip(-600., 600.)
                    bb = 2. * nc.hh * self.grid.freq**3 / nc.cc**2 / (np.exp(np.float64(dum)) - 1.0)
                    fnuaccdisk[ix, :] = bb * np.pi * (self.grid.xi[ix + 1]**2 - self.grid.xi[ix]**2) / nc.pc**2
                self.fnuaccdisk = fnuaccdisk.sum(0) * np.cos(incl/180.*np.pi)
            else:
                self.fnuaccdisk = np.zeros([self.grid.nwav], dtype=float)
        else:
            self.fnuaccdisk = np.zeros([self.grid.nwav], dtype=float)

    def getAccdiskStellarTemplates(self, ppar=None, grid=None):
        """Calculates the stellar template for continuous starlike sources for modeling a viscous accretion disk.


        Parameters
        ----------

        ppar : dictionary
               Dictionary containing all input parameters keys should include:
               * mstar   : stellar mass
               * rstar   : stellar radius
               * accrate : accretion rate

               NOTE, that for the calculation of the effective disk temperature only the first
               star is used if more than one values are given in mstar and rstar. 

        grid : radmc3dGrid, optional 
               An instance of a radmc3dGrid class containing the spatial and wavelength grid

        """

        if self.incl_accretion:
            self.getAccdiskTemperature(ppar=ppar, grid=grid)
            self.cstemptype = 1
            self.cststar = -self.tacc
            self.csmstar = self.cststar * 0. + 1.
            self.csrstar = self.cststar * 0. + 1.
            self.csntemplate = self.grid.nx
        else:

            self.cstemptype = 1
            self.cststar = np.zeros(self.grid.nx)
            self.csmstar = self.cststar * 0.
            self.csrstar = self.cststar * 0.
            self.csntemplate = self.grid.nx

    def getAccdiskStellarDensity(self, grid=None):
        """Calculates the stellar density for continuous starlike sources for modeling a viscous accretion disk.

        Parameters
        ----------
        grid : radmc3dGrid, optional 
               An instance of a radmc3dGrid class containing the spatial and wavelength grid

        """

        if grid is not None:
            self.grid = grid

        self.csdens = np.zeros([self.grid.nx, self.grid.ny, self.grid.nz, self.grid.nx], dtype=float)
        vol = self.grid.getCellVolume()

        if self.incl_accretion:
            if self.grid.crd_sys != 'sph':
                raise RuntimeError('Viscous accretion is currently available only in spherical coordinate system')
            else:
                if abs(self.grid.yi[self.grid.ny] - np.pi) < 1e-8:

                    for ix in range(self.grid.nx):
                        dA = 2.0 * (self.grid.xi[ix + 1]**2 - self.grid.xi[ix]**2) * np.pi * (
                            self.grid.zi[1:] - self.grid.zi[:-1]) / (2. * np.pi)
                        dV = vol[ix, self.grid.ny / 2 - 1, :] + vol[ix, self.grid.ny / 2, :]
                        self.csdens[ix, self.grid.ny / 2 - 1, :, ix] = dA / dV / (4. * np.pi)
                        self.csdens[ix, self.grid.ny / 2, :, ix] = dA / dV / (4. * np.pi)

                elif abs(self.grid.yi[self.grid.ny] - np.pi / 2.) < 1e-8:
                    for ix in range(self.grid.nx):
                        dA = 2.0 * (self.grid.xi[ix + 1]**2 - self.grid.xi[ix]**2) * np.pi * (
                            self.grid.zi[1:] - self.grid.zi[:-1]) / (2. * np.pi)
                        dV = vol[ix, self.grid.ny - 1, :] * 2.
                        self.csdens[ix, self.grid.ny - 1, :, ix] = dA / dV / (4. * np.pi)

        return True

    def readStellarsrcTemplates(self, fname='stellarsrc_templates.inp'):
        """Reads the stellar template of a continuous starlike source.

        Parameters
        ----------

        fname : str, optional 
                Name of the file from which the stellar templates will be read. If omitted the default
                'stellarsrc_templates.inp' will be used.
        """

        with open(fname, 'r') as rfile:

            self.grid = readGrid()

            hdr = np.fromfile(fname, count=3, sep="\n", dtype=int)

            if self.grid.nwav != hdr[2]:
                raise ValueError('Number of grid points in wavelength_micron.inp is not equal to that in ' + fname)
            else:
                self.csntemplate = hdr[1]
                dum = ''

                # Read the header
                for i in range(3):
                    dum = rfile.readline()

                # Read the frequency grid
                for i in range(hdr[2]):
                    dummy = float(rfile.readline())
                    if np.abs((nc.cc / dummy * 1e4 - self.grid.wav[i]) / self.grid.wav[i]) > 1e-4:
                        raise ValueError('The wavelength grid in wavelength_micron.inp is different '
                                         + 'from that in ' + fname)

                dum = rfile.readline()
                if float(dum) > 0.0:
                    self.cstemp = np.zeros([self.csntemplate, self.grid.nwav], dtype=float)
                    self.cststar = []
                    self.csrstar = []
                    self.csmstar = []

                    self.cstemp[0, 0] = float(dum)
                    for inu in range(1, self.grid.nwav):
                        dum = rfile.readline()
                        self.cstemp[0, inu] = float(dum)

                    for itemp in range(1, self.csntemplate):
                        for inu in range(self.grid.nwav):
                            dum = rfile.readline()
                            self.cstemp[itemp, inu] = float(dum)

                else:
                    self.cstemp = []
                    self.cststar = np.zeros(self.csntemplate, dtype=float)
                    self.csrstar = np.zeros(self.csntemplate, dtype=float)
                    self.csmstar = np.zeros(self.csntemplate, dtype=float)

                    self.cststar[0] = float(dum)
                    dum = rfile.readline()
                    self.csrstar[0] = float(dum)
                    dum = rfile.readline()
                    self.csmstar[0] = float(dum)

                    for i in range(1, self.csntemplate):
                        dum = rfile.readline()
                        self.cststar[i] = float(dum)
                        dum = rfile.readline()
                        self.csrstar[i] = float(dum)
                        dum = rfile.readline()
                        self.csmstar[i] = float(dum)

    def writeStellarsrcTemplates(self, fname='stellarsrc_templates.inp'):
        """Writes the stellar template of a continuous starlike source.

        Parameters
        ----------

        fname : str, optional
                Name of the file into which the stellar templates will be written. If omitted the default
                'stellarsrc_templates.inp' will be used.
        """
        # First check if we'd need to write anything at al

        if len(self.cststar) == 0:
            if len(self.cstemp) == 0:
                if os.path.exists('stellarsrc_templates.inp'):
                    print('The continuous starlike source seems to be inactive (zero input luminosity) '
                          + 'still stellarsrc_templates.inp file is present in the current working directory.')
                    dum = input('Can it be deleted (yes/no)')
                    if dum.strip().lower()[0] == 'y':
                        os.system('rm stellarsrc_templates.inp')
                    return
                return
            else:
                if np.abs(self.cstemp).max() == 0.:
                    if os.path.exists('stellarsrc_templates.inp'):
                        print('The continuous starlike source seems to be inactive (zero input luminosity)'
                              + ' still stellarsrc_templates.inp file is present in the current working directory.')
                        dum = input('Can it be deleted (yes/no)')
                        if dum.strip().lower()[0] == 'y':
                            os.system('rm stellarsrc_templates.inp')
                        return
                    return

        else:
            if np.abs(self.cststar).max() == 0.:
                if os.path.exists('stellarsrc_templates.inp'):
                    print('The continuous starlike source seems to be inactive (zero input luminosity) '
                          + ' still stellarsrc_templates.inp file is present in the current working directory.')
                    dum = input('Can it be deleted (yes/no)')
                    if dum.strip().lower()[0] == 'y':
                        os.system('rm stellarsrc_templates.inp')
                    return
                return

        print('Writing ' + fname)
        with open(fname, 'w') as wfile:
            # Format number
            wfile.write("%d\n" % 1)
            # Nr of templates
            wfile.write("%d\n" % self.csntemplate)
            # Nr of wavelengths
            wfile.write("%d\n" % self.grid.nwav)
            # Write the wavelength grid (in micron!)
            for ilam in range(self.grid.nwav):
                wfile.write("%.9e\n" % self.grid.freq[ilam])

            # Now write the templates
            # Similar to the discrete stellar imput if the first number is negative it means that
            #  instead of a full frequency-dependent spectrum only the blackbody temperature is given.
            #  Thus I'd only need to give the temperatures as negative numbers and radmc-3d will take care
            #  of calculating the Planck-function. This will save some harddisk space.

            if self.cstemptype == 1:
                for itemp in range(self.csntemplate):
                    # Effective temperature
                    if self.cststar[itemp] > 0:
                        wfile.write("%.9e\n" % (-self.cststar[itemp]))
                    else:
                        wfile.write("%.9e\n" % self.cststar[itemp])
                    # "Stellar radius"
                    wfile.write("%.9e\n" % self.csrstar[itemp])
                    # "Stellar mass"
                    wfile.write("%.9e\n" % self.csmstar[itemp])
            elif self.cstemptype == 2:
                for itemp in range(self.csntemplate):
                    for inu in range(self.grid.nwav):
                        wfile.write("%.9e\n" % self.cstemp[itemp, inu])
            else:
                raise ValueError('Unknown cstemptype for the continuous starlike source')

    def readStellarsrcDensity(self, fname=None, binary=False):
        """Reads the stellar density of a continuous starlike source.

        Parameters
        ----------

        fname  : str, optional
                 Name of the file from which the stellar templates will be read. If omitted the default
                 'stellarsrc_templates.inp' will be used.

        binary : bool, optional
                 If True the file should contain a C-style binary stream, if False the file should be 
                 written as formatted ASCII
        """

        if fname is None:
            if binary:
                fname = 'stellarsrc_density.binp'
            else:
                fname = 'stellarsrc_density.inp'

        self.grid = readGrid()
        self.csdens = None

        if binary:
            hdr = np.fromfile(fname, count=4, dtype=int)

            if hdr[2] != (self.grid.nx * self.grid.ny * self.grid.nz):
                raise ValueError('Number of grid points in ' + fname + ' is different from that in amr_grid.inp\n',
                                 (self.grid.nx * self.grid.ny * self.grid.nz), hdr[2])

            if hdr[1] == 8:
                data = np.fromfile(fname, count=-1, dtype=np.float64)
            elif hdr[1] == 4:
                data = np.fromfile(fname, count=-1, dtype=float)
            else:
                raise TypeError('Unknown datatype in ' + fname + ' radmc3dPy supports only 4 byte float or '
                                + ' 8 byte double types.')

            data = np.reshape(data[4:], [hdr[3], self.grid.nz, self.grid.ny, self.grid.nx])
            data = np.swapaxes(data, 0, 3)
            data = np.swapaxes(data, 1, 2)
        else:
            hdr = np.fromfile(fname, count=3, sep="\n", dtype=int)

            if (self.grid.nx * self.grid.ny * self.grid.nz) != hdr[1]:
                raise ValueError('Number of grid points in amr_grid.inp is not equal to that in ' + fname,
                                 (self.grid.nx * self.grid.ny * self.grid.nz), hdr[1])
            else:

                data = np.fromfile(fname, count=-1, sep="\n", dtype=np.float64)
                data = np.reshape(data[3:], [hdr[2], self.grid.nz, self.grid.ny, self.grid.nx])
                # We need to change the axis orders as Numpy always reads  in C-order while RADMC-3D
                # uses Fortran-order
                data = np.swapaxes(data, 0, 3)
                data = np.swapaxes(data, 1, 2)

                self.csdens = data

    def writeStellarsrcDensity(self, fname='', binary=False):
        """Writes the stellar density of a continuous starlike source.

        Parameters 
        ----------

        fname : str, optional
                Name of the file into which the stellar templates will be written. If omitted the default
                'stellarsrc_templates.inp' will be used.

        binary : bool, optional
                If True the output will be written in a C-style binary stream, if False the output will be 
                formatted ASCII
        """

        # First check if we'd need to write anything at al
        # Both the stellar temperature and the stellar template spectrum arrays are empty
        # So the continuous starlike source needs to be deactivated
        # Let's check if there are any files, and if so ask them to be deleted
        if len(self.cststar) == 0:
            if len(self.cstemp) == 0:
                if (os.path.exists('stellarsrc_density.inp')) | (os.path.exists('stellarsrc_density.binp')):
                    print('The continuous starlike source seems to be inactive (zero input luminosity still '
                          + 'stellarsrc_density.inp/stellarsrc_density.binp file is present in the current '
                          + 'working directory.')
                    dum = input('Can it be deleted (yes/no)')
                    if dum.strip().lower()[0] == 'y':
                        os.system('rm stellarsrc_density.inp')
                        os.system('rm stellarsrc_density.binp')
                    return
                return
            else:
                if np.abs(self.cstemp).max() == 0.:
                    if (os.path.exists('stellarsrc_density.inp')) | (os.path.exists('stellarsrc_density.binp')):
                        print('The continuous starlike source seems to be inactive (zero input luminosity)'
                              + ' still stellarsrc_density.inp/stellarsrc_density.binp file is present in the current '
                              + ' working directory.')
                        dum = input('Can it be deleted (yes/no)')
                        if dum.strip().lower()[0] == 'y':
                            os.system('rm stellarsrc_density.inp')
                            os.system('rm stellarsrc_density.binp')
                        return
                    return

        else:
            if np.abs(self.cststar).max() == 0.:
                if (os.path.exists('stellarsrc_density.inp')) | (os.path.exists('stellarsrc_density.binp')):
                    print('The continuous starlike source seems to be inactive (zero input luminosity)'
                          + ' still stellarsrc_density.inp/stellarsrc_density.binp file is present in the current '
                          + ' working directory.')
                    dum = input('Can it be deleted (yes/no)')
                    if dum.strip().lower()[0] == 'y':
                        os.system('rm stellarsrc_density.inp')
                        os.system('rm stellarsrc_density.binp')
                    return
                return

        if binary:
            if fname.strip() == '':
                fname = 'stellarsrc_density.binp'
            print('Writing ' + fname)
            wfile = open(fname, 'w')
            hdr = np.array([1, 8, self.grid.nx * self.grid.ny * self.grid.nz, self.csntemplate], dtype=int)
            hdr.tofile(wfile)
            data = np.array(self.csdens)
            data = np.swapaxes(data, 0, 3)
            data = np.swapaxes(data, 1, 2)
            data.tofile(wfile)

        else:
            if fname.strip() == '':
                fname = 'stellarsrc_density.inp'

            with open(fname, 'w') as wfile:
                hdr = np.array([1, self.grid.nx * self.grid.ny * self.grid.nz, self.csntemplate], dtype=int)
                hdr.tofile(wfile, sep=" ", format="%d\n")
                data = np.array(self.csdens)
                data = np.swapaxes(data, 0, 3)
                data = np.swapaxes(data, 1, 2)
                data.tofile(wfile, sep=" ", format="%.9e\n")


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


def readOpac(ext=None, idust=None, scatmat=None, old=False):
    """Reads the dust opacity files.
    This function is an interface to radmc3dDustOpac.readOpac()

    Parameters
    ----------
    ext   : list
            Each element of the list is be a string, the file name extension 
            (file names should look like 'dustkappa_ext.inp')

    idust : list
            Each element of the list is an integer, the index of the dust species in the master opacity file 
            (dustopac.inp')

    scatmat: list
            If specified, its elements should be booleans indicating whether the opacity file 
            contains also the full scattering matrix (True) or only dust opacities (False)

    old   : bool, optional
            If set to True the file format of the previous, 2D version of radmc will be used

    Returns
    -------
        Returns an instance of the radmc3dDustOpac class 
    """

    res = radmc3dDustOpac()
    res.readOpac(ext=ext, idust=idust, scatmat=scatmat, old=old)

    return res


def readData(ddens=False, dtemp=False, gdens=False, gtemp=False, gvel=False, ispec=None, vturb=False, grid=None,
             binary=True, old=False, octree=False):
    """Reads the physical variables of the model (e.g. density, velocity, temperature).

    Parameters
    ----------

    ddens : bool
            If True dust density will be read (all dust species and grain sizes)

    dtemp : bool
            If True dust temperature will be read (all dust species and grain sizes)

    gdens : bool
            If True gas density will be read (NOTE: the gas density will be number density in 1/cm^3)

    gtemp : bool
            If True gas temperature will be read (all dust species and grain sizes)

    gvel  : bool
            If True the velocity field will be read

    ispec : str
            Name of the molecule in the 'molecule_ispec.inp' filename

    vturb : bool
            If True the microturbulent velocity field will be read

    grid  : radmc3dGrid
            An instance of radmc3dGrid containing the spatial and frequency grid of the model. If the grid
            is passed to the function it will not be read again from file. This can be useful for octree
            models to save time. 

    old   : bool, optional
            If set to True the file format of the previous, 2D version of radmc will be used

    binary: bool
            Set it to True for C-style binary and False for formatted ASCII files

    octree: bool
            True for models with octree AMR and False for models with regular grid

    Returns
    -------
    Returns an instance of the radmc3dData class 
    """

    if grid is not None:
        res = radmc3dData(grid=grid)
    else:
        res = radmc3dData()

        if octree:
            res.grid = radmc3dOctree()
            res.grid.readSpatialGrid()
        else:
            res.grid = radmc3dGrid()
            res.grid.readSpatialGrid(old=old)

    if ddens:
        res.readDustDens(binary=binary, old=old, octree=octree)
    if dtemp:
        res.readDustTemp(binary=binary, old=old, octree=octree)
    if gvel:
        res.readGasVel(binary=binary, octree=octree)
    if gtemp:
        res.readGasTemp(binary=binary, octree=octree)
    if vturb:
        res.readVTurb(binary=binary, octree=octree)
    if gdens:
        if not ispec:
            raise ValueError('Unknown ispec.\n'
                             + 'No gas species is specified!\n'
                             + 'The ispec input keyword should be set to the name of the gas species as it appears in\n'
                             + 'numberdens_gasspecname.inp')
        else:
            res.readGasDens(ispec=ispec, binary=binary, octree=octree)

    return res

def readStars(fname=''):
    """
    Reads the data (mass, radius, temperature, spectrum) of discrete stellar sources

    Parameters
    ----------
    fname       : str
                  Name of the file to be read (if omitted the default value is stars.inp)

    Returns
    -------
    An instance of radmc3dRadSources containing the stellar data
    """
    if fname == '':
        fname = 'stars.inp'

    res = radmc3dRadSources()
    res.readStarsinp(fname=fname)

    return res

def readGrid(sgrid=True, wgrid=True, sgrid_fname='amr_grid.inp', wgrid_fname='wavelength_micron.inp'):
    """Reads the spatial and frequency grid.
    This function is an interface to radmc3dGrid.readGrid().

    Parameters
    ----------

    sgrid       : bool
                  If True the spatial grid will be read

    wgrid       : bool
                  If True the wavelength grid will be read
    
    sgrid_fname : str
                  File containing the spatial grid (default: amr_grid.inp)
    
    wgrid_fname : str
                  File containing the wavelength grid (default: wavelength_micron.inp)

    Returns
    -------

    Returns an instance of the radmc3dGrid (for regular grid) or radmc3dOctree (for octree AMR) class 
    """

    grid = None

    if wgrid:
        grid = radmc3dGrid()
        grid.readWavelengthGrid(fname=wgrid_fname)

    if sgrid:
        if grid is None:
            grid = radmc3dGrid()
        #
        # Check the grid type
        #
        hdr = np.fromfile(sgrid_fname, count=7, sep="\n", dtype=np.int)
        if hdr[1] == 0:
            grid = radmc3dGrid()
        elif hdr[1] == 1:
            grid = radmc3dOctree()
        else:
            raise ValueError('Unsupported amr_style' + ("%d" % hdr[1]) + '\n '
                             + 'Only regular (0) or octree-like (1) AMR styles are supported')
        grid.readSpatialGrid(fname=sgrid_fname)

    return grid


def readParams():
    """Reads the problem_params.inp file.
    This function is an interface to radmc3dPar.readPar().

    Returns
    -------
    Returns an instance of the radmc3dPar class 

    """

    dum = radmc3dPar()
    dum.readPar()
    return dum


def writeDefaultParfile(model='', fname=''):
    """Writes a parameter file (problem_params.inp) with default parameters for a given model.

    Parameters
    ----------

    model : str
            Name of the model whose parameter should be written to the file

    fname : str, optional
            Name of the parameter file to be written (if omitted problem_params.inp will be used)
    """

    if model == '':
        raise ValueError('Unknown model. \n No model name is given. ')

    dum = radmc3dPar()
    dum.loadDefaults(model=model)
    dum.writeParfile(fname=fname)


def readSpectrum(fname='', old=False):
    """Reads the spectrum / SED


    Parameters
    -----------
    fname : str, optional
            Name of the file to be read

    old   : bool, optional
            If set to True the file format of the previous, 2D version of radmc will be used


    Returns
    -------

        Returns an ndarray with [Nwavelength, 2] dimensions 
        [Nwavelength,0] is the wavelength / velocity and
        [Nwavelength,1] is the flux density

    """

    if not old:
        if fname.strip() == '':
            fname = 'spectrum.out'

        with open(fname, 'r') as rfile:
            # Read the format number
            dum = rfile.readline()
            # Read the number of wavelengths
            nwav = int(rfile.readline())
            # Read a blank line
            dum = rfile.readline()

            res = np.zeros([nwav, 2], dtype=np.float64)
            for iwav in range(nwav):
                dum = rfile.readline().split()
                res[iwav, 0] = float(dum[0])
                res[iwav, 1] = float(dum[1])

    else:
        if fname.strip() == '':
            fname = 'spectrum.dat'

        with open(fname, 'r') as rfile:
            # Read the number of wavelengths
            nwav = int(rfile.readline())
            rfile.readline()
            cc = 29979245800.

            res = np.zeros([nwav, 2], dtype=float)
            for iwav in range(nwav):
                dum = rfile.readline().split()
                res[iwav, 0] = cc / float(dum[0]) * 1e4
                res[iwav, 1] = float(dum[1])

    return res


def getDensVstruct(data=None, vmean_temp=False, ispec_tgas=0, gsize=None, idust=None, mstar=None):
    """Calculates the vertical hydrostatic equilibrium

    Parameters
    ----------
    data        : radmc3dData
                  An instance of the radmc3DData class containing the density structure of the model

    vmean_temp  : bool
                  If True (T(z) = T(-z) = 0.5*(T(z) + T(-z))) if False (T(z)!=T(-z)) 

    idust       : list
                  List of dust indices whose structure must be calculated

    mstar       : float
                  Stellar mass

    ispec_tgas  : int
                  Index of dust species whose temperature is taken to be the gas temperature

    gsize       : ndarray
                  Dust grain sizes - If specified, the gas temperature is calculated as the average temperature
                  of all dust grains in the grid cell weighted by the total surface area of dust grains with given
                  size - NOTE: this approach assumes that all dust grains of a given size have the same bulk density

    Returns
    -------
    Returns an ndarray with the dust density
    """

    # Fix the mean molecular weight to 2.3
    mu = 2.3

    # Pre-calculate some constants
    A = mu * nc.mp * nc.gg * mstar / nc.kk
    cost = np.cos(data.grid.y)
    costi = np.cos(data.grid.yi)

    if mstar is None:
        raise ValueError('Unkonwn mstar. \n The stellar mass is required to calculate the '
                         + ' vertical structure of the disk')

    if idust is None:
        print(' Unknown idust. No dust index was given for which the vertical structure should be calculated, '
              + ' so we do for all dust species')
        idust = range(data.rhodust.shape[3])
    else:
        if isinstance(idust, int) | isinstance(idust, float):
            idust = [int(idust)]

    # To improve the smoothness of the temperature structure, if the density structure is
    #  symmetric to the disk midplane we use T_new(theta) = T_new(pi-theta) = 0.5 * (T(theta) + T(pi-theta))
    if vmean_temp:
        if abs(data.grid.yi[data.grid.nyi - 1] - np.pi / 2.) < 1e-8:
            raise RuntimeError("Cannot average temperature in the vertical direction if theta mirroring is active")
        else:
            print(' Smoothing the vertical temperature structure by averaging the temperature of the two half \n'
                  ' planes above and below the disk midplane')
            dusttemp_dummy = np.zeros(data.dusttemp.shape, dtype=np.float64)
            for iy in range(int(data.grid.ny / 2)):
                print(iy)
                dusttemp_dummy[:, iy, :, :] = 0.5 * (data.dusttemp[:, iy, :, :]
                                                     + data.dusttemp[:, data.grid.ny - 1 - iy, :, :])
                dusttemp_dummy[:, data.grid.ny - 1 - iy, :, :] = dusttemp_dummy[:, iy, :, :]
    # Calculate the vertical hydrostatic equilibrium for the two half space (z<0, z>0) separately
    else:
        dusttemp_dummy = data.dusttemp[:, :, :, ispec_tgas]

    # rho_new = np.zeros(data.rhodust.shape, dtype=np.float64)
    rho_new = np.array(data.rhodust)

    if len(gsize) != 0:
        dusttemp = np.zeros([data.grid.nx, data.grid.ny, data.grid.nz], dtype=np.float64)
        w = np.zeros(data.rhodust.shape, dtype=np.float64)
        for ispec in idust:
            w[:, :, :, ispec] = gsize[ispec] ** 2 * (data.rhodust[:, :, :, ispec] / gsize[ispec] ** 3)

        wnorm = w.sum(3)
        for ispec in idust:
            w[:, :, :, ispec] = w[:, :, :, ispec] / wnorm

        for ispec in idust:
            dusttemp = dusttemp + data.dusttemp[:, :, :, ispec] * w[:, :, :, ispec]
    else:
        dusttemp = np.array(dusttemp_dummy)

    # Loop over all dust species where we should calculate the vertical structure
    for ispec in idust:
        rho_new[:, :, :, ispec] = 0.
        for ir in range(data.grid.nx):
            print(ir, data.grid.nx - 1)
            r = data.grid.x[ir]
            z = r * cost
            zi = r * costi
            dz = z[:-1] - z[1:]
            const = A / r ** 3

            # Do we have theta mirroring active?
            if abs(data.grid.yi[data.grid.nyi - 1] - np.pi / 2.) < 1e-8:
                for ip in range(data.grid.nz):
                    # dlgrho = np.log(data.rhodust[ir, 1:, ip, ispec]) - np.log(data.rhodust[ir, :-1, ip, ispec])
                    temp = dusttemp[ir, :, ip]

                    it = data.grid.ny - 1
                    temp[it] = 0.5 * (temp[it] + temp[it - 1])

                    dlgtemp = np.log(temp[1:]) - np.log(temp[:-1])
                    zpt = z / temp
                    zpt = 0.5 * (zpt[1:] + zpt[:-1])

                    # Calculate the normalized (rho[z=0] = 1.0) density
                    rho_new[ir, data.grid.ny - 1, ip, ispec] = 1.0

                    for it in range(data.grid.ny - 1, 0, -1):
                        rho_new[ir, it, ip, ispec] = rho_new[ir, it + 1, ip, ispec] * np.exp(
                            -(const * zpt[it] + dlgtemp[it] / dz[it]) * dz[it])

                    rho_new = rho_new.clip(1e-90, 1e90)

                    # Now re-normalize the surface density to the input value
                    sigma = (data.rhodust[ir, :, ip, ispec] * (zi[1:] - zi[:-1])).sum()
                    sigma_new = (rho_new[ir, :, ip, ispec] * (zi[1:] - zi[:-1])).sum()

                    rho_new[ir, :, ip, ispec] = rho_new[ir, :, ip, ispec] * sigma / sigma_new

            else:
                for ip in range(data.grid.nz):
                    # dlgrho = np.log(data.rhodust[ir, 1:, ip, ispec]) - np.log(data.rhodust[ir, :-1, ip, ispec])
                    # if len(ispec_weights) != 0:
                    #     temp = (dusttemp[ir, :, ip, ispec] * ispec_weights).sum()
                    # else:
                    #     temp = dusttemp[ir, :, ip, ispec]

                    temp = dusttemp[ir, :, ip]
                    dlgtemp = np.log(temp[1:]) - np.log(temp[:-1])
                    zpt = z / temp
                    zpt = 0.5 * (zpt[1:] + zpt[:-1])

                    # Calculate the normalized (rho[z=0] = 1.0) density
                    rho_new[ir, data.grid.ny / 2 - 1, ip, ispec] = 1.0
                    rho_new[ir, data.grid.ny / 2, ip, ispec] = 1.0

                    for it in range(int(data.grid.ny / 2), 0, -1):
                        rho_new[ir, it - 1, ip, ispec] = rho_new[ir, it, ip, ispec] \
                                                         * np.exp(-(const * zpt[it] + dlgtemp[it] / dz[it]) * dz[it])
                    for it in range(int(data.grid.ny / 2), data.grid.ny - 1, -1):
                        rho_new[ir, it, ip, ispec] = rho_new[ir, it - 1, ip, ispec] \
                                                     * np.exp((const * zpt[it - 1] + dlgtemp[it - 1]
                                                               / dz[it - 1]) * dz[it - 1])

                    rho_new = rho_new.clip(1e-90, 1e90)

                    # Now re-normalize the surface density to the input value
                    sigma = (data.rhodust[ir, :, ip, ispec] * (zi[1:] - zi[:-1])).sum()
                    sigma_new = (rho_new[ir, :, ip, ispec] * (zi[1:] - zi[:-1])).sum()

                    rho_new[ir, :, ip, ispec] = rho_new[ir, :, ip, ispec] * sigma / sigma_new

            print(rho_new[ir, int(data.grid.ny / 2) - 1, ip, ispec])

    return rho_new


def readMol(mol='', fname=''):
    """ Wrapper around the radmc3dMolecule.read() method

       Parameters
       ----------
       mol             : str
                        molecule name (e.g. 'co') if the file name is in the form of 'molecule_<mol>.inp'

       fname           : str
                        full file name
    """

    m = radmc3dMolecule()
    if m.read(mol=mol, fname=fname) is True:
        return m
    else:
        return


def plotSpectrum(a, ev=False, kev=False, micron=False, jy=False, lsun=False,
                 lnu=False, nulnu=False, fnu=False, nufnu=False, dpc=1.e0,
                 oplot=False, xlg=False, ylg=False, obs=False,
                 mol=None, ilin=None):
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

    fnu             : bool
                     True --> F_nu in units of erg/s/cm^2/Hz(default L_nu)

    nufnu           : bool
                     True --> nu*F_nu in units of erg/s/cm^2 (default L_nu)

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

    mol             : radmc3dMolecule
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
    lam = a[:, 0]
    fluxnu = a[:, 1]
    #
    # Calculate frequency in Hz
    #
    freq = 1e4 * nc.cc / lam
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
    sed = True
    ylum = False
    # The flags:
    if jy:
        sed = False
    if fnu:
        sed = False
        ylum = False
    if lnu:
        sed = False
        ylum = True
    if nulnu:
        sed = True
        ylum = True
    if fnu:
        sed = False
        ylum = False
    if nufnu:
        sed = True
        ylum = False
    if jy:
        ylum = False
    if lsun:
        ylum = True
        sed = True
    #
    # If ilin is set, then override the above and use instead the line
    # as a reference and use km/s as x-axis
    #
    if ilin is not None:
        if mol is None:
            raise ValueError("Unknown mol. If ilin is set, the molecular data should also be provided as mol=...")
        else:
            freq0 = mol.freq[ilin - 1]
            xcoord = nc.cc * (freq0 - freq) / freq0 / 1.e5
            xtitle = '$\Delta v [\mathrm{km/h}]$'
    #
    # Which plot to make? Lum or flux?
    #
    if not ylum:
        #
        # Plot spectrum as flux at a certain distance
        #
        if not obs:
            distfact = 1.0 / (dpc ** 2)
        else:
            distfact = 1.0
        #
        # Set the vertical axis name
        #
        if not jy:
            if not sed:
                lumfact = 1.0
                ytitle = '$F_{\\nu}\; [\mathrm{erg}\,\mathrm{cm}^{-2}\,\mathrm{Hz}^{-1}\, \mathrm{s}^{-1}]$'
            else:
                lumfact = 1.0 * freq
                ytitle = '$\\nu F_{\\nu}\; [\mathrm{erg}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}]$'
        else:
            if not sed:
                lumfact = 1e+23
                ytitle = '$F_{\\nu} [Jy]$'
            else:
                lumfact = 1e+23 * freq
                ytitle = '$\\nu F_{\\nu} [JyHz]$'
    else:
        #
        # Plot spectrum as luminosity
        #
        if not obs:
            distfact = 1.1965280793e38  # = 4*pi*(1 parsec)^2 = 1.19d38 cm^2
        else:
            distfact = dpc ** 2 * 1.1965280793e38

        if not sed:
            lumfact = 1.e0
            ytitle = 'L_{\\nu}\; [\mathrm{erg}\,\mathrm{Hz}^{-1}\, \mathrm{s}^{-1}]'
        else:
            if not lsun:
                lumfact = 1.0 * freq
                ytitle = '\\nu L_{\\nu}\; [\mathrm{erg}\, \mathrm{s}^{-1}]'
            else:
                lumfact = freq * 2.5956986e-34
                ytitle = '\\nu L_{\\nu}\; [L_{\odot}]'

    #
    # The data on the y axis
    #
    ycoord = distfact * lumfact * fluxnu
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
    plb.plot(xcoord, ycoord)


def gmass(x=None, y=None, z=None, dx=None, dy=None, dz=None, model=None, ppar=None, **kwargs):
    """
    Example function to be used as decision function for resolving cells in tree building. It calculates the gas density
    at a random sample of coordinates within a given cell than take the ratio of the max/min density. If it is larger
    than a certain threshold value it will return True (i.e. the cell should be resolved) if the density variation is 
    less than the threshold it returns False (i.e. the cell should not be resolved)

    Parameters
    ----------

    x       : ndarray
              Cell centre coordinates of the cells in the first dimension

    y       : ndarray
              Cell centre coordinates of the cells in the second dimension

    z       : ndarray
              Cell centre coordinates of the cells in the third dimension

    dx      : ndarray
              Half size of the cells in the first dimension

    dy      : ndarray
              Half size of the cells in the second dimension

    dz      : ndarray
              Half size of the cells in the third dimension

    model   : object
              A radmc3dPy model (must contain a getGasDensity() function) 

    ppar    : dictionary
              All parameters of the problem (from the problem_params.inp file). It is not used here, but must be present 
              for compatibility reasons.

    **kwargs: dictionary
              Parameters used to decide whether the cell should be resolved. It should contain the following keywords; 
              'nsample', which sets the number of random points the gas desity is sampled at within the cell and 
              'threshold' that sets the threshold value for max(gasdens)/min(gasdens) above which the cell should 
              be resolved.
    """

    ncell = x.shape[0]
    rho = np.zeros([ncell, kwargs['nsample']], dtype=np.float64)

    for isample in range(int(kwargs['nsample'])):
        xoffset = (np.random.random_sample(ncell) - 0.5) * dx * 4.0
        yoffset = (np.random.random_sample(ncell) - 0.5) * dy * 4.0
        zoffset = (np.random.random_sample(ncell) - 0.5) * dz * 4.0
        rho[:, isample] = model.getGasDensity(x + xoffset, y + yoffset, z + zoffset, ppar=ppar)

    mass = rho.max(1) * dx * dy * dz * 8.0
    jj = (mass > ppar['threshold'])

    decision = np.zeros(ncell, dtype=bool)
    if True in jj:
        decision[jj] = True

    return decision


def gdensMinMax(x=None, y=None, z=None, dx=None, dy=None, dz=None, model=None, ppar=None, **kwargs):
    """
    Example function to be used as decision function for resolving cells in tree building. It calculates the gas density
    at a random sample of coordinates within a given cell than take the ratio of the max/min density. If it is larger
    than a certain threshold value it will return True (i.e. the cell should be resolved) if the density variation is 
    less than the threshold it returns False (i.e. the cell should not be resolved)

    Parameters
    ----------

    x       : ndarray
              Cell centre coordinates of the cells in the first dimension

    y       : ndarray
              Cell centre coordinates of the cells in the second dimension

    z       : ndarray
              Cell centre coordinates of the cells in the third dimension

    dx      : ndarray
              Half size of the cells in the first dimension

    dy      : ndarray
              Half size of the cells in the second dimension

    dz      : ndarray
              Half size of the cells in the third dimension

    model   : object
              A radmc3dPy model (must contain a getGasDensity() function) 

    ppar    : dictionary
              All parameters of the problem (from the problem_params.inp file). It is not used here, but must be present 
              for compatibility reasons.

    **kwargs: dictionary
              Parameters used to decide whether the cell should be resolved. It should the following keywords; 
              'nsample', which sets the number of random points the gas desity is sampled at within the cell and 
              'threshold' that sets the threshold value for max(gasdens)/min(gasdens) above which the cell should 
              be resolved.
    """

    ncell = x.shape[0]
    rho = np.zeros([ncell, kwargs['nsample']], dtype=np.float64)

    for isample in range(kwargs['nsample']):
        xoffset = (np.random.random_sample(ncell) - 0.5) * dx * 4.0
        yoffset = (np.random.random_sample(ncell) - 0.5) * dy * 4.0
        zoffset = (np.random.random_sample(ncell) - 0.5) * dz * 4.0
        rho[:, isample] = model.getGasDensity(x + xoffset, y + yoffset, z + zoffset, ppar=ppar)

    rho_max = rho.max(axis=1)
    rho_min = rho.min(axis=1)
    # jj      = ((rho_max/rho_min)>ppar['threshold'])
    jj = ((rho_max - rho_min) / rho_max > ppar['threshold'])

    decision = np.zeros(ncell, dtype=bool)
    if True in jj:
        decision[jj] = True
    return decision


def findContainerLeafID(cellCRD=None, cellHW=None, xi=None, yi=None, zi=None, childID=None, isLeaf=None, nChild=None,
                        crd=None):
    """
    Function to find the tree index of a leaf cell containing a given point in space, i.e. if the following is true : 
    xcell - dxcell <= xpoint < xcell + dxcell for each dimension. This function is to be used in multiprocessing.

    Parameters
    ----------

    cellCRD         : ndarray
                      Array with dimensions [ncell, 3] containing the cell centre coordiantes of the tree

    cellHW          : ndarray
                      Array with dimensions [ncell, 3] containing the half width of cells in the tree

    xi              : ndarray
                      Array of cell interface indices in the base grid in the first dimension

    yi              : ndarray
                      Array of cell interface indices in the base grid in the second dimension

    zi              : ndarray
                      Array of cell interface indices in the base grid in the third dimension

    childID         : ndarray
                      Child index array

    isLeaf          : ndarray
                      Boolean array containing the node type for each cell (True - leaf, False - branch) 

    nChild          : int
                      Number of children (8,4,2 for 3,2,1 active dimensions)

    crd             : ndarray
                      Array of length 3 containing the coordinates of the point whose container
                      leaf is to be found

    Returns
    -------
    Tree index of the container leaf if it is found. If the point is outside of the base grid -1 is returned.

    """

    leafID = -1

    if (crd[0] < xi[0]) | (crd[0] > xi[-1]):
        return leafID
    if (crd[1] < yi[0]) | (crd[1] > yi[-1]):
        return leafID
    if (crd[2] < zi[0]) | (crd[2] > zi[-1]):
        return leafID

    ix = np.searchsorted(xi, crd[0])
    iy = np.searchsorted(yi, crd[1])
    iz = np.searchsorted(zi, crd[2])

    if xi[ix] != crd[0]:
        ix -= 1
    if yi[iy] != crd[1]:
        iy -= 1
    if zi[iz] != crd[2]:
        iz -= 1

    nxRoot = xi.shape[0] - 1
    nyRoot = yi.shape[0] - 1
    nzRoot = zi.shape[0] - 1

    if crd[0] == xi[-1]:
        ix = nxRoot - 1
    if crd[1] == yi[-1]:
        iy = nyRoot - 1
    if crd[2] == zi[-1]:
        iz = nzRoot - 1

    ind = iz * nyRoot * nxRoot + iy * nxRoot + ix
    dum = findContainerLeafIDRec(cellCRD[:, 0], cellCRD[:, 1], cellCRD[:, 2], cellHW[:, 0], cellHW[:, 1], cellHW[:, 2],
                                 childID, isLeaf, nChild, crd, ind)
    if dum is None:
        leafID = -1
    else:
        leafID = dum

    return leafID


def findContainerLeafIDRec(x=None, y=None, z=None, dx=None, dy=None, dz=None, childID=None, isLeaf=None, nChild=None,
                           crd=(), cellID=None):
    """
    Recursive function to find the leaf cell in the tree that contains a given point in space

    Parameters
    ----------

    x                 : ndarray
                        Tree cell center array in the first dimension

    y                 : ndarray
                        Tree cell center array in the second dimension

    z                 : ndarray
                        Tree cell center array in the tird dimension

    dx                : ndarray
                        Tree cell halfwidth array in the first dimension

    dy                : ndarray
                        Tree cell halfwidth array in the second dimension

    dz                : ndarray
                        Tree cell halfwidth array in the third dimension

    childID           : list
                        List of children indices. Each list element is an ndarray with nChild elements containing 
                        the child indices

    isLeaf            : ndarray
                        Boolean array for the cell type (True - leaf, False - branch)

    nChild            : int
                        Nr of children (i.e. 8, 4, or 2 for 3, 2, 1 active dimensions, respectively)

    crd               : ndarray
                        Three element list/tuple/array containing the point coordinates

    cellID            : int
                        Index of cell to be tested
    """

    xmin = x[cellID] - dx[cellID]
    xmax = x[cellID] + dx[cellID]
    ymin = y[cellID] - dy[cellID]
    ymax = y[cellID] + dy[cellID]
    zmin = z[cellID] - dz[cellID]
    zmax = z[cellID] + dz[cellID]

    if isLeaf[cellID]:
        if (((crd[0] >= xmin) & (crd[0] < xmax)) &
                ((crd[1] >= ymin) & (crd[1] < ymax)) &
                ((crd[2] >= zmin) & (crd[2] < zmax))):
            return cellID
        else:
            return None

    else:
        dum = None
        for i in range(nChild):
            dum = findContainerLeafIDRec(x, y, z, dx, dy, dz, childID, isLeaf, nChild, crd, childID[cellID][i])
            if dum is not None:
                break

        return dum


def interpolateOctree(data=None, x=None, y=None, z=None, var=None, nproc=1):
    """
    Nearest neighbour inteprolation on an octree

    data        : radmc3dData
                  Data container

    x           : ndarray
                  Coordiantes of the point to be interpolated on in the first dimension

    y           : ndarray
                  Coordiantes of the point to be interpolated on in the second dimension

    z           : ndarray
                  Coordiantes of the point to be interpolated on in the third dimension

    var         : list
                  Name of the variables to be interpolated, supported names are:
                  ddens, dtemp, gdens, ndens, gtemp, gvel, vturb

    nproc       : int
                  Number of processes to be used (for parallel computing) 

    Returns:
    --------
    A dictionary with the interpolated fields 

    """

    if var is None:
        var = ['ddens']

    if nproc == 1:
        print("Nearest neighbour interpolation using " + ("%d" % nproc) + ' process')
    else:
        print("Nearest neighbour interpolation using " + ("%d" % nproc) + ' processes')
    nx = x.shape[0]
    ny = y.shape[0]
    nz = z.shape[0]

    npoint = nx * ny * nz

    if nproc > 1:
        cellCRD = np.zeros([data.grid.nCell, 3], dtype=np.float64)
        cellHW = np.zeros([data.grid.nCell, 3], dtype=np.float64)

        cellCRD[:, 0] = data.grid.x
        cellCRD[:, 1] = data.grid.y
        cellCRD[:, 2] = data.grid.z
        cellHW[:, 0] = data.grid.dx
        cellHW[:, 1] = data.grid.dy
        cellHW[:, 2] = data.grid.dz
        childID = np.array(data.grid.childID)

        chunkSize = int(np.ceil(float(npoint) / nproc))

        crdList = np.zeros([npoint, 3], dtype=np.float64)
        for ind in range(npoint):
            ix = int(np.floor(ind / ny / nz))
            iy = int(np.floor((ind - ix * ny * nz) / nz))
            iz = int(ind - ix * ny * nz - iy * nz)

            crdList[ind, 0] = x[ix]
            crdList[ind, 1] = y[iy]
            crdList[ind, 2] = z[iz]

        pool = Pool(processes=nproc)
        target = partial(findContainerLeafID, cellCRD, cellHW, data.grid.xi, data.grid.yi, data.grid.zi, childID,
                         data.grid.isLeaf, data.grid.nChild)
        res = pool.map(target, crdList, chunksize=chunkSize)
        res = np.array(res)
        pool.close()

        idata = {'cellID': None, 'rhodust': None, 'dusttemp': None, 'rhogas': None, 'ndens_mol': None, 'gtemp': None,
                 'vturb': None, 'gvel': None}

        idata['cellID'] = res

        if 'ddens' in var:
            ndust = data.rhodust.shape[1]
            idata['rhodust'] = np.zeros([nx, ny, nz, ndust], dtype=np.float64)
            for ind in range(npoint):
                ix = int(np.floor(ind / ny / nz))
                iy = int(np.floor((ind - ix * ny * nz) / nz))
                iz = int(ind - ix * ny * nz - iy * nz)

                if res[ind] >= 0:
                    idata['rhodust'][ix, iy, iz, :] = data.rhodust[data.grid.leafID[res[ind]], :]
                else:
                    idata['rhodust'][ix, iy, iz, :] = 0

        if 'dtemp' in var:
            ndust = data.dusttemp.shape[1]
            idata['dusttemp'] = np.zeros([nx, ny, nz, ndust], dtype=np.float64)
            for ind in range(npoint):
                ix = int(np.floor(ind / ny / nz))
                iy = int(np.floor((ind - ix * ny * nz) / nz))
                iz = int(ind - ix * ny * nz - iy * nz)

                if res[ind] is not None:
                    idata['dusttemp'][ix, iy, iz, :] = data.dusttemp[data.grid.leafID[res[ind]], :]
                else:
                    idata['dusttemp'][ix, iy, iz, :] = 0

        if 'gdens' in var:
            idata['rhogas'] = np.zeros([nx, ny, nz], dtype=np.float64)
            for ind in range(npoint):
                ix = int(np.floor(ind / ny / nz))
                iy = int(np.floor((ind - ix * ny * nz) / nz))
                iz = int(ind - ix * ny * nz - iy * nz)

                if res[ind] is not None:
                    idata['rhogas'][ix, iy, iz] = data.rhogas[data.grid.leafID[res[ind]]]
                else:
                    idata['rhogas'][ix, iy, iz] = 0

        if 'ndens' in var:
            idata['ndens_mol'] = np.zeros([nx, ny, nz], dtype=np.float64)
            for ind in range(npoint):
                ix = int(np.floor(ind / ny / nz))
                iy = int(np.floor((ind - ix * ny * nz) / nz))
                iz = int(ind - ix * ny * nz - iy * nz)

                if res[ind] is not None:
                    idata['ndens_mol'][ix, iy, iz] = data.ndens_mol[data.grid.leafID[res[ind]]]
                else:
                    idata['ndens_mol'][ix, iy, iz] = 0

        if 'gtemp' in var:
            idata['gastemp'] = np.zeros([nx, ny, nz], dtype=np.float64)
            for ind in range(npoint):
                ix = int(np.floor(ind / ny / nz))
                iy = int(np.floor((ind - ix * ny * nz) / nz))
                iz = int(ind - ix * ny * nz - iy * nz)

                if res[ind] is not None:
                    idata['gastemp'][ix, iy, iz] = data.gastemp[data.grid.leafID[res[ind]]]
                else:
                    idata['gastemp'][ix, iy, iz] = 0

        if 'vturb' in var:
            idata['vturb'] = np.zeros([nx, ny, nz], dtype=np.float64)
            for ind in range(npoint):
                ix = int(np.floor(ind / ny / nz))
                iy = int(np.floor((ind - ix * ny * nz) / nz))
                iz = int(ind - ix * ny * nz - iy * nz)

                if res[ind] is not None:
                    idata['vturb'][ix, iy, iz] = data.vturb[data.grid.leafID[res[ind]]]
                else:
                    idata['vturb'][ix, iy, iz] = 0

        if 'gvel' in var:
            idata['gvel'] = np.zeros([nx, ny, nz, 3], dtype=np.float64)
            for ind in range(npoint):
                ix = int(np.floor(ind / ny / nz))
                iy = int(np.floor((ind - ix * ny * nz) / nz))
                iz = int(ind - ix * ny * nz - iy * nz)

                if res[ind] is not None:
                    idata['gasvel'][ix, iy, iz, :] = data.gasvel[data.grid.leafID[res[ind]], :]
                else:
                    idata['gasvel'][ix, iy, iz, :] = 0
    else:

        idata = {}
        if 'ddens' in var:
            ndust = data.rhodust.shape[1]
            idata['rhodust'] = np.zeros([nx, ny, nz, ndust], dtype=np.float64)
        if 'dtemp' in var:
            ndust = data.dusttemp.shape[1]
            idata['dusttemp'] = np.zeros([nx, ny, nz, ndust], dtype=np.float64)
        if 'gdens' in var:
            idata['rhogas'] = np.zeros([nx, ny, nz], dtype=np.float64)
        if 'ndens' in var:
            idata['ndens_mol'] = np.zeros([nx, ny, nz], dtype=np.float64)
        if 'gtemp' in var:
            idata['gastemp'] = np.zeros([nx, ny, nz], dtype=np.float64)
        if 'vturb' in var:
            idata['vturb'] = np.zeros([nx, ny, nz], dtype=np.float64)
        if 'gvel' in var:
            idata['gvel'] = np.zeros([nx, ny, nz, 3], dtype=np.float64)

        idata['cellID'] = np.zeros(npoint, dtype=np.int)

        for ind in range(npoint):

            ix = int(np.floor(ind / ny / nz))
            iy = int(np.floor((ind - ix * ny * nz) / nz))
            iz = int(ind - ix * ny * nz - iy * nz)

            cellID = data.grid.getContainerLeafID((x[ix], y[iy], z[iz]))

            if cellID is not None:
                idata['cellID'][ind] = cellID
                if 'ddens' in var:
                    if cellID >= 0:
                        idata['rhodust'][ix, iy, iz, :] = data.rhodust[data.grid.leafID[cellID]]
                    else:
                        idata['rhodust'][ix, iy, iz, :] = 0.0
                if 'dtemp' in var:
                    if cellID >= 0:
                        idata['dusttemp'][ix, iy, iz, :] = data.dusttemp[data.grid.leafID[cellID]]
                    else:
                        idata['dusttemp'][ix, iy, iz, :] = 0.0
                if 'gdens' in var:
                    if cellID >= 0:
                        idata['rhogas'][ix, iy, iz] = data.ndens_mol[data.grid.leafID[cellID]]
                    else:
                        idata['rhogas'][ix, iy, iz] = 0.0
                if 'ndens' in var:
                    if cellID >= 0:
                        idata['ndens_mol'][ix, iy, iz] = data.ndens_mol[data.grid.leafID[cellID]]
                    else:
                        idata['ndens_mol'][ix, iy, iz] = 0.0
                if 'gtemp' in var:
                    if cellID >= 0:
                        idata['gastemp'][ix, iy, iz] = data.gastemp[data.grid.leafID[cellID]]
                    else:
                        idata['gastemp'][ix, iy, iz] = 0.0
                if 'vturb' in var:
                    if cellID >= 0:
                        idata['vturb'][ix, iy, iz] = data.vturb[data.grid.leafID[cellID]]
                    else:
                        idata['vturb'][ix, iy, iz] = 0.0
                if 'gvel' in var:
                    if cellID >= 0:
                        idata['gasvel'][ix, iy, iz, :] = data.gasvel[data.grid.leafID[cellID], :]
                    else:
                        idata['gasvel'][ix, iy, iz, :] = 0.0

    return idata


def plotSlice2D(data=None, var='ddens', plane='xy', crd3=0.0, icrd3=None, ispec=-1, xlim=(), ylim=(), log=False,
                linunit='cm', angunit='rad',
                nx=100, ny=100, showgrid=False, gridcolor='k', gridalpha=1.0, nproc=1,
                contours=False, clev=None, clmin=None, clmax=None, ncl=None, cllog=False, clcol='k', clcmap=None,
                cllabel=False, cllabel_fontsize=10, cllabel_fmt="%.1f", clalpha=1.0, ax=None, lattitude=True,
                **kwargs):
    """
    Function to plot an axis-aligned 2D slice of the variables in the model. Any additional keyword
    argument above the listed ones will be passed on to matplotlib.pylab.pcolormesh(). For an octree grid the variables 
    are interpolated onto a regular grid using nearest neighbour interpolation before plotting. 
    The size and resolution of the regular image grid can be set at input. 


    Parameters
    ----------

    data            : radmc3dData
                      Instance of radmc3dData containing the field variable to be displayed

    var             : {'ddens', 'dtemp', 'gdens', 'ndens', 'gtemp', 'vturb', 'vx', 'vy', 'vz', 'taux', 'tauy'}
                      Variable to be displayed

    plane           : {'xy', 'xz', 'yz', 'yx, 'zx', 'yz'}
                      Plane to be displayed           

    crd3            : float
                      Coordinate of the third dimension (i.e. when plotting a slice in the x-y plane, crd3 is the 
                      z-coordinate)

    icrd3           : int
                      Index of the third coordinate in the grid (only for regular grid!)

    ispec           : int
                      Dust species index. If negative dust densities will be summed up and the total cumulative density 
                      will be displayed

    xlim            : tuple
                      Coordinate boundaries in the first dimension of the plot (also the coordinate boundary of the 
                      regular grid data on
                      AMR grids are interpolated to)

    ylim            : tuple
                      Coordinate boundaries in the second dimension of the plot (also the coordinate boundary of the 
                      regular grid data on
                      AMR grids are interpolated to)

    log             : bool
                      If True the contour/image will be displayed on a logarithmic stretch

    linunit         : {'cm', 'au', 'pc', 'rs'}
                      Unit selection for linear image coordinate axes.

    nx              : int
                      Number of horizontal pixels in the interpolated image if the data is defined in an Octree

    ny              : int
                      Number of vertical pixels in the interpolated image if the data is defined in an Octree

    showgrid        : bool
                      If True the spatial grid will be overlayed 

    gridcolor       : str
                      Color of the spatial grid overlay

    gridalpha       : float
                      Opacity of the lines in the spatial grid overlay (0.0 - fully transparent, 1.0 - fully opaque)

    angunit         : {'rad', 'deg'}
                      Unit selection for angular image coordinate axes (only if spherical coordinate system is used).

    nproc           : int
                      Number of parallel processes to be used for interpolation. 

    contours        : bool
                      If True contour lines are plotted, if False a colorscale plot will be created

    clev            : ndarray  
                      A numpy ndarray containing the levels to be displayed with contour lines. If clev is set
                      then clmin, clmax and ncl are omitted

    clmin           : float
                      Min. contour level (for setting auto-contours between clmin and clmax at ncl values)

    clmax           : float
                      Max. contour level (for setting auto-contours between clmin and clmax at ncl values)

    ncl             : float
                      Number of contour levels (for setting auto-contours between clmin and clmax at ncl values)

    cllog           : bool
                      If clmin, clmax and ncl are used to generate the contour levels, then if cllog is True
                      the contours will be log-scaled

    clcol           : str
                      Color-code for the contour lines for single color contours

    clcmap          : matplotlib colormap
                      Colormap used for the contour lines

    cllabel         : bool
                      If True the contour line values will be displayed, if False only the contour lines will be
                      displayed (default = False)

    cllabel_fontsize: int
                      Size of the font used to displaye the contour line values

    cllabel_fmt     : str
                      Format of the contour line labels (default "%.1f")

    clalpha         : float
                      Transparency of the contour lines (1.0 fully opaque, 0.0 fully transparent)

    lattitude       : bool
                      If the coordinate sytem used in RADMC-3D is spherical, then the 2nd coordiante is the 
                      co-lattitude. If lattitude is set to True then the 2nd coordinate in the RADMC-3D grid will be 
                      transformet to true lattitude (i.e. pi/2.-colattitude). If set to false the original co-lattitude
                      will be used. 

    ax              : matplotlib.axes.Axes
                      Matplotlib axis to plot to


    Keyword Arguments :
                      All other keyword arugments will be passed to pcolormesh()

    """

    octree = False
    #
    # Check the input consistency
    #
    if data is None:
        raise ValueError('Unkonwn data. Data to be plotted is not specified.')

    if data.grid is None:
        raise AttributeError('Missing grid information in data. Plots cannot be made without a spatial grid.')
    else:
        if isinstance(data.grid, radmc3dOctree):
            octree = True

    if not octree:
        if icrd3 is None:
            if crd3 is None:
                raise ValueError('Unknown coordinate for the third dimension (icrd3/crd3)')

    var = var.strip().lower()
    varFound = False
    if var == 'ddens':
        varFound = True
    if var == 'dtemp':
        varFound = True
    if var == 'gdens':
        varFound = True
    if var == 'ndens':
        varFound = True
    if var == 'gtemp':
        varFound = True
    if var == 'vx':
        varFound = True
    if var == 'vy':
        varFound = True
    if var == 'vz':
        varFound = True
    if var == 'vturb':
        varFound = True
    if var == 'taux':
        varFound = True
        if octree:
            raise RuntimeError('Optical depth calculation has not yet been implemented for octrees')
    if var == 'tauy':
        varFound = True
        if octree:
            raise RuntimeError('Optical depth calculation has not yet been implemented for octrees')

    if not varFound:
        raise ValueError('Unknown variable to be plotted : ', var, '\n Allowed variable names are : ddens, dtemp, '
                         + 'gdens, ndens, dtemp, vx, vy, vz, vturb, taux, '
                         + 'tauy')

    #
    # Get the units
    #
    if linunit.strip().lower() == 'cm':
        linunit_label = '[cm]'
        linunit_norm = 1.
    elif linunit.strip().lower() == 'au':
        linunit_label = '[au]'
        linunit_norm = 1. / nc.au
    elif linunit.strip().lower() == 'pc':
        linunit_label = '[pc]'
        linunit_norm = 1. / nc.pc
    elif linunit.strip().lower() == 'rs':
        linunit_label = r'[R$_\odot$]'
        linunit_norm = 1. / nc.rs
    else:
        raise ValueError('Unknown linunit ', linunit, '\nSupported units are : cm, au, pc, rs')

    if angunit.strip().lower() == 'rad':
        angunit_label = '[rad]'
        angunit_norm = 1.0
    elif angunit.strip().lower() == 'deg':
        angunit_label = '[deg]'
        angunit_norm = np.pi / 180.
    else:
        raise ValueError('Unknown angunit ', angunit, '\nSupported units are : deg, rad')

    #
    # Now check which plane to be plotted
    #
    swapDim = False
    xnorm = 1.0
    ynorm = 1.0
    znorm = 1.0
    iplane = -1

    if octree:
        plot_x = xlim[0] + (xlim[1] - xlim[0]) * np.arange(nx, dtype=float) / float(nx - 1)
        plot_y = ylim[0] + (ylim[1] - ylim[0]) * np.arange(ny, dtype=float) / float(ny - 1)
        plot_z = np.array([crd3])
    else:
        plot_x = None
        plot_y = None
        plot_z = None

    if 'x' in plane:
        # xy plane
        if data.grid.nx <= 1:
            msg = 'The x dimension is switched off or has only a single grid cell, thus a 2D slice plot cannot ' \
                  'be done in the '+plane+' plane.'
            raise ValueError(msg)
        if 'y' in plane:
            if data.grid.ny <= 1:
                msg = 'The y dimension is switched off or has only a single grid cell, thus a 2D slice plot cannot ' \
                      'be done in the '+plane+' plane.'
                raise ValueError(msg)

            iplane = 2
            if not octree:
                plot_x = np.copy(data.grid.x)
                plot_y = np.copy(data.grid.y)
                if lattitude:
                    plot_y = np.pi / 2.0 - data.grid.y

                if icrd3 is None:
                    icrd3 = np.abs(data.grid.z - crd3).argmin()
            if data.grid.crd_sys == 'car':
                xlabel = r'x ' + linunit_label
                ylabel = r'y ' + linunit_label
                xnorm = linunit_norm
                ynorm = linunit_norm
            else:
                xlabel = 'r ' + linunit_label

                if lattitude:
                    ylabel = '$\\pi/2-\\theta$ ' + angunit_label
                else:
                    ylabel = '$\\theta$ ' + angunit_label

                xnorm = linunit_norm
                ynorm = angunit_norm

            if plane == 'yx':
                swapDim = True
                if data.grid.crd_sys == 'car':
                    xlabel = r'y ' + linunit_label
                    ylabel = r'x ' + linunit_label
                    xnorm = linunit_norm
                    ynorm = linunit_norm
                else:
                    xlabel = '$\\theta$ ' + angunit_label
                    ylabel = 'r ' + linunit_label
                    xnorm = angunit_norm
                    ynorm = linunit_norm

            if octree:
                idata = interpolateOctree(data, x=plot_x / xnorm, y=plot_y / ynorm, z=plot_z, var=var, nproc=nproc)
        # xz plane
        elif 'z' in plane:
            if data.grid.nz <= 1:
                msg = 'The z dimension is switched off or has only a single grid cell, thus a 2D slice plot cannot ' \
                      'be done in the '+plane+' plane.'
                raise ValueError(msg)

            iplane = 1
            if not octree:
                plot_x = np.copy(data.grid.x)
                plot_y = np.copy(data.grid.z)
                if icrd3 is None:
                    icrd3 = np.abs(data.grid.y - crd3).argmin()
                    if lattitude:
                        icrd3 = np.abs((np.pi / 2. - data.grid.y) - crd3).argmin()

            if data.grid.crd_sys == 'car':
                xlabel = r'x ' + linunit_label
                ylabel = r'z ' + linunit_label
                xnorm = linunit_norm
                ynorm = linunit_norm
            else:
                xlabel = 'r ' + linunit_label
                ylabel = '$\phi$ ' + angunit_label
                xnorm = linunit_norm
                ynorm = angunit_norm
            if plane == 'zx':
                swapDim = True
                if data.grid.crd_sys == 'car':
                    xlabel = r'z ' + linunit_label
                    ylabel = r'x ' + linunit_label
                    xnorm = linunit_norm
                    ynorm = linunit_norm
                else:
                    xlabel = '$\phi$ ' + angunit_label
                    ylabel = 'r ' + linunit_label
                    xnorm = angunit_norm
                    ynorm = linunit_norm

            if octree:
                idata = interpolateOctree(data, x=plot_x / xnorm, y=plot_z, z=plot_y / ynorm, var=var, nproc=nproc)
    # yz plane
    else:
        if data.grid.ny <= 1:
            msg = 'The y dimension is switched off or has only a single grid cell, thus a 2D slice plot cannot ' \
                  'be done in the '+plane+' plane.'
            raise ValueError(msg)

        if data.grid.nz <= 1:
            msg = 'The z dimension is switched off or has only a single grid cell, thus a 2D slice plot cannot ' \
                  'be done in the '+plane+' plane.'
            raise ValueError(msg)
        iplane = 0
        if not octree:
            plot_x = np.copy(data.grid.y)
            plot_y = np.copy(data.grid.z)
            if lattitude:
                plot_x = (np.pi / 2.0 - data.grid.y)
            if icrd3 is None:
                icrd3 = np.abs(data.grid.x - crd3).argmin()

        if data.grid.crd_sys == 'car':
            xlabel = r'y ' + linunit_label
            ylabel = r'z ' + linunit_label
            xnorm = linunit_norm
            ynorm = linunit_norm

        else:
            if lattitude:
                xlabel = '$\\pi/2-\\theta$ ' + angunit_label
            else:
                xlabel = '$\\theta$ ' + angunit_label
            ylabel = '$\phi$ ' + angunit_label
            xnorm = angunit_norm
            ynorm = angunit_norm
        if plane == 'zy':
            swapDim = True
            if data.grid.crd_sys == 'car':
                xlabel = r'z ' + linunit_label
                ylabel = r'y ' + linunit_label
                xnorm = linunit_norm
                ynorm = linunit_norm
            else:
                xlabel = '$\phi$ ' + angunit_label
                if lattitude:
                    ylabel = '$\\pi/2-\\theta$ ' + angunit_label
                else:
                    ylabel = '$\\theta$ ' + angunit_label
                xnorm = angunit_norm
                ynorm = angunit_norm

        if octree:
            idata = interpolateOctree(data, x=plot_z, y=plot_x / xnorm, z=plot_y / ynorm, var=var, nproc=nproc)

    #
    # Get the variable to be plotted
    #
    if var == 'ddens':
        if isinstance(data.rhodust, int):
            raise ValueError('Dust density is not present in the passed radmc3dData instance')
        else:
            if octree:
                if ispec >= 0:
                    pdata = np.squeeze(idata['rhodust'][:, :, :, ispec])
                else:
                    pdata = np.squeeze(idata['rhodust'].sum(3))

            else:
                if iplane == 0:
                    if ispec >= 0:
                        pdata = data.rhodust[icrd3, :, :, ispec]
                    else:
                        pdata = data.rhodust[icrd3, :, :, :].sum(2)
                elif iplane == 1:
                    if ispec >= 0:
                        pdata = data.rhodust[:, icrd3, :, ispec]
                    else:
                        pdata = data.rhodust[:, icrd3, :, :].sum(2)
                elif iplane == 2:
                    if ispec >= 0:
                        pdata = data.rhodust[:, :, icrd3, ispec]
                    else:
                        pdata = data.rhodust[:, :, icrd3, :].sum(2)

            cblabel = r'$\rho_{\rm dust}$ [g/cm$^3$]'

    elif var == 'dtemp':
        if isinstance(data.dusttemp, int):
            raise ValueError('Dust temperature is not present in the passed radmc3dData instance')
        else:
            if octree:
                if ispec >= 0:
                    pdata = np.squeeze(idata['dusttemp'][:, :, :, ispec])
                else:
                    raise IndexError('Negative dust species index.')

            else:
                if iplane == 0:
                    if ispec >= 0:
                        pdata = data.dusttemp[icrd3, :, :, ispec]
                    else:
                        raise IndexError('Negative dust species index.')
                elif iplane == 1:
                    if ispec >= 0:
                        pdata = data.dusttemp[:, icrd3, :, ispec]
                    else:
                        raise IndexError('Negative dust species index.')
                elif iplane == 2:
                    if ispec >= 0:
                        pdata = data.dusttemp[:, :, icrd3, ispec]
                    else:
                        raise IndexError('Negative dust species index : ', ispec)
            cblabel = r'$T_{\rm dust}$ [K]'

    elif var == 'gdens':
        if isinstance(data.rhogas, int):
            raise ValueError('Gas density is not present in the passed radmc3dData instance')
        else:
            if octree:
                pdata = np.squeeze(idata['rhogas'])
            else:
                if iplane == 0:
                    pdata = data.rhogas[icrd3, :, :]
                elif iplane == 1:
                    pdata = data.rhogas[:, icrd3, :]
                elif iplane == 2:
                    pdata = data.rhogas[:, :, icrd3]
            cblabel = r'$\rho_{\rm gas}$ [g/cm$^3$]'

    elif var == 'ndens':
        if isinstance(data.ndens_mol, int):
            raise ValueError('Gas number density is not present in the passed radmc3dData instance')
        else:
            if octree:
                pdata = np.squeeze(idata['ndens_mol'])
            else:
                if iplane == 0:
                    pdata = data.ndens_mol[icrd3, :, :]
                elif iplane == 1:
                    pdata = data.ndens_mol[:, icrd3, :]
                elif iplane == 2:
                    pdata = data.ndens_mol[:, :, icrd3]
            cblabel = r'$n_{\rm gas}$ [molecule/cm$^3$]'

    elif var == 'gtemp':
        if isinstance(data.gastemp, int):
            raise ValueError('Gas temperature is not present in the passed radmc3dData instance')
        else:
            if octree:
                pdata = np.squeeze(idata['gastemp'])
            else:
                if iplane == 0:
                    pdata = data.gastemp[icrd3, :, :]
                elif iplane == 1:
                    pdata = data.gastemp[:, icrd3, :]
                elif iplane == 2:
                    pdata = data.gastemp[:, :, icrd3]
            cblabel = r'$T_{\rm gas}$ [K]'

    elif var == 'vx':
        if isinstance(data.gasvel, int):
            raise ValueError('Gas velocity is not present in the passed radmc3dData instance')
        else:
            if octree:
                pdata = np.squeeze(idata['gasvel'][:, :, :, 0])
            else:
                if iplane == 0:
                    pdata = data.gasvel[icrd3, :, :, 0]
                elif iplane == 1:
                    pdata = data.gasvel[:, icrd3, :, 0]
                elif iplane == 2:
                    pdata = data.gasvel[:, :, icrd3, 0]
            cblabel = r'$v_{\rm x}$ [cm/s]'

    elif var == 'vy':
        if isinstance(data.gasvel, int):
            raise ValueError('Gas velocity is not present in the passed radmc3dData instance')
        else:
            if octree:
                pdata = np.squeeze(idata['gasvel'][:, :, :, 1])
            else:
                if iplane == 0:
                    pdata = data.gasvel[icrd3, :, :, 1]
                elif iplane == 1:
                    pdata = data.gasvel[:, icrd3, :, 1]
                elif iplane == 2:
                    pdata = data.gasvel[:, :, icrd3, 1]
            cblabel = r'$v_{\rm y}$ [cm/s]'
    elif var == 'vz':
        if isinstance(data.gasvel, int):
            raise ValueError('Gas velocity is not present in the passed radmc3dData instance')
        else:
            if octree:
                pdata = np.squeeze(idata['gasvel'][:, :, :, 2])
            else:
                if iplane == 0:
                    pdata = data.gasvel[icrd3, :, :, 2]
                elif iplane == 1:
                    pdata = data.gasvel[:, icrd3, :, 2]
                elif iplane == 2:
                    pdata = data.gasvel[:, :, icrd3, 2]
            cblabel = r'$v_{\rm z}$ [cm/s]'

    elif var == 'vturb':
        if isinstance(data.vturb, int):
            raise ValueError('Microturbulent velocity is not present in the passed radmc3dData instance')
        else:
            if octree:
                pdata = np.squeeze(idata['vturb'])
            else:
                if iplane == 0:
                    pdata = data.vturb[icrd3, :, :]
                elif iplane == 1:
                    pdata = data.vturb[:, icrd3, :]
                elif iplane == 2:
                    pdata = data.vturb[:, :, icrd3]
            cblabel = r'$v_{\rm turb}$ [cm/s]'

    if var == 'taux':
        if isinstance(data.taux, int):
            raise ValueError('Optical depth is not present in the passed radmc3dData instance')
        else:
            if octree:
                raise RuntimeError('Optical depth calculation has not yet been implemented for octrees.')
            else:
                if (data.taux.shape[0]==0):
                    raise ValueError('Optical depth has not been calculated yet. Run radmc3dData.getTau(wav=X) to '
                                     'calculate the optical depth before calling plotSlice2D')
                if iplane == 0:
                    pdata = data.taux[icrd3, :, :]
                elif iplane == 1:
                    pdata = data.taux[:, icrd3, :]
                elif iplane == 2:
                    pdata = data.taux[:, :, icrd3]
            cblabel = r'$\tau_{\rm r}$'

    if var == 'tauy':
        if isinstance(data.tauy, int):
            raise ValueError('Optical depth is not present in the passed radmc3dData instance')
        else:
            if octree:
                raise RuntimeError('Optical depth calculation has not yet been implemented for octrees.')
            else:
                if iplane == 0:
                    pdata = data.tauy[icrd3, :, :]
                elif iplane == 1:
                    pdata = data.tauy[:, icrd3, :]
                elif iplane == 2:
                    pdata = data.tauy[:, :, icrd3]
            cblabel = r'$\tau_{\rm \theta}$'

    #
    # Apparently there is some inconsistency in the dimensionality of the gas arrays. I.e. when data is read from file
    #  there is always four dimension with the last dimension, for dust used for the dust species index, set to one.
    #  I guess it was meant to be if we have multiple gas species. However, if data is created by the model functions
    #  directly the gas density is usually returned as a three dimensional array. So the quick and dirty fix is to
    #  check for the dimensionality here and decrease it by one stripping the last dimension either by a specified
    #  'ispec' index
    #

    if len(pdata.shape) == 3:
        if pdata.shape[2] == 1:
            pdata = pdata[:, :, 0]
        else:
            raise ValueError('The plotted data has four dimension (3 spatial + 1 species) but the species index is '
                             + ' not set. Specify ispec keyword which of the dimensinos should be plotted')

    #
    # If the dimensions are flipped in the plotted plane (i.e. yx, zy, or xz) flip the array dimensions
    #
    if swapDim:
        pdata = pdata.swapaxes(0, 1)
        dum = np.array(plot_x)
        plot_y = np.array(plot_x)
        plot_x = dum

    if not octree:
        plot_x *= xnorm
        plot_y *= ynorm
    #
    # Get the data min/max values
    #

    if ax is None:
        ax = plb.gca()

    #
    # Now do the colorscale plotting but only if contours is set to False
    #
    if not contours:
        if log:
            if 'vmin' not in kwargs:
                vmin = pdata.min()
            else:
                vmin = kwargs['vmin']

            if 'vmax' not in kwargs:
                vmax = pdata.max()
            else:
                vmax = kwargs['vmax']

            if pdata.min() <= 0:
                pdata = pdata.clip(1e-90)
            p = ax.pcolormesh(plot_x, plot_y, pdata.T, norm=LogNorm(vmin, vmax), **kwargs)
            # p = ax.imshow(pdata.T, origin='lower', extent=(plot_x[0]-dx*0.5, plot_x[-1]+dx*0.5,
            # plot_y[0]-dy*0.5, plot_y[-1]+dy*0.5),
            # norm=LogNorm(vmin, vmax), interpolation='nearest', aspect='auto', **kwargs)
            # Enable rasterization to enable easy save to file
            p.set_rasterized(True)
        else:
            p = ax.pcolormesh(plot_x, plot_y, pdata.T, **kwargs)
            # p = ax.imshow(pdata.T, origin='lower', extent=(plot_x[0]-dx*0.5, plot_x[-1]+dx*0.5,
            # plot_y[0]-dy*0.5, plot_y[-1]+dy*0.5),
            # vmin=vmin, vmax=vmax, interpolation='nearest', aspect='auto', **kwargs)
            # Enable rasterization to enable easy save to file
            p.set_rasterized(True)

        #
        # Generate the colorbar
        #
        cb = plb.colorbar(p)
        cb.set_label(cblabel)

    #
    # Show the grid (currently only for otctree AMR)
    #
    if showgrid:
        if octree:
            ind = 0
            plottedInd = np.zeros(idata['cellID'].shape[0], dtype=np.int) - 1

            if plane.strip().lower() == 'xy':
                for i in idata['cellID']:
                    if i not in plottedInd:
                        ind += 1
                        bottomleft = ((data.grid.x[i] - data.grid.dx[i]) * xnorm,
                                      (data.grid.y[i] - data.grid.dy[i]) * ynorm)
                        ax.add_patch(patches.Rectangle(bottomleft, data.grid.dx[i] * 2 * xnorm,
                                                       data.grid.dy[i] * 2 * ynorm, fill=False, edgecolor=gridcolor,
                                                       alpha=gridalpha))
                        plottedInd[ind] = i

            elif plane.strip().lower() == 'yx':

                for i in idata['cellID']:
                    if i not in plottedInd:
                        ind += 1
                        bottomleft = ((data.grid.y[i] - data.grid.dy[i]) * xnorm,
                                      (data.grid.x[i] - data.grid.dx[i]) * ynorm)
                        ax.add_patch(patches.Rectangle(bottomleft, data.grid.dy[i] * 2 * xnorm,
                                                       data.grid.dx[i] * 2 * ynorm, fill=False, edgecolor=gridcolor,
                                                       alpha=gridalpha))
                        plottedInd[ind] = i

            elif plane.strip().lower() == 'xz':

                for i in idata['cellID']:
                    if i not in plottedInd:
                        ind += 1
                        bottomleft = ((data.grid.x[i] - data.grid.dx[i]) * xnorm,
                                      (data.grid.z[i] - data.grid.dz[i]) * ynorm)
                        ax.add_patch(patches.Rectangle(bottomleft, data.grid.dx[i] * 2 * xnorm,
                                                       data.grid.dz[i] * 2 * ynorm, fill=False, edgecolor=gridcolor,
                                                       alpha=gridalpha))
                        plottedInd[ind] = i

            elif plane.strip().lower() == 'zx':

                for i in idata['cellID']:
                    if i not in plottedInd:
                        ind += 1
                        bottomleft = ((data.grid.z[i] - data.grid.dz[i]) * xnorm,
                                      (data.grid.x[i] - data.grid.dx[i]) * ynorm)
                        ax.add_patch(patches.Rectangle(bottomleft, data.grid.dz[i] * 2 * xnorm,
                                                       data.grid.dx[i] * 2 * ynorm, fill=False, edgecolor=gridcolor,
                                                       alpha=gridalpha))
                        plottedInd[ind] = i

            elif plane.strip().lower() == 'yz':

                for i in idata['cellID']:
                    if i not in plottedInd:
                        ind += 1
                        bottomleft = ((data.grid.y[i] - data.grid.y[i]) * xnorm,
                                      (data.grid.z[i] - data.grid.dz[i]) * ynorm)
                        ax.add_patch(patches.Rectangle(bottomleft, data.grid.dy[i] * 2 * xnorm,
                                                       data.grid.dz[i] * 2 * ynorm, fill=False, edgecolor=gridcolor,
                                                       alpha=gridalpha))
                        plottedInd[ind] = i

            elif plane.strip().lower() == 'zy':

                for i in idata['cellID']:
                    if i not in plottedInd:
                        ind += 1
                        bottomleft = ((data.grid.z[i] - data.grid.dz[i]) * xnorm,
                                      (data.grid.y[i] - data.grid.dy[i]) * ynorm)
                        ax.add_patch(patches.Rectangle(bottomleft, data.grid.dz[i] * 2 * xnorm,
                                                       data.grid.dy[i] * 2 * ynorm, fill=False, edgecolor=gridcolor,
                                                       alpha=gridalpha))
                        plottedInd[ind] = i

        else:

            if plane.strip().lower() == 'xy':
                px = data.grid.xi * xnorm
                if (data.grid.crd_sys == 'sph') & lattitude:
                    py = (np.pi / 2. - data.grid.yi) * ynorm
                else:
                    py = data.grid.yi * ynorm

            elif plane.strip().lower() == 'yx':
                py = data.grid.xi * xnorm
                if (data.grid.crd_sys == 'sph') & lattitude:
                    px = (np.pi / 2. - data.grid.yi) * ynorm
                else:
                    px = data.grid.yi * ynorm

            elif plane.strip().lower() == 'xz':
                px = data.grid.xi * xnorm
                py = data.grid.zi * znorm

            elif plane.strip().lower() == 'zx':
                py = data.grid.xi * xnorm
                px = data.grid.zi * znorm

            elif plane.strip().lower() == 'yz':
                if (data.grid.crd_sys == 'sph') & lattitude:
                    px = (np.pi / 2. - data.grid.yi) * ynorm
                else:
                    px = data.grid.yi * ynorm
                py = data.grid.zi * znorm

            elif plane.strip().lower() == 'zy':
                if (data.grid.crd_sys == 'sph') & lattitude:
                    py = (np.pi / 2. - data.grid.yi) * ynorm
                else:
                    py = data.grid.yi * ynorm
                px = data.grid.zi * znorm

            for ix in range(data.grid.nxi):
                ax.add_line(ml.Line2D((px[ix], px[ix]), (py[0], py[-1]), color=gridcolor, alpha=gridalpha))
            for iy in range(data.grid.nyi):
                ax.add_line(ml.Line2D((px[0], px[-1]), (py[iy], py[iy]), color=gridcolor, alpha=gridalpha))

    if contours is True:
        if clmin is None:
            clmin = pdata.min()
        if clmax is None:
            clmax = pdata.max()
        if ncl is None:
            ncl = 20

        # Generate the contour levels
        if clev is None:
            if cllog is True:
                clev = clmin * (clmax / clmin) ** (np.arange(ncl, dtype=float) / float(ncl - 1))
            else:
                clev = clmin + (clmax - clmin) * (np.arange(ncl, dtype=float) / float(ncl - 1))
        if (clcol == 'none') | (clcol is None):
            if 'cmap' in kwargs:
                c = ax.contour(plot_x, plot_y, pdata, clev, kwargs['cmap'], alpha=clalpha)
            else:
                c = ax.contour(plot_x, plot_y, pdata, clev, alpha=clalpha)
        else:
            c = ax.contour(plot_x, plot_y, pdata.T, clev, colors=clcol, alpha=clalpha)

        if cllabel:
            plb.clabel(c, inline=1, fontsize=cllabel_fontsize, fmt=cllabel_fmt)
    #
    # Set the axis limits
    #
    if not octree:
        if len(xlim) == 0:
            xlim = (plot_x[0], plot_x[-1])
        if len(ylim) == 0:
            ylim = (plot_y[0], plot_y[-1])

    plb.xlim(xlim[0], xlim[1])
    plb.ylim(ylim[0], ylim[1])
    #
    # Set the axis labels
    #
    plb.xlabel(xlabel)
    plb.ylabel(ylabel)
