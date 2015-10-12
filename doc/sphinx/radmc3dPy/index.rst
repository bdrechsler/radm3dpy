.. radmc3dPy documentation master file, created by
   sphinx-quickstart on Fri Jun  6 11:10:00 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to radmc3dPy's documentation!
=====================================

radmc3dPy is a python package/user-interface to the 3D radiative transfer code `RADMC-3D <http://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/>`_.
As of now radmc3dPy provides functionality to

* set up continuum and gas models
* handle a library of models (which are easy to extend user-defined models)
* read and write data files in ASCII and C-style binary format 
* provide interface to the Mie-scattering code shipped with RADMC-3D to calculate dust opacities
* set up spherical or cartesian spatial grids with "separable grid refinement"  (No AMR yet)
* write legacy VTK format (spherical grids only)
* set up discrete and continuous starlike radiation sources
* read and display images (single or multiple frequency, i.e. continuum or line channel maps, total intensity or full 4D polarization, ASCII or C-style binary)
* write fits files with CASA-compatible fits headers
* do simple image manipulations
    * calculation of moment maps for line channel maps
    * convolve image with a 2D Gaussian or an Airy-type psf
    * add coronographic masks
* calculate interferometric visibilities of images for given projected baselines and position angles


Copyright
---------

radmc3dPy was developed by Attila Juhasz at the Leiden Observatory between 2011 and 2014 and from 2014 at the Institute of Astronomy in Cambridge. 

Disclaimer
----------
I/We reject all responsibility for the use of this package. The package is provided as-is, and we are not responsible for any damage to hardware or software, 
nor for incorrect results that may result from the software. The package is still in development and it may contain bugs. In case of any bug, please, contact
the author of the package (Attila Juhasz <juhasz@ast.cam.ac.uk>).  

Requirements
------------

radmc3dPy requires the following packages to be installed:

* `Numpy <http://www.numpy.org>`_ v1.6.2 or later
* `SciPy <http://www.scipy.org>`_ 0.11.0 or later
* `matplotlib <http://matplotlib.org>`_ 1.2.0 or later
* `AstroPy <http://www.astropy.org>`_ v0.3 or later or `Pyfits <http://www.stsci.edu/institute/software_hardware/pyfits>`_ v3.0.7 

Download
--------

The latest version of radmc3dPy (0.28) can be downloaded from :download:`here <radmc3dPy-0.28.tar.gz>`


Contents
--------

.. toctree::
   :maxdepth: 2

   tutorial
   models
   parfile

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

