.. _changes:

*******
Changes
*******

**Planned for v0.31**

* Nested mesh support
* Data inspector GUI
* Line RT analysis tools 

**v0.30**

* Python3 compatibility
* Integration of Kees' Mie-scattering code (:mod:`radmc3dPy.miescat`, meth:`radmc3dPy.dustopac.radmc3dDustOpac.computeDustOpacMie`) for dust opacity calculation
* New package (`sphtool <http://www.ast.cam.ac.uk/~juhasz/sphtoolDoc/index.html>`_) for re-gridding 3D Smoothed Particle Hydrodynamics (SPH) simulations to regular or AMR grids, 
  including some basic visualisation of the SPH simulations
* Possibility for modular module setup (i.e. generation of grids, variables, radiation sources, etc independently of each other) using the new 
  class :class:`~radmc3dPy.setup.radmc3dModel`.
* Functions to plot dust opacities (:func:`radmc3dPy.analyze.plotDustOpac`), and scattering matrix elements (:func:`radmc3dPy.analyze.plotScatmat`)
* New function (:meth:`radmc3dPy.dustopac.radmc3dDustOpac.writeOpac`) to write dust opacities to file
* Internal structural changes (splitting the analyze module to smaller modules, for easier maintenance, but all previous functionality 
  of the analyze module is kept intact)
* Numpy's fromfile() function is now used for reading dust opacities in :class:`radmc3dPy.dustopac.radmc3dDustOpac`, which makes the file reading less sensitive to the presence of empty lines in the file dust opacity files

**v0.29**

* Octree AMR support (read and write grid and data fields, generate models with an octree grid, see Octree tutorial)
* Plot 2D contour maps of data fields along axis aligned slices showing also the cell boundaries for AMR models (:func:`radmc3dPy.analyze.plotSlice2D`)
* Plot the polarisation direction for full stokes images (:func:`radmc3dPy.image.plotImage` thanks Kees!)
* Read molecular data files (:func:`radmc3dPy.analyze.readMol` thanks Kees!)
* Plot spectrum (:func:`radmc3dPy.analyze.plotSpectrum` thanks Kees!)
* Functions to calculate dust and gas mass in a model (:meth:`radmc3dPy.data.radmc3dData.getDustMass`, :meth:`radmc3dPy.data.radmc3dData.getGasMass`)

* Updated tutorials for models using an octree AMR grid



