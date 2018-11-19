.. _downloads:

*********
Downloads
*********

You can find below links to download the various versions of radmc3dPy

**Latest stable version**

* :download:`v0.30.2 <radmc3dPy-0.30.2.tar.gz>`
    
    *Important bugfixes*

        * When scalarfields (e.g. dust density) was read from a data file in ascii format, the file header was read wrong, 
          which may resulted in a RuntimeException due to inconsistency in data size between the mesh and the scalar
          data.
        * When scalarfields defined on an octree AMR mesh were read from a file in ascii format the header had an indexing
          error, which may resulted in a RuntimeException due to inconsistency in data size between the mesh and the scalar
          data.
        * Image convolution always returned positive pixel values. While images calculated in total intensity (i.e. Stokes I),
          should always be positive, this is not true for the other Stokes parameters. This might have resulted in wrong
          intensity values when images with full stokes parameters were convolved. 
        * Axis labels in the plotSpectrum() function were drawn wrong.

**Previous versions**

    * :download:`v0.30.1 <radmc3dPy-0.30.1.tar.gz>`
        
        *New features*

            * New Fortran implementation of the Mie-code with Python interface (f2Py) replacing the pure Python implementation. 

    * :download:`v0.30 <radmc3dPy-0.30.tar.gz>`
        
        *New features*

            * Python 3 compatibility
            * Pure python implementation of a Mie-scattering code (Kees)
            * New package (sphtool) for re-gridding 3D Smoothed Particle Hydrodynamics (SPH) simulations
            * Modular model setup and a new model class (radmc3dModel)
            * Functions to plot dust opacities and scattering matrix elements
            * Function to write dust opacities to file
            * Internal structural changes (splitting the analyze modules to smaller modules)



    * :download:`v0.29 <radmc3dPy-0.29.tar.gz>`
    
        *New features*

            * Octree AMR support (read and write grid and data fields, generate models with an octree grid)
            * Plot 2D contour maps of data fields along axis aligned slices (show the cell boundaries for AMR models) 
            * Plot the polarisation direction for full stokes images (thanks Kees!)
            * Read molecular data files (thanks Kees!)
            * Plot spectrum (thanks Kees!)
            * Functions to calculate dust and gas mass in a model
            * Updated tutorials for models using an octree AMR grid

    * :download:`v0.28.2 <radmc3dPy-0.28.2.tar.gz>`


