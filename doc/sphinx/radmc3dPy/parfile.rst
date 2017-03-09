.. _parfile:

**********************
Parameter file
**********************

.. _parameters-parfile-syntax:

Parameter file syntax
=====================

The parameter file is split up to blocks of parameters. Each parameter block begins with a line ::

    # Block: BLOCKNAME

Anything in the parameter file followed by a ``#`` sign is taken to be a comment, except for the block 
name definitions. Block name definitions begin with a ``#`` followed by the word ``Block``, a colon and then the
actual name of the block. 

Block names are followed by the lines containing the parameters and their values. Each parameter definition contains
three parts, name of the parameter, value and some explanatory comment::

    PARAMETER_NAME = PARAMETER_VALUE  # PARAMETER_DESCRIPTION

Long expressions and comments can be broken into multiple lines with a backslash (\\)
at the end of of the broken line. 

In reality a parameter block for the radiation sources can looks like this::

    # -----------------------------------------------------------------------------
    # Block: Radiation sources
    # -----------------------------------------------------------------------------
    incl_cont_stellarsrc      = False  # # Switches on (True) or off (False) continuous stellar sources )
    incl_disc_stellarsrc      = True  # # Switches on (True) or off (False) discrete stellar sources)
    mstar                     = [1.0*ms]  # # Mass of the star(s)
    pstar                     = [0.0, 0.0, 0.0]  # # Position of the star(s) (cartesian coordinates)
    rstar                     = [2.0*rs]  # # Radius of the star(s)
    tstar                     = [4000.0]  # # Effective temperature of the star(s) [K]

List of blocknames

    * ``Radiation sources``
    * ``Grid parameters``
    * ``Dust opacity``
    * ``Gas line RT``
    * ``Code parameters``
    * ``Model``

These block-names should not be modified as the reader function will look for these strings. 

.. _parameters-radiation-sources:

Radiation sources
=================

    * |  **incl_disc_stellarsrc :** list 
      |  Switches on (True) or off (False) discrete stellar radiation sources.
    * |  **mstar :** list 
      |  Mass of the star. Each element of the list contains the mass of an individual star as a float.        
    * |  **pstar :** list     
      |  Coordinates of the star. Each element of the list contains a three element vector containing the 3D cartesian coordinates of each individual star.
    * |  **rstar :** list    
      |  Stellar radius. Each element of the list contains the radius of an individual star as a float.        
    * |  **tstar :** list     
      |  Effective temperature. Each element of the list contains the effective temperature of an individual star as a float.        
    * |  **incl_cont_stellarsrc :** list 
      |  Switches on (True) or off (False) continuous stellar radiation sources
      |  NOTE, the model should have the appropriate functions (:ref:`getStellarsrcDensity() <models-model-functions-stellarsrc>`, :ref:`getStellarsrcTemplates() <models-model-functions-stellarsrc>`)


.. _parameters-grid:

Grid
====

    
    * |  **crd_sys :** {'sph', 'car'}      
      |  Coordinate system type
    * |  **nw :** list     
      |  Number of wavelength points in the wavelength grid, ``nw[i]`` sets the number of grid points in the ``[wbound[i], wbound[i+1])`` interval (or on the ``[wbound[-2], wbound[-1]]`` for the last interval)
    * |  **nx :** list     
      |  Number of grid cells in the first spatial coordinate, ``nx[i]`` sets the number of grid points in the ``[xbound[i], xbound[i+1])`` interval (or on the ``[xbound[-2], xbound[-1]]`` for the last interval)
    * |  **ny :** list     
      |  Number of grid cells in the second spatial coordinate, ``ny[i]`` sets the number of grid points in the ``[ybound[i], ybound[i+1])`` interval (or on the ``[ybound[-2], ybound[-1]]`` for the last interval)
    * |  **nz :** list     
      |  Number of grid cells in the third spatial coordinate, ``nz[i]`` sets the number of grid points in the ``[zbound[i], zbound[i+1])`` interval (or on the ``[zbound[-2], zbound[-1]]`` for the last interval)
    * |  **wbound :** list     
      |  Boundaries of the wavelength grid
    * |  **xbound :** list     
      |  Boundaries of the spatial grid in the first coordinate, ``nx[i]`` sets the number of grid points in the ``[xbound[i], xbound[i+1])`` interval (or on the ``[xbound[-2], xbound[-1]]`` for the last interval)
    * |  **xres_nlev :** float    
      |  Radial grid refinement in spherical coordinate system
    * |  **xres_nspan :** float    
      |  Radial grid refinement in spherical coordinate system
    * |  **xres_nstep :** int      
      |  Radial grid refinement in spherical coordinate system
    * |  **ybound :** list     
      |  Boundaries of the spatial grid in the second coordinate, ``ny[i]`` sets the number of grid points in the ``[ybound[i], ybound[i+1])`` interval (or on the ``[ybound[-2], ybound[-1]]`` for the last interval)
    * |  **zbound :** list     
      |  Boundaries of the spatial grid in the third coordinate, ``nz[i]`` sets the number of grid points in the ``[zbound[i], zbound[i+1])`` interval (or on the ``[zbound[-2], zbound[-1]]`` for the last interval)

.. _parameters-grid-separable-refinement:

Separable mesh refinement
-------------------------
    Spatial and wavelength grid definitions allow 'separable mesh refinement', i.e. refinement of the wavelength or the spatial mesh along individual
    axes. Let us take now the wavelength grid for an example. If we wish to cover the :math:`10^{-2}-10^4\mu{\rm m}` interval with 200 grid points 
    we should set ``wbound = [0.01, 1e4]`` and ``nw = [200]``. This results in a logarithmic wavelength grid between ``wbound[0]`` and ``wbound[1]``
    containing ``nw[0]`` grid points. This wavelenght grid might be fine enough to sample the radiation field of the sources and the thermal emission
    of the dust in the model, but too coarse to study e.g. the shape of the silicate features in the mid-infrared. If we are interested in the 
    silicate feature around :math:`10\mu{\rm m}` only, we can refine this region in the wavelength grid by setting ``wbound = [0.01, 7.5, 13.5, 1e4]`` and
    ``nw = [50,100,50]``.  This grid setup will result in 50, 100, 50 grid points in the :math:`[0.01\mu{\rm m},7.5\mu{\rm m})`,
    :math:`[7.5\mu{\rm m},13.5\mu{\rm m})` , :math:`[13.5\mu{\rm m},10^4\mu{\rm m}]` intervals, respectively.
    
    **Note**, the number of grid points are defined always on a right-open interval, except in the last, rightmost interval, where the interval is closed.


.. _parameters-grid-xrefinement:

Grid refinement at the inner boundary
-------------------------------------

    Even for logarithmic radial grids the innermost parts of the model can still be optically thick for centrally concentrated density distributions.
    With the use of the ``xres_nlev``, ``xres_nspan`` and ``xres_nstep`` parameters we can introduce additional grid refinement in the radial grid of
    a spherical coordinate system. The grid refinement is done in the following way. First a logarithmic radial grid is set up on the basis of the
    ``nx`` and ``xbound`` parameters. Then the interval between the innermost cell interface and the ``xres_nspan`` th 
    cell interface (i.e. ``xres_nspan``-1 grid cell) will be taken and split into ``xres_nlev`` grid cells. Then the innermost grid cell will be
    split into ``xres_nlev`` 'new' grid cells, then again the innermost, refined grid cell will be taken and split into ``xres_nlev`` cells.
    The splitting of the innermost cell will be done ``xres_nstep`` times. 

.. _parameters-dust-opacity:

Dust opacity
============

    * |  **dustkappa_ext :** str
      |  File name tag in the dust opacity file. Dust opacity files should have names like e.g., ``dustkappa_EXT.inp``, where the ``dustkappa_ext`` parameter should contain the 'EXT' tags from the file name (e.g. for ``dustkappa_ext = 'silicate'`` the dust opacity file should be ``dustkappa_silicate.inp``.
    * |  **gdens :** float
      |  Bulk density of the material
    * |  **gsdist_powex :** float
      |  Grain size distribution power exponent 
    * |  **gsmax :** float
      |  Maximum grain size in the distribution 
    * |  **gsmin :** float
      |  Minimum grain size in the distribution 
    * |  **lnk_fname :** list
      |  File name list (including full path) containing optical constants (NOTE, the file should contain three columns: wavelength [micron], n, k)
    * |  **mixabun :** list
      |  If multiple species specified their mass absorption coefficients can be mixed according to the mixing ratios (mass fractions) in mixabun. 
    * |  **ngs :** float
      |  Number of grain sizes in the grain size distribution

.. _parameters-gas-lines:

Gas lines
=========
    
    * |  **gasspec_colpart_abun :** float
      |  Abundance of the collisional partner
    * |  **gasspec_colpart_name :** float
      |  Name of the collisional partner 
    * |  **gasspec_mol_abun :** float
      |  Molecular abundance 
    * |  **gasspec_mol_dbase_type :** {'leiden', 'linelist'}
      |  Database type of the molecular data (see the `RADMC-3D manual <http://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/>`_ for the definitions of various formats).
    * |  **gasspec_mol_name :** str
      |  Name of the molecular species whose lines should be calculated

.. _parameters-code:

Code
====
    * |  **istar_sphere :** int
      |  If 0 discrete stars are taken to be point-like, if 1 the finite extent of the star is taken into account
    * |  **itemdecoup :** int
      |  Allows (0) or prevents (1) the decoupling of the temperature of different dust species 
    * |  **lines_mode :** int
      |  Line mode (for the definitions of the individual line modes see the `RADMC-3D manual <http://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/>`_):
            
            * 1 - LTE
            * 2 - User-defined populations I
            * 3 - LVG populations
            * 4 - Optically thin NLTE level populations method
            * 5 - User-defined populations II 

    * |  **nphot :** int
      |  Number of photons in the thermal Monte Carlo simulations
    * |  **nphot_scat :** int
      |  Number of photons used for the scattering Monte Carlo simulations when images are calculated 
    * |  **nphot_spec :** int
      |  Number of photons used for the scattering Monte Carlo simultaions when SEDs/spectra are calculated 
    * |  **rto_style :** int
      |  Output format: 1 - Formatted ASCII, 3 - C-style binary
    * |  **scattering_mode_max :** int
      |  Scattering mode :

            * 0 - Scattering is switched off 
            * 1 - Isotropic scattering
            * 2 - Anysotropic scattering with Henyey-Greenstein phase function 
            * 3 - Anysotropic scattering with tabulated phase function
            * 4 - Anysotropic scattering with polarization but the full scattering matrix is only used for the last scattering
            * 5 - Anysotropic scattering with scattering matrix, full treatment

    * |  **tgas_eq_tdust :** int
      |  Dust temperature is taken to be the gas kinetic temperature 
    * |  **modified_random_walk :** int
      |  Switches on (1) and off (0) modified random walk 
   
