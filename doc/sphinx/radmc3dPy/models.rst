.. _models:

**********************
Model source files
**********************

.. _models-file-names-and-locations:

File names and location
=======================

Each model in radmc3dPy should be named as `modelname.py`. The file should be located either in the
`radmc3dPy/models` directory or in the current workging directory. When a model is set up the
:meth:`~radmc3dPy.setup.problemSetupDust` and :meth:`~radmc3dPy.setup.problemSetupGas` methods
first try to import the model from the current working directory. If the model can't be imported
from the current directory, then the radmc3dPy model library will be checked for the model.

.. _models-add-model:

Add a model to the library
==========================

To add a model to the radmc3dPy model library, we need write permissions to the directory where
radmc3dPy was installed. First we need to copy the model file to the `radmc3dPy/models` directory.
Then the :meth:`~radmc3dPy.models.updateModelList` function in the :class:`~radmc3dPy.models` class
needs to be called. It will update the list of available models in the package. Note, that for the
update to take effect we need to reaload the imported modules / package or we need to quit the 
current python session and restart it again. 

However, it is not necessary to add the model to the library to be able to use it. The setup functions
look for the model source file in the current working directory first and if such file is present it 
will use it otherwise they will look for the model in the library. 
**Note**, if the current working directory contains a file with the same name as one of
the models in the model library the file in the current working directory will be used for the setup!

 
.. _models-model-functions:

Model functions
===============

Each model source file should contain functions to provide information on the current model and provide
the parameters of the model with default values. Furthermore each model should contain functions 
to describe the dust and/or gas structure in the model and optionally to describe the distribution
and spectra of continuous star-like radiation sources. The names of these functions are listed below.
Apart from these functions any number and type of user-defined auxiliary functions can be implemented
in the model source files. 

.. _models-model-functions-init:

Model init functions
--------------------

Every model source file should contain the following functions:
    
    * **getModelDesc()** - Returns the brief description of the model.

    * **getDefaultParams()** - Returns the default parameters of the model. This function should return a list with default parameter values. Each element
      of the list should contain another at least three element list: 1) name of the parameter, 2) parameter value, 3) Explanatory comment, 4) Block name (to which this parameter should be placed in the master parameter file). All parameter values should be given as strings as they should be written to the master parameter file. 

.. _models-model-functions-dust:

Dust structure functions
------------------------

This function must be implemented for all dust continuum models to describe the dust density structure:

    * **getDustDensity(ppar=None, grid=None)** - Calculates the dust density. It takes two parameters :
        * **ppar** - Dictionary containing all parameters of the models with the parameter names as keys
        * **grid** - An instance of the :class:`~radmc3dPy.analyze.radmc3dGrid` class containing the spatial 
                    and wavelength/frequency grid of the model
        **Returns** an ndarray with [nx,ny,nz,ndust] dimensions

.. _models-model-functions-lines:

Gas structure functions
-----------------------
   
The following functions must be implemented for a gas line model:

    * **getGasDensity(ppar=None, grid=None)** - Calculates the gas density. It takes two parameters :
        * **ppar** - Dictionary containing all parameters of the models with the parameter names as keys
        * **grid** - An instance of the :class:`~radmc3dPy.analyze.radmc3dGrid` class containing the spatial 
                    and wavelength/frequency grid of the model
        **Returns** an ndarray with [nx,ny,nz] dimensions
 
    * **getGasAbundance(ppar=None, grid=None)** - Calculates the gas abundance. It takes two parameters :
        * **ppar** - Dictionary containing all parameters of the models with the parameter names as keys
        * **grid** - An instance of the :class:`~radmc3dPy.analyze.radmc3dGrid` class containing the spatial 
                    and wavelength/frequency grid of the model
        **Returns** an ndarray with [nx,ny,nz] dimensions

    * **getVTurb(ppar=None, grid=None)** - Calculates the microturbulent velocity field. It takes two parameters :
        * **ppar** - Dictionary containing all parameters of the models with the parameter names as keys
        * **grid** - An instance of the :class:`~radmc3dPy.analyze.radmc3dGrid` class containing the spatial 
                    and wavelength/frequency grid of the model
        **Returns** an ndarray with [nx,ny,nz] dimensions

    * **getVelocity(ppar=None, grid=None)** - Calculates the gas velocity. It takes two parameters :
        * **ppar** - Dictionary containing all parameters of the models with the parameter names as keys
        * **grid** - An instance of the :class:`~radmc3dPy.analyze.radmc3dGrid` class containing the spatial 
                    and wavelength/frequency grid of the model
        **Returns** an ndarray with [nx,ny,nz,3] dimensions, with [i,j,k,:] containing the three velocity component (vx,vy,vz) for cartesian and (vr,vtheta,vphi) for spherical coordinate system.
 

.. _models-model-functions-stellarsrc:

Radiation source functions
--------------------------

These functions should be implemented if a continuous star-like radiation source is added to the model. 
These functions describe the "density" distribution of the continuous starlike sources as well as their
spectra. 

    * **getStellarsrcDensity(ppar=None, grid=None)** - Calculates the stellar density for continuous star-like radiation sources 
        * **ppar** - Dictionary containing all parameters of the models with the parameter names as keys
        * **grid** - An instance of the :class:`~radmc3dPy.analyze.radmc3dGrid` class containing the spatial 
                    and wavelength/frequency grid of the model
        **Returns** an ndarray with [nx,ny,nz] dimensions

    * **getStellarsrcTemplates()** - Calculates the stellar templates for continuous star-like radiation sources 
        * **ppar** - Dictionary containing all parameters of the models with the parameter names as keys
        * **grid** - An instance of the :class:`~radmc3dPy.analyze.radmc3dGrid` class containing the spatial 
                    and wavelength/frequency grid of the model
        **Returns** an ndarray with [ntemplate,3] or [ntemplate,nwavelength] dimensions. If the returned array has [ntemplate,3] dimension, the [:,0], [:,1], [:,2] elements of the array should contain the stellar effective temperature (as negative numbers!), stellar radius and mass, respectively. If the returned dimension is [ntemplate, nwavelenght] the array should contain the frequency dependent spectrum of each stellar template. The setup functions will check the [0,0] element of the array. If it is negative, it assumes the template is defined as stellar temperature, radius and mass, if the [0,0] element is positive it assumes the array contains the full frequency dependent spectrum of each template.  

