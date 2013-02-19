"""
PYTHON module for RADMC3D 
(c) Attila Juhasz 2011,2012,2013

This sub-module functions to set up a RADMC3D model for dust and/or line simulations
For help on the syntax or functionality of each function see the help of the individual functions

FUNCTIONS:
----------
    get_model_names()    - Returns the list of available models 
    get_model_desc()     - Returns the brief description of a model (if the model file contains a get_desc() function)
    problem_setup_dust() - Function to set up a dust model
    problem_setup_gas()  - Function to set up a line simulation
    write_lines_inp()    - Writes the lines.inp master command file for line simulations
    write_radmc3d_inp()  - Writes the radmc3d.inp master command file required for all RADMC3D runs

"""

try:
    from numpy import *
except:
    print 'ERROR'
    print ' Numpy cannot be imported '
    print ' To use the python module of RADMC-3D you need to install Numpy'


from radmc3dPy.natconst import *
try:
    from matplotlib.pylab import *
except:
    print ' WARNING'
    print ' matploblib.pylab cannot be imported ' 
    print ' To used the visualization functionality of the python module of RADMC-3D you need to install matplotlib'
    print ' Without matplotlib you can use the python module to set up a model but you will not be able to plot things or'
    print ' display images'

from radmc3dPy.crd_trans import vrot
from radmc3dPy.analyze import *
from subprocess import Popen, PIPE
import os, sys, copy



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def get_model_names():
    """
    Returns the name of the available models
    """

    mod_names = []
   
    # Get the name of all model files in the module directory
    import radmc3dPy
    mod_path = radmc3dPy.__file__.strip()[:-12]
    dum = Popen(['ls -1 '+mod_path+'/model_*.py'], shell=True, stdout=PIPE, stderr=PIPE).communicate()[0].split()

    for i in range(len(dum)):
        id1 = dum[i].index('model_')+6
        id2 = dum[i].index('.py')
        mod_names.append(dum[i][id1:id2])

    # Get the name of all model files in the current working directory
    if os.getcwd().strip()!=mod_path:
        mod_path = os.getcwd()
        dum = Popen(['ls -1 '+mod_path+'/model_*.py'], shell=True, stdout=PIPE, stderr=PIPE).communicate()[0].split()

        for i in range(len(dum)):
            id1 = dum[i].index('model_')+6
            id2 = dum[i].index('.py')
            mod_names.append(dum[i][id1:id2])

    return mod_names
    

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def get_model_desc(model=''):
    """
    Returns the brief description of the selected model
    """

    if model.strip()=='':
        print 'ERROR'
        print 'No model name is given'
        return

    try:
        mdl = __import__('model_'+model)
    except:
        try:
            mdl  = __import__('radmc3dPy.model_'+model, fromlist=['']) 
        except:
            print 'ERROR'
            print ' model_'+model+'.py could not be imported'
            print ' The model files should either be in the current working directory or'
            print ' in the radmc3d python module directory'
            return -1

    
    if callable(getattr(mdl, 'get_desc')):
        return mdl.get_desc()
    else:
        print 'ERROR'
        print 'model_'+model+'.py does not contain a get_desc() function.'
        return





# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def problem_setup_dust(model='', **kwargs):
    """
    Function to set up a dust model for RADMC3D 
    
    INPUT:
    ------
        model : Name of the model that should be used to create the density structure.
                The file should be in a directory from where it can directly be imported 
                (i.e. the directory should be in the PYTHON_PATH environment variable or
                it should be in the current working directory)
                and the file name should be 'model_xxx.py', where xxx stands for the string
                that should be specified in this variable

    OPTIONS:
    --------
        Any varible name in problem_params.inp can be used as a keyword argument.
        At first all variables are read from problem_params.in to a dictionary called ppar. Then 
        if there is any keyword argument set in the call of problem_setup_dust the ppar dictionary 
        is searched for this key. If found the value belonging to that key in the ppar dictionary 
        is changed to the value of the keyword argument. If no such key is found then the dictionary 
        is simply extended by the keyword argument. 
       
    FILES WRITTEN DURING THE SETUP:
    -------------------------------
        amr_grid.inp - spatial grid
        wavelength_micron.inp - wavelength grid
        dust_density.inp - dust density distribution
        dustopac.inp - dust opacity master file
        stars.inp - input radiation field 
        radmc3d.inp - parameters for RADMC3D (e.g. Nr of photons to be used, scattering type, etc)
    """
  
    # Read the parameters from the problem_params.inp file 
    dum = readparams()
    ppar = dum.ppar


    if model=='':
        print 'ERROR'
        print 'No model name is given'
        return

    if not ppar:
        print 'ERROR'
        print 'problem_params.inp was not found'
        return 
# --------------------------------------------------------------------------------------------
# If there is any additional keyword argument (**kwargs) then check
#   if there is such key in the ppar dictionary and if is change its value that of
#   the keyword argument. If there is no such key in the ppar dictionary then add the keyword
#   to the dictionary
# --------------------------------------------------------------------------------------------
    if kwargs:
        for keys in kwargs.keys():
            ppar[keys] = kwargs[keys]

# --------------------------------------------------------------------------------------------
# Create the grid
# --------------------------------------------------------------------------------------------
    grid = radmc3dGrid()
    
    # Wavelength grid
    grid.make_wav_grid(ppar=ppar)

    # Spatial grid
    grid.make_spatial_grid(ppar=ppar)

# --------------------------------------------------------------------------------------------
# Dust opacity
# --------------------------------------------------------------------------------------------
    if ppar.has_key('dustkappa_ext'):
        opac=radmc3dDustOpac()
        #Master dust opacity file
        opac.write_masteropac(ext=ppar['dustkappa_ext'])
    else:
        opac=radmc3dDustOpac()
        # Calculate the opacities and write the master opacity file
        opac.makeopac(ppar=ppar)
# --------------------------------------------------------------------------------------------
# Create the input radiation field (stars at this point) 
# --------------------------------------------------------------------------------------------

    stars = radmc3dStars(ppar=ppar)
    stars.get_stellar_spectrum(tstar=ppar['tstar'], rstar=ppar['rstar'], wav=grid.wav)

# --------------------------------------------------------------------------------------------
# Create the dust density distribution 
# --------------------------------------------------------------------------------------------

    try:
        mdl = __import__('model_'+model)
    except:
        try:
            mdl  = __import__('radmc3dPy.model_'+model, fromlist=['']) 
        except:
            print 'ERROR'
            print ' model_'+model+'.py could not be imported'
            print ' The model files should either be in the current working directory or'
            print ' in the radmc3d python module directory'
            return

    data = radmc3dData(grid)
    if dir(mdl).__contains__('get_dust_density'):
        if callable(getattr(mdl, 'get_dust_density')):
            data.rhodust = mdl.get_dust_density(grid=grid, ppar=ppar)
        else:
            print 'WARNING'
            print 'model_'+model+'.py does not contain a get_dust_density() function, therefore, '
            print ' dust_density.inp cannot be written'
            return 
    else:
        print 'WARNING'
        print 'model_'+model+'.py does not contain a get_dust_density() function, therefore, '
        print ' dust_density.inp cannot be written'
        return 
    #data.rhodust = mdl.get_density(grid=grid, ppar=ppar) * ppar['dusttogas']
# --------------------------------------------------------------------------------------------
# Now write out everything 
# --------------------------------------------------------------------------------------------

    #Frequency grid
    grid.write_wav_grid()
    #Spatial grid
    grid.write_spatial_grid()
    #Input radiation field
    stars.write_starsinp(wav=grid.wav, pstar=ppar['pstar'], tstar=ppar['tstar'])
    #Dust density distribution
    data.write_dustdens()
    #radmc3d.inp
    write_radmc3d_inp(ppar=ppar)

# --------------------------------------------------------------------------------------------
# Calculate optical depth for diagnostics purposes
# --------------------------------------------------------------------------------------------
    
    stars = radmc3dStars()
    stars.read_starsinp()
    pwav = stars.find_peak_starspec()[0]
    #data.get_tau(wav=pwav, usedkappa=False)
    #print 'Radial optical depth at '+("%.2f"%pwav)+'um : ', data.taux.max()

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def problem_setup_gas(model='', fullsetup=False, **kwargs):
    """
    Function to set up a gas model for RADMC3D 
    
    INPUT:
    ------
        model : Name of the model that should be used to create the density structure
                the file should be in a directory from where it can directly be imported 
                (i.e. the directory should be in the PYTHON_PATH environment variable, or
                it should be the current working directory)
                and the file name should be 'model_xxx.py', where xxx stands for the string
                that should be specified in this variable
    
        fullsetup : if False only the files related to the gas simulation is written out
                    (i.e. no grid and radmc3d master command file is written)
                    if True the spatial and wavelength grid as well as the input radiation field
                    and the radmc3d master command file will be (over)written. 

    OPTIONS:
    --------
        Any varible name in problem_params.inp can be used as a keyword argument.
        At first all variables are read from problem_params.in to a dictionary called ppar. Then 
        if there is any keyword argument set in the call of problem_setup_dust the ppar dictionary 
        is searched for such key. If found the value belonging to that key in the ppar dictionary 
        is changed to the value of the keyword argument. If no such key is found then the dictionary 
        is simply extended by the keyword argument. 

       
    FILES WRITTEN DURING THE SETUP:
    -------------------------------
        fullsetup = True
            amr_grid.inp - spatial grid
            wavelength_micron.inp - wavelength grid
            stars.inp - input radiation field 
            radmc3d.inp - parameters for RADMC3D (e.g. Nr of photons to be used, scattering type, etc)
            lines.inp -  line mode master command file
            numberdens_xxx.inp - number density of molecule/atomic species 'xxx'
            gas_velocity.inp - Gas velocity
            microturbulence.inp - The standard deviation of the Gaussian line profile caused by turbulent 
                                  broadening (doublecheck if it is really the standard deviation or a factor
                                  of sqrt(2) less than that!)
            gas_temperature.inp - Gas temperature (which may be different from the dust temperature). If
                                  tgas_eq_tdust is set to zero in radmc3d.inp the gas temperature in this
                                  file will be used instead of the dust temperature. 
        fullsetup = False
            lines.inp -  line mode master command file
            numberdens_xxx.inp - number density of molecule/atomic species 'xxx'
            gas_velocity.inp - Gas velocity
            microturbulence.inp - The standard deviation of the Gaussian line profile caused by turbulent 
                                  broadening (doublecheck if it is really the standard deviation or a factor
                                  of sqrt(2) less than that!)
            gas_temperature.inp - Gas temperature (which may be different from the dust temperature). If
                                  tgas_eq_tdust is set to zero in radmc3d.inp the gas temperature in this
                                  file will be used instead of the dust temperature. 

    """
    
    dum = readparams()
    ppar = dum.ppar

    if not ppar:
        print 'ERROR'
        print 'problem_params.inp was not found'
        return
    
# --------------------------------------------------------------------------------------------
# If there is any additional keyword argument (**kwargs) then check
#   if there is such key in the ppar dictionary and if is change its value that of
#   the keyword argument. If there is no such key in the ppar dictionary then add the keyword
#   to the dictionary
# --------------------------------------------------------------------------------------------
    if kwargs:
        for keys in kwargs.keys():
            ppar[keys] = kwargs[keys]

# --------------------------------------------------------------------------------------------
# If the current working directory is empty (i.e. no dust setup is present) then
#   we must make a complete setup and dump the spatial and wavelength grids as well
#   as the parameters in the radmc3d.inp file
# --------------------------------------------------------------------------------------------
    if fullsetup:

# --------------------------------------------------------------------------------------------
# Create the grid
# --------------------------------------------------------------------------------------------
        grid = radmc3dGrid()
    
        # Wavelength grid
        grid.make_wav_grid(ppar=ppar)

        # Spatial grid
        grid.make_spatial_grid(ppar=ppar)

# --------------------------------------------------------------------------------------------
# Create the input radiation field (stars at this point) 
# --------------------------------------------------------------------------------------------

        stars = radmc3dStars(ppar=ppar)
        stars.get_stellar_spectrum(tstar=ppar['tstar'], rstar=ppar['rstar'], wav=grid.wav)

# --------------------------------------------------------------------------------------------
# Now write out everything 
# --------------------------------------------------------------------------------------------

        #Frequency grid
        grid.write_wav_grid()
        #Spatial grid
        grid.write_spatial_grid()
        #Input radiation field
        stars.write_starsinp(wav=grid.wav, pstar=ppar['pstar'], tstar=ppar['tstar'])
        #radmc3d.inp
        write_radmc3d_inp(ppar=ppar)
# --------------------------------------------------------------------------------------------
# If the current working directory contains already a dust setup then we can use the
#   already existing grid files 
# --------------------------------------------------------------------------------------------
    else:
        grid=read_grid()
# --------------------------------------------------------------------------------------------
# Create the gas density distribution 
# --------------------------------------------------------------------------------------------
    try:
        import os
        imp_path = os.getcwd()
        mdl = __import__('model_'+model)
    except:
        try:
            mdl  = __import__('radmc3dPy.model_'+model, fromlist=['']) 
        except:
            print 'ERROR'
            print ' model_'+model+'.py could not be imported'
            print ' The model files should either be in the current working directory or'
            print ' in the radmc3d python module directory'
            return 

    # Create the data structure
    data = radmc3dData(grid)
    # Calculate the gas density and velocity
    # NOTE: the density function in the model sub-modules should provide the gas volume density
    #       in g/cm^3 but RADMC3D needs the number density in 1/cm^3 so we should convert the
    #       output of the get_density() function to number density using ppar['gasspec_abun']
    #       which is the abundance of the gas species with respect to hydrogen divided by the
    #       mean molecular weight
    if dir(mdl).__contains__('get_gas_density'):
        if callable(getattr(mdl, 'get_gas_density')):
            data.rhogas = mdl.get_gas_density(grid=grid, ppar=ppar) / (2.4*mp)
    else:
        print 'WARNING'
        print 'model_'+model+'.py does not contain a get_gas_density() function, therefore, '
        print ' numberdens_***.inp cannot be written'
        return 

    if ppar['gasspec_abun']:
        data.rhogas = data.rhogas * ppar['gasspec_abun']
    else:
        if dir(mdl).__contains__('get_gas_abundance'):
            if callable(getattr(mdl, 'get_gas_abundance')):
                data.gasabun = mdl.get_gas_abundance(grid=grid, ppar=ppar)
        else:
            print 'WARNING'
            print 'model_'+model+'.py does not contain a get_gas_abundance() function, and no "gasspec_abun" '
            print ' parameter is found in the problem_setup.inp file. numberdens_***.inp cannot be written'
            return



    #data.rhogas = mdl.get_density(grid=grid, ppar=ppar) * ppar['gasspec_abun']
    data.gasvel = mdl.get_velocity(grid=grid, ppar=ppar)
    # Write the gas density
    data.write_gasdens(ispec=ppar['gasspec_name'])
    # Write the gas velocity
    data.write_gasvel()
    # Write the gas temperature if specified 
    if ppar['write_gastemp']:
        if dir(mdl).__contains__('get_gastemp'):
            if callable(getattr(mdl, 'get_gastemp')):
                data.gastemp = mdl.get_gastemp(grid=grid, ppar=ppar)
                # Write the gas temperature
                data.write_gastemp() 
        else:
            print 'WARNING'
            print 'model_'+model+'.py does not contain a get_gastemp() function, therefore, '
            print ' gas_temperature.inp cannot be written'
            return

    # Write the turbulent velocity field (important for the width of the spectral lines)
    # If the ppar dictionary does not have a key 'gasspec_vturb' then check if the model
    # sub-modul contains a fuction 'vturb()' that provides the turbulent velocity field 
    if not ppar.has_key('gasspec_vturb'):
        try:
            if callable(getattr(mdl, 'get_vturb')):
                data.vturb = mdl.get_vturb(grid=grid, ppar=ppar)
                # Write the turbulent velocity field
                data.write_vturb() 
        except:
            print 'WARNING'
            print 'model_'+model+'.py does not contain a get_vturb() function, nor does problem_params.inp contain '
            print " a gasspec_vturb variable therefore, microturbulence.inp cannot be written"
            return

    else:
        # If there is a 'gasspec_vturb' key in the ppar dictionary but it is set to a negative value
        # then again look at the model sub-module if it contains a provider function 'vturb()'
        if (ppar['gasspec_vturb']<0.):
            if dir(mdl).__contains__('get_vturb'):
                if callable(getattr(mdl, 'get_vturb')):
                    data.vturb = mdl.get_vturb(grid=grid, ppar=ppar)
                    # Write the turbulent velocity field
                    data.write_vturb() 
            else:
                data.vturb = zeros([grid.nx, grid.ny, grid.nz], dtype=float64)
                data.vturb[:,:,:] = 0.
                data.write_vturb()
                return
        else:
            # If ppar['gasspec_vturb'] is set to a positive value then it is assumed that
            # the turbulent velocity field is homogeneous and its value is ppar['gasspec_vturb']  
            data.vturb = zeros([grid.nx, grid.ny, grid.nz], dtype=float64)
            data.vturb[:,:,:] = ppar['gasspec_vturb']
            data.write_vturb()

    # Write the lines.inp the main control file for the line RT
    write_lines_inp(ppar=ppar)
# --------------------------------------------------------------------------------------------------
def write_radmc3d_inp(ppar=None):
    """
    Function to write the radmc3d.inp master command file for RADMC3D

    INPUT:
    ------
        ppar : dictionary containing all parameters of a RADMC3D run 

    """

    print 'Writing radmc3d.inp'

    wfile = open('radmc3d.inp', 'w')
    radmc3d_runkeys = ['nphot', 'nphot_scat', 'scattering_mode_max', 'lines_mode', 'istar_sphere', 'itempdecoup', 'tgas_eq_tdust']
    radmc3d_runkeys.sort()

    for key in radmc3d_runkeys:
        if ppar.has_key(key):
            wfile.write('%s %d\n'%(key+' =',ppar[key]))

    wfile.close()

# --------------------------------------------------------------------------------------------------
def write_lines_inp(ppar=None):
    """
    Function to write the lines.inp master command file for line simulation in RADMC3D

    INPUT:
    ------
        ppar : dictionary containing all parameters of a RADMC3D run 

    """

    print 'Writing lines.inp'
    wfile = open('lines.inp', 'w')
    # File format
    wfile.write("%d\n"%(1))
    # Nr of gas species
    wfile.write("%d\n"%(1))
    # Gas species name and database type
    wfile.write("%s %s %d %d\n"%(ppar['gasspec_name'], ppar['gasspec_dbase_type'], 0, 0))
    wfile.close()

# --------------------------------------------------------------------------------------------------
