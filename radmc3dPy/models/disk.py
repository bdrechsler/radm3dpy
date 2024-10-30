"""This is a radmc3dPy model template 

This template is an empty model, i.e. all model functions return zeros in the appropriate arrays and dimensions. 
The purpose of this model is to illustrate the names and syntax of the model functions. Hence, this file can
be a starting point for implementing new models in the library. 

A radmc3dPy model file can contain any / all of the functions below

    * getDefaultParams()
    * getModelDesc()
    * getDustDensity()
    * getDustTemperature()
    * getGasAbundance()
    * getGasDensity()
    * getGasTemperature()
    * getVelocity()
    * getVTurb()

The description of the individual functions can be found in the docstrings below the function name.
If a model does not provide a variable or the variable should be calculated by RADMC-3D 
(e.g. dust temperature) the corresponding function (e.g. get_dust_temperature) should be removed from
or commented out in the model file. 

NOTE: When using this template it is strongly advised to rename the template model (to e.g. mydisk.py)
as the get_model_names() function in the setup module removes the name 'template' from the list of available
models. 

"""

from __future__ import absolute_import
from __future__ import print_function
import traceback

try:
    import numpy as np
except ImportError:
    np = None
    print(" Numpy cannot be imported ")
    print(" To use the python module of RADMC-3D you need to install Numpy")
    print(traceback.format_exc())

from ..natconst import *
from scipy.integrate import trapezoid


def getModelDesc():
    """Provides a brief description of the model"""

    return "basic disk model"


def getDefaultParams():
    """Provides default parameter values


    Returns a list whose elements are also lists with three elements:
    1) parameter name, 2) parameter value, 3) parameter description
    All three elements should be strings. The string of the parameter
    value will be directly written out to the parameter file if requested,
    and the value of the string expression will be evaluated and be put
    to radmc3dData.ppar. The third element contains the description of the
    parameter which will be written in the comment field of the line when
    a parameter file is written.
    """

    # defpar = [
    #     ["mstar", "1.0*ms", "Mass of the star(s)"],
    #     ["rin", "0.1*au", "Inner disk radius"],
    #     ["rdisk", "100*au", "Outter disk radius"],
    #     ["rhopl", "2.37", "Density power law exponent"],
    #     ["h0", "0.1*au", "Scale height at 1 au"],
    #     ["hpl", "1.29", "Scale height power law exponent"],
    #     ["t0", "315.", "Temperature at 1 au"],
    #     ["tpl", "0.25", "Temperature power law exponent"],
    #     ["gasspec_vturb", "1e5", "Microturbulent line width"],
    # ]
    defpar = [
        ["xres_nlev", "3", "Number of refinement levels"],
        ["xres_nspan", "3", "Number of the original grid cells to refine"],
        ["xres_nstep", "3", "Number of grid cells to create in a refinement level"],
        ["nx", "[30,100]", "Number of grid points in the first dimension"],
        ["xbound", "[0.01*au,1.05*au, 300.0*au]", "Number of radial grid points"],
        ["ny", "50", "Number of grid points in the first dimension"],
        ["ybound", "[0.01, pi / 2]", "Number of radial grid points"],
        ["nz", "1", "Number of grid points in the first dimension"],
        ["zbound", "[0., 2.0*pi]", "Number of radial grid points"],
        ["gasspec_mol_name", "['co']", ""],
        ["gasspec_mol_abun", "[1e-4]", ""],
        ["gasspec_mol_dbase_type", "['leiden']", ""],
        # ['gasspec_mol_dissoc_taulim', '[1.0]', 'Continuum optical depth limit below which all molecules dissociate'],
        # ['gasspec_mol_freezeout_temp', '[19.0]', 'Freeze-out temperature of the molecules in Kelvin'],
        # ['gasspec_mol_freezeout_dfact', '[1e-3]',
        # 'Factor by which the molecular abundance should be decreased in the frezze-out zone'],
        ["gasspec_vturb", "0.2e5", "Microturbulent line width"],
        ["rin", "0.1*au", " Inner radius of the disk"],
        ["rdisk", "100.0*au", " Outer radius of the disk"],
        # ['hrdisk', '0.1', ' Ratio of the pressure scale height over radius at hrpivot'],
        # ['hrpivot', "100.0*au", ' Reference radius at which Hp/R is taken'],
        ["h0", "0.1*au", "Scale height at 1 au"],
        ["hpl", "1.29", "Scale height power law exponent"],
        # ['plsig1', '-1.0', ' Power exponent of the surface density distribution as a function of radius'],
        # ['sig0', '0.0', ' Surface density at rdisk'],
        [
            "mdisk",
            "1e-3*ms",
            " Mass of the disk (either sig0 or mdisk should be set to zero or commented out)",
        ],
        ["t0", "315.", "Temperature at 1 au"],
        ["tpl", "0.25", "Temperature power law exponent"],
        ["rhopl", "2.37", "Density power law exponent"],
        # ['bgdens', '1e-30', ' Background density (g/cm^3)'],
        # ['srim_rout', '0.0', 'Outer boundary of the smoothing in the inner rim in terms of rin'],
        # ['srim_plsig', '0.0', 'Power exponent of the density reduction inside of srim_rout*rin'],
        # ['prim_rout', '0.0', 'Outer boundary of the puffed-up inner rim in terms of rin'],
        # ['hpr_prim_rout', '0.0', 'Pressure scale height at rin'],
        # ['gap_rin', '[0e0*au]', ' Inner radius of the gap'],
        # ['gap_rout', '[0e0*au]', ' Outer radius of the gap'],
        # ['gap_drfact', '[0e0]', ' Density reduction factor in the gap'],
        # ['sigma_type', '0',
        #  ' Surface density type (0 - polynomial, 1 - exponential outer edge (viscous self-similar solution)'],
        ["dusttogas", "0.01", " Dust-to-gas mass ratio"],
    ]

    return defpar


def getGasTemperature(grid=None, ppar=None):
    """Calculates/sets the gas temperature

    Parameters
    ----------
    grid : radmc3dGrid
            An instance of the radmc3dGrid class containing the spatial and wavelength grid

    ppar : dictionary
            Dictionary containing all parameters of the model

    Returns
    -------
    Returns the gas temperature in K
    """
    # get coordinates (assume spherical)
    r_sph, th, ph = np.meshgrid(grid.x, grid.y, grid.z, indexing="ij")
    r_cyl = r_sph * np.sin(th)  # cylindrical radius
    zz = r_sph * np.cos(th)  # z coordinate

    # load in relavant parameters
    t0 = ppar["t0"]  # temp at 1 au
    tpl = ppar["tpl"]  # temp power law
    rdisk = ppar["rdisk"]  # outer disk radius
    rin = ppar["rin"]  # inner disk radius

    # initialize temperature array
    tgas = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64)

    # calculate temperature everywhere in the grid
    t = t0 * (r_cyl / au) ** (-tpl)  # create temperature array
    t[(r_cyl >= rdisk) ^ (r_cyl <= rin)] = 0e0  # set temperature outside of disk to 0
    t[t > 10000.0] = 10000.0  # enforce 10,000K as max temperature

    # populate the tgas array
    tgas[:, :, :] = t

    return tgas


def getDustTemperature(grid=None, ppar=None):
    """Calculates/sets the dust temperature

    Parameters
    ----------
    grid : radmc3dGrid
            An instance of the radmc3dGrid class containing the spatial and wavelength grid

    ppar : dictionary
            Dictionary containing all parameters of the model

    Returns
    -------
    Returns the dust temperature in K
    """

    # assume tdust = tgas
    # initialize tdust array
    tdust = np.zeros([grid.nx, grid.ny, grid.nz, 1], dtype=np.float64)

    # get gas temperature
    t = getGasTemperature(grid=grid, ppar=ppar)
    # I think the extra index is to distinguish multiple dust species
    # start with one dust species for simplicity
    tdust[:, :, :, 0] = t
    return tdust


def getGasAbundance(grid=None, ppar=None, ispec=""):
    """Calculates/sets the molecular abundance of species ispec
    The number density of a molecule is rhogas * abun

    Parameters
    ----------
    grid  : radmc3dGrid
            An instance of the radmc3dGrid class containing the spatial and wavelength grid

    ppar  : dictionary
            Dictionary containing all parameters of the model

    ispec : str
            The name of the gas species whose abundance should be calculated

    Returns
    -------
    Returns the abundance as an ndarray
    """

    gasabun = -1
    if ppar["gasspec_mol_name"].__contains__(ispec):
        ind = ppar["gasspec_mol_name"].index(ispec)
        gasabun[:, :, :] = ppar["gasspec_mol_abun"][ind]

    elif ppar["gasspec_colpart_name"].__contains__(ispec):
        gasabun = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64)
        ind = ppar["gasspec_colpart_name"].index(ispec)
        gasabun[:, :, :] = ppar["gasspec_colpart_abun"][ind]
    else:
        raise ValueError(
            ' The abundance of "' + ispec + '" is not specified in the parameter file'
        )

    return gasabun


def getGasDensity(grid=None, ppar=None):
    """Calculates the total gas density distribution

    Parameters
    ----------
    grid : radmc3dGrid
            An instance of the radmc3dGrid class containing the spatial and wavelength grid

    ppar : dictionary
            Dictionary containing all parameters of the model

    Returns
    -------
    Returns the gas volume density in g/cm^3
    """

    # set up the coordinates (assume spherical)
    r_sph, th, ph = np.meshgrid(grid.x, grid.y, grid.z, indexing="ij")
    r_cyl = r_sph * np.sin(th)  # cylindrical radius coordinate
    zz = r_sph * np.cos(th)  # z coordinate

    # get relavant parameters
    rin = ppar["rin"]
    rdisk = ppar["rdisk"]
    mdisk = ppar["mdisk"]
    h0 = ppar["h0"]
    hpl = ppar["hpl"]
    gamma = ppar["rhopl"] - ppar["hpl"]

    # function to calculate the surface density
    def surface_density(r, normalize=True):

        # Sigma0 = (2-gamma)*mdisk/(2*np.pi*au**(gamma)) / (rdisk**(-gamma+2) - rin**(-gamma+2))
        # Sigma0 = (2 - gamma) * mdisk / (2 * np.pi * rdisk**2)
        # Sigma = Sigma0 * (r / rdisk) ** (-gamma) * np.exp(-((r / rdisk) ** (2 - gamma)))
        Sigma0 = (
            (2 - gamma)
            * mdisk
            / (2 * np.pi * (1 * au) ** (gamma))
            / (rdisk ** (-gamma + 2) - rin ** (-gamma + 2))
        )

        Sigma = Sigma0 * (r / au) ** (-gamma)

        # set surface density of disk to 0 outside of the disk
        Sigma[(r >= rdisk) ^ (r <= rin)] = 0e0

        # deal with singularity at r==0
        # dr = r[r > 0].min()  # get the smallest (positive) value of r
        # Sigma[r == 0] = (
        #     Sigma0
        #     * (0.7 * dr / rdisk) ** (-gamma)
        #     * np.exp(-((0.7 * dr / rdisk) ** (2 - gamma)))
        # )

        if normalize:
            r_high = np.logspace(np.log10(rin), np.log10(rdisk), 1000)
            Sigma_high = surface_density(r_high, normalize=False)

            scale = mdisk / (2 * np.pi * trapezoid(r_high * Sigma_high, r_high))

            Sigma *= scale

        return Sigma

    # function to calculate the scale height
    def scale_height(r):
        return h0 * (r / au) ** hpl

    # initialize density array
    rhogas = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64)
    # get surface density
    Sigma = surface_density(r_cyl)
    # get scale height
    h = scale_height(r_cyl)
    # calculate density
    rho = Sigma / (np.sqrt(2 * np.pi) * h) * np.exp(-0.5 * (zz / h) ** 2)

    rhogas[:, :, :] = rho

    return rhogas


def getDustDensity(grid=None, ppar=None):
    """Calculates the dust density distribution

    Parameters
    ----------
    grid : radmc3dGrid
            An instance of the radmc3dGrid class containing the spatial and wavelength grid

    ppar : dictionary
            Dictionary containing all parameters of the model

    Returns
    -------
    Returns the dust volume density in g/cm^3
    """
    # get the gas density
    rhogas = getGasDensity(grid=grid, ppar=ppar)
    # initialize the dust density grid (assume one dust species)
    rhodust = np.zeros([grid.nx, grid.ny, grid.nz, 1], dtype=np.float64)
    # convert gas to dust density with dust to gas ratio
    rhodust[:, :, :, 0] = rhogas * ppar["dusttogas"]
    return rhodust


def getVTurb(grid=None, ppar=None):
    """Calculates/sets the turbulent velocity field

    Parameters
    ----------
    grid : radmc3dGrid
            An instance of the radmc3dGrid class containing the spatial and wavelength grid

    ppar : dictionary
            Dictionary containing all parameters of the model

    Returns
    -------
    Returns the turbulent velocity in cm/s
    """

    vturb = (
        np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64) + ppar["gasspec_vturb"]
    )
    return vturb


def getVelocity(grid=None, ppar=None):
    """Calculates/sets the gas velocity field

    Parameters
    ----------
    grid : radmc3dGrid
            An instance of the radmc3dGrid class containing the spatial and wavelength grid

    ppar : dictionary
            Dictionary containing all parameters of the model

    Returns
    -------
    Returns the turbulent velocity in cm/s
    """

    # set up the coordinates (assume spherical)
    r_sph, th, ph = np.meshgrid(grid.x, grid.y, grid.z, indexing="ij")
    r_cyl = r_sph * np.sin(th)  # cylindrical radius coordinate
    zz = r_sph * np.cos(th)  # z coordinate

    # get relavant parameters
    mstar = ppar["mstar"]

    # set azimuthal velocity to keplerian
    v_kep = np.sqrt(gg * mstar * r_cyl**2 / r_sph**3)

    # initialize velocity array
    vel = np.zeros([grid.nx, grid.ny, grid.nz, 3], dtype=np.float64)

    # set phi velocity to keplerian velocity
    vel[:, :, :, 2] = v_kep
    return vel
