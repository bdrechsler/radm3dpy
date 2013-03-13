# Import the radmc3dPy module
import radmc3dPy

# Write the parameter file with the default parameters
radmc3dPy.analyze.write_default_parfile('simple_1')

# Dust model setup with ascii input files
radmc3dPy.setup.problem_setup_dust('simple_1', binary=False)
