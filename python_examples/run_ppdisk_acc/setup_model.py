# Import the radmc3dPy module
import radmc3dPy
import os

# Write the parameter file with the default parameters
radmc3dPy.analyze.writeDefaultParfile('ppdisk_acc')

# Dust model setup with ascii input files
radmc3dPy.setup.problemSetupDust('ppdisk_acc', binary=False)

# Copy the dust opacity 
os.system('cp -v ../datafiles/dustkappa_silicate.inp .')

# Calculate the dust temperature

os.system('radmc3d mctherm')

