#
# The sdist option will create a source distribution packge in the 'dist' directory
# The files included in the source distribution are specified in MANIFEST.in 
# The MANIFEST file is created when a package is made and contains the actually
# packaged files. The version number of the package should be present in the setup.py file 
#  
python setup.py sdist

# INSTALL
# 
# Untar the tarball and run
# $>python setup install 
# This is going to build the package according to the specification in the setup.py
# and also installs it to the 'site-packages' directory 
# 
# To only build the package without installation run
# $>python setup build
