.. _installation:

************
Installation
************

The first step to install radmc3dPy is to download the source of the latest version from :ref:`downloads`.

Once you obtained the source you need to extract the tarball::

    $>tar -xzvf radmc3dPy-0.29.0.tar.gz

There are two ways of installing radmc3dPy using distutils, depending on the
write access of Python's global site-packages directory. 

1. Simplest install (write access of /usr or /opt is needed)
   Go to the radmc3dPy directory and execute the following command::
   
   $>python setup.py install

This will install radmc3dPy to Python's global site-packages directory.

2. Manual install
Add the directory where this file is located to the PYTHONPATH environment variable, and
that is it. The best way to do so is in the .bashrc or .csrc files so it will automatically
be loaded at every shell startup. I.e. for bash add the following line to the 
.bashrc file::

    export PYTHONPATH=${PYTHONPATH}:ABSOLUTE-PATH-TO-THIS-DIRECTORY

or for c-shell add the following line to the .csrc file::

    setenv PYTHONPATH ${PYTHONPATH}:ABSOLUTE-PATH-TO-THIS-DIRECTORY

