.. _installation:

************
Installation
************

The first step to install radmc3dPy is to download the source of the latest version from :ref:`downloads`.

Once you obtained the source you need to extract the tarball::

    $>tar -xzvf radmc3dPy-0.30.2.tar.gz

There are two ways of installing radmc3dPy using distutils, depending on the
write access of Python's global site-packages directory. 

.. _installation-global:

Global installation
===================

Simplest install (write access of /usr or /opt is needed). Go to the radmc3dPy directory and execute the following command::
   
   $>python setup.py install

This will install radmc3dPy to Python's global site-packages directory. 

.. _installation-local:

Local installation
==================

In case you do not have writing permissions to this directory it is possible to use the following alternative::

   $>python setup.py install --user

This will install radmc3dPy to $HOME/.local/lib/python3.X/site-packages or  $HOME/.local/lib/python2.X/site-packages
depending on the python version used. 

Alternatively it is also possible to install radmc3dPy to an arbitrary directory by using the ``--prefix`` option::

   $>python setup.py install --prefix=DIR

where ``DIR`` is the full path of the directory one wishes to install radmc3dPy. Note, that in this case radmc3dPy will be installed
to DIR/lib/python2.X/site-packages/radmc3dPy or DIR/lib/python3.X/site-packages/radmc3dPy depending on the python version. Therefore, 
it is important to make sure that DIR/lib/python2.X/site-packages or DIR/lib/python3.X/site-packages is in the ``PYTHONPATH`` environment
variable. This can be done by e.g. adding the following lines to the .bashrc or .cscrc files, depending on the shell type used.

For bash add the following line to the .bashrc file::
    
    export PYTHONPATH=${PYTHONPATH}:DIR/lib/python3.X/site-packages

for c-shell add the following line to the .cshrc file::

    setenv PYTHONPATH ${PYTHONPATH}:DIR/lib/python3.X/site-packages

Again, replace the ``3.x`` tag with the proper python version and ``DIR`` with the path you used in the ``--prefix`` option. 

