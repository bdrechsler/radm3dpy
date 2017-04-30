rm -rf build/*
rm radmc3dPy/modules.rst
rm radmc3dPy/radmc3dPy.rst
rm radmc3dPy/radmc3dPy.models.rst
sphinx-apidoc-3.5 -o radmc3dPy radmc3dPy
sphinx-build-3.5 radmc3dPy build
