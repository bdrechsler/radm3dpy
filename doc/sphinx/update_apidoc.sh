rm -rf build/*
rm radmc3dPy/modules.rst
rm radmc3dPy/radmc3dPy.rst
rm radmc3dPy/radmc3dPy.models.rst
sphinx-apidoc-2.7 -o radmc3dPy radmc3dPy
sphinx-build-2.7 radmc3dPy build
