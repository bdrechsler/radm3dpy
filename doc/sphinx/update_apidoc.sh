rm -rf build/*
rm radmc3dPy/modules.rst
rm radmc3dPy/radmc3dPy.rst
rm radmc3dPy/radmc3dPy.models.rst
sphinx-apidoc -o radmc3dPy radmc3dPy
sphinx-build radmc3dPy build

#cp -rv build/* ../html/
