.. _tutorial-octree:

***********************
Octree AMR in radmc3dPy
***********************

Apart from regular grids `RADMC-3D <http://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/>`_  supports two types of mesh refinements,
octree and layered adaptive mesh refinement (AMR, see Chapter 10 in the RADMC-3D manual 
(click `here <http://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/radmc-3d_v0.35.pdf>`_ for the manual of RADMC-3D v0.35)
both in caresian and spherical coordinate systems. As of v0.29 radmc3dPy supports regular grid and octree AMR (no layered mesh refinement / nested grid yet).

.. _tutorial-octree-grid-building:

Grid building 
=============

Octree AMR grid in `RADMC-3D <http://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/>`_ starts with a uniform, regular base grid.  The base grid must
have a fixed cell size, which means that we cannot use logarithmic grid for the radial coordinate in spherical coordinate system. For cartesian grids the cell size should
be the same in all active dimensions (i.e. cubical cells). After the base grid has been created, each cell
is tested for refinement based on user defined criteria. If the cell is refined, the cell is halved in each spatial dimension producing 2,4,8 sub-cells in 1,2,3 dimensions,
respectively. Each sub-cell is then tested for refinement, then each sub-sub-cell, etc. Thus cells are refined recursively building a tree. radmc3dPy provides a new class 
:class:`~radmc3dPy.analyze.radmc3dOctree` for handling octree AMR grid. To understand the model setup and the building of an octree mesh with radmc3dPy let us have a look 
at the grid structure and how they are stored in :class:`~radmc3dPy.analyze.radmc3dOctree`.

.. _tutorial-octree-grid-building-radmc3d:

RADMC-3D
--------
Octree AMR mesh in `RADMC-3D <http://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/>`_ starts from a regular base grid. For simplicity let us take a 2D grid (in which
case the grid is called a quadtree). We loop through all cells of the base grid and check if a cell is refined. If it is, we follow the refinement of that cell to the highest 
level before moving on to the next base grid cell. If we refine the bottom left and top right cells of the 2x2 base grid, then refining again the top right and bottom left of the
refined 2x2 cell structure, the buiding of the grid/tree in `RADMC-3D <http://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/>`_ proceeds this way:

.. image:: screenshots/build_radmc3d_allstages_grid_struct.png
    :align: center
    :width: 900px

.. _tutorial-octree-grid-building-radmc3dPy:

radmc3dPy
---------
Building the mesh in radmc3dPy this way can be very slow depending on the refinement criteria. If a base grid cell is resolved we need to loop over the
child cells of that base grid cell to test whether any of them need to be further refined. If so, we need to loop over the children of that child cell, etc. Thus building 
the grid consist of looping over every parent cell at every level and loop over their children to test whether they need to be further refined. 
This involves lot of nested loops which is usually very inefficient. It is much more efficient to do fewer but larger loops and to use numPy array operations. 
radmc3dPy does just that. First all cells at the highest actual refinement level are selected and tested for refinement. During the test and refinement all cells
processed together as arrays. Thus instead of looping over every parent as in case of RADMC-3D we loop over refinement level and process all cells simultaneously
at a given level. This way we have only one loop, over the refinement level. The steps of refinement for the example 2x2 base grid looks like this:

.. image:: screenshots/build_radmc3dPy_allstages_grid_struct.png
    :align: center
    :width: 900px

.. _tutorial-octree-array-layout:

Array layout
============

The grid as well as any physical variable (e.g. density, temperature, velocity, etc) is stored in linear arrays. The cell layout within the arrays of grid variables
(e.g. cell centers) is different depending on whether the grid is built on the fly by radmc3dPy or read from file (`amr_grid.inp`). Physical variables (e.g. dust density,
dust temperature, etc.) are stored in arrays of different length and layout compared to grid arrays. 

.. _tutorial-octree-array-layout-grid:

Grid array layout 
-----------------
Arrays containing information on the grid (i.e. cell centre coordinates, cell widths, boolean arrays to indicate a node is a leaf or a branch) etc, always contain the full
tree including both branch and leaf nodes. When the grid is read from file the array layout is the same as the order of cells in `amr_grid.inp`. The array starts with the 
list of the regular base grid cells. In the rest of the array we loop through the base grid cells again, but we add the higher resolution sub-trees for each branch node.
When a base grid cell is resolved, first we add all higher refinement levels before we would move on to the next base grid cell, i.e. we follow the tree immediately. 
This is the same layout as cell information is written to `amr_grid.inp`. The layout of such array looks like this:

.. image:: screenshots/build_radmc3d_allstages_arrayindex.png
    :align: center
    :width: 900px

Blue and red numbers indicate the index of the array element in the *full tree*, blue for leaf and red for branch nodes. The green numbers show the refinement level of the 
node, starting with zero for the regular base grid. Black numbers indicate the cell/leaf indices, i.e. the index of the cell/leaf node in an array containing only leaf nodes 
("true" cells). The arrows indicate parent/child relationships between nodes. In this layout all nodes with the same ancestor base grid cell are located next to each other. 


The order of the cells within an array is somewhat different if we build the grid with radmc3dPy. The reason is that the grid is built in a different way in radmc3dPy than
in RADMC-3D. As discussed above, radmc3dPy tests the refinement criteria for all leaf nodes at the actual highest refinement level and adds one level of refinement to *all* 
desired cells. This translates to the following layout of cells within the array. Just like in RADMC-3D style layout the array starts with the list of all base grid cells, but
it is followed by the list of *all* cells at the first level of refinement, then comes *all* cells at the second level, etc.. Thus the resulting array layout looks like this:

.. image:: screenshots/build_radmc3dPy_allstages_arrayindex.png
    :align: center
    :width: 900px

In this case all grid cells at the same refinement level are packed together within the array.  


.. _tutorial-octree-array-layout-data:

Data array layout 
-----------------
The layout of the data arrays of :class:`~radmc3dPy.analyze.radmc3dData` is different from that of the grid arrays, since data (i.e. 
physical variables) are only stored in leaf nodes, which are the "true" cells. Branch nodes are only used in the grid to navigate within the tree between leaf nodes. 
Therefore data arrays are smaller in size than grid arrays and they contain the leaf cells in increasing leaf/cell index order:

.. image:: screenshots/build_radmc3dPy_allstages_leafindex.png 
    :align: center
    :width: 900px

