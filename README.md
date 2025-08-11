# Cryspy: Computational Crystallography in Python: *A Python Toolbox for (some, specialized) EBSD Data Analyses*

Cryspy is a toolbox for computational crystallography in Python.
It intends to provide a set of open source tools for
visualization, analysis, and postprocessing of EBSD data.

## COMPONENTS

* io   : modules for import and export of EBSD data
* rot  : modules for rotation representation, manipulation, and conversion
* vis  : modules for data visualization
* ebsd : modules for EBSD data analysis and manipulation
* xtal : modules for crystallography and symmetry calculations
* util : utilities

## CONTENT HIGHLIGHTS

io
------------------------------------------------------------------------------
* loadang:: support for loading TSL *.ang files

rot
------------------------------------------------------------------------------
* angax  : class for the angle/axis convention
* bunge  : class for Bunge Euler angles
* quat   : class for quaternion representations
* rodri  : class for Rodrigues vector representations
* tvec   : class for transformation vectors (reshaped rotation matrices)

xtal
------------------------------------------------------------------------------
* lattvec : class for representing lattice vectors
* miller  : class for representing lattice plane normals
* unitcell: class for representing unit cells
* rotsymm : class for rotational symmetry objects

vis
------------------------------------------------------------------------------
* stereoproj : class for plotting stereographic projections
* eaproj     : class for plotting equal area projections

util
------------------------------------------------------------------------------
* rationalize : function for turning fractional vectors into rational indices
* sigdec      : function for rounding to n significant digits
