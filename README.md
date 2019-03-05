# BlockLab
Blocky landscape evolution model, intended to simulate canyon development into layered rock landscapes.

Developed by Rachel Glade and Charles Shobe at the University of Colorado Department of Geological Sciences.

License: MIT.

Paper: Glade, R.C.\*, Shobe, C.M.\*, Anderson, R.S., and Tucker, G.E. (in revision) _Canyon evolution governed by channel-hillslope feedbacks_. (\* equal contributions)

A quick guide:

* brakeLandlab: The Blocky River and Knickpoint Evolution (BRaKE) model (Shobe et al., 2016; 2018) implemented in Landlab as a Python package.
* hogLab: The Hogback 2-D hillslope evolution model (Glade et al., 2017; Glade and Anderson, 2018) implemented in Landlab as a Python package.
* `cfuncs` files: Cython implementations of some of the slower parts of the model. Compile with `python setup.py build_ext --inplace`.
* `setup.py`: Run to compile the Cython functionality.
