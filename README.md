# BlockLab
Blocky landscape evolution model, intended to simulate canyon development into layered rock landscapes.

Developed by Rachel Glade and Charles Shobe at the University of Colorado Department of Geological Sciences.

License: MIT.

[![DOI](https://zenodo.org/badge/151611886.svg)](https://zenodo.org/badge/latestdoi/151611886)

Paper: Glade, R.C.\*, Shobe, C.M.\*, Anderson, R.S., and Tucker, G.E. (_in revision_) Canyon evolution governed by channel-hillslope feedbacks. (\* equal contributions)

A quick guide:

* brakeLandlab: The Blocky River and Knickpoint Evolution (BRaKE) model (Shobe et al., 2016; 2018) implemented in Landlab as a Python package.
* hogLab: The Hogback 2-D hillslope evolution model (Glade et al., 2017; Glade and Anderson, 2018) implemented in Landlab as a Python package.
* driver scripts: Each one of these corresponds to a model realization presented in Glade, Shobe, et al (_in revision_). See the paper for a description of model parameters and how they vary between driver scripts.
* `cfuncs` files: Cython implementations of some of the slower parts of the model. Compile with `python setup.py build_ext --inplace`.
* `setup.py`: Run to compile the Cython functionality.
