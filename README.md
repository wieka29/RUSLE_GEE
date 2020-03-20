# RUSLE_GEE
The RUSLE model on Google Earth Engine

This model is entirely written in Python code that needs to be transformed to make the model run on GEE with the help of Python API for GEE. The model consists out of 3 scripts: 
  Input_File.py : derives input data from GEE and exports them to tif files which are then read and output is transformed to arrays. The needed processing is done (reprojection, selection of test region etc.)
  Adj_RUSLE.py: This is the main model, where the RUSLE factors are calculated based on different methods from literature
  Output_File.py: The run script

In addition the erosion model has 2 functions that need to be imported to make the model run:
  erosivity.so and erosivity.f90
  scaling_func.py
