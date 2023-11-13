# DENV_In_Vitro

This repository contains the data and source code files used to obtain the results described in the paper "Modelling the impact of JNJ-1802, a first-in-class dengue inhibitor blocking the NS3-NS4B interaction, on in-vitro DENV-2 dynamics" (Plos Computational Biology, 2023).  

The code was developed using Visual Studio 2019 (Version 16.11.25), and all output from the model is saved in .csv files.

Fixed parameter values can be inputted by reading in values from a file such as Fixed_Params.txt, and parameters to be estimated are specified in a file such as Fit_Params.txt. The 'Data' folder contains the following files:
- Data_Raw.txt (raw data)
- Data_LOD.txt (data used for model fitting, with left censoring at limit of detection during model fitting)
- Data_LOQ.txt (data used for model fitting, with left censoring at limit of quantification during model fitting)
