# Variable Selection with Missing Data
This is the instruction file to run the simulations for 'Variable Selection in the Missing data after Multiple Imputation with Adaptive Weights and Beyond'.

## Overview of the R codes
* ***degree_freedom.R*** : function for calculating degrees of freedom for each penalty under weighted least squares
* ***degree_freedom_noweight.R*** : function for calculating degrees of freedom for each penalty without any weights
* ***parallel_compute.R*** : code for running the simulation with multiple imputed data in parallel
* ***parallel_compute_fulldata.R*** : code for running the simulation with full data in parallel
* ***sim_aoas_data_generate.R*** : function for generating the data for simulation
* ***sim_imput_combined.R*** : function for obtaining results of the simulated multiple imputed data for all four data weighing scheme
* ***sim_imput_combined.R*** : function for obtaining results of the simulated full data for all four data weighing scheme
