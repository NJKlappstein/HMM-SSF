# HMM-SSF
Fit state-switching step selection functions to animal movement data

Natasha Klappstein

Working directory note: All scripts were written with the root directory of the repository as the working directory. Sourcing the necessary files must be done from this directory. 

This repository contains the code to reproduce the illustrative example of Klappstein et al. ([https://doi.org/10.1101/2022.11.30.518554](https://doi.org/10.1101/2022.11.30.518554)). The habitat raster and zebra tracking data were provided by Simon Chamaill\'e-Jammes. This README describes 3 folders (`code`, `data`, and `documents`). Note, this repo is a work-in-progress and will be updated periodically.

## Code folder
There are 2 main folders within the code folder. The main functions for data processing and model fitting are found in `functions`. Code to fit the HMM-SSF to the zebra data is found `illustration`, and this folder includes all code to produce the plots in the manuscript. There are also two RMarkdown files found in `documents` that describe the steps of the data processing (`HMM-SSF_data.Rmd`) and model fitting/plotting/interpretation (`HMM-SSF_fitting.Rmd`). 

### HMM-SSF functions
All functions for data processing and model fitting are found within the `functions` folder. 

#### Model fitting functions
- *"fitHMMSSF.R"*: main function to fit HMM-SSF (calls all other sub-functions)
- *"nllk.R"*: function to calculate the negative log-likelihood with the forward algorithm 
- *"state_dens_Rcpp.cpp"*: C++ function to calculate the state-dependent densities
- *"sampling_dens.R"*: function to calculate the control sampling densities (currently only supports gamma-distributed steps and uniform sampling)
- *"hessian_CI.R"*: function to calculate confidence intervals from the hessian matrix returned by the numerical optimisation
- *"format_par.R"*: function to format parameters from a vector (format for optim) to a matrix (format of import/export)
-*"viterbi_decoding.R"*: function to use the viterbi algorithm for global decoding (adpated from moveHMM)
-*"local_decoding.R"*: function to calculate local state probabilities via the forward-backward algorithm
-*"predict_tpm.R"*: function to predict the time-varying transition probability matrix from the fitted model
-*"predict_delta.R"*: function to predict stationary state probabilities with standard errors (via the delta method)


#### Data processing functions
- *"sim_controls.R"*: function to simulate control locations from either a uniform spatial distribution or a gamma distribution of step lengths
- *"get_step_and_angle.R"*: function to calculate the step length and turning angle from each case location to the next case and all corresponding controls

#### other functions
- *"beta_to_mean.R"*: calculates the mean/sd of gamma distribution from estimated betas
-*"get_shape_scale.R"*: calculates the shape and scale of the gamma distribution from the estimated betas
- *"plot_step_angle.R"*: plots a histogram of data with the overlaid estimated state distributions


### Zebra example scripts
The folder `illustration.R` includes the code needed to reproduce the zebra example. This includes the following scripts:

- *"1_zebra_data_processing.R"*: script to process the zebra data (generating control locations, obtaining covariates, etc.)
- *"2_initial_par.R"*: script to set multiple sets of initial values (and save as `initial_par.RData`)
- *"3_zebra_fit.R"*: script to fit the model to the zebra data (using all sets of starting parameters and picking the best fitted model, and saving as `best_fit.RData`
- *"plot_hb.R"*: plot the movement data and habitat raster
- *"zebra_plot_decode.R"*: plot the viterbi sequence and local state probabilities
- *"zebra_plot_estimates.R"*: plot the estimated movement distributions and habitat selection estimates
- *"zebra_plot_tp.R"*: plot the estimated transition probabilities and stationary state probabilities based on the fitted model

## Data folder
This folder contains the data for the zebra example:

- *"zebra.csv"*: raw location data obtained from Michelot et al. (2020)
- *"zebra_controls.RData"*: location data with control locations (no covariates)
- *"zebra_processed.RData"*: zebra data with habitat covariates interpolated and formatted
- *"vegetation2.grd"* and *"vegetation2.gri"*: habitat raster (both files needed to load as raster)




