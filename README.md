# HMM-SSF
Fit state-switching step selection functions to animal movement data

Natasha Klappstein

Working directory note: All scripts were written with the root directory of the repository as the working directory. Sourcing the necessary files must be done from this directory.  

This repository is composed of 2 main folders (`code`, `data`), described here. The other folders (`results`, `writing`) just contain exported data/figures and the dissertation files. 

## Code folder
There are 3 main folders within the code folder. The main functions for data processing and model fitting are found in `functions`. Code to run simulations and fit the model to the zebra data are found in `simulations` and `illustration`, respectively.

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

#### Data processing and simulation functions
- *"sim_controls.R"*: function to simulate control locations from either a uniform spatial distribution or a gamma distribution of step lengths
-*"get_step_and_angle.R"*: function to calculate the step length and turning angle from each case location to the next case and all corresponding controls
- *"simRaster.R"*: function to simulate a covariate raster (adapted from the localGibbs package)
- *"simHMMSSF.R"*: function to simulate a movement track from the HMM-SSF
- *"cov_df.R"*: function to extract covariates from raster list in simulation function

#### other functions
- *"beta_to_mean.R"*: calculates the mean/sd of gamma distribution from estimated betas
-*"get_shape_scale.R"*: calculates the shape and scale of the gamma distribution from the estimated betas
- *"plot_step_angle.R"*: plots a histogram of data with the overlaid estimated state distributions

### Simulation scripts
Within the folder `simulations`, these scripts can reproduce the simulations in Chapter 3. 
- `simulations/importance_sampling_IS/IS_sim.R` to reproduce the importance sampling simulation, run the following file. All other scripts in this folder are to plot the results.
- `simulations/trans_prob_TPM/TPM_sim.R` to reproduce the transition probability simulation, run the following file. All other scripts in this folder are to plot the results.
- the scripts titled *"3_state.R"* and *"overlap_sim.R"* run a 3 state model and a model with more state overlap (not in thesis). 

### Zebra example scripts
The folder `illustration.R` includes the code needed to reproduce the zebra example. This includes the following scripts:

- *"starting_val_zebra.R": script to process the zebra data and fit model with several starting values and pick the best based on the lowest negative log likelihood
- *"zebra_fit.R"*: script to fit the model to the zebra data (using best starting values)
- *"zebra_plot.R"*: plot the zebra results
- *"plot_hb.R"*: plot the movement data and habitat raster

## Data folder
This folder contains the data for the zebra example:

- *"zebra.csv"*: raw location data obtained from Michelot et al. (2020)
- *"zebra_processed.RData"*: zebra data with habitat covariates interpolated and formatted
- *"vegetation2.grd"*: habitat raster




