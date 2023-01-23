source("code/functions/fitHMMSSF.R")
source("code/functions/nllk.R")
source("code/functions/format_par.R")
source("code/functions/hessian_CI.R")
source("code/functions/sampling_dens.R")
source("code/functions/viterbi_decoding.R")
source("code/functions/local_decoding.R")
library(Rcpp)
sourceCpp("code/functions/state_dens_Rcpp.cpp")
library(MASS)
library(ggplot2)
theme_set(theme_bw())

###################
## Model fitting ##
###################
# set model formulation
ssf_formula <- ~ step + log(step) + cos(angle) + veg
tpm_formula <- ~ cos(2*pi*tod/24) + sin(2*pi*tod/24)
n_states <- 2

# load starting values
initial_par <- readRDS("code/illustration/initial_par.RData")

#load movement data
data <- readRDS("data/zebra_processed.RData")
data$tod <- as.numeric(data$tod)

fits <- list()
times <- matrix(NA, ncol = 1, nrow = length(initial_par))

 for(i in 1:length(initial_par)) {
    start_t <- Sys.time()
    fits[[i]] <- fitHMMSSF(ssf_formula = ssf_formula, 
                           tpm_formula = tpm_formula,
                           data = data, 
                           par0 = initial_par[[i]], 
                           n_states = n_states, 
                           dist = "gamma", 
                           optim_opts = list(trace = 1,
                                             maxit = 1e4))
    end_t <- Sys.time()
    times[i,] <- as.numeric(end_t - start_t)
}


saveRDS(fits, "results/illustration/fits.RData")
saveRDS(times, "results/illustration/times.RData")

nllk_out <- c(fits[[1]]$nllk, fits[[2]]$nllk, fits[[3]]$nllk, 
          fits[[4]]$nllk, fits[[5]]$nllk, fits[[6]]$nllk)
which(min(nllk_out) == nllk_out)

saveRDS(fits[[3]], "results/illustration/best_fit.RData")

####################
## State decoding ##
####################
# Viterbi
decode <- viterbi_decoding(ssf_formula = ssf_formula, 
                           tpm_formula = tpm_formula, 
                           data = data, 
                           fit = fits[[3]], 
                           n_states = n_states)
saveRDS(decode, file = "results/illustration/zebra_viterbi.RData")

# local decoding
sp <- local_decoding(ssf_formula = ssf_formula, 
                     tpm_formula = tpm_formula, 
                     data = data, 
                     fit = fits[[3]], 
                     n_states = n_states)

saveRDS(sp, file = "results/illustration/zebra_sp.RData")

ggplot(subset(data, obs == 1), aes(x, y, col = sp[,2])) +
    geom_point(size = 0.5) +
    geom_path(size = 0.5) +
    scale_color_gradient(low = "firebrick", high = "royalblue") +
    coord_equal()







