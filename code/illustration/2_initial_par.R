library(moveHMM)
library(fitdistrplus)
library(circular)

# load movement track
track <- read.csv("data/zebra.csv")
track_processed <- readRDS("data/zebra_processed.RData")

# fit regular HMM
track <- prepData(track, type = "UTM")
fitHMM <- fitHMM(track, 
                 nbStates = 2, 
                 stepPar0 = c(0.1, 1, 0.1, 1, 0.25, 0.1), 
                 anglePar0 = c(0, 0, 0.5, 5), 
                 stepDist = "gamma", 
                 angleDist = "vm")


# get starting values from HMM
shapeHMM <- fitHMM$mle$stepPar[c(1, 4)]^2 / fitHMM$mle$stepPar[c(2, 5)]^2
scaleHMM <- fitHMM$mle$stepPar[c(2, 5)]^2 / fitHMM$mle$stepPar[c(1, 4)]
b1HMM <- -1 / scaleHMM
b2HMM <- shapeHMM - 2
angleHMM <- fitHMM$mle$anglePar[c(2, 4)]

# get starting values based on pooled movement
obs <- subset(track_processed, obs == 1 & !is.na(step))
gamma_par <- fitdist(obs$step, distr = "gamma")
b1 <- as.vector(rep(-1 / gamma_par$estimate[2], times = 2))
b2 <- as.vector(rep(gamma_par$estimate[1] - 2, times = 2))
vm_par <- mle.vonmises(obs$angle)
angle <- rep(vm_par$kappa, times = 2)

# set initial values of alphas
alphas0 <- matrix(c(-2, -2, 
                    0, 0, 
                    0, 0), 
                  ncol = 2, 
                  byrow = T)

# set scenarios with movement based states

# no selection
par1 <- list(betas = matrix(c(b1HMM,
                              b2HMM,  
                              angleHMM, 
                              0, 0, 
                              0, 0, 
                              0, 0),
                            ncol = 2,
                            byrow = TRUE), 
             alphas = alphas0)

# pos, neg
par2 <- list(betas = matrix(c(b1HMM,
                              b2HMM,  
                              angleHMM, 
                              2, -2, 
                              2, -2, 
                              2, -2),
                            ncol = 2,
                            byrow = TRUE), 
             alphas = alphas0)

# neg, pos
par3 <- list(betas = matrix(c(b1HMM,
                              b2HMM,  
                              angleHMM, 
                              -2, 2, 
                              -2, 2, 
                              -2, 2),
                            ncol = 2,
                            byrow = TRUE), 
             alphas = alphas0)

# pos, pos
par4 <- list(betas = matrix(c(b1HMM,
                              b2HMM,  
                              angleHMM, 
                              2, 2, 
                              2, 2, 
                              2, 2),
                            ncol = 2,
                            byrow = TRUE), 
             alphas = alphas0)


# states driven by hs
par5 <- list(betas = matrix(c(b1,
                              b2,  
                              angle, 
                              2, -2, 
                              2, -2, 
                              2, -2),
                            ncol = 2,
                            byrow = TRUE), 
             alphas = alphas0)

par6 <- list(betas = matrix(c(b1,
                              b2,  
                              angle, 
                              -2, 2, 
                              -2, 2, 
                              -2, 2),
                            ncol = 2,
                            byrow = TRUE), 
             alphas = alphas0)


par <- list(par1, par2, par3, par4, par5, par6)


saveRDS(par, "code/illustration/initial_par.RData")



