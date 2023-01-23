source("code/functions/sampling_dens.R")
source("code/functions/viterbi_decoding.R")
source("code/functions/local_decoding.R")
library(Rcpp)
sourceCpp("code/functions/state_dens_Rcpp.cpp")
library(ggplot2)
theme_set(theme_light(13))
library(raster)
library(wesanderson)

# model formulation
ssf_formula <- ~ step + log(step) + cos(angle) + veg
tpm_formula <- ~ cos(2*pi*tod/24) + sin(2*pi*tod/24)
n_states <- 2

# load zebra data and results
data <- readRDS("data/zebra_processed.RData")
fit <- readRDS("best_fit.RData")

# load habitat raster and format for plot
hb <- raster("data/vegetation2.grd")
pal <- c("#9BA6A0", "#5F7A6B", "#365543", "#1C332D")
lab <- c("grassland", "bushed\ngrassland", "bushland", "woodland")
covmap <- data.frame(coordinates(hb), val = values(hb))

# Viterbi
global <- viterbi_decoding(ssf_formula = ssf_formula, 
                           tpm_formula = tpm_formula, 
                           data = data, 
                           fit = fit, 
                           n_states = n_states)

# local decoding
local <- local_decoding(ssf_formula = ssf_formula, 
                        tpm_formula = tpm_formula, 
                        data = data, 
                        fit = fit, 
                        n_states = n_states)

obs <- subset(data, obs == 1)
obs$global <- paste0("state ", global)
obs$sp2 <- local[,2]
obs$global_name <- ifelse(obs$global == "state 1", "encamped", "exploratory")

##################################################
#### separate plots for global and local #########
##################################################
pal_points <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

# global
p1 <- ggplot(covmap, aes(x, y)) + 
    geom_raster(aes(fill = factor(val)), alpha = 0.6) +
    coord_equal() + 
    labs(x = "Easting (km)", y = "Northing (km)") +
    scale_fill_manual(values = pal, guide = "none") +
    geom_point(data = obs, aes(color = global),
               size = 0.8, alpha = 0.6) +
    scale_color_manual(values = pal_points[c(5,6)], guide = "none") +
    facet_wrap("global_name") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(strip.background = element_blank(), 
          strip.text = element_text(colour = "black", size = 13))

# local
p2 <- ggplot(covmap, aes(x, y)) + 
    geom_raster(aes(fill = factor(val)), alpha = 0.6) +
    coord_equal() + 
    labs(x = "Easting (km)", y = "Northing (km)") +
    scale_fill_manual(values = pal, guide = "none") +
    scale_color_gradient(low = pal_points[5], high = pal_points[6], 
                         name = expression(Pr(S[t] == 2))) +
    geom_point(aes(col = sp2), data = obs, size = 0.8, alpha = 0.5) +
    theme(legend.key.size = unit(1.2, "lines")) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0))


################################################
#### Comparison of Viterbi and state probs #####
################################################
obs$global_num <- ifelse(obs$global == "state 1", 1, 2)
p3 <- ggplot(obs[100:250,], aes(time, sp2)) +
    geom_point(aes(y = global_num - 1), color = "grey65") +
    geom_line(aes(y = global_num - 1), color = "grey65") +
    geom_line(color = "black") +
    xlab("") +
    scale_y_continuous(
        breaks = c(0, 0.5, 1),
        expression(Pr(S[t] == 2)), 
        sec.axis = sec_axis(~ . + 1, name = "Viterbi state", breaks = c(1, 2))) +
    theme(axis.title.y.right = element_text(colour = "grey65"), 
          axis.title.y.left = element_text(color = "black"))

pdf("documents/figures/vit_sp_seq.pdf", width = 9, height = 3)
p3
dev.off()

# to create plots to go in Klappstein et al. 
pdf("writing/figures/vit_map.pdf", width = 6, height = 3)
cowplot::plot_grid(p1, labels = "a")
dev.off()
pdf("writing/figures/vit_sp_seq.pdf", width = 9, height = 3)
cowplot::plot_grid(p3, labels = "b")
dev.off()



