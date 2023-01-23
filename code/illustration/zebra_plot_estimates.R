source("code/functions/plot_step_angle.R")
source("code/functions/beta_to_mean.R")
library(dplyr)
library(CircStats)
pal <- c("#0072B2", "#D55E00")

# load results
fit <- readRDS("results/illustration/best_fit.RData")
decode <- readRDS("results/illustration/zebra_viterbi.RData")
data <- readRDS("data/zebra_processed.RData")

n_states <- 2
p_list <- plot_step_angle(fit, decode, data, n_states, 
                          pal = pal[c(1,2)], as_list = TRUE)
p1 <- p_list$step + theme_light(13) + 
  scale_color_manual("state", values = pal[c(1,2)], 
                     labels = c("encamped", "exploratory")) +  
  theme(legend.position = c(0.7, 0.80), 
        legend.title = element_blank())
p2 <- p_list$angle + theme_light(13) + theme(legend.position = "none")

ssf_df <- data.frame(state = rep(c(1, 2), each = 3), 
                     covariate = rep(c("bushed \n grassland", 
                                       "bushland", 
                                       "woodland"),
                                     times = 2), 
                     estimate = c(fit$betas$estimate[c(4:6, 10:12)]), 
                     upper = c(fit$betas$upper[c(4:6, 10:12)]),
                     lower = c(fit$betas$lower[c(4:6, 10:12)]))

pd <- position_dodge(0.25)
p3 <- ggplot(ssf_df, aes(y = estimate, 
                         x = covariate, 
                         group = state, 
                         color = as.factor(state))) + 
  geom_point(position = pd, size = 2) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                position = pd,
                size = 0.5,
                width = 0.15) + 
  xlab("habitat type") + 
  scale_color_manual("state", values = pal[c(1:2)]) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  theme_light(13) + theme(legend.position = "none")

p4 <- plot_grid(p1, p2, p3, labels = "auto", nrow = 1)

pdf("writing/figures/zebra_estimates.pdf", width = 10, height = 2.8)
plot(p4)
dev.off()









