library(ggplot2)
library(cowplot)
theme_set(theme_light(13))
source("code/functions/predict_tpm.R")
source("code/functions/predict_delta.R")

# load results
fit <- readRDS("results/illustration/best_fit.RData")
new_data <- data.frame(tod = seq(0, 24, length = 100))
tpm_formula <- ~ cos(2*pi*tod/24) + sin(2*pi*tod/24)

##############################
## Transition probabilities ##
##############################
tpm <- predict_tpm(tpm_formula = tpm_formula, 
                   new_data = new_data, 
                   fit = fit, 
                   n_states = 2, 
                   return_CI = TRUE)

df <- as.data.frame.table(tpm$mle)
colnames(df) <- c("from", "to", "tod", "value")
df$tod <- rep(new_data$tod, each = 4)
df$low <- as.data.frame.table(tpm$lower)[,4]
df$upp <- as.data.frame.table(tpm$upper)[,4]
df$from <- ifelse(df$from == 1, "encamped", "exploratory")
df$to <- ifelse(df$to == 1, "encamped", "exploratory")
df$prob <- paste0(df$from, " - ", df$to)

p1 <- ggplot(subset(df, from != to), aes(tod, value, linetype = prob)) + 
  geom_ribbon(aes(ymin = low, ymax = upp), alpha = 0.3) +
  geom_line() +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "time of day", y = "transition probability") +
  scale_x_continuous(breaks = seq(0, 24, by = 4)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.5), labels = c("0", "0.5", "1")) +
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 10.25))
p1 
####################################
## Stationary state probabilities ##
####################################
delta <- predict_delta(tpm_formula = tpm_formula, 
                       new_data = new_data, 
                       fit = fit, 
                       n_states = 2, 
                       return_CI = TRUE)

df <- as.data.frame.table(delta$mle)
colnames(df) <- c("tod", "state", "value")
df$tod <- new_data$tod
df$low <- as.data.frame.table(delta$lower)[,3]
df$upp <- as.data.frame.table(delta$upper)[,3]
df$state <- ifelse(df$state == 1, "encamped", "exploratory")
pal <- c("#0072B2", "#D55E00")


p2 <- ggplot(df, aes(tod, value, group = state, col = factor(state), fill = factor(state))) +
    geom_ribbon(aes(ymin = low, ymax = upp), col = NA, alpha = 0.3) +
    geom_line() +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "time of day", y = "stationary probability") +
    scale_color_manual(values = pal, name = NULL) +
    scale_fill_manual(values = pal, name = NULL) +
    scale_x_continuous(breaks = seq(0, 24, by = 4)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.5), labels = c("0", "0.5", "1")) +
    # theme(legend.position = c(0.1, 0.8))
    theme(legend.position = "bottom", 
          legend.text = element_text(size = 10.5))
p2

pdf(file = "writing/figures/zebra_tp.pdf", height = 4, width = 9)
cowplot::plot_grid(p1, p2, nrow = 1, labels = "auto")
dev.off()






