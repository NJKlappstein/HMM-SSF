library(raster)
source("code/functions/sim_controls.R")
source("code/functions/get_step_and_angle.R")

# load movement data
load("data/zebHMM.RData")
track <- subset(zebHMM, id == "Apero")
track <- data.frame(ID = 1, 
                    x = track$x, 
                    y = track$y,
                    time = track$expectTime)
track$x <- track$x / 1000
track$y <- track$y / 1000
write.csv(track, "data/zebra.csv", row.names = F)

# load habitat data
hb <- raster("data/vegetation2.grd")

set.seed(250)
data <- sim_controls(obs = track, 
                     n_controls = 25, 
                     step_dist = "gamma")
saveRDS(data, "data/zebra_controls.RData")

# get covariates
data <- get_step_and_angle(data)
notNA <-  which(!is.na(data$x))
data$veg[notNA] <- extract(hb, data[notNA, c("x", "y")], method = "simple")
data$veg <- factor(data$veg)
levels(data$veg) <- c("grassland", "bushed grassland", "bushland", "woodland")
data$tod <- as.numeric(format(data$time, "%H")) + 
  as.numeric(format(data$time, "%M"))/60
data$step[which(data$step == 0)] <- data$step[which(data$step == 0)] +
  runif(length(data$step[which(data$step == 0)]), 0, min(data$step[which(data$step > 0)]))

# write data for fitting later
saveRDS(data, file = "data/zebra_processed.RData")


 


