library(raster)
library(ggplot2)
library(dplyr)

# load movement data
track <- read.csv("data/zebra.csv")
names(track) <- c("ID", "x", "y", "time")
hb <- raster("data/vegetation2.grd")

# set color palettes and labels
pal <- c("#9BA6A0", "#5F7A6B", "#365543", "#1C332D")
lab <- c("grassland", "bushed\ngrassland", "bushland", "woodland")
covmap <- data.frame(coordinates(hb), val = values(hb))

# plot habitat map with track
# pdf(file = "writing/figures/hb_plot.pdf", width = 5.25, height = 3.75)
ggplot(covmap, aes(x, y)) + 
  geom_raster(aes(fill = factor(val))) +
  coord_equal() + 
  xlab("Easting (km)") + 
  ylab("Northing (km)") +
  scale_fill_manual(values = pal, name = NULL, labels = lab) +
  geom_path(aes(x, y), track, size = 0.35, color = "burlywood1", alpha = 0.6) +
  geom_point(aes(x, y), track, size = 0.2, color = "burlywood1", alpha = 0.6) +
  theme(legend.key.size = unit(1.4, "lines")) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))
# dev.off()


