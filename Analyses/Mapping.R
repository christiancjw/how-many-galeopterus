# Assuming your original data is called cologu. Note you don't need the PCs or anything, just the data for how you want to split/colour the points and the lat longs
## Load libraries
library(tidyverse)
library(sf)
library(rgeos)
library(rgdal)
library(Hmisc)
install.packages("sf")
install.packages("rgeos")
install.packages("rgdal")
install.packages("Hmsic")

# Helper functions for plotting
remove_y <- 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
remove_x <-   
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


##----------------------------------------------

# First convert the lat long into "points" data for plotting
colugo <-
  pc_data %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)


##----------------------------------------------
# Make a base map of the land
##----------------------------------------------
baseMap <- 
  rnaturalearth::ne_countries(returnclass = 'sf') %>%
  st_union()

##-------------------------------------------------
## Choose coordinates to limit to Asia (may need to modify to zoom in or out more)
##-------------------------------------------------
# puts the coords into the order expected down in ggmap coords
asia_bbox <- c(68, -10, 
               130, 40)
xlim_as <- c(asia_bbox[1], asia_bbox[3])
ylim_as <- c(asia_bbox[2], asia_bbox[4])
##-------------------------------------------------
# Make map - only the first bits are special
# changing colours etc is same as normal ggplot code
##-------------------------------------------------
ggplot(baseMap) +
  geom_sf() +
  # Add points
  geom_sf(alpha = 0.5, aes(colour = Region, fill = Region),
          data = specs, show.legend = TRUE, size = 0.5) +
  # restrict map to just Asia
  coord_sf(xlim = xlim_as, ylim = ylim_as, expand = TRUE) +
  theme_bw() +
  remove_x +
  remove_y
# Save plots  
ggsave("figures/colugo-maps.png")