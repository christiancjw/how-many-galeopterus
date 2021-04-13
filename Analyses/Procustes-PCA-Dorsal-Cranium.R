# SETUP ---- 
## Geometric morphometric analyses of Galeopterus Skulls 
## Set working directory to folder for project, or open project at new directory. 
## Version control repository established with GitHub

# Load libraries:
## tidyverse (Use of data viewing functions), geomorph (For landmark analyses),
## ggplot2 (For creating figures), dplyr (for data manipulation)
library(geomorph)
library(tidyverse)
library(ggplot2)
library(dplyr)

# Inputting of coordinate data ----
# TPS file created from scaled landmarks measured in ImageJ
# Read in landmark coordinates
lands <- readland.tps(file = "Rawdata/dorsal.crania.tps", specID = "ID")
lands

# Plot y against x (Without procrustes - plotted points represent landmarks prior
## to rotation and scaling)
plot(lands[,2,]~lands[,1,])

# Normalisation of data through Generalised Procustes Analysis ----
## Procrustes analyses - uses reference specimen to line up all specimens 
gpa.lands <- gpagen(lands)

# Visualising procrustes analysis
# Black points are mean position of coordinates with grey points showing variation
plot(gpa.lands)
# Shows average landmark coordinates
gpa.lands
# Shows the scaled coordinates of all specimens
gpa.lands$coords
# Expanded information of gpa.lands. (Str:structure)
str(gpa.lands)

# Reading in of metadata ----
metadata <- read.csv("Rawdata/dermopteradata.csv")
glimpse(metadata)

## PCA Analysis ----
pca.landmarks <- gm.prcomp(gpa.lands$coords)
pca.landmarks

# How many PCs to include up to 95%? 
summary(pca.landmarks)
## Can be seen that 16 principal components explains for 95.7% of variation

# Merge PC scores and metadata ----
# Extract PC scores and make ID name from TPS into a taxon column for 
## combining with the metadata
pc_scores <- data.frame(specimenID = rownames(pca.landmarks$x), 
                        pca.landmarks$x)
pc_scores

# Make column names begin with PC not Comp
colnames(pc_scores) <- gsub("Comp", "PC", colnames(pc_scores))

# Merge with metadata ----
pc_data <- full_join(pc_scores, metadata, by = c("specimenID" = "DorsalID"))

# Making a new dataset of just galeopterus (excluding volans) ----
# The ! is used to reverse arguements
g_data <- pc_data
gv_data <- filter(g_data, !Region == "Philippines")

## Plotting dorsal landmarks without volans. Only  principle component 1
ggplot(gv_data, aes(x=PC1, y=PC2, colour=Region)) + 
  geom_point()

# Write to file ----
write_csv(pc_data, path = "~/Users/christianching/Desktop/
          Galeopterus/RProject/galeopterus.skulls/Data/colugo-pca-data-dorsal.csv") 

## Which landmarks are contributing towards each principle component
pca.landmarks$rotation

## Plot a reference shape
ref <- mshape(gpa.lands$coords)

#plotting the minimum and maximum of x & y values 
plot(ref)
plotRefToTarget(ref, pca.landmarks$shapes$shapes.comp1$min)
plotRefToTarget(ref, pca.landmarks$shapes$shapes.comp1$max)
plotRefToTarget(ref, pca.landmarks$shapes$shapes.comp2$min)
plotRefToTarget(ref, pca.landmarks$shapes$shapes.comp2$max)

# Plotting PC Data. ----
## Aesthetics(X var, Y Var, shape based variable, colour based variable)
## geom_label(aes(label=XX))
ggplot(pc_data, aes(x=PC1, y=PC2, shape=CurrentSp, colour=Region)) + 
  geom_point()
  
  geom_label(aes(label=SpecimenID))

print(pc_data)
options(max.print = 100000)

