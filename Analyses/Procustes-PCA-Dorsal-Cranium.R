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
library(patchwork)

# Inputting of coordinate data ----

# TPS file created from scaled landmarks measured in ImageJ
# Read in landmark coordinates
lands <- readland.tps(file = "Rawdata/dorsal.crania.tps", specID = "ID")
lands

# Input of only G. variegatus tps data
sundalands <- readland.tps(file = "Rawdata/sunda.dorsal.crania.tps", specID = "ID")
sundalands

# Plot raw y against raw x (Without procrustes - plotted points represent landmarks prior
## to rotation and scaling)
plot(lands[,2,]~lands[,1,])

# Plot G. Varigatus raw points
plot(sundalands[,2,]~sundalands[,1,])

# Normalisation of data through Generalised Procustes Analysis ----
## Procrustes analyses - uses reference specimen to line up all specimens 
gpa.lands <- gpagen(lands)

# Same with G. variegatus data
gpa.sunda <- gpagen(sundalands)

# Visualising procrustes analysis ----

# Black points are mean position of coordinates with grey points showing variation
plot(gpa.lands)
plot(gpa.sunda)

# Shows average landmark coordinates
gpa.lands
gpa.sunda

# Shows the scaled coordinates of all specimens
gpa.lands$coords
gpa.sunda$coords

# Expanded information of gpa.lands. (Str:structure)
str(gpa.lands)
str(gpa.sunda)

# Reading in of metadata ----
metadata <- read.csv("Rawdata/dermopteradata.csv")
glimpse(metadata)

# PCA Analysis ----

# PCA Analysis (Geomorph Principal and phylogenetically-aligned 
## components analysis of shape data)
pca.landmarks <- gm.prcomp(gpa.lands$coords)

# Ditto for G.variegatus data
sunda.pca.landmarks <- gm.prcomp(gpa.sunda$coords)

# How many PCs to include up to 95%? 
summary(pca.landmarks)
summary(sunda.pca.landmarks)

## Which landmarks are contributing towards each principle component
pca.landmarks$rotation
sunda.pca.landmarks$rotation

## Can be seen that 16 principal components explains for 95.7% of variation (within full
## data set including volans and variegatus)

# Merge PC scores and metadata ----

# Extract PC scores and make ID name from TPS into a taxon column for 
## combining with the metadata
pc_scores <- data.frame(specimenID = rownames(pca.landmarks$x), 
                        pca.landmarks$x)
pc_scores

# same with variegatus
sunda_scores <- data.frame(specimenID = rownames(sunda.pca.landmarks$x), 
                           sunda.pca.landmarks$x)
sunda_scores

# Make column names begin with PC not Comp
colnames(pc_scores) <- gsub("Comp", "PC", colnames(pc_scores))
colnames(sunda_scores) <- gsub("Comp", "PC", colnames(sunda_scores))

# Merge with metadata ----
pc_data <- full_join(pc_scores, metadata, by = c("specimenID" = "Dorsal.ID"))
sunda_pc_data <- full_join(sunda_scores, metadata, by = c("specimenID" = "Dorsal.ID"))

# Write PC scores to new csv files ----
write_csv(pc_data, file = "Rawdata/colugo-pca-data-dorsal.csv") 
write_csv(sunda_pc_data, file = "Rawdata/variegatus-pca-data-dorsal.csv") 

# Creation of wireframes ----

## Plot a reference shape
ref <- mshape(gpa.lands$coords)
sundaref <- mshape(gpa.sunda$coords)

# Plotting the minimum and maximum of x & y values 
plot(ref)
plotRefToTarget(ref, pca.landmarks$shapes$shapes.comp1$min)
plotRefToTarget(ref, pca.landmarks$shapes$shapes.comp1$max)
plotRefToTarget(ref, pca.landmarks$shapes$shapes.comp2$min)
plotRefToTarget(ref, pca.landmarks$shapes$shapes.comp2$max)
plotRefToTarget(ref, pca.landmarks$shapes$shapes.comp3$min)
plotRefToTarget(ref, pca.landmarks$shapes$shapes.comp3$max)

# Same for G.variegatus
plot(sundaref)
plotRefToTarget(sundaref, sunda.pca.landmarks$shapes$shapes.comp1$min)
plotRefToTarget(sundaref, sunda.pca.landmarks$shapes$shapes.comp1$max)
plotRefToTarget(sundaref, sunda.pca.landmarks$shapes$shapes.comp2$min)
plotRefToTarget(sundaref, sunda.pca.landmarks$shapes$shapes.comp2$max)
plotRefToTarget(sundaref, sunda.pca.landmarks$shapes$shapes.comp3$min)
plotRefToTarget(sundaref, sunda.pca.landmarks$shapes$shapes.comp3$max)

# Plotting PC Data. ----

## Plotting dorsal landmarks
derm.PC1 <- ggplot(pc_data, aes(x=PC1, y=PC2, shape=CurrentSp, colour=Region)) + 
  geom_point() + 
  theme_bw(base_size = 10) +
  labs(shape = "Current species")
derm.PC2 <- ggplot(pc_data, aes(x=PC1, y=PC3, shape=CurrentSp, colour=Region)) + 
  geom_point() + 
  theme_bw(base_size = 10) +
  labs(shape = "Current species")
ggplot(pc_data, aes(x=PC1, y=PC4, shape=CurrentSp, colour=Region)) + 
  geom_point() + 
  theme_bw(base_size = 10) +
  labs(shape = "Current species")
ggplot(pc_data, aes(x=PC2, y=PC3, shape=CurrentSp, colour=Region)) + 
  geom_point() + 
  theme_bw(base_size = 10) +
  labs(shape = "Current species")
ggplot(pc_data, aes(x=PC2, y=PC4, shape=CurrentSp, colour=Region)) + 
  geom_point() + 
  theme_bw(base_size = 10) +
  labs(shape = "Current species")
ggplot(pc_data, aes(x=PC3, y=PC4, shape=CurrentSp, colour=Region)) + 
  geom_point() + 
  theme_bw(base_size = 10) +
  labs(shape = "Current species")

derm.PC1 + derm.PC2
# Plots for G.variegatus principal component analysis 
ggplot(sunda_pc_data, aes(x=PC1, y=PC2, colour=Region)) + 
  geom_point() + 
  theme_bw(base_size = 10)
ggplot(sunda_pc_data, aes(x=PC1, y=PC3, colour=Region)) + 
  geom_point() + 
  theme_bw(base_size = 10)
ggplot(sunda_pc_data, aes(x=PC1, y=PC4, colour=Region)) + 
  geom_point() + 
  theme_bw(base_size = 10)
ggplot(sunda_pc_data, aes(x=PC2, y=PC3, colour=Region)) + 
  geom_point() + 
  theme_bw(base_size = 10)
ggplot(sunda_pc_data, aes(x=PC2, y=PC4, colour=Region)) + 
  geom_point() + 
  theme_bw(base_size = 10)
ggplot(sunda_pc_data, aes(x=PC3, y=PC4, colour=Region)) + 
  geom_point() + 
  theme_bw(base_size = 10)

# Misc Code ----
  geom_label(aes(label=SpecimenID))

print(pc_data)
options(max.print = 20)

## Aesthetics(X var, Y Var, shape based variable, colour based variable)
## geom_label(aes(label=XX))
