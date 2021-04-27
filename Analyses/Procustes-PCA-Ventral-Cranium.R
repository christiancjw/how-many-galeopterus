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
library(grid)
library(gridExtra)
library(cowplot)

# Inputting of coordinate data ----

# TPS file created from scaled landmarks measured in ImageJ
# Read in landmark coordinates
ventral.lands <- readland.tps(file = "Rawdata/ventral.crania.tps", specID = "ID")
ventral.lands

# Input of only G. variegatus tps data
sunda.ventral.lands <- readland.tps(file = "Rawdata/sunda.ventral.crania.tps", specID = "ID")
sunda.ventral.lands

# Plot raw y against raw x (Without procrustes - plotted points represent landmarks prior
## to rotation and scaling)
plot(ventral.lands[,2,]~ventral.lands[,1,])

# Plot G. Varigatus raw points
plot(sunda.ventral.lands[,2,]~sunda.ventral.lands[,1,])

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
sunda_pc_data_unfixed <- full_join(sunda_scores, metadata, by = c("specimenID" = "Dorsal.ID"))

sunda_pc_data <- sunda_pc_data_unfixed %>% 
  filter(!Region == "Philippines")

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

derm.plot.1 <- ggplot(pc_data, aes(x=PC1, y=PC2, shape=CurrentSp, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  theme_bw(base_size = 7) +
  labs(shape = "Current species") + 
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8)

derm.plot.2 <- ggplot(pc_data, aes(x=PC1, y=PC3, shape=CurrentSp, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  theme_bw(base_size = 7) +
  labs(shape = "Current species") +
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8)

derm.plot.3 <- ggplot(pc_data, aes(x=PC1, y=PC4, shape=CurrentSp, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  theme_bw(base_size = 7) +
  labs(shape = "Current species") +
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8)

derm.plot.4 <- ggplot(pc_data, aes(x=PC2, y=PC3, shape=CurrentSp, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  theme_bw(base_size = 7) +
  labs(shape = "Current species") +
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8)

derm.plot.5 <- ggplot(pc_data, aes(x=PC2, y=PC4, shape=CurrentSp, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  theme_bw(base_size = 7) +
  labs(shape = "Current species") +
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8)

derm.plot.6 <- ggplot(pc_data, aes(x=PC3, y=PC4, shape=CurrentSp, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  theme_bw(base_size = 7) +
  labs(shape = "Current species") + 
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8)

# Isolation of legend through creation of new plot and cowplot pkg
derml <- ggplot(pc_data, aes(x=PC3, y=PC4, shape=CurrentSp, colour=Region)) + 
  geom_point() + 
  theme_bw(base_size = 10) +
  labs(shape = "Current species") 

dermlegend <- cowplot::get_legend(derml)
grid.newpage()
grid.draw(legend)

plot(legend)

layout <- "
AABBCC###
AABBCC###
AABBCCGGG
DDEEFFGGG
DDEEFF###
DDEEFF###
"

# Using patchwork to make composite figure
derm.pc.plots <- derm.plot.1 + derm.plot.2 + derm.plot.3 + derm.plot.4 + derm.plot.5 + 
  derm.plot.6 + dermlegend + plot_layout(design = layout)

derm.pc.plots
ggsave("dermpcplots.png")

# Plots for G.variegatus principal component analysis ----

sunda.plot.1 <- ggplot(sunda_pc_data, aes(x=PC1, y=PC2, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  theme_bw(base_size = 7) +
theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8)

sunda.plot.2 <- ggplot(sunda_pc_data, aes(x=PC1, y=PC3, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  theme_bw(base_size = 7) +
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8)

sunda.plot.3 <- ggplot(sunda_pc_data, aes(x=PC1, y=PC4, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  theme_bw(base_size = 7) +
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8)

sunda.plot.4 <- ggplot(sunda_pc_data, aes(x=PC2, y=PC3, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  theme_bw(base_size = 7) +
  scale_x_continuous(breaks = c(-0.02, 0, 0.02, 0.04)) +
  theme(legend.position = "NONE")+
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8)

sunda.plot.5 <- ggplot(sunda_pc_data, aes(x=PC2, y=PC4, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  theme_bw(base_size = 7) +
  scale_x_continuous(breaks = c(-0.02, 0, 0.02, 0.04)) +
  theme(legend.position = "NONE")+
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8)

sunda.plot.6 <- ggplot(sunda_pc_data, aes(x=PC3, y=PC4, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  theme_bw(base_size = 7) +
  scale_x_continuous(breaks = c(-0.02, 0, 0.02)) +
  theme(legend.position = "NONE")+
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8)

sundal <- ggplot(sunda_pc_data, aes(x=PC3, y=PC4, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  theme_bw(base_size = 10) 

sundalegend <- cowplot::get_legend(sundal)
grid.newpage()
grid.draw(sundalegend)
plot(sundalegend)

sunda.pc.plots <- sunda.plot.1 + sunda.plot.2 + sunda.plot.3 + sunda.plot.4 + sunda.plot.5 +
  sunda.plot.6 + sundalegend + plot_layout(design = layout)

sunda.pc.plots
ggsave("sundapcplots.png")

# Composite PCA plots ----
composite.pc.plot.data <- pc_data %>%
  dplyr::select(CurrentSp, Region, PC1:PC19) %>%
  pivot_longer(PC1:PC18, "PC", "value") %>%
  mutate(PC = factor(PC, levels = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6",
                                    "PC7", "PC8", "PC9", "PC10", "PC11", "PC12",
                                    "PC13", "PC14", "PC15", "PC16", "PC17", "PC18"))) 

composite.pc.plot.data


megaplot_plot <- 
  ggplot(composite.pc.plot.data, aes(x = PC, y = value, colour = Region, shape = CurrentSp)) +
  scale_fill_manual(values=c("#083D77", "#65532F", "#63B0CD", "#8AB17D", "#E9C46A", "#E94F37", "#0CF574"))+
  geom_boxplot() +
  geom_jitter(alpha = 0.7) +
  theme_bw(base_size = 14) +
  scale_color_manual(values = c("#083D77", "#65532F", "#63B0CD", "#8AB17D", "#E9C46A", "#E94F37", "#0CF574")) +
  ylim(-0.2,0.2) +
  coord_flip() +
  xlab("PC axis") +
  ylab("PC score") +
  labs(shape = "Current species") +
  xlim(-0.1,0.1)

xlim(-0.1,0.1)
  
megaplot_plot

# Misc Code ----
  geom_label(aes(label=SpecimenID))

print(pc_data)
options(max.print = 20)

## Aesthetics(X var, Y Var, shape based variable, colour based variable)
## geom_label(aes(label=XX))
