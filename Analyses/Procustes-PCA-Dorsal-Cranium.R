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
dorsal.lands <- readland.tps(file = "Rawdata/dorsal.crania.tps", specID = "ID")
dorsal.lands

# Input of only G. variegatus tps data
dorsal.gv.lands <- readland.tps(file = "Rawdata/sunda.dorsal.crania.tps", specID = "ID")
dorsal.gv.lands

# Plot raw y against raw x (Without procrustes - plotted points represent landmarks prior
## to rotation and scaling)
plot(dorsal.lands[,2,]~dorsal.lands[,1,])

# Plot G. Varigatus raw points
plot(dorsal.gv.lands[,2,]~dorsal.gv.lands[,1,])

# Normalisation of data through Generalised Procustes Analysis ----
## Procrustes analyses - uses reference specimen to line up all specimens 
gpa.dorsal.lands <- gpagen(dorsal.lands)

# Same with G. variegatus data
gpa.dorsal.gv.lands <- gpagen(dorsal.gv.lands)

# Visualising procrustes analysis ----

# Black points are mean position of coordinates with grey points showing variation
plot(gpa.dorsal.lands)
plot(gpa.dorsal.gv.lands)

# Shows average landmark coordinates
gpa.dorsal.lands
gpa.dorsal.gv.lands

# Shows the scaled coordinates of all specimens
gpa.dorsal.lands$coords
gpa.dorsal.gv.lands$coords

# Expanded information of gpa.lands. (Str:structure)
str(gpa.dorsal.lands)
str(gpa.dorsal.gv.lands)

# Reading in of metadata ----
metadata <- read.csv("Rawdata/dermopteradata.csv")
glimpse(metadata)

# PCA Analysis ----

# PCA Analysis (Geomorph Principal and phylogenetically-aligned 
## components analysis of shape data)
dorsal.pca.landmarks <- gm.prcomp(gpa.dorsal.lands$coords)

# Ditto for G.variegatus data
dorsal.pca.gv.landmarks <- gm.prcomp(gpa.dorsal.gv.lands$coords)

# How many PCs to include up to 95%? 
summary(dorsal.pca.landmarks)
summary(dorsal.pca.gv.landmarks)

## Which landmarks are contributing towards each principle component
dorsal.pca.landmarks$rotation
dorsal.pca.gv.landmarks$rotation

## Can be seen that 16 principal components explains for 95.7% of variation (within full
## data set including volans and variegatus)

# Read in metadata and prep PC data for merging ----

# Extract PC scores and make ID name from TPS into a taxon column for 
## combining with the metadata
dorsal.pc.scores <- data.frame(specimenID = rownames(dorsal.pca.landmarks$x), 
                        dorsal.pca.landmarks$x)
dorsal.pc.scores

# same with variegatus
dorsal.gv.pc.scores <- data.frame(specimenID = rownames(dorsal.pca.gv.landmarks$x), 
                                  dorsal.pca.gv.landmarks$x)
dorsal.gv.pc.scores

# Make column names begin with PC not Comp
colnames(dorsal.pc.scores) <- gsub("Comp", "PC", colnames(dorsal.pc.scores))
colnames(dorsal.gv.pc.scores) <- gsub("Comp", "PC", colnames(dorsal.gv.pc.scores))

# Merge with metadata ----
dorsal.pc.data <- full_join(dorsal.pc.scores, metadata, by = c("specimenID" = "Dorsal.ID"))
dorsal.gv.pc.data.unfixed <- full_join(dorsal.gv.pc.scores, metadata, by = c("specimenID" = "Dorsal.ID"))

dorsal.gv.pc.data <- dorsal.gv.pc.data.unfixed %>% 
  filter(!Region == "Philippines")

# Write PC scores to new csv files ----
write_csv(dorsal.pc.data, file = "Rawdata/dorsal_dermoptera_pca_data.csv") 
write_csv(dorsal.gv.pc.data, file = "Rawdata/dorsal_variegatus_pca_data.csv") 

# Creation of wireframes ----

## Plot a reference shape
dorsal.ref <- mshape(gpa.dorsal.lands$coords)
dorsal.gv.ref <- mshape(gpa.dorsal.gv.lands$coords)

# Plotting the minimum and maximum of x & y values 
plot(dorsal.ref)
plotRefToTarget(dorsal.ref, dorsal.pca.landmarks$shapes$shapes.comp1$min)
plotRefToTarget(dorsal.ref, dorsal.pca.landmarks$shapes$shapes.comp1$max)
plotRefToTarget(dorsal.ref, dorsal.pca.landmarks$shapes$shapes.comp2$min)
plotRefToTarget(dorsal.ref, dorsal.pca.landmarks$shapes$shapes.comp2$max)
plotRefToTarget(dorsal.ref, dorsal.pca.landmarks$shapes$shapes.comp3$min)
plotRefToTarget(dorsal.ref, dorsal.pca.landmarks$shapes$shapes.comp3$max)

# Same for G.variegatus
plot(dorsal.gv.ref)
plotRefToTarget(dorsal.gv.ref, dorsal.pca.gv.landmarks$shapes$shapes.comp1$min)
plotRefToTarget(dorsal.gv.ref, dorsal.pca.gv.landmarks$shapes$shapes.comp1$max)
plotRefToTarget(dorsal.gv.ref, dorsal.pca.gv.landmarks$shapes$shapes.comp2$min)
plotRefToTarget(dorsal.gv.ref, dorsal.pca.gv.landmarks$shapes$shapes.comp2$max)
plotRefToTarget(dorsal.gv.ref, dorsal.pca.gv.landmarks$shapes$shapes.comp3$min)
plotRefToTarget(dorsal.gv.ref, dorsal.pca.gv.landmarks$shapes$shapes.comp3$max)

# Plots for dorsal principal component analysis ----

## Plotting dorsal landmarks

dorsal.plot1 <- ggplot(dorsal.pc.data, aes(x=PC1, y=PC2, shape=CurrentSp, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  theme_bw(base_size = 7) +
  labs(shape = "Current species") + 
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.05,0.10), ylim=c(-0.05,0.05)) 

dorsal.plot2 <- ggplot(dorsal.pc.data, aes(x=PC1, y=PC3, shape=CurrentSp, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  theme_bw(base_size = 7) +
  labs(shape = "Current species") +
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.05,0.10), ylim=c(-0.05,0.05))

dorsal.plot3 <- ggplot(dorsal.pc.data, aes(x=PC1, y=PC4, shape=CurrentSp, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  theme_bw(base_size = 7) +
  labs(shape = "Current species") +
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.05,0.10), ylim=c(-0.05,0.05))

dorsal.plot4 <- ggplot(dorsal.pc.data, aes(x=PC2, y=PC3, shape=CurrentSp, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  theme_bw(base_size = 7) +
  labs(shape = "Current species") +
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.05,0.05), ylim=c(-0.04,0.04))

dorsal.plot5 <- ggplot(dorsal.pc.data, aes(x=PC2, y=PC4, shape=CurrentSp, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  theme_bw(base_size = 7) +
  labs(shape = "Current species") +
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.05,0.05), ylim=c(-0.04,0.04))

dorsal.plot6 <- ggplot(dorsal.pc.data, aes(x=PC3, y=PC4, shape=CurrentSp, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  theme_bw(base_size = 7) +
  labs(shape = "Current species") + 
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.05,0.05), ylim=c(-0.04,0.04))

# Isolation of legend through creation of new plot and cowplot pkg
dorsal.legend.ref <- ggplot(dorsal.pc.data, aes(x=PC3, y=PC4, shape=CurrentSp, colour=Region)) + 
  geom_point() + 
  theme_bw(base_size = 10) +
  labs(shape = "Current species") 

dorsal.legend <- cowplot::get_legend(dorsal.legend.ref)
grid.newpage()
grid.draw(dorsal.legend)

plot(dorsal.legend)

layout <- "
AABBCC###
AABBCC###
AABBCCGGG
DDEEFFGGG
DDEEFF###
DDEEFF###
"

# Using patchwork to make composite figure
dorsal.pc.plots <- dorsal.plot1 + dorsal.plot2 + dorsal.plot3 + dorsal.plot4 + dorsal.plot5 + 
  dorsal.plot6 + dorsal.legend + plot_layout(design = layout)

dorsal.pc.plots
ggsave("dorsal_dermoptera_pc_plots.png")

# Plots for dorsal G.variegatus principal component analysis ----

gv.dorsal.plot1 <- ggplot(dorsal.gv.pc.data, aes(x=PC1, y=PC2, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  theme_bw(base_size = 7) +
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.04,0.04), ylim=c(-0.05,0.05))

gv.dorsal.plot2 <- ggplot(dorsal.gv.pc.data, aes(x=PC1, y=PC3, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  theme_bw(base_size = 7) +
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.04,0.04), ylim=c(-0.04,0.04))

gv.dorsal.plot3 <- ggplot(dorsal.gv.pc.data, aes(x=PC1, y=PC4, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  theme_bw(base_size = 7) +
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.04,0.04), ylim=c(-0.04,0.04))

gv.dorsal.plot4 <- ggplot(dorsal.gv.pc.data, aes(x=PC2, y=PC3, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  theme_bw(base_size = 7) +
  scale_x_continuous(breaks = c(-0.02, 0, 0.02, 0.04)) +
  theme(legend.position = "NONE")+
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.04,0.04), ylim=c(-0.04,0.04))

gv.dorsal.plot5 <- ggplot(dorsal.gv.pc.data, aes(x=PC2, y=PC4, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  theme_bw(base_size = 7) +
  scale_x_continuous(breaks = c(-0.02, 0, 0.02, 0.04)) +
  theme(legend.position = "NONE")+
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.04,0.04), ylim=c(-0.04,0.04))

gv.dorsal.plot6 <- ggplot(dorsal.gv.pc.data, aes(x=PC3, y=PC4, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  theme_bw(base_size = 7) +
  scale_x_continuous(breaks = c(-0.02, 0, 0.02)) +
  theme(legend.position = "NONE")+
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.04,0.04), ylim=c(-0.04,0.04))

gvdl <- ggplot(dorsal.gv.pc.data, aes(x=PC3, y=PC4, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  theme_bw(base_size = 10) 

gvdorsallegend <- cowplot::get_legend(gvdl)
grid.newpage()
grid.draw(gvdorsallegend)
plot(gvdorsallegend)

dorsal.gv.pc.plots <- gv.dorsal.plot1 + gv.dorsal.plot2 + gv.dorsal.plot3 + 
  gv.dorsal.plot4 + gv.dorsal.plot5 + gv.dorsal.plot6 + gvdorsallegend + plot_layout(design = layout)

dorsal.gv.pc.plots
ggsave("dorsal_gv_pc_plots.png")

# Composite PCA plots ----
composite.pc.plot.data <- dorsal.pc.data %>%
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
  ylim(-0.05,0.1) +
  coord_flip() +
  xlab("PC axis") +
  ylab("PC score") +
  labs(shape = "Current species")

  
megaplot_plot

# Misc Code ----
  geom_label(aes(label=SpecimenID))

print(pc_data)
options(max.print = 20)

## Aesthetics(X var, Y Var, shape based variable, colour based variable)
## geom_label(aes(label=XX))
