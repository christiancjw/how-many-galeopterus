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
d.lands.derm <- readland.tps(file = "Rawdata/dorsal.derm.crania.tps", specID = "ID")
d.lands.derm

# Input of only G. variegatus tps data
d.lands.gv <- readland.tps(file = "Rawdata/dorsal.sunda.crania.tps", specID = "ID")
d.lands.gv

# Imput of mainland G. variegatus data
d.lands.mainl <- readland.tps(file = "Rawdata/dorsal.mainland.crania.tps", specID = "ID")
d.lands.mainl

# Plot raw y against raw x (Without procrustes - plotted points represent landmarks prior
## to rotation and scaling)
plot(d.lands.derm[,2,]~d.lands.derm[,1,])
plot(d.lands.gv[,2,]~d.lands.gv[,1,])
plot(d.lands.mainl[,2,]~d.lands.mainl[,1,])

# Normalisation of data through Generalised Procustes Analysis ----
## Procrustes analyses - uses reference specimen to line up all specimens 
# gpal - Generalised Procrustes Analysis Landmarks
d.gpal.derm <- gpagen(d.lands.derm)
d.gpal.gv <- gpagen(d.lands.gv)
d.gpal.mainl <- gpagen(d.lands.mainl)

# Visualising procrustes analysis ----

# Black points are mean position of coordinates with grey points showing variation
plot(d.gpal.derm)
plot(d.gpal.gv)
plot(d.gpal.mainl)

# Shows average landmark coordinates
d.gpal.derm
d.gpal.gv
d.gpal.mainl

# Shows the scaled coordinates of all specimens
d.gpal.derm$coords
d.gpal.gv$coords
d.gpal.mainl$coords

# Expanded information of gpa.lands. (Str:structure)
str(d.gpal.gv)
str(d.gpal.gv)
str(d.gpal.mainl)

# Reading in of metadata ----
metadata <- read.csv("Rawdata/dermopteradata.csv")
glimpse(metadata)

# PCA Analysis ----

# PCA Analysis (Geomorph Principal and phylogenetically-aligned 
## components analysis of shape data). pcal - Principal component analysis landmarks
d.pcal.derm <- gm.prcomp(d.gpal.derm$coords)

# Ditto for G.variegatus data
d.pcal.gv <- gm.prcomp(d.gpal.gv$coords)

# Mainland data
d.pcal.mainl <- gm.prcomp(d.gpal.mainl$coords)

# How many PCs to include up to 95%? 
summary(d.pcal.derm)
summary(d.pcal.gv)
summary(d.pcal.mainl)

## Which landmarks are contributing towards each principle component
d.pcal.derm$rotation
d.pcal.gv$rotation
d.pcal.mainl$rotation

## Can be seen that 16 principal components explains for 95.7% of variation (within full
## data set including volans and variegatus)

# Read in metadata and prep PC data for merging ----

# Extract PC scores and make ID name from TPS into a taxon column for 
## combining with the metadata. pcs - PC score
d.pcs.derm <- data.frame(specimenID = rownames(d.pcal.derm$x), 
                         d.pcal.derm$x)
d.pcs.derm

# same with variegatus
d.pcs.gv <- data.frame(specimenID = rownames(d.pcal.gv$x), 
                       d.pcal.gv$x)
d.pcs.gv

# Mainland data
d.pcs.mainl <-data.frame(specimenID = rownames(d.pcal.mainl$x), 
                         d.pcal.mainl$x)
d.pcs.mainl

# Make column names begin with PC not Comp
colnames(d.pcs.derm) <- gsub("Comp", "PC", colnames(d.pcs.derm))
colnames(d.pcs.gv) <- gsub("Comp", "PC", colnames(d.pcs.gv))
colnames(d.pcs.mainl) <- gsub("Comp", "PC", colnames(d.pcs.mainl))

# Merge with metadata ----
d.pcdata.derm <- full_join(d.pcs.derm, metadata, by = c("specimenID" = "Dorsal.ID"))
d.pcdata.gv <- full_join(d.pcs.gv, metadata, by = c("specimenID" = "Dorsal.ID")) %>% 
  filter(!Region == "Philippines")
d.pcdata.mainl <- full_join(d.pcs.mainl, metadata, by = c("specimenID" = "Dorsal.ID")) %>% 
  filter(!Landmass == "NA", !Region == "Philippines")

# Write PC scores to new csv files ----
write_csv(d.pcdata.derm, file = "Rawdata/csvfiles/dorsal_dermoptera_pca_data.csv") 
write_csv(d.pcdata.gv, file = "Rawdata/csvfiles/dorsal_variegatus_pca_data.csv") 
write_csv(d.pcdata.mainl, file = "Rawdata/csvfiles/dorsal_mainland_pca_data.csv") 


# Creation of wireframes ----

## Plot a reference shape
d.ref.derm <- mshape(d.gpal.derm$coords)
d.ref.gv <- mshape(d.gpal.gv$coords)

# Plotting the minimum and maximum of x & y values 
plot(d.ref.derm)
plotRefToTarget(d.ref.derm, d.pcal.derm$shapes$shapes.comp1$min)
plotRefToTarget(d.ref.derm, d.pcal.derm$shapes$shapes.comp1$max)
plotRefToTarget(d.ref.derm, d.pcal.derm$shapes$shapes.comp2$min)
plotRefToTarget(d.ref.derm, d.pcal.derm$shapes$shapes.comp2$max)
plotRefToTarget(d.ref.derm, d.pcal.derm$shapes$shapes.comp3$min)
plotRefToTarget(d.ref.derm, d.pcal.derm$shapes$shapes.comp3$max)

# Same for G.variegatus
plot(d.ref.gv)
plotRefToTarget(d.ref.gv, d.pcal.gv$shapes$shapes.comp1$min)
plotRefToTarget(d.ref.gv, d.pcal.gv$shapes$shapes.comp1$max)
plotRefToTarget(d.ref.gv, d.pcal.gv$shapes$shapes.comp2$min)
plotRefToTarget(d.ref.gv, d.pcal.gv$shapes$shapes.comp2$max)
plotRefToTarget(d.ref.gv, d.pcal.gv$shapes$shapes.comp3$min)
plotRefToTarget(d.ref.gv, d.pcal.gv$shapes$shapes.comp3$max)

# Plots for dorsal Dermoptera principal component analysis ----

## Plotting dorsal landmarks
d.pcplot.derm1 <- ggplot(d.pcdata.derm, aes(x=PC1, y=PC2, shape=CurrentSp, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  scale_fill_manual(values=c("#FF155B", "#00CD6C", "#009ADE", "#AF58BA", "#FFC61E")) +
  scale_color_manual(values = c("#FF155B", "#00CD6C", "#009ADE", "#AF58BA", "#FFC61E")) +
  theme_bw(base_size = 7) +
  labs(shape = "Current species") + 
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.05,0.10), ylim=c(-0.05,0.05)) 

d.pcplot.derm2 <- ggplot(d.pcdata.derm, aes(x=PC1, y=PC3, shape=CurrentSp, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  scale_fill_manual(values=c("#FF155B", "#00CD6C", "#009ADE", "#AF58BA", "#FFC61E")) +
  scale_color_manual(values = c("#FF155B", "#00CD6C", "#009ADE", "#AF58BA", "#FFC61E")) + 
  theme_bw(base_size = 7) +
  labs(shape = "Current species") +
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.05,0.10), ylim=c(-0.05,0.05))

d.pcplot.derm3 <- ggplot(d.pcdata.derm, aes(x=PC1, y=PC4, shape=CurrentSp, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  scale_fill_manual(values=c("#FF155B", "#00CD6C", "#009ADE", "#AF58BA", "#FFC61E")) +
  scale_color_manual(values = c("#FF155B", "#00CD6C", "#009ADE", "#AF58BA", "#FFC61E")) +
  theme_bw(base_size = 7) +
  labs(shape = "Current species") +
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.05,0.10), ylim=c(-0.05,0.05))

d.pcplot.derm4 <- ggplot(d.pcdata.derm, aes(x=PC2, y=PC3, shape=CurrentSp, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  scale_fill_manual(values=c("#FF155B", "#00CD6C", "#009ADE", "#AF58BA", "#FFC61E")) +
  scale_color_manual(values = c("#FF155B", "#00CD6C", "#009ADE", "#AF58BA", "#FFC61E")) +
  theme_bw(base_size = 7) +
  labs(shape = "Current species") +
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.05,0.05), ylim=c(-0.04,0.04))

d.pcplot.derm5 <- ggplot(d.pcdata.derm, aes(x=PC2, y=PC4, shape=CurrentSp, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  scale_fill_manual(values=c("#FF155B", "#00CD6C", "#009ADE", "#AF58BA", "#FFC61E")) +
  scale_color_manual(values = c("#FF155B", "#00CD6C", "#009ADE", "#AF58BA", "#FFC61E")) +
  theme_bw(base_size = 7) +
  labs(shape = "Current species") +
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.05,0.05), ylim=c(-0.04,0.04))

d.pcplot.derm6 <- ggplot(d.pcdata.derm, aes(x=PC3, y=PC4, shape=CurrentSp, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  scale_fill_manual(values=c("#FF155B", "#00CD6C", "#009ADE", "#AF58BA", "#FFC61E")) +
  scale_color_manual(values = c("#FF155B", "#00CD6C", "#009ADE", "#AF58BA", "#FFC61E")) +
  theme_bw(base_size = 7) +
  labs(shape = "Current species") + 
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.05,0.05), ylim=c(-0.04,0.04))

# Isolation of legend through creation of new plot and cowplot pkg
d.pcplot.leg.d <- ggplot(d.pcdata.derm, aes(x=PC3, y=PC4, shape=CurrentSp, colour=Region)) + 
  geom_point() + 
  theme_bw(base_size = 10) +
  labs(shape = "Current species") +
  scale_fill_manual(values=c("#FF155B", "#00CD6C", "#009ADE", "#AF58BA", "#FFC61E")) +
  scale_color_manual(values = c("#FF155B", "#00CD6C", "#009ADE", "#AF58BA", "#FFC61E"))

d.pc.legend <- cowplot::get_legend(d.pcplot.leg.d)
grid.newpage()
grid.draw(d.pc.legend)

plot(d.pc.legend)

layout <- "
AABBCC###
AABBCC###
AABBCCGGG
DDEEFFGGG
DDEEFF###
DDEEFF###
"

# Using patchwork to make composite figure
d.pc.plots.derm <- d.pcplot.derm1 + d.pcplot.derm2 + d.pcplot.derm3 + d.pcplot.derm4 + d.pcplot.derm5 + 
  d.pcplot.derm6 + d.pc.legend + plot_layout(design = layout)

d.pc.plots.derm
ggsave(file = "dorsal_dermoptera_pc_plots.png", height = 4.95, width = 11.475, dpi = 900)

# Plots for dorsal G.variegatus principal component analysis ----

d.pcplot.gv1 <- ggplot(d.pcdata.gv, aes(x=PC1, y=PC2, colour=Region)) + 
  scale_fill_manual(values=c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) +
  scale_color_manual(values = c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) +
  geom_point(alpha = 0.75) + 
  theme_bw(base_size = 7) +
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.04,0.04), ylim=c(-0.05,0.05))

d.pcplot.gv2 <- ggplot(d.pcdata.gv, aes(x=PC1, y=PC3, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  scale_fill_manual(values=c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) +
  scale_color_manual(values = c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) +
  theme_bw(base_size = 7) +
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.04,0.04), ylim=c(-0.04,0.04))

d.pcplot.gv3 <- ggplot(d.pcdata.gv, aes(x=PC1, y=PC4, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  scale_fill_manual(values=c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) +
  scale_color_manual(values = c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) +
  theme_bw(base_size = 7) +
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.04,0.04), ylim=c(-0.04,0.04))

d.pcplot.gv4 <- ggplot(d.pcdata.gv, aes(x=PC2, y=PC3, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  scale_fill_manual(values=c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) +
  scale_color_manual(values = c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) +
  theme_bw(base_size = 7) +
  scale_x_continuous(breaks = c(-0.02, 0, 0.02, 0.04)) +
  theme(legend.position = "NONE")+
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.04,0.04), ylim=c(-0.04,0.04))

d.pcplot.gv5 <- ggplot(d.pcdata.gv, aes(x=PC2, y=PC4, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  scale_fill_manual(values=c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) +
  scale_color_manual(values = c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) +
  theme_bw(base_size = 7) +
  scale_x_continuous(breaks = c(-0.02, 0, 0.02, 0.04)) +
  theme(legend.position = "NONE")+
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.04,0.04), ylim=c(-0.04,0.04))

d.pcplot.gv6 <- ggplot(d.pcdata.gv, aes(x=PC3, y=PC4, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  scale_fill_manual(values=c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) +
  scale_color_manual(values = c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) +
  theme_bw(base_size = 7) +
  scale_x_continuous(breaks = c(-0.02, 0, 0.02)) +
  theme(legend.position = "NONE")+
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.04,0.04), ylim=c(-0.04,0.04))

d.pcplot.leg.d2 <- ggplot(d.pcdata.gv, aes(x=PC3, y=PC4, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  theme_bw(base_size = 10) + 
  scale_fill_manual(values=c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) +
  scale_color_manual(values = c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) 

d.pcplot.leg.gv <- cowplot::get_legend(d.pcplot.leg.d2)
grid.newpage()
grid.draw(d.pcplot.leg.gv)
plot(d.pcplot.leg.gv)

d.pc.plots.gv <- d.pcplot.gv1 + d.pcplot.gv2 + d.pcplot.gv3 + 
  d.pcplot.gv4 + d.pcplot.gv5 + d.pcplot.gv6 + d.pcplot.leg.gv + plot_layout(design = layout)

d.pc.plots.gv
ggsave(file = "dorsal_gv_pc_plots.png", height = 4.95, width = 11.475, dpi = 900)


# Plots for dorsal Mainland principal component analysis ----

d.pcplot.mainl1 <- ggplot(d.pcdata.mainl, aes(x=PC1, y=PC2, colour=Region)) + 
  scale_fill_manual(values=c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) +
  scale_color_manual(values = c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) +
  geom_point(alpha = 0.75) + 
  theme_bw(base_size = 7) +
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.04,0.04), ylim=c(-0.05,0.05))

d.pcplot.mainl2 <- ggplot(d.pcdata.mainl, aes(x=PC1, y=PC3, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  scale_fill_manual(values=c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) +
  scale_color_manual(values = c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) +
  theme_bw(base_size = 7) +
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.04,0.04), ylim=c(-0.04,0.04))

d.pcplot.mainl3 <- ggplot(d.pcdata.mainl, aes(x=PC1, y=PC4, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  scale_fill_manual(values=c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) +
  scale_color_manual(values = c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) +
  theme_bw(base_size = 7) +
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.04,0.04), ylim=c(-0.04,0.04))

d.pcplot.mainl4 <- ggplot(d.pcdata.mainl, aes(x=PC2, y=PC3, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  scale_fill_manual(values=c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) +
  scale_color_manual(values = c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) +
  theme_bw(base_size = 7) +
  scale_x_continuous(breaks = c(-0.02, 0, 0.02, 0.04)) +
  theme(legend.position = "NONE")+
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.04,0.04), ylim=c(-0.04,0.04))

d.pcplot.mainl5 <- ggplot(d.pcdata.mainl, aes(x=PC2, y=PC4, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  scale_fill_manual(values=c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) +
  scale_color_manual(values = c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) +
  theme_bw(base_size = 7) +
  scale_x_continuous(breaks = c(-0.02, 0, 0.02, 0.04)) +
  theme(legend.position = "NONE")+
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.04,0.04), ylim=c(-0.04,0.04))

d.pcplot.mainl6 <- ggplot(d.pcdata.mainl, aes(x=PC3, y=PC4, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  scale_fill_manual(values=c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) +
  scale_color_manual(values = c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) +
  theme_bw(base_size = 7) +
  scale_x_continuous(breaks = c(-0.02, 0, 0.02)) +
  theme(legend.position = "NONE")+
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.04,0.04), ylim=c(-0.04,0.04))

d.pcplot.leg.d3 <- ggplot(d.pcdata.mainl, aes(x=PC3, y=PC4, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  theme_bw(base_size = 10) + 
  scale_fill_manual(values=c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) +
  scale_color_manual(values = c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) 

d.pcplot.leg.mail <- cowplot::get_legend(d.pcplot.leg.d3)
grid.newpage()
grid.draw(d.pcplot.leg.mail)
plot(d.pcplot.leg.mail)

d.pc.plots.mainl <- d.pcplot.mainl1 + d.pcplot.mainl2 + d.pcplot.mainl3 + 
  d.pcplot.mainl4 + d.pcplot.mainl5+ d.pcplot.mainl6 + d.pcplot.leg.mail + plot_layout(design = layout)

d.pc.plots.mainl
ggsave(file = "dorsal_mainl_pc_plots.png", height = 4.95, width = 11.475, dpi = 900)



# Composite PCA plots ----

# Dermoptera composite plot data processing
d.compdata.derm <- d.pcdata.derm %>%
  dplyr::select(CurrentSp, Region, PC1:PC19) %>%
  pivot_longer(PC1:PC18, "PC", "value") %>%
  mutate(PC = factor(PC, levels = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6",
                                    "PC7", "PC8", "PC9", "PC10", "PC11", "PC12",
                                    "PC13", "PC14", "PC15", "PC16", "PC17", "PC18"))) 
# Processed Data
d.compdata.derm
# Plotting
d.compplot.derm <- 
  ggplot(d.compdata.derm, aes(x = PC, y = value, colour = Region, shape = CurrentSp)) +
  scale_fill_manual(values=c("#FF155B", "#00CD6C", "#009ADE", "#AF58BA", "#FFC61E"))+
  geom_boxplot() +
  geom_jitter(alpha = 0.5) +
  theme_bw(base_size = 14) +
  scale_color_manual(values = c("#FF155B", "#00CD6C", "#009ADE", "#AF58BA", "#FFC61E")) +
  ylim(-0.05,0.1) +
  theme(legend.position = "NONE") +
  coord_flip() +
  xlab("PC axis") +
  ylab("PC score") +
  labs(shape = "Current species")
#Check
Dorsal.composite.pc.plot
ggsave(file = "figures/dorsal_derm_composite_plot.png", height = 4.95, width = 11.475, dpi = 900)


# G. variegatus composite plot data processing
d.compdata.gv <- d.pcdata.gv %>%
  dplyr::select(Region, PC1:PC20) %>%
  pivot_longer(PC1:PC19, "PC", "value") %>%
  mutate(PC = factor(PC, levels = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6",
                                    "PC7", "PC8", "PC9", "PC10", "PC11", "PC12",
                                    "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19"))) 
# Processed Data
d.compdata.gv

# Plotting
d.compplot.derm <- 
  ggplot(dorsal.gv.pc.plot.data, aes(x = PC, y = value, colour = Region)) +
  scale_fill_manual(values=c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E"))+
  geom_boxplot() +
  geom_jitter(alpha = 0.5) +
  theme_bw(base_size = 14) +
  scale_color_manual(values = c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) +
  ylim(-0.05,0.05) +
  theme(legend.position = "NONE")+
  coord_flip() +
  xlab("PC axis") +
  ylab("PC score") 

#Check
d.compplot.derm
ggsave(file = "figures/dorsal_gv_composite_plot.png", height = 4.95, width = 11.475, dpi = 900)

# Mainland Composite plot data processing
d.compdata.mainl <- d.pcdata.mainl %>%
  dplyr::select(Landmass, PC1:PC17) %>%
  pivot_longer(PC1:PC16, "PC", "value") %>%
  mutate(PC = factor(PC, levels = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6",
                                    "PC7", "PC8", "PC9", "PC10", "PC11", "PC12",
                                    "PC13", "PC14", "PC15", "PC16"))) 

# Processed Data
d.compdata.mainl

# Plotting
d.compplot.mainl <- 
  ggplot(d.compdata.mainl, aes(x = PC, y = value, colour = Landmass)) +
  scale_fill_manual(values=c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E"))+
  geom_boxplot() +
  geom_jitter(alpha = 0.5) +
  theme_bw(base_size = 14) +
  scale_color_manual(values = c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) +
  ylim(-0.05,0.05) +
  coord_flip() +
  xlab("PC axis") +
  theme(legend.position = "NONE")+
  ylab("PC score") 

#Check
d.compplot.mainl

ggsave(file = "figures/dorsal_mainland_composite_plot.png", width = 5, height = 8, dpi = 900)

#
ggplot(dorsal.mainland.pc.data, aes(x=PC1, y=PC2, colour=Landmass)) + 
  geom_point(alpha = 0.75) + 
  scale_fill_manual(values=c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) +
  scale_color_manual(values = c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) +
  theme_bw(base_size = 7) +
  scale_x_continuous(breaks = c(-0.04, 0, 0.02)) +
  theme(legend.position = "NONE")+
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.1,0.04), ylim=c(-0.04,0.04)) +
  geom_label(aes(label=specimenID)) 


 # Misc Code ----
  geom_label(aes(label=SpecimenID))

print(pc_data)
options(max.print = 20)

## Aesthetics(X var, Y Var, shape based variable, colour based variable)
## geom_label(aes(label=XX))
