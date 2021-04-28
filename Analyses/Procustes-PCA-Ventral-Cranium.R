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
ventral.gv.lands <- readland.tps(file = "Rawdata/sunda.ventral.crania.tps", specID = "ID")
ventral.gv.lands

# Plot raw y against raw x (Without procrustes - plotted points represent landmarks prior
## to rotation and scaling)
plot(ventral.lands[,2,]~ventral.lands[,1,])

# Plot G. Varigatus raw points
plot(ventral.gv.lands[,2,]~ventral.gv.lands[,1,])

# Normalisation of data through Generalised Procustes Analysis ----
## Procrustes analyses - uses reference specimen to line up all specimens 
gpa.ventral.lands <- gpagen(ventral.lands)

# Same with G. variegatus data
gpa.ventral.gv.lands <- gpagen(ventral.gv.lands)

# Visualising procrustes analysis ----

# Black points are mean position of coordinates with grey points showing variation
plot(gpa.ventral.lands)
plot(gpa.ventral.gv.lands)

# Shows average landmark coordinates
gpa.ventral.lands
gpa.ventral.gv.lands

# Shows the scaled coordinates of all specimens
gpa.ventral.lands$coords
gpa.ventral.gv.lands$coords

# Expanded information of gpa.lands. (Str:structure)
str(gpa.ventral.lands)
str(gpa.ventral.gv.lands)

# Reading in of metadata ----
metadata <- read.csv("Rawdata/dermopteradata.csv")
glimpse(metadata)

# PCA Analysis ----

# PCA Analysis (Geomorph Principal and phylogenetically-aligned 
## components analysis of shape data)
ventral.pca.landmarks <- gm.prcomp(gpa.ventral.lands$coords)

# Ditto for G.variegatus data
ventral.pca.gv.landmarks <- gm.prcomp(gpa.ventral.gv.lands$coords)

# How many PCs to include up to 95%? 
summary(ventral.pca.landmarks)
summary(ventral.pca.gv.landmarks)

## Which landmarks are contributing towards each principle component
ventral.pca.landmarks$rotation
ventral.pca.gv.landmarks$rotation

## Can be seen that 16 principal components explains for 95.7% of variation (within full
## data set including volans and variegatus)

# Read in metadata and prep PC data for merging  ----

# Extract PC scores and make ID name from TPS into a taxon column for 
## combining with the metadata
ventral.pc.scores <- data.frame(specimenID = rownames(ventral.pca.landmarks$x), 
                             ventral.pca.landmarks$x)
ventral.pc.scores

# same with variegatus
ventral.sunda.scores <- data.frame(specimenID = rownames(ventral.pca.gv.landmarks$x), 
                             ventral.pca.gv.landmarks$x)
ventral.sunda.scores

# Make column names begin with PC not Comp
colnames(ventral.pc.scores) <- gsub("Comp", "PC", colnames(ventral.pc.scores))
colnames(ventral.sunda.scores) <- gsub("Comp", "PC", colnames(ventral.sunda.scores))

# Merge with metadata ----
ventral.pc.data <- full_join(ventral.pc.scores, metadata, by = c("specimenID" = "Ventral.ID"))
ventral.gv.pc.data.unfixed <- full_join(ventral.sunda.scores, metadata, by = c("specimenID" = "Ventral.ID"))

ventral.gv.pc.data <- ventral.gv.pc.data.unfixed %>% 
  filter(!Region == "Philippines")


# Write PC scores to new csv files ----
write_csv(dorsal.pc.data, file = "Rawdata/ventral_dermoptera_pca_data.csv") 
write_csv(gv.pc.data, file = "Rawdata/ventral_variegatus_pca_data.csv") 

# Creation of wireframes ----

## Plot a reference shape
ventral.ref <- mshape(gpa.ventral.lands$coords)
sunda.ventral.ref <- mshape(gpa.ventral.gv.lands$coords)

# Plotting the minimum and maximum of x & y values 
plot(ventral.ref)
plotRefToTarget(ventral.ref, ventral.pca.landmarks$shapes$shapes.comp1$min)
plotRefToTarget(ventral.ref, ventral.pca.landmarks$shapes$shapes.comp1$max)
plotRefToTarget(ventral.ref, ventral.pca.landmarks$shapes$shapes.comp2$min)
plotRefToTarget(ventral.ref, ventral.pca.landmarks$shapes$shapes.comp2$max)
plotRefToTarget(ventral.ref, ventral.pca.landmarks$shapes$shapes.comp3$min)
plotRefToTarget(ventral.ref, ventral.pca.landmarks$shapes$shapes.comp3$max)

# Same for G.variegatus
plot(sunda.ventral.ref)
plotRefToTarget(sunda.ventral.ref, ventral.pca.gv.landmarks$shapes$shapes.comp1$min)
plotRefToTarget(sunda.ventral.ref, ventral.pca.gv.landmarks$shapes$shapes.comp1$max)
plotRefToTarget(sunda.ventral.ref, ventral.pca.gv.landmarks$shapes$shapes.comp2$min)
plotRefToTarget(sunda.ventral.ref, ventral.pca.gv.landmarks$shapes$shapes.comp2$max)
plotRefToTarget(sunda.ventral.ref, ventral.pca.gv.landmarks$shapes$shapes.comp3$min)
plotRefToTarget(sunda.ventral.ref, ventral.pca.gv.landmarks$shapes$shapes.comp3$max)

# Plots for ventral principal component analysis ----

ventral.plot1 <- ggplot(ventral.pc.data, aes(x=PC1, y=PC2, shape=CurrentSp, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  theme_bw(base_size = 7) +
  labs(shape = "Current species") + 
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.07,0.07), ylim=c(-0.04,0.04))

ventral.plot2 <- ggplot(ventral.pc.data, aes(x=PC1, y=PC3, shape=CurrentSp, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  theme_bw(base_size = 7) +
  labs(shape = "Current species") +
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.07,0.07), ylim=c(-0.04,0.04))

ventral.plot3 <- ggplot(ventral.pc.data, aes(x=PC1, y=PC4, shape=CurrentSp, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  theme_bw(base_size = 7) +
  labs(shape = "Current species") +
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +  
  coord_cartesian(xlim=c(-0.07,0.07), ylim=c(-0.04,0.04))

ventral.plot4 <- ggplot(ventral.pc.data, aes(x=PC2, y=PC3, shape=CurrentSp, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  theme_bw(base_size = 7) +
  labs(shape = "Current species") +
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.04,0.04), ylim=c(-0.04,0.04))

ventral.plot5 <- ggplot(ventral.pc.data, aes(x=PC2, y=PC4, shape=CurrentSp, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  theme_bw(base_size = 7) +
  labs(shape = "Current species") +
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.04,0.04), ylim=c(-0.04,0.04))

ventral.plot6 <- ggplot(ventral.pc.data, aes(x=PC3, y=PC4, shape=CurrentSp, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  theme_bw(base_size = 7) +
  labs(shape = "Current species") + 
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.04,0.04), ylim=c(-0.04,0.04))

# Isolation of legend through creation of new plot and cowplot pkg
gvvl <- ggplot(ventral.pc.data, aes(x=PC3, y=PC4, shape=CurrentSp, colour=Region)) + 
  geom_point() + 
  theme_bw(base_size = 10) +
  labs(shape = "Current species") 

gvventrallegend <- cowplot::get_legend(gvvl)
grid.newpage()
grid.draw(gvventrallegend)

plot(gvventrallegend)

layout <- "
AABBCC###
AABBCC###
AABBCCGGG
DDEEFFGGG
DDEEFF###
DDEEFF###
"

# Using patchwork to make composite figure
vderm.pc.plots <- ventral.plot1 + ventral.plot2 + ventral.plot3 + ventral.plot4 + ventral.plot5 + 
  ventral.plot6 + gvventrallegend + plot_layout(design = layout)

vderm.pc.plots
ggsave("ventral_dermoptera_pc_plots.png")

# Plots for ventral G.variegatus principal component analysis ----

gv.ventral.plot1 <- ggplot(ventral.gv.pc.data, aes(x=PC1, y=PC2, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  theme_bw(base_size = 7) +
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.04,0.04), ylim=c(-0.04,0.04))


gv.ventral.plot2 <- ggplot(ventral.gv.pc.data, aes(x=PC1, y=PC3, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  theme_bw(base_size = 7) +
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.04,0.04), ylim=c(-0.04,0.04))


gv.ventral.plot3 <- ggplot(ventral.gv.pc.data, aes(x=PC1, y=PC4, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  theme_bw(base_size = 7) +
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.04,0.04), ylim=c(-0.04,0.04))

gv.ventral.plot4 <- ggplot(ventral.gv.pc.data, aes(x=PC2, y=PC3, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  theme_bw(base_size = 7) +
  scale_x_continuous(breaks = c(-0.02, 0, 0.02, 0.04)) +
  theme(legend.position = "NONE")+
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.04,0.04), ylim=c(-0.04,0.04))

gv.ventral.plot5 <- ggplot(ventral.gv.pc.data, aes(x=PC2, y=PC4, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  theme_bw(base_size = 7) +
  scale_x_continuous(breaks = c(-0.02, 0, 0.02, 0.04)) +
  theme(legend.position = "NONE")+
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.04,0.04), ylim=c(-0.04,0.04))


gv.ventral.plot6 <- ggplot(ventral.gv.pc.data, aes(x=PC3, y=PC4, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  theme_bw(base_size = 7) +
  scale_x_continuous(breaks = c(-0.02, 0, 0.02)) +
  theme(legend.position = "NONE")+
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.04,0.04), ylim=c(-0.04,0.04))

gvdorsallegend <- ggplot(ventral.gv.pc.data, aes(x=PC3, y=PC4, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  theme_bw(base_size = 10) 

gv.dorsal.legend <- cowplot::get_legend(gvdorsallegend)
grid.newpage()
grid.draw(gv.dorsal.legend)
plot(gv.dorsal.legend)

ventral.gv.pc.plots <- gv.ventral.plot1 + gv.ventral.plot2 + gv.ventral.plot3 + gv.ventral.plot4 +
  gv.ventral.plot5 + gv.ventral.plot6 + gv.dorsal.legend + plot_layout(design = layout)

ventral.gv.pc.plots
ggsave("ventral_gv_pc_plots.png")

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
