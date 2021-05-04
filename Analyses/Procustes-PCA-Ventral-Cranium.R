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
v.lands.derm <- readland.tps(file = "Rawdata/ventral.derm.crania.tps", specID = "ID")
v.lands.derm

# Input of only G. variegatus tps data
v.lands.gv <- readland.tps(file = "Rawdata/ventral.sunda.crania.tps", specID = "ID")
v.lands.gv

#input of G.v. mainland data
v.lands.mainl <- readland.tps(file = "Rawdata/ventral.mainland.crania.tps", specID = "ID")
v.lands.mainl


# Plot raw y against raw x (Without procrustes - plotted points represent landmarks prior
## to rotation and scaling)
plot(v.lands.derm[,2,]~v.lands.derm[,1,])
plot(v.lands.gv[,2,]~v.lands.gv[,1,])
plot(v.lands.mainl[,2,]~v.lands.mainl[,1,])

# Normalisation of data through Generalised Procustes Analysis ----
## Procrustes analyses - uses reference specimen to line up all specimens 
v.gpal.derm <- gpagen(v.lands.derm)
v.gpal.gv <- gpagen(v.lands.gv)
v.gpal.mainl <- gpagen(v.lands.mainl)

# Visualising procrustes analysis ----

# Black points are mean position of coordinates with grey points showing variation
plot(v.gpal.derm)
plot(v.gpal.gv)
plot(v.gpal.mainl)

# Shows average landmark coordinates
v.gpal.derm
v.gpal.gv
v.gpal.mainl

# Shows the scaled coordinates of all specimens
v.gpal.derm$coords
v.gpal.gv$coords
v.gpal.mainl$coords

# Expanded information of gpa.lands. (Str:structure)
str(v.gpal.derm)
str(v.gpal.gv)
str(v.gpal.mainl)

# Reading in of metadata ----
metadata <- read.csv("Rawdata/dermopteradata.csv")
glimpse(metadata)

# PCA Analysis ----

# PCA Analysis (Geomorph Principal and phylogenetically-aligned 
## components analysis of shape data)
v.pcal.derm <- gm.prcomp(v.gpal.derm$coords)
v.pcal.gv <- gm.prcomp(v.gpal.gv$coords)
v.pcal.mainl  <- gm.prcomp(v.gpal.mainl$coords)

# How many PCs to include up to 95%? 
summary(v.pcal.derm)
summary(v.pcal.gv)
summary(v.pcal.mainl)

## Which landmarks are contributing towards each principle component
v.pcal.derm$rotation
v.pcal.gv$rotation
v.pcal.mainl$rotation

## Can be seen that 16 principal components explains for 95.7% of variation (within full
## data set including volans and variegatus)

# Read in metadata and prep PC data for merging  ----

# Extract PC scores and make ID name from TPS into a taxon column for 
## combining with the metadata. pcs - PC score
v.pcs.derm <- data.frame(specimenID = rownames(v.pcal.derm$x), 
                        v.pcal.derm$x)
v.pcs.derm

# same with variegatus
v.pcs.gv <- data.frame(specimenID = rownames(v.pcal.gv$x), 
                        v.pcal.gv$x)
v.pcs.gv

# Mainland data
v.pcs.mainl <- data.frame(specimenID = rownames(v.pcal.mainl$x), 
                        v.pcal.mainl$x)
v.pcs.mainl

# Make column names begin with PC not Comp
colnames(v.pcs.derm) <- gsub("Comp", "PC", colnames(v.pcs.derm))
colnames(v.pcs.gv) <- gsub("Comp", "PC", colnames(v.pcs.gv))
colnames(v.pcs.mainl) <- gsub("Comp", "PC", colnames(v.pcs.mainl))

# Merge with metadata ----
v.pcdata.derm <- full_join(v.pcs.derm, metadata, by = c("specimenID" = "Ventral.ID"))
v.pcdata.gv <- full_join(v.pcs.gv, metadata, by = c("specimenID" = "Ventral.ID")) %>% 
  filter(!Region == "Philippines")
v.pcdata.mainl <- full_join(v.pcs.mainl, metadata, by = c("specimenID" = "Ventral.ID"))  %>% 
  filter(!Landmass == "NA", !Region == "Philippines")

# Write PC scores to new csv files ----
write_csv(v.pcdata.derm, file = "Rawdata/ventral_dermoptera_pca_data.csv") 
write_csv(v.pcdata.gv, file = "Rawdata/ventral_variegatus_pca_data.csv") 
write_csv(v.pcdata.mainl, file = "Rawdata/ventral_mainland_pca_data.csv") 

# Creation of wireframes ----

## Plot a reference shape
v.ref.derm <- mshape(v.gpal.derm$coords)
v.ref.gv <- mshape(v.gpal.gv$coords)

# Plotting the minimum and maximum of x & y values 
plot(v.ref.derm)
plotRefToTarget(v.ref.derm, v.pcal.derm$shapes$shapes.comp1$min)
plotRefToTarget(v.ref.derm, v.pcal.derm$shapes$shapes.comp1$max)
plotRefToTarget(v.ref.derm, v.pcal.derm$shapes$shapes.comp2$min)
plotRefToTarget(v.ref.derm, v.pcal.derm$shapes$shapes.comp2$max)
plotRefToTarget(v.ref.derm, v.pcal.derm$shapes$shapes.comp3$min)
plotRefToTarget(v.ref.derm, v.pcal.derm$shapes$shapes.comp3$max)

# Same for G.variegatus
plot(v.ref.gv)
plotRefToTarget(v.ref.gv, v.pcal.gv$shapes$shapes.comp1$min)
plotRefToTarget(v.ref.gv, v.pcal.gv$shapes$shapes.comp1$max)
plotRefToTarget(v.ref.gv, v.pcal.gv$shapes$shapes.comp2$min)
plotRefToTarget(v.ref.gv, v.pcal.gv$shapes$shapes.comp2$max)
plotRefToTarget(v.ref.gv, v.pcal.gv$shapes$shapes.comp3$min)
plotRefToTarget(v.ref.gv, v.pcal.gv$shapes$shapes.comp3$max)

# Plots for ventral principal component analysis ----

v.pcplot.derm1 <- ggplot(v.pcdata.derm, aes(x=PC1, y=PC2, shape=CurrentSp, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  scale_fill_manual(values=c("#FF155B", "#00CD6C", "#009ADE", "#AF58BA", "#FFC61E")) +
  scale_color_manual(values = c("#FF155B", "#00CD6C", "#009ADE", "#AF58BA", "#FFC61E")) +
  theme_bw(base_size = 7) +
  labs(shape = "Current species") + 
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.07,0.07), ylim=c(-0.04,0.04))

v.pcplot.derm2 <- ggplot(v.pcdata.derm, aes(x=PC1, y=PC3, shape=CurrentSp, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  scale_fill_manual(values=c("#FF155B", "#00CD6C", "#009ADE", "#AF58BA", "#FFC61E")) +
  scale_color_manual(values = c("#FF155B", "#00CD6C", "#009ADE", "#AF58BA", "#FFC61E")) +
  theme_bw(base_size = 7) +
  labs(shape = "Current species") +
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.07,0.07), ylim=c(-0.04,0.04))

v.pcplot.derm3 <- ggplot(v.pcdata.derm, aes(x=PC1, y=PC4, shape=CurrentSp, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  scale_fill_manual(values=c("#FF155B", "#00CD6C", "#009ADE", "#AF58BA", "#FFC61E")) +
  scale_color_manual(values = c("#FF155B", "#00CD6C", "#009ADE", "#AF58BA", "#FFC61E")) +
  theme_bw(base_size = 7) +
  labs(shape = "Current species") +
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +  
  coord_cartesian(xlim=c(-0.07,0.07), ylim=c(-0.04,0.04))

v.pcplot.derm4 <- ggplot(v.pcdata.derm, aes(x=PC2, y=PC3, shape=CurrentSp, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  scale_fill_manual(values=c("#FF155B", "#00CD6C", "#009ADE", "#AF58BA", "#FFC61E")) +
  scale_color_manual(values = c("#FF155B", "#00CD6C", "#009ADE", "#AF58BA", "#FFC61E")) +
  theme_bw(base_size = 7) +
  labs(shape = "Current species") +
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.04,0.04), ylim=c(-0.04,0.04))

v.pcplot.derm5 <- ggplot(v.pcdata.derm, aes(x=PC2, y=PC4, shape=CurrentSp, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  scale_fill_manual(values=c("#FF155B", "#00CD6C", "#009ADE", "#AF58BA", "#FFC61E")) +
  scale_color_manual(values = c("#FF155B", "#00CD6C", "#009ADE", "#AF58BA", "#FFC61E")) +
  theme_bw(base_size = 7) +
  labs(shape = "Current species") +
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.04,0.04), ylim=c(-0.04,0.04))

v.pcplot.derm6 <- ggplot(v.pcdata.derm, aes(x=PC3, y=PC4, shape=CurrentSp, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  scale_fill_manual(values=c("#FF155B", "#00CD6C", "#009ADE", "#AF58BA", "#FFC61E")) +
  scale_color_manual(values = c("#FF155B", "#00CD6C", "#009ADE", "#AF58BA", "#FFC61E")) +
  theme_bw(base_size = 7) +
  labs(shape = "Current species") + 
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.04,0.04), ylim=c(-0.04,0.04))

# Isolation of legend through creation of new plot and cowplot pkg
v.pcplot.leg.d <- ggplot(v.pcdata.derm, aes(x=PC3, y=PC4, shape=CurrentSp, colour=Region)) + 
  geom_point() + 
  theme_bw(base_size = 10) +
  labs(shape = "Current species") +  
  scale_fill_manual(values=c("#FF155B", "#00CD6C", "#009ADE", "#AF58BA", "#FFC61E")) +
  scale_color_manual(values = c("#FF155B", "#00CD6C", "#009ADE", "#AF58BA", "#FFC61E")) 

v.pcplot.leg <- cowplot::get_legend(v.pcplot.leg.d)
grid.newpage()
grid.draw(v.pcplot.leg)

plot(v.pcplot.leg)

layout <- "
AABBCC###
AABBCC###
AABBCCGGG
DDEEFFGGG
DDEEFF###
DDEEFF###
"

# Using patchwork to make composite figure
v.pc.plots.derm <- v.pcplot.derm1 + v.pcplot.derm2 + v.pcplot.derm3 + v.pcplot.derm4 + v.pcplot.derm5 + 
  v.pcplot.derm6 + v.pcplot.leg + plot_layout(design = layout)

v.pc.plots.derm
ggsave("ventral_dermoptera_pc_plots.png", dpi = 900)

# Plots for ventral G.variegatus principal component analysis ----

v.pcplot.gv1 <- ggplot(v.pcdata.gv, aes(x=PC1, y=PC2, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  scale_fill_manual(values=c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) +
  scale_color_manual(values = c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) +
  theme_bw(base_size = 7) +
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.04,0.04), ylim=c(-0.04,0.04))


v.pcplot.gv2 <- ggplot(v.pcdata.gv, aes(x=PC1, y=PC3, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  scale_fill_manual(values=c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) +
  scale_color_manual(values = c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) +
  theme_bw(base_size = 7) +
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.04,0.04), ylim=c(-0.04,0.04))


v.pcplot.gv3 <- ggplot(v.pcdata.gv, aes(x=PC1, y=PC4, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  scale_fill_manual(values=c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) +
  scale_color_manual(values = c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) +
  theme_bw(base_size = 7) +
  theme(legend.position = "NONE") +
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.04,0.04), ylim=c(-0.04,0.04))

v.pcplot.gv4 <- ggplot(v.pcdata.gv, aes(x=PC2, y=PC3, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  scale_fill_manual(values=c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) +
  scale_color_manual(values = c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) +
  theme_bw(base_size = 7) +
  scale_x_continuous(breaks = c(-0.02, 0, 0.02, 0.04)) +
  theme(legend.position = "NONE")+
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.04,0.04), ylim=c(-0.04,0.04))

v.pcplot.gv5 <- ggplot(v.pcdata.gv, aes(x=PC2, y=PC4, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  scale_fill_manual(values=c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) +
  scale_color_manual(values = c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) +
  theme_bw(base_size = 7) +
  scale_x_continuous(breaks = c(-0.02, 0, 0.02, 0.04)) +
  theme(legend.position = "NONE")+
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.04,0.04), ylim=c(-0.04,0.04))


v.pcplot.gv6 <- ggplot(v.pcdata.gv, aes(x=PC3, y=PC4, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  scale_fill_manual(values=c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) +
  scale_color_manual(values = c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) +
  theme_bw(base_size = 7) +
  scale_x_continuous(breaks = c(-0.02, 0, 0.02)) +
  theme(legend.position = "NONE")+
  geom_vline(xintercept = 0, linetype='dotted', alpha = 0.8) +
  geom_hline(yintercept = 0, linetype='dotted', alpha = 0.8) +
  coord_cartesian(xlim=c(-0.04,0.04), ylim=c(-0.04,0.04))

v.pcplot.leg.d2 <- ggplot(v.pcdata.gv, aes(x=PC3, y=PC4, colour=Region)) + 
  geom_point(alpha = 0.75) + 
  theme_bw(base_size = 10) +
  scale_fill_manual(values=c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) +
  scale_color_manual(values = c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E"))

v.pcplot.leg.gv <- cowplot::get_legend(v.pcplot.leg.d2)
grid.newpage()
grid.draw(v.pcplot.leg.gv)
plot(v.pcplot.leg.gv)

v.pc.plots.gv <- v.pcplot.gv1 + v.pcplot.gv2 + v.pcplot.gv3 + v.pcplot.gv4 +
  v.pcplot.gv5 + v.pcplot.gv6 + v.pcplot.leg.gv + plot_layout(design = layout)

v.pc.plots.gv
ggsave("ventral_gv_pc_plots.png")


# Dermoptera composite PCA plot ---------

# Dermoptera composite plot data processing (PCs up to 85%)
v.compdata.derm <- v.pcdata.derm %>%
  dplyr::select(CurrentSp, Region, PC1:PC19) %>%
  pivot_longer(PC1:PC18, "PC", "value") %>%
  mutate(PC = factor(PC, levels = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6",
                                    "PC7", "PC8", "PC9", "PC10", "PC11", "PC12",
                                    "PC13", "PC14", "PC15", "PC16", "PC17", "PC18"))) 

# Processed Data
v.compdata.derm

# Plotting
v.compplot.derm <- 
  ggplot(v.compdata.derm, aes(x = PC, y = value, colour = Region, shape = CurrentSp)) +
  scale_fill_manual(values=c("#FF155B", "#00CD6C", "#009ADE", "#AF58BA", "#FFC61E"))+
  geom_boxplot() +
  geom_jitter(alpha = 0.5) +
  theme_bw(base_size = 14) +
  scale_color_manual(values = c("#FF155B", "#00CD6C", "#009ADE", "#AF58BA", "#FFC61E")) +
  ylim(-0.05,0.1) +
  coord_flip() +
  xlab("PC axis") +
  ylab("PC score") +
  labs(shape = "Current species") + 
  theme(legend.position = "NONE") 

#Check
v.compplot.derm

ggsave(file = "figures/ventral_derm_composite_plot.png", width = 5, height = 8, dpi = 900)


# G. variegatus composite PCA plot

# Plot data processing
v.compdata.gv <- v.pcdata.gv %>%
  dplyr::select(Region, PC1:PC20) %>%
  pivot_longer(PC1:PC19, "PC", "value") %>%
  mutate(PC = factor(PC, levels = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6",
                                    "PC7", "PC8", "PC9", "PC10", "PC11", "PC12",
                                    "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19"))) 

# Processed Data
v.compdata.gv

# Plotting
v.compplot.gv <- 
  ggplot(v.compdata.gv, aes(x = PC, y = value, colour = Region)) +
  scale_fill_manual(values=c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E"))+
  geom_boxplot() +
  geom_jitter(alpha = 0.5) +
  theme_bw(base_size = 14) +
  scale_color_manual(values = c("#FF155B", "#00CD6C", "#009ADE", "#FFC61E")) +
  ylim(-0.05,0.05) +
  coord_flip() +
  xlab("PC axis") +
  ylab("PC score") + 
  theme(legend.position = "NONE")

#Check
v.compplot.gv
# Save
ggsave(file = "figures/ventral_gv_composite_plot.png", width = 5, height = 8, dpi = 900)


# Mainland Composite plot data processing
v.compdata.mainl <- v.pcdata.mainl %>%
  dplyr::select(Landmass, PC1:PC18) %>%
  pivot_longer(PC1:PC17, "PC", "value") %>%
  mutate(PC = factor(PC, levels = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6",
                                    "PC7", "PC8", "PC9", "PC10", "PC11", "PC12",
                                    "PC13", "PC14", "PC15", "PC16", "PC17"))) 

# Processed Data
v.compdata.mainl

# Plotting
v.compplot.mainl <- 
  ggplot(v.compdata.mainl, aes(x = PC, y = value, colour = Landmass)) +
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

# Check
v.compplot.mainl

# Save
ggsave(file = "figures/ventral_mainland_composite_plot.png", width = 5, height = 8, dpi = 900)

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
