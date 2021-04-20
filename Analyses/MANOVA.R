# Input  of data ----

## Reading in data for dermoptera + G.variegatus groups
sunda_pc_data <- read_csv(file = "Rawdata/variegatus-pca-data-dorsal.csv")
pc_data <- read_csv(file = "Rawdata/colugo-pca-data-dorsal.csv")

# Multivariate analysis of variance for Dermoptera Data ----

# MANOVA of Dermoptera PCA data by species (Philippine vs Sunda)
derm.sp.manova <- manova(cbind(pc_data$PC1, pc_data$PC2) ~ CurrentSp,
                   data = pc_data)
summary(derm.sp.manova)

# MANOVA of Dermopetra PCA data by region
derm.r.manova <- manova(cbind(pc_data$PC1, pc_data$PC2) ~ Region,
                      data = pc_data)
summary(derm.r.manova)

# MANOVA of Dermoptera PCA data by sex
derm.s.manova <- manova(cbind(pc_data$PC1, pc_data$PC2) ~ Sex,
                        data = pc_data)
summary(derm.s.manova)

# MANOVA of Dermopetra PCA data by date photographed
derm.d.manova <- manova(cbind(pc_data$PC1, pc_data$PC2) ~ DatePhotographed,
                        data = pc_data)
summary(derm.d.manova)

# Multivariate analysis of variance for Sunda Data ----

# MANOVA of G.variegatus PCA data by region
sunda.manova <- manova(cbind(sunda_pc_data$PC1, sunda_pc_data$PC2) ~ Region,
                       data = sunda_pc_data)
summary(sunda.manova)
anova(sunda.manova)

# MANOVA of G.variegatus PCA data by sex
sunda.s.manova <- manova(cbind(pc_data$PC1, pc_data$PC2) ~ Sex,
                        data = sunda_pc_data)
summary(derm.s.manova)

# MANOVA of G.variegatus PCA data by date photographed
sunda.d.manova <- manova(cbind(pc_data$PC1, pc_data$PC2) ~ DatePhotographed,
                        data = pc_data)
summary(derm.d.manova)

# Misc Notes ----

# Column bind – merge function that can combine two data frames with the 
# same number of multiple rows into a single data frame.