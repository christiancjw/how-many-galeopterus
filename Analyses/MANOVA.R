# Input  of data ----

## Reading in data for dermoptera + G.variegatus groups
sunda_pc_data <- read_csv(file = "Rawdata/variegatus-pca-data-dorsal.csv")
pc_data <- read_csv(file = "Rawdata/colugo-pca-data-dorsal.csv")

# Multivariate analysis of variance for Dermoptera Data ----

# MANOVA of Dermoptera PCA data by species (Philippine vs Sunda)
derm.manova <- manova(cbind(pc_data$PC1, pc_data$PC2) ~ CurrentSp,
                   data = pc_data)
anova(derm.manova)

# Multivariate analysis of variance for Sunda Data ----

# MANOVA of G.variegatus PCA data by region
sunda.manova <- manova(cbind(sunda_pc_data$PC1, sunda_pc_data$PC2) ~ Region,
                       data = sunda_pc_data)

anova(sunda.manova)

# Misc Notes ----

# Column bind – merge function that can combine two data frames with the 
# same number of multiple rows into a single data frame.