# Input  of data ----

## Reading in data for dermoptera + G.variegatus groups
sunda_pc_data <- read_csv(file = "Rawdata/variegatus-pca-data-dorsal.csv")
pc_data <- read_csv(file = "Rawdata/colugo-pca-data-dorsal.csv")

# Multivariate analysis of variance for Dermoptera Data ----

# MANOVA of Dermoptera PCA data by species (Philippine vs Sunda)
# As matrix allows the use of multiple PC columns
derm.species.manova <- manova(as.matrix(pc_data[,2:19]) ~ CurrentSp,
                   data = pc_data)
summary(derm.species.manova)

# MANOVA of Dermopetra PCA data by region
derm.region.manova <- manova(as.matrix(pc_data[,2:19]) ~ Region,
                      data = pc_data)
summary(derm.region.manova)

# MANOVA of Dermoptera PCA data by sex
derm.sex.manova <- manova(as.matrix(pc_data[,2:19]) ~ Sex,
                        data = pc_data)
summary(derm.sex.manova)

# MANOVA of Dermopetra PCA data by date photographed
derm.date.manova <- manova(as.matrix(pc_data[,2:19]) ~ Date.Photographed,
                        data = pc_data)
summary(derm.date.manova)

# Multivariate analysis of variance for Sunda Data ----

# MANOVA of G.variegatus PCA data by region
sunda.region.manova <- manova(as.matrix(sunda_pc_data[,2:20]) ~ Region,
                       data = sunda_pc_data)
summary(sunda.region.manova)

# MANOVA of G.variegatus PCA data by sex
sunda.sex.manova <- manova(as.matrix(sunda_pc_data[,2:20]) ~ Sex,
                        data = sunda_pc_data)
summary(sunda.sex.manova)

# MANOVA of G.variegatus PCA data by date photographed
sunda.date.manova <- manova(as.matrix(sunda_pc_data[,2:20]) ~ Date.Photographed,
                        data = sunda_pc_data)
summary(sunda.date.manova)

# MANOVA of G.variegatus PCA data by sex and region
sunda.sr.manova <- manova(as.matrix(sunda_pc_data[,2:20]) ~ Sex * Region * Date.Photographed,
                          data = sunda_pc_data)
summary(sunda.sr.manova)

# Misc Notes ----

# Column bind – merge function that can combine two data frames with the 
# same number of multiple rows into a single data frame.