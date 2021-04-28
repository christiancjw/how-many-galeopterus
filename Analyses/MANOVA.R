# Input of data ----

## Reading in of dorsal data
dorsal.pc.data <- read_csv(file = "Rawdata/dorsal_dermoptera_pca_data.csv")
dorsal.gv.pc.data <- read_csv(file = "Rawdata/dorsal_variegatus_pca_data.csv")

## Reading in of ventral data
ventral.pc.data <- read_csv(file = "Rawdata/ventral_dermoptera_pca_data.csv")
ventral.gv.pc.data <- read_csv(file = "Rawdata/ventral_variegatus_pca_data.csv")

# Multivariate analysis of variance for dorsal Dermoptera data ----

# MANOVA of Dermoptera PCA data by species (Philippine vs Sunda)
# As matrix allows the use of multiple PC columns
# DD = Dorsal Dermoptera
DD.species.manova <- manova(as.matrix(dorsal.pc.data[,2:19]) ~ CurrentSp,
                   data = dorsal.pc.data)
summary(DD.species.manova)

# MANOVA of Dermopetra PCA data by region
DD.region.manova <- manova(as.matrix(dorsal.pc.data[,2:19]) ~ Region,
                      data = dorsal.pc.data)
summary(DD.region.manova)

# MANOVA of Dermoptera PCA data by sex
DD.sex.manova <- manova(as.matrix(dorsal.pc.data[,2:19]) ~ Sex,
                        data = dorsal.pc.data)
summary(DD.sex.manova)

# MANOVA of Dermopetra PCA data by date photographed 
## DV = Dorsal Vareigatus
DD.date.manova <- manova(as.matrix(dorsal.pc.data[,2:19]) ~ Date.Photographed,
                        data = dorsal.pc.data)
summary(DD.date.manova)

# Multivariate analysis of variance for dorsal G.variegatus data ----

# MANOVA of G.variegatus PCA data by region
DV.region.manova <- manova(as.matrix(dorsal.gv.pc.data[,2:20]) ~ Region,
                       data = dorsal.gv.pc.data)
summary(DV.region.manova)

# MANOVA of G.variegatus PCA data by sex
DV.sex.manova <- manova(as.matrix(dorsal.gv.pc.data[,2:20]) ~ Sex,
                        data = dorsal.gv.pc.data)
summary(DV.sex.manova)

# MANOVA of G.variegatus PCA data by date photographed
DV.date.manova <- manova(as.matrix(dorsal.gv.pc.data[,2:20]) ~ Date.Photographed,
                        data = dorsal.gv.pc.data)
summary(DV.date.manova)

# MANOVA of G.variegatus PCA data by sex and region
DV.sr.manova <- manova(as.matrix(dorsal.gv.pc.data[,2:20]) ~ Sex * Region * Date.Photographed,
                          data = dorsal.gv.pc.data)
summary(DV.sr.manova)

# Multivariate analysis of variance for ventral Dermoptera data ----
# As matrix allows the use of multiple PC columns
# VD = Ventral Dermoptera

# MANOVA of Dermoptera PCA data by species (Philippine vs Sunda)
VD.species.manova <- manova(as.matrix(ventral.pc.data[,2:27]) ~ CurrentSp,
                            data = ventral.pc.data)
summary(VD.species.manova)

# MANOVA of Ventral Dermopetra PCA data by region
VD.region.manova <- manova(as.matrix(ventral.pc.data[,2:27]) ~ Region,
                           data = ventral.pc.data)
summary(VD.region.manova)

# MANOVA of Ventral Dermoptera PCA data by sex
VD.sex.manova <- manova(as.matrix(ventral.pc.data[,2:27]) ~ Sex,
                        data = ventral.pc.data)
summary(VD.sex.manova)

# MANOVA of Ventral Dermopetra PCA data by date photographed 
VD.date.manova <- manova(as.matrix(ventral.pc.data[,2:27]) ~ Date.Photographed,
                         data = ventral.pc.data)
summary(VD.date.manova)

# Multivariate analysis of variance for ventral G.variegatus data ----
# VV = Ventral variegatus

# MANOVA of G.variegatus PCA data by region
VV.region.manova <- manova(as.matrix(ventral.gv.pc.data[,2:20]) ~ Region,
                           data = ventral.gv.pc.data)
summary(VV.region.manova)

# MANOVA of G.variegatus PCA data by sex
VV.sex.manova <- manova(as.matrix(ventral.gv.pc.data[,2:20]) ~ Sex,
                        data = ventral.gv.pc.data)
summary(VV.sex.manova)

# MANOVA of G.variegatus PCA data by date photographed
VV.date.manova <- manova(as.matrix(ventral.gv.pc.data[,2:20]) ~ Date.Photographed,
                         data = ventral.gv.pc.data)
summary(VV.date.manova)

# MANOVA of G.variegatus PCA data by sex and region
VV.sr.manova <- manova(as.matrix(ventral.gv.pc.data[,2:20]) ~ Sex * Region,
                       data = ventral.gv.pc.data)
summary(VV.sr.manova)

# Misc Notes ----

# Column bind – merge function that can combine two data frames with the 
# same number of multiple rows into a single data frame.