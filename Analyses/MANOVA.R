# Input of data ----
library(tidyverse)

## Reading in of dorsal data
d.pcdata.derm<- read_csv(file = "Rawdata/csvfiles/dorsal_dermoptera_pca_data.csv")
d.pcdata.gv <- read_csv(file = "Rawdata/csvfiles/dorsal_variegatus_pca_data.csv")
d.pcdata.mainl <- read_csv(file = "Rawdata/csvfiles/dorsal_mainland_pca_data.csv")

## Reading in of ventral data
v.pcdata.derm <- read_csv(file = "Rawdata/csvfiles/ventral_dermoptera_pca_data.csv")
v.pcdata.gv <- read_csv(file = "Rawdata/csvfiles/ventral_variegatus_pca_data.csv")
v.pcdata.mainl <- read_csv(file = "Rawdata/csvfiles/ventral_mainland_pca_data.csv")

# Multivariate analysis of variance for dorsal Dermoptera data ----

# MANOVA of dorsal Dermoptera PCA data by species (Philippine vs Sunda)
# As matrix allows the use of multiple PC columns
derm.manova.dorsal.species <- manova(as.matrix(d.pcdata.derm[,2:19]) ~ CurrentSp,
                   data = d.pcdata.derm)
summary(derm.manova.dorsal.species)

# MANOVA of dorsal Dermopetra PCA data by region
derm.manova.dorsal.region <- manova(as.matrix(d.pcdata.derm[,2:19]) ~ Region,
                      data = d.pcdata.derm)
summary(derm.manova.dorsal.region)

# MANOVA of dorsal Dermoptera PCA data by sex
derm.manova.dorsal.sex <- manova(as.matrix(d.pcdata.derm[,2:19]) ~ Sex,
                        data = d.pcdata.derm)
summary(derm.manova.dorsal.sex)

# MANOVA of dorsal Dermopetra PCA data by date photographed 
## DV = Dorsal Vareigatus
derm.manova.dorsal.date <- manova(as.matrix(d.pcdata.derm[,2:19]) ~ Date.Photographed,
                        data = d.pcdata.derm)
summary(derm.manova.dorsal.date)

# Multivariate analysis of variance for dorsal G.variegatus data ----

# MANOVA of dorsal G.variegatus PCA data by region
gv.manova.dorsal.region <- manova(as.matrix(d.pcdata.gv[,2:20]) ~ Region,
                       data = d.pcdata.gv)
summary(gv.manova.dorsal.region)

# MANOVA of dorsal G.variegatus PCA data by sex
gv.manova.dorsal.sex <- manova(as.matrix(d.pcdata.gv[,2:20]) ~ Sex,
                        data = d.pcdata.gv)
summary(gv.manova.dorsal.sex)

# MANOVA of dorsal G.variegatus PCA data by date photographed
gv.manova.dorsal.date <- manova(as.matrix(d.pcdata.gv[,2:20]) ~ Date.Photographed,
                        data = d.pcdata.gv)
summary(gv.manova.dorsal.date)

# MANOVA of dorsal G.variegatus PCA data by insular group
gv.manova.dorsal.is <- manova(as.matrix(d.pcdata.gv[,2:20]) ~ Insular.gp,
                                data = d.pcdata.gv)
summary(gv.manova.dorsal.is)


# MANOVA of dorsal G.variegatus PCA data by sex and region
gv.manova.dorsal.sex.region <- manova(as.matrix(d.pcdata.gv[,2:20]) ~ Sex * Region * Date.Photographed,
                          data = d.pcdata.gv)
summary(gv.manova.dorsal.sex.region)

# Multivariate analysis of variance for dorsal mainland G.v. data ----

# MANOVA of dorsal mainland G.v. PCA data by region
mainl.manova.dorsal.landmass <- manova(as.matrix(d.pcdata.mainl[,2:17]) ~ Landmass,
                           data = d.pcdata.mainl)
summary(mainl.manova.dorsal.landmass)

# MANOVA of dorsal mainland G.v. PCA data by sex
mainl.manova.dorsal.sex <- manova(as.matrix(d.pcdata.mainl[,2:17]) ~ Sex,
                        data = d.pcdata.mainl)
summary(mainl.manova.dorsal.sex)

# MANOVA of dorsal mainland G.v. PCA data by date photographed
mainl.manova.dorsal.date <- manova(as.matrix(d.pcdata.mainl[,2:17]) ~ Date.Photographed,
                         data = d.pcdata.mainl)
summary(mainl.manova.dorsal.date)

# MANOVA of dorsal mainland G.v. PCA data by sex and region
mainl.manova.dorsal.sex.region <- manova(as.matrix(d.pcdata.mainl[,2:17]) ~ Sex * Region * Date.Photographed,
                       data = d.pcdata.mainl)
summary(mainl.manova.dorsal.sex.region)


# Multivariate analysis of variance for ventral Dermoptera data ----
# As matrix allows the use of multiple PC columns
# VD = Ventral Dermoptera

# MANOVA of ventral Dermoptera PCA data by species (Philippine vs Sunda)
derm.manova.ventral.species <- manova(as.matrix(v.pcdata.derm[,2:27]) ~ CurrentSp,
                            data = v.pcdata.derm)
summary(derm.manova.ventral.species)

# MANOVA of Ventral Dermopetra PCA data by region
derm.manova.ventral.region <- manova(as.matrix(v.pcdata.derm[,2:27]) ~ Region,
                           data = v.pcdata.derm)
summary(derm.manova.ventral.region)

# MANOVA of Ventral Dermoptera PCA data by sex
derm.manova.ventral.sex <- manova(as.matrix(v.pcdata.derm[,2:27]) ~ Sex,
                        data = v.pcdata.derm)
summary(derm.manova.ventral.sex)

# MANOVA of Ventral Dermopetra PCA data by date
derm.manova.ventral.dateregion <- manova(as.matrix(v.pcdata.derm[,2:27]) ~ Date.Photographed,
                                         data = v.pcdata.derm)
summary(derm.manova.ventral.date)


# MANOVA of Ventral Dermopetra PCA data for errors using date
derm.manova.ventral.dateregion <- manova(as.matrix(v.pcdata.derm[,2:27]) ~ Date.Photographed * Region,
                         data = v.pcdata.derm)
summary(derm.manova.ventral.dateregion)

derm.manova.ventral.regiondate <- manova(as.matrix(v.pcdata.derm[,2:27]) ~ Region * Date.Photographed,
                                         data = v.pcdata.derm)
summary(derm.manova.ventral.dateregion)

# Multivariate analysis of variance for ventral G.variegatus data ----

# MANOVA of G.variegatus PCA data by region
gv.manova.ventral.region <- manova(as.matrix(v.pcdata.gv[,2:28]) ~ Region,
                           data = v.pcdata.gv)
summary(gv.manova.ventral.region)

# MANOVA of G.variegatus PCA data by sex
gv.manova.ventral.sex <- manova(as.matrix(v.pcdata.gv[,2:28]) ~ Sex,
                        data = v.pcdata.gv)
summary(gv.manova.ventral.sex)

# MANOVA of ventral G.variegatus PCA data by date photographed
gv.manova.ventral.date <- manova(as.matrix(v.pcdata.gv[,2:28]) ~ Date.Photographed,
                         data = v.pcdata.gv)
summary(gv.manova.ventral.date)

# MANOVA of ventral G.variegatus PCA data by insular group
gv.manova.ventral.is <- manova(as.matrix(v.pcdata.gv[,2:28]) ~ Insular.gp,
                              data = v.pcdata.gv)
summary(gv.manova.ventral.is)


# MANOVA of ventral G.variegatus PCA data by sex and region
gv.manova.ventral.sex.region <- manova(as.matrix(v.pcdata.gv[,2:28]) ~ Sex * Region,
                       data = v.pcdata.gv)
summary(gv.manova.ventral.sex.region)

# Multivariate analysis of variance for ventral  mainland G.v. PCA data ----

# MANOVA of mainland G.v. PCA data by landmass
mainl.manova.ventral.landmass <- manova(as.matrix(v.pcdata.mainl[,2:23]) ~ Landmass,
                           data = v.pcdata.mainl)
summary(mainl.manova.ventral.landmass)

# MANOVA of mainland G.v. PCA data by sex
DV.sex.manova <- manova(as.matrix(v.pcdata.mainl[,2:23]) ~ Sex,
                        data = v.pcdata.mainl)
summary(DV.sex.manova)

# MANOVA of mainland G.v. PCA data by date photographed
DV.date.manova <- manova(as.matrix(v.pcdata.mainl[,2:23]) ~ Date.Photographed,
                         data = v.pcdata.mainl)
summary(DV.date.manova)

# MANOVA of mainland G.v. PCA data by sex and region
DV.sr.manova <- manova(as.matrix(v.pcdata.mainl[,2:20]) ~ Sex * Region * Date.Photographed,
                       data = v.pcdata.mainl)
summary(DV.sr.manova)


# Misc Notes ----

# Column bind – merge function that can combine two data frames with the 
# same number of multiple rows into a single data frame.
