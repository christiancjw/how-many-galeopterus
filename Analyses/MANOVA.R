## Reading in data for dermoptera + just variegatus group
sunda_pc_data <- read_csv(file = "Rawdata/variegatus-pca-data-dorsal.csv")
pc_data <- read_csv(file = "Rawdata/colugo-pca-data-dorsal.csv")

## MANOVA
## Comma 
sunda.manova <- manova(cbind(sunda_pc_data$PC1, sunda_pc_data$PC2) ~ Region,
                       data = sunda_pc_data)

anova(sunda.manova)

derm.manova <- manova(cbind(pc_data$PC1, pc_data$PC2) ~ CurrentSp,
                   data = pc_data)
anova(derm.manova)
