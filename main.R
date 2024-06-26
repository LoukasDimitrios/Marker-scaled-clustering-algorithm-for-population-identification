#Libraries
source('Scripts/Libraries.R')

#Load data
load("Scripts/data_count.Rdata")

# Create a Subset
pops <- expr_count$population
subpops <- expr_count$subpopulation
lineage_channels <- c("CD57", "CD19", "CD4", "CD8", "IgD", "CD11c", "CD16", "CD3", "CD38", "CD27", "CD14", "CXCR5", "CCR7", "CD45RA", "CD20", "CD127", "CD33", "CD28", "CD161", "TCRgd", "CD123", "CD56", "HLADR", "CD25")
expr <- expr_count[,lineage_channels]


# Subset data
set.seed(42)
subset <- sample(nrow(expr), 10000)
expr_sub <- expr[subset,]
pops_sub <- pops[subset]
subpops_sub <- subpops[subset]

# Load marker scaling groups defined by experts immunologists
source('Scripts/Groups.R')

# Load Normalization chunk for marker scaling groups
source('Scripts/Normalization.R')

#Examine Populations using PCA
source('Scripts/Populations.R')

#Subpopulations
source('Scripts/Subpopulations.R')

#Perform GridSearch on the Unscaled data
source('Scripts/GridSearch_unscaled.R')

#Perform GridSearch on the Scaled data with groups defined by expert immunologists
source('Scripts/GridSearch_scaled.R')

#Perform GridSearch on the Scaled data with groups based on Log Fold Change
source('Scripts/Log-Fold-Change_groups.R')
source('Scripts/GridSearch_scaled_lfc.R')


# Perform Clustering on a bigger subset

# Load marker scaling groups defined by experts immunologists
source('Scripts/Groups.R')

# Load Normalization chunk for marker scaling groups
source('Scripts/Normalization.R')

# Create a bigger Subset data
set.seed(42)
subset <- sample(nrow(expr), 500000)
expr_sub <- expr[subset,]
pops_sub <- pops[subset]
subpops_sub <- subpops[subset]

#Seurat
source('Scripts/Seurat_unscaled_clustering')
source('Scripts/Seurat_scaled_clustering')

#Phenograph
source('Scripts/Phenograph_clustering_unscaled.R')
source('Scripts/Phenograph_clustering_scaled.R')
