#Libraries
source('Libraries.R')

#Load data
load("data_count.Rdata")

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
source('Groups.R')

# Load Normalization chunk for marker scaling groups
source('Normalization.R')

#Examine Populations using PCA
source('Populations.R')

#Subpopulations
source('Subpopulations.R')

#Perform GridSearch on the Unscaled data
source('GridSearch_unscaled.R')

#Perform GridSearch on the Scaled data with groups defined by expert immunologists
source('GridSearch_scaled.R')

#Perform GridSearch on the Scaled data with groups based on Log Fold Change
source('Log-Fold-Change_groups.R')
source('GridSearch_scaled_lfc.R')


# Perform Clustering on a bigger subset

# Load marker scaling groups defined by experts immunologists
source('Groups.R')

# Load Normalization chunk for marker scaling groups
source('Normalization.R')

# Create a bigger Subset data
set.seed(42)
subset <- sample(nrow(expr), 500000)
expr_sub <- expr[subset,]
pops_sub <- pops[subset]
subpops_sub <- subpops[subset]

#Seurat
source('Seurat_unscaled_clustering')
source('Seurat_scaled_clustering')

#Phenograph
source('Phenograph_clustering_unscaled.R')
source('Phenograph_clustering_scaled.R')
