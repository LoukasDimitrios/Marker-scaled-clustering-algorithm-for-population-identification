source("Libraries.R")
source("Normalization.R")
source("Groups.R")


# Subset count
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


expr_sub_exp <- normalization(expr_sub, groups)


# PCA

# Define a list of subpopulations and their corresponding indices
subpopulations <- list(
  "T cells" = c("CD4pos", "CD8pos", "NKT", "TCRgd"),
  "B cells" = c("B cells"),
  "Monocytes" = c("Monocytes")
)

# Iterate over each subpopulation
for (subpop_name in names(subpopulations)) {
  subpop_indices <- pops_sub %in% subpopulations[[subpop_name]]
  expr_sub_subpop <- expr_sub[subpop_indices, ]
  pops_sub_subpop <- pops_sub[subpop_indices]
  subpops_sub_subpop <- subpops_sub[subpop_indices]
  
  # Perform PCA for the current subpopulation
  pca_result <- prcomp(expr_sub_subpop, scale. = TRUE)
  
  # Create a scree plot to visualize the variance explained by each PC
  scree_plot <- fviz_eig(pca_result, addlabels = TRUE)
  
  # Create a PCA plot with subpopulation labels
  pca_plot <- fviz_pca_ind(pca_result, 
                           geom.ind = "point", 
                           col.ind = as.factor(subpops_sub_subpop),  # Convert to factor
                           addEllipses = TRUE, 
                           ellipse.type = "confidence",
                           legend.title = "Subpopulation",
                           pointshape = 19) +
    labs(title = paste("Unscaled - PCA for", subpop_name))
  
  # Normalize the data with division by PC1^2
  expr_sub_exp <- NULL
  for (i in 1:length(groups)) {
    pca <- prcomp(expr_sub_subpop[, groups[[i]]], scale. = TRUE)
    expr_sub_exp <- cbind(expr_sub_exp, as.matrix(expr_sub_subpop[, groups[[i]]] / pca$sdev[1]^2))
  }
  
  # Perform PCA for the scaled data
  pca_result_scaled <- prcomp(expr_sub_exp, scale. = TRUE)
  
  # Create a scree plot for the scaled data
  scree_plot_scaled <- fviz_eig(pca_result_scaled, addlabels = TRUE)
  
  # Create a PCA plot for the scaled data with subpopulation labels
  pca_plot_scaled <- fviz_pca_ind(pca_result_scaled, 
                                  geom.ind = "point", 
                                  col.ind = as.factor(subpops_sub_subpop),  # Convert to factor
                                  addEllipses = TRUE, 
                                  ellipse.type = "confidence",
                                  legend.title = "Subpopulation",
                                  pointshape = 19) +
    labs(title = paste("Scaled - PCA for", subpop_name))
  
  # Perform UMAP for the unscaled data
  umap_result <- umap(expr_sub_subpop)
  
  # Create a data frame from UMAP results for the unscaled data
  umap_data <- data.frame(V1 = umap_result[, 1], V2 = umap_result[, 2])
  
  # Add subpopulation information to the data frame for the unscaled data
  umap_data$subpop <- subpops_sub_subpop
  
  # UMAP plot for the unscaled data with consistent colors
  umap_plot <- ggplot(umap_data, aes(x = V1, y = V2)) +
    geom_point(aes(color = subpops_sub_subpop)) +  # Color points by subpopulation
    labs(x = "UMAP 1", y = "UMAP 2") +  # Customize axis labels
    ggtitle(paste("Unscaled - UMAP Visualization of", subpop_name)) +  # Add a title
    guides(color = guide_legend(override.aes = list(size = 4)))  # Adjust legend point size
  
  # Perform UMAP for the scaled data
  umap_result_scaled <- umap(expr_sub_exp)
  
  # Create a data frame from UMAP results for the scaled data
  umap_data_scaled <- data.frame(V1 = umap_result_scaled[, 1], V2 = umap_result_scaled[, 2])
  
  # Add subpopulation information to the data frame for the scaled data
  umap_data_scaled$subpop <- subpops_sub_subpop
  
  # UMAP plot for the scaled data with consistent colors
  umap_plot_scaled <- ggplot(umap_data_scaled, aes(x = V1, y = V2)) +
    geom_point(aes(color = subpops_sub_subpop)) +  # Color points by subpopulation
    labs(x = "UMAP 1", y = "UMAP 2") +  # Customize axis labels
    ggtitle(paste("Scaled - UMAP Visualization of", subpop_name)) +  # Add a title
    guides(color = guide_legend(override.aes = list(size = 4)))  # Adjust legend point size
  
  # Print the plots
  print(scree_plot)
  print(scree_plot_scaled)
  print(pca_plot)
  print(pca_plot_scaled)
  print(umap_plot)
  print(umap_plot_scaled)
}
