source("Rphenograph_algorithm.R")

### PCA - UNSCALED vs UNSCALED


# PCA - POPULATION without normalization (without division by PC1^2)
pca_result <- prcomp(expr_sub, scale. = TRUE)  # Standardize the data

# Create a scree plot to visualize the variance explained by each PC
scree_plot <- fviz_eig(pca_result,
                       addlabels = TRUE,
                              title = "Scree plot UNSCALED")

# PCA plot with population names in legends using the default color palette
pca_plot <- fviz_pca_ind(pca_result, 
                         geom.ind = "point", 
                         col.ind = pops_sub,  # Color points by population
                         addEllipses = TRUE, 
                         ellipse.type = "confidence",
                         legend.title = "Population",
                         pointshape = 19,
                         fill = pops_sub,
                         title = "PCA visualization UNSCALED")  # Set fill aesthetic to population names

# Display the scree plot
print(scree_plot)

# Display the PCA plot
print(pca_plot)


result <- Rphenograph2(pca_result$x[, 1:best_p_unscaled], k = best_k_neighbors_unscaled, resolution = best_resolution_unscaled)

# Get cluster assignments from Phenograph
cluster_assignments <- result[[2]]$membership

# Perform UMAP on the first 5 PCs
umap_result <- umap(pca_result$x[, ])
umap(pca_result$x[, ], batch = TRUE, n_threads = 32, verbose = FALSE, n_neighbors = 15, n_sgd_threads = "auto")

# Create a data frame with UMAP results and cluster assignments
umap_data <- data.frame(V1 = umap_result[, 1], V2 = umap_result[, 2], Cluster = factor(cluster_assignments))

# Plot UMAP with cluster colors
umap_plot <- ggplot(umap_data, aes(x = V1, y = V2, color = Cluster)) +
  geom_point() +
  labs(title = "UMAP Visualization with Phenograph Clusters (Unscaled)") +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme_minimal()

# Plot UMAP with cluster colors
umap_plot_pops <- ggplot(umap_data, aes(x = V1, y = V2, color = pops_sub)) +
  geom_point() +
  labs(title = "UMAP Visualization with Populations (Unscaled)") +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme_minimal()

# Display the UMAP plot
print(umap_plot_pops)
print(umap_plot)
