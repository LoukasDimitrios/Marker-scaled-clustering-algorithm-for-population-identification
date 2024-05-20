### PCA - UNSCALED vs UNSCALED


# PCA - POPULATION without normalization (without division by PC1^2)
pca_result_scaled <- prcomp(expr_sub_exp, scale. = TRUE)  # Standardize the data

# Create a scree plot to visualize the variance explained by each PC
scree_plot <- fviz_eig(pca_result_scaled,
                       addlabels = TRUE,
                              title = "Scree plot SCALED")

# PCA plot with population names in legends using the default color palette
pca_plot <- fviz_pca_ind(pca_result_scaled, 
                         geom.ind = "point", 
                         col.ind = pops_sub,  # Color points by population
                         addEllipses = TRUE, 
                         ellipse.type = "confidence",
                         legend.title = "Population",
                         pointshape = 19,
                         fill = pops_sub,
                         title = "PCA visualization SCALED")  # Set fill aesthetic to population names

# Display the scree plot
print(scree_plot)

# Display the PCA plot
print(pca_plot)


result_scaled <- Rphenograph2(pca_result_scaled$x[, 1:best_p_scaled], k = best_k_neighbors_scaled, resolution = best_resolution_scaled)

# Get cluster assignments from Phenograph
cluster_assignments_scaled <- result_scaled[[2]]$membership

# Perform UMAP on the first 5 PCs
#umap_result <- umap(pca_result$x[, ])
umap_result_scaled <- umap(pca_result_scaled$x[, ], batch = TRUE, n_threads = 32, verbose = TRUE, n_neighbors = 15, n_sgd_threads = "auto")

# Create a data frame with UMAP results and cluster assignments
umap_data_scaled <- data.frame(V1 = umap_result_scaled[, 1], V2 = umap_result_scaled[, 2], Cluster = factor(cluster_assignments_scaled))

# Plot UMAP with cluster colors
umap_plot_scaled <- ggplot(umap_data_scaled, aes(x = V1, y = V2, color = Cluster)) +
  geom_point() +
  labs(title = "UMAP Visualization with Phenograph Clusters (scaled)") +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme_minimal()

# Plot UMAP with cluster colors
umap_plot_scaled_pops <- ggplot(umap_data_scaled, aes(x = V1, y = V2, color = pops_sub)) +
  geom_point() +
  labs(title = "UMAP Visualization with Populations (scaled)") +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme_minimal()

# Display the UMAP plot
print(umap_plot_scaled_pops)
print(umap_plot_scaled)
