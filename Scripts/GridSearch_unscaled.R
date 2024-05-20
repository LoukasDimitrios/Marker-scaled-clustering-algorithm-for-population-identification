# Function to perform clustering using Seurat and calculate ARI
perform_clustering_seurat <- function(expr_data, resolution, algorithm, pops_sub, k_neighbors, p) {
  
  matrix_data <- as.matrix(expr_data)
  matrix_data <- t(matrix_data)
  # Convert the matrix to a dgCMatrix
  dgCMatrix_data <- as(matrix_data, "dgCMatrix")
  
  # Create Seurat object
  data <- CreateSeuratObject(counts = dgCMatrix_data)
  
  # Preprocessing steps
  data <- NormalizeData(data, verbose = FALSE)
  data <- ScaleData(data, verbose = FALSE)
  data <- FindVariableFeatures(data, nfeatures = 24, verbose = FALSE)
  data <- RunPCA(data, npcs = 24, verbose = FALSE)
  data <- FindNeighbors(data, compute.SNN = TRUE, dims = 1:p, k.param = k_neighbors, verbose = FALSE)

  # Perform clustering based on the specified algorithm
  if (algorithm == "louvain") {
    data <- FindClusters(data, resolution = resolution, algorithm = 1, reduction.method = "pca", verbose = FALSE)
  } else if (algorithm == "leiden") {
    data <- FindClusters(data, resolution = resolution, algorithm = 4, reduction.method = "pca", verbose = FALSE)
  } else if (algorithm == "slm") {
    data <- FindClusters(data, resolution = resolution, algorithm = 3, reduction.method = "pca", verbose = FALSE)
  } else {
    stop("Invalid algorithm. Choose 'louvain', 'leiden', or 'slm'")
  }

  # Get cluster assignments from Seurat
  cluster_assignments <- data@meta.data$seurat_clusters

  pops_sub_factor <- as.factor(pops_sub)
  # Use cluster.stats to calculate adjusted Rand index
  ari <- adj.rand.index(cluster_assignments, pops_sub_factor)

  return(list(data, cluster_assignments, ari))
}

# Initialize variables to store best parameters and ARI for unscaled
best_resolution_unscaled <- 0
best_algorithm_unscaled <- ""
best_ari_unscaled <- -Inf
best_cluster_assignments_unscaled <- NULL
best_k_neighbors_unscaled <- 0
best_p_unscaled <- 0

# Define parameter grids
param_grid_unscaled <- expand.grid(
  algorithm = c("louvain", "leiden", "slm"),
  resolution = seq(0.1, 0.1, by = 0.1),
  # k_neighbors = c(5, 10, 15, 20, 25),
  # p = c(5, 10, 15, 20, 23),
  k_neighbors = c(5),
  p = c(5)
)

# Grid search loop for unscaled
iteration_counter_unscaled <- 0
for (i in seq_len(nrow(param_grid_unscaled))) {
  iteration_counter_unscaled <- iteration_counter_unscaled + 1
  cat("Iteration (Unscaled):", iteration_counter_unscaled, "/", nrow(param_grid_unscaled), "\n")
  
  # Extract parameters for the current iteration
  params_unscaled <- param_grid_unscaled[i, ]

  # Perform clustering using Seurat and calculate ARI for unscaled
  result_unscaled <- perform_clustering_seurat(expr_sub, params_unscaled$resolution, params_unscaled$algorithm, pops_sub, params_unscaled$k_neighbors, params_unscaled$p)
  data_unscaled <- result_unscaled[[1]]
  cluster_assignments_unscaled <- result_unscaled[[2]]
  ari_unscaled <- result_unscaled[[3]]

  # Update best parameters if ARI is improved for unscaled
  if (ari_unscaled > best_ari_unscaled) {
    best_resolution_unscaled <- params_unscaled$resolution
    best_algorithm_unscaled <- params_unscaled$algorithm
    best_ari_unscaled <- ari_unscaled
    best_cluster_assignments_unscaled <- cluster_assignments_unscaled
    best_k_neighbors_unscaled <- params_unscaled$k_neighbors
    best_p_unscaled <- params_unscaled$p
  }
}

# Print the best parameters and ARI for unscaled
cat("Best Algorithm (Unscaled):", best_algorithm_unscaled, "\n")
cat("Best Resolution (Unscaled):", best_resolution_unscaled, "\n")
cat("Best k-Neighbors (Unscaled):", best_k_neighbors_unscaled, "\n")
cat("Best Number of PCs (p) (Unscaled):", best_p_unscaled, "\n")
cat("Best Adjusted Rand Index (Unscaled):", best_ari_unscaled, "\n")

# PCA - POPULATION without normalization (without division by PC1^2)
pca_result_unscaled <- prcomp(expr_sub, scale. = TRUE)  # Standardize the data

# Get PCA from Seurat object for unscaled
selected_pcs_unscaled <- pca_result_unscaled$x[, 1:best_p_unscaled]

# UMAP for unscaled
umap_result_unscaled <- umap(selected_pcs_unscaled, batch = TRUE, n_threads = 32, verbose = FALSE, n_neighbors = best_k_neighbors_unscaled, n_sgd_threads = "auto")

# SEURAT CLUSTERS - PLOT for unscaled
umap_data_clusters_unscaled <- data.frame(V1 = umap_result_unscaled[, 1], V2 = umap_result_unscaled[, 2], Cluster = factor(best_cluster_assignments_unscaled))
umap_plot_clusters_unscaled <- ggplot(umap_data_clusters_unscaled, aes(x = V1, y = V2, color = Cluster)) +
  geom_point() +
  labs(title = paste("UMAP - Seurat Clusters (Unscaled) - Algorithm:", best_algorithm_unscaled, ", Resolution:", best_resolution_unscaled)) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme_minimal()

# POPULATION - PLOT for unscaled
umap_data_population_unscaled <- data.frame(V1 = umap_result_unscaled[, 1], V2 = umap_result_unscaled[, 2])
umap_data_population_unscaled$pops_sub <- pops_sub
umap_plot_population_unscaled <- ggplot(umap_data_population_unscaled, aes(x = V1, y = V2, color = pops_sub)) +
  geom_point() +
  labs(title = paste("UMAP - Populations (Unscaled) - Algorithm:", best_algorithm_unscaled, ", Resolution:", best_resolution_unscaled)) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme_minimal()

# Display plots for unscaled
print(umap_plot_clusters_unscaled)
print(umap_plot_population_unscaled)


# Save each variable in an RData file
save(
  best_algorithm_unscaled,
  best_resolution_unscaled,
  best_k_neighbors_unscaled,
  best_p_unscaled,
  best_ari_unscaled,
  file = "best_params_unscaled.RData"
)
