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

# Initialize variables to store best parameters and ARI
best_resolution <- 0
best_algorithm <- ""
best_ari <- -Inf
best_cluster_assignments <- NULL
best_k_neighbors <- 0
best_p <- 0

# Define parameter grids
param_grid <- expand.grid(
  algorithm = c("louvain", "leiden", "slm"),
  resolution = seq(0.1, 1, by = 0.1),
  k_neighbors = c(5, 10, 15, 20, 25),
  p = c(5, 10, 15, 20, 23)
)

# Grid search loop
iteration_counter <- 0
for (i in seq_len(nrow(param_grid))) {
  iteration_counter <- iteration_counter + 1
  cat("Iteration:", iteration_counter, "/", nrow(param_grid), "\n")
  
  # Extract parameters for the current iteration
  params <- param_grid[i, ]

  # Perform clustering using Seurat and calculate ARI
  result <- perform_clustering_seurat(expr_sub, params$resolution, params$algorithm, pops_sub, params$k_neighbors, params$p)
  data <- result[[1]]
  cluster_assignments <- result[[2]]
  ari <- result[[3]]

  # Update best parameters if ARI is improved
  if (ari > best_ari) {
    best_resolution <- params$resolution
    best_algorithm <- params$algorithm
    best_ari <- ari
    best_cluster_assignments <- cluster_assignments
    best_k_neighbors <- params$k_neighbors
    best_p <- params$p
  }
}

# Print the best parameters and ARI
cat("Best Algorithm:", best_algorithm, "\n")
cat("Best Resolution:", best_resolution, "\n")
cat("Best k-Neighbors:", best_k_neighbors, "\n")
cat("Best Number of PCs (p):", best_p, "\n")
cat("Best Adjusted Rand Index:", best_ari, "\n")

# PCA - POPULATION without normalization (without division by PC1^2)
pca_result <- prcomp(expr_sub, scale. = TRUE)  # Standardize the data

# Get PCA from Seurat object
selected_pcs <- pca_result$x[, 1:best_p]

# UMAP
umap_result <- umap(selected_pcs, batch = TRUE, n_threads = 32, verbose = FALSE, n_neighbors = best_k_neighbors, n_sgd_threads = "auto")

# SEURAT CLUSTERS - PLOT
umap_data_clusters <- data.frame(V1 = umap_result[, 1], V2 = umap_result[, 2], Cluster = factor(best_cluster_assignments))
umap_plot_clusters <- ggplot(umap_data_clusters, aes(x = V1, y = V2, color = Cluster)) +
  geom_point() +
  labs(title = paste("UMAP - Seurat Clusters (Unscaled) - Algorithm:", best_algorithm, ", Resolution:", best_resolution)) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme_minimal()

# POPULATION - PLOT
umap_data_population <- data.frame(V1 = umap_result[, 1], V2 = umap_result[, 2])
umap_data_population$pops_sub <- pops_sub
umap_plot_population <- ggplot(umap_data_population, aes(x = V1, y = V2, color = pops_sub)) +
  geom_point() +
  labs(title = paste("UMAP - Populations (Unscaled) - Algorithm:", best_algorithm, ", Resolution:", best_resolution)) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme_minimal()

# Display plots
print(umap_plot_clusters)
print(umap_plot_population)
