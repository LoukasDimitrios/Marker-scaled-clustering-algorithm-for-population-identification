# Function to perform clustering using Seurat and calculate ARI
perform_clustering_seurat <- function(expr_data, resolution, algorithm, pops_sub, k_neighbors, p) {
  
  # Identify columns with empty names
  empty_col_indices <- which(colnames(expr_data) == "")
  
  # Rename empty columns
  colnames(expr_data)[empty_col_indices] <- c("V22", "V23", "V30")
  
  # Make the column names unique
  colnames(expr_data) <- make.unique(colnames(expr_data), sep = ".")
  
  
  matrix_data <- as.matrix(expr_data)
  matrix_data <- t(matrix_data)
  # Convert the matrix to a dgCMatrix
  dgCMatrix_data <- as(matrix_data, "dgCMatrix")
  
  # Create Seurat object
  data <- CreateSeuratObject(counts = dgCMatrix_data)
  
  # Preprocessing steps
  data <- NormalizeData(data, verbose = FALSE)
  data <- ScaleData(data, verbose = FALSE)
  data <- FindVariableFeatures(data, nfeatures = 55, verbose = FALSE)
  data <- RunPCA(data, npcs = 55, verbose = FALSE)
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


# Initialize variables to store best parameters and ARI for scaled
best_resolution_scaled <- 0
best_algorithm_scaled <- ""
best_ari_scaled <- -Inf
best_cluster_assignments_scaled <- NULL
best_k_neighbors_scaled <- 0
best_p_scaled <- 0

# Define parameter grids
param_grid_scaled <- expand.grid(
  algorithm = c("louvain", "leiden", "slm"),
  resolution = seq(0.1, 0.1, by = 0.1),
  k_neighbors = c(5),
  p = c(5)
  # k_neighbors = c(5, 10, 15, 20, 25),
  # p = c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 54)
)

# Grid search loop for scaled
iteration_counter_scaled <- 0
for (i in seq_len(nrow(param_grid_scaled))) {
  iteration_counter_scaled <- iteration_counter_scaled + 1
  cat("Iteration (Scaled):", iteration_counter_scaled, "/", nrow(param_grid_scaled), "\n")
  
  # Extract parameters for the current iteration
  params_scaled <- param_grid_scaled[i, ]

  # Perform clustering using Seurat and calculate ARI for Scaled
  result_scaled <- perform_clustering_seurat(expr_sub_exp, params_scaled$resolution, params_scaled$algorithm, pops_sub, params_scaled$k_neighbors, params_scaled$p)
  data_scaled <- result_scaled[[1]]
  cluster_assignments_scaled <- result_scaled[[2]]
  ari_scaled <- result_scaled[[3]]

  # Update best parameters if ARI is improved for Scaled
  if (ari_scaled > best_ari_scaled) {
    best_resolution_scaled <- params_scaled$resolution
    best_algorithm_scaled <- params_scaled$algorithm
    best_ari_scaled <- ari_scaled
    best_cluster_assignments_scaled <- cluster_assignments_scaled
    best_k_neighbors_scaled <- params_scaled$k_neighbors
    best_p_scaled <- params_scaled$p
  }
}

# Print the best parameters and ARI for Scaled
cat("Best Algorithm (Scaled):", best_algorithm_scaled, "\n")
cat("Best Resolution (Scaled):", best_resolution_scaled, "\n")
cat("Best k-Neighbors (Scaled):", best_k_neighbors_scaled, "\n")
cat("Best Number of PCs (p) (Scaled):", best_p_scaled, "\n")
cat("Best Adjusted Rand Index (Scaled):", best_ari_scaled, "\n")

# PCA - POPULATION without normalization (without division by PC1^2)
pca_result_scaled <- prcomp(expr_sub_exp, scale. = TRUE)  # Standardize the data

# Get PCA from Seurat object for Scaled
selected_pcs_scaled <- pca_result_scaled$x[, 1:best_p_scaled]

# UMAP for Scaled
umap_result_scaled <- umap(selected_pcs_scaled, batch = TRUE, n_threads = 32, verbose = FALSE, n_neighbors = best_k_neighbors_scaled, n_sgd_threads = "auto")

# SEURAT CLUSTERS - PLOT for Scaled
umap_data_clusters_scaled <- data.frame(V1 = umap_result_scaled[, 1], V2 = umap_result_scaled[, 2], Cluster = factor(best_cluster_assignments_scaled))
umap_plot_clusters_scaled <- ggplot(umap_data_clusters_scaled, aes(x = V1, y = V2, color = Cluster)) +
  geom_point() +
  labs(title = paste("UMAP - Seurat Clusters (Scaled) - Algorithm:", best_algorithm_scaled, ", Resolution:", best_resolution_scaled, ", k:", best_k_neighbors_scaled, ", P:", best_p_scaled)) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme_minimal()

# POPULATION - PLOT for Scaled
umap_data_population_scaled <- data.frame(V1 = umap_result_scaled[, 1], V2 = umap_result_scaled[, 2])
umap_data_population_scaled$pops_sub <- pops_sub
umap_plot_population_scaled <- ggplot(umap_data_population_scaled, aes(x = V1, y = V2, color = pops_sub)) +
  geom_point() +
  labs(title = paste("UMAP - Populations (Scaled) - Algorithm:", best_algorithm_scaled, ", Resolution:", best_resolution_scaled, ", k:", best_k_neighbors_scaled, ", P:", best_p_scaled)) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme_minimal()

# Display plots for Scaled
print(umap_plot_clusters_scaled)
print(umap_plot_population_scaled)
