# Run PCA for Log Fold Change groups
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

# Save the PCA plot as a PNG file
#ggsave("pca_plot_unscaled_COUNTS.png", pca_plot, width = 8, height = 6, units = "in", dpi = 300)



# Load required libraries
library(Seurat)
library(uwot)
library(factoextra)
library(ggplot2)
library(fpc)  # Ensure the fpc package is installed
library(pdfCluster)

set.seed(42)

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
best_resolution_scaled_fc <- 0
best_algorithm_scaled_fc <- ""
best_ari_scaled_fc <- -Inf
best_cluster_assignments_scaled_fc <- NULL
best_k_neighbors_scaled_fc <- 0
best_p_scaled_fc <- 0

# Define parameter grids
param_grid_scaled_fc <- expand.grid(
  algorithm = c("louvain", "leiden", "slm"),
  resolution = seq(0.1, 1, by = 0.1),
  k_neighbors = c(5, 10, 15, 20, 25),
  p = c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 54)
)

# Grid search loop for scaled
iteration_counter_scaled_fc <- 0
for (i in seq_len(nrow(param_grid_scaled_fc))) {
  iteration_counter_scaled_fc <- iteration_counter_scaled_fc + 1
  cat("Iteration (Scaled):", iteration_counter_scaled_fc, "/", nrow(param_grid_scaled_fc), "\n")
  
  # Extract parameters for the current iteration
  params_scaled_fc <- param_grid_scaled_fc[i, ]

  # Perform clustering using Seurat and calculate ARI for Scaled
  result_scaled_fc <- perform_clustering_seurat(expr_sub_exp, params_scaled_fc$resolution, params_scaled_fc$algorithm, pops_sub, params_scaled_fc$k_neighbors, params_scaled_fc$p)
  data_scaled_fc <- result_scaled_fc[[1]]
  cluster_assignments_scaled_fc <- result_scaled_fc[[2]]
  ari_scaled_fc <- result_scaled_fc[[3]]

  # Update best parameters if ARI is improved for Scaled
  if (ari_scaled_fc > best_ari_scaled_fc) {
    best_resolution_scaled_fc <- params_scaled_fc$resolution
    best_algorithm_scaled_fc <- params_scaled_fc$algorithm
    best_ari_scaled_fc <- ari_scaled_fc
    best_cluster_assignments_scaled_fc <- cluster_assignments_scaled_fc
    best_k_neighbors_scaled_fc <- params_scaled_fc$k_neighbors
    best_p_scaled_fc <- params_scaled_fc$p
  }
}

# Print the best parameters and ARI for Scaled
cat("Best Algorithm (Scaled):", best_algorithm_scaled_fc, "\n")
cat("Best Resolution (Scaled):", best_resolution_scaled_fc, "\n")
cat("Best k-Neighbors (Scaled):", best_k_neighbors_scaled_fc, "\n")
cat("Best Number of PCs (p) (Scaled):", best_p_scaled_fc, "\n")
cat("Best Adjusted Rand Index (Scaled):", best_ari_scaled_fc, "\n")

# PCA - POPULATION without normalization (without division by PC1^2)
pca_result_scaled_fc <- prcomp(expr_sub_exp, scale. = TRUE)  # Standardize the data

# Get PCA from Seurat object for Scaled
selected_pcs_scaled_fc <- pca_result_scaled_fc$x[, 1:best_p_scaled_fc]

# UMAP for Scaled
umap_result_scaled_fc <- umap(selected_pcs_scaled_fc, batch = TRUE, n_threads = 32, verbose = FALSE, n_neighbors = best_k_neighbors_scaled_fc, n_sgd_threads = "auto")

# SEURAT CLUSTERS - PLOT for Scaled
umap_data_clusters_scaled_fc <- data.frame(V1 = umap_result_scaled_fc[, 1], V2 = umap_result_scaled_fc[, 2], Cluster = factor(best_cluster_assignments_scaled_fc))
umap_plot_clusters_scaled_fc <- ggplot(umap_data_clusters_scaled_fc, aes(x = V1, y = V2, color = Cluster)) +
  geom_point() +
  labs(title = paste("UMAP - Seurat Clusters (Scaled) - Algorithm:", best_algorithm_scaled_fc, ", Resolution:", best_resolution_scaled_fc)) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme_minimal()

# POPULATION - PLOT for Scaled
umap_data_population_scaled_fc <- data.frame(V1 = umap_result_scaled_fc[, 1], V2 = umap_result_scaled_fc[, 2])
umap_data_population_scaled_fc$pops_sub <- pops_sub
umap_plot_population_scaled_fc <- ggplot(umap_data_population_scaled_fc, aes(x = V1, y = V2, color = pops_sub)) +
  geom_point() +
  labs(title = paste("UMAP - Populations (Scaled) - Algorithm:", best_algorithm_scaled_fc, ", Resolution:", best_resolution_scaled_fc)) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme_minimal()


# Save each variable in an RData file
save(
  best_algorithm_scaled_fc,
  best_resolution_scaled_fc,
  best_k_neighbors_scaled_fc,
  best_p_scaled_fc,
  best_ari_scaled_fc,
  file = "best_params_scaled_fc.RData"
)

# Display plots for Scaled
print(umap_plot_clusters_scaled_fc)
print(umap_plot_population_scaled_fc)
