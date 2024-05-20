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

pca_plot + theme_minimal()
# Normalization and PCA
expr_sub_exp <- NULL
for(i in 1:length(groups)) {
  pca <- prcomp(expr_sub[,groups[[i]]], scale. = TRUE)
  expr_sub_exp <- cbind(expr_sub_exp, as.matrix(expr_sub[,groups[[i]]]/pca$sdev[1]^2))
}

# PCA scaled
pca_result_scaled <- prcomp(expr_sub_exp, scale. = TRUE)

scree_plot_scaled <- fviz_eig(pca_result_scaled,
                              addlabels = TRUE,
                              title = "Scree plot SCALED")

# PCA plot with population names in legends using the default color palette
pca_plot_scaled <- fviz_pca_ind(pca_result_scaled, 
                                geom.ind = "point", 
                                col.ind = pops_sub,  # Color points by population
                                addEllipses = TRUE, 
                                ellipse.type = "confidence",
                                legend.title = "Population",
                                pointshape = 19,
                                fill = pops_sub,
                                title = "PCA visualization SCALED")  # Set fill aesthetic to population names


# UMAP - POPULATION plot UNSCALED
#umap_result <- umap(expr_sub)

# UMAP - POPULATION plot UNSCALED
umap_result <- umap(expr_sub, batch = TRUE, n_threads = 32, verbose = TRUE, n_neighbors = 20, n_sgd_threads = "auto")

# Create a data frame from UMAP results
umap_data <- data.frame(V1 = umap_result[, 1], V2 = umap_result[, 2])

# Add population information to the data frame
umap_data$pops_sub <- pops_sub

# UMAP plot with consistent colors and aesthetics
umap_plot <- ggplot(umap_data, aes(x = V1, y = V2)) +
  geom_point(aes(color = pops_sub)) +  # Color points by population
  labs(x = "UMAP 1", y = "UMAP 2", color = "Populations") +  # Customize axis labels and legend title
  #ggtitle("UMAP Visualization UNSCALED") +  # Add a title
  guides(color = guide_legend(override.aes = list(size = 4))) +  # Adjust legend point size
  theme_minimal()  # Set a minimal theme to match PCA plot style


# UMAP - POPULATION plot SCALED
#umap_result <- umap(expr_sub_exp)

# UMAP - POPULATION plot UNSCALED
umap_result <- umap(expr_sub_exp, batch = TRUE, n_threads = 32, verbose = TRUE, n_neighbors = 15, n_sgd_threads = "auto")

# Create a data frame from UMAP results
umap_data <- data.frame(V1 = umap_result[, 1], V2 = umap_result[, 2])

# Add population information to the data frame
umap_data$pops_sub <- pops_sub

# UMAP plot with consistent colors
umap_plot_scaled <- ggplot(umap_data, aes(x = V1, y = V2)) +
  geom_point(aes(color = pops_sub)) +  # Color points by population
  labs(x = "UMAP 1", y = "UMAP 2") +  # Customize axis labels
  #ggtitle("UMAP Visualization SCALED") +  # Add a title
  guides(color = guide_legend(override.aes = list(size = 4)))  +  # Adjust legend point size
  theme_minimal()  # Set a minimal theme to match PCA plot style


# Display the scree plot
print(scree_plot)

# Display the scree plot
print(scree_plot_scaled)

# Display the PCA plot
print(pca_plot)

# Display the PCA plot
print(pca_plot_scaled) 

# Display the UMAP plot
print(umap_plot)

# Display the UMAP plot
print(umap_plot_scaled) 
