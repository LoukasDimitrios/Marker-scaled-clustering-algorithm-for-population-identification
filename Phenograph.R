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
