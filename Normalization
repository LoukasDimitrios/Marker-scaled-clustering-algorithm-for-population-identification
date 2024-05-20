normalization <- function(expr_sub, groups) {
  expr_sub_exp <- NULL
  
  # Loop through each group
  for(i in 1:length(groups)) {
    # Perform PCA
    pca <- prcomp(expr_sub[, groups[[i]]], scale. = TRUE)
    
    # Normalize and add to expr_sub_exp
    normalized_data <- as.matrix(expr_sub[, groups[[i]]]) / pca$sdev[1]^2
    expr_sub_exp <- cbind(expr_sub_exp, normalized_data)
  }
  
  return(expr_sub_exp)
}
