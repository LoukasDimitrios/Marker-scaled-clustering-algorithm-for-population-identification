groups <- list()
for(p in unique(pops)) {
  fc <- c()
  for(i in 1:ncol(expr)) {
    fc <- c(fc,log2(mean(expr[pops == p, i])/mean(expr[!pops == p, i])))
  }
  names(fc) <- colnames(expr)
  groups[[p]] <- names(fc[which(fc > 0)])
}

# Define the groups
groups <- list()

for (p in unique(pops_sub)) {
  fc <- c()
  for (i in 1:ncol(expr_sub)) {
    mean_p <- mean(expr_sub[pops_sub == p, i])
    mean_not_p <- mean(expr_sub[!pops_sub == p, i])
    
    # print(paste("Population:", p, "Gene Index:", i))
    # print(paste("Mean Expression in Population:", mean_p))
    # print(paste("Mean Expression in Other Population:", mean_not_p))
    
    fc <- c(fc, log2(mean_p / mean_not_p))
  }
  
  names(fc) <- colnames(expr_sub)
  groups[[p]] <- names(fc[which(fc > 0)])
}

# Print the resulting groups
print(groups)

expr_sub_exp <- NULL
for(i in 1:length(groups)) {
  pca <- prcomp(expr_sub[,groups[[i]]], scale. = TRUE)
  expr_sub_exp <- cbind(expr_sub_exp, as.matrix(expr_sub[,groups[[i]]]/pca$sdev[1]^2))
}
