Rphenograph2 <- function(data, k=20, resolution = 10){
    if(is.data.frame(data))
        data <- as.matrix(data)
    
    if(!is.matrix(data))
        stop("Wrong input data, should be a data frame of matrix!")
    
    if(k<1){
        stop("k must be a positive integer!")
    }else if (k > nrow(data)-2){
        stop("k must be smaller than the total number of points!")
    }
    
    message("Run Rphenograph starts:","\n", 
        "  -Input data of ", nrow(data)," rows and ", ncol(data), " columns","\n",
        "  -k is set to ", k)
    
    cat("  Finding nearest neighbors...")
    t1 <- system.time(neighborMatrix <- find_neighbors(data, k=k+1)[,-1])
    cat("DONE ~",t1[3],"s\n", " Compute jaccard coefficient between nearest-neighbor sets...")
    t2 <- system.time(links <- jaccard_coeff(neighborMatrix))

    cat("DONE ~",t2[3],"s\n", " Build undirected graph from the weighted links...")
    links <- links[links[,1]>0, ]
    relations <- as.data.frame(links)
    colnames(relations)<- c("from","to","weight")
    t3 <- system.time(g <- graph.data.frame(relations, directed=FALSE))
    
    # Other community detection algorithms: 
    #    cluster_walktrap, cluster_spinglass, 
    #    cluster_leading_eigen, cluster_edge_betweenness, 
    #    cluster_fast_greedy, cluster_label_prop  
    cat(" Run louvain clustering on the graph with resolution =", resolution, "...")
    t4 <- system.time(community <- cluster_louvain(g, resolution = resolution))
    cat("DONE ~", t4[3], "s\n")
    
    message("Run Rphenograph DONE, totally takes ", sum(c(t1[3],t2[3],t3[3],t4[3])), "s.")
    cat("  Return a community class\n  -Modularity value:", modularity(community),"\n")
    cat("  -Number of clusters:", length(unique(membership(community))))
    
    return(list(g, community))
}



find_neighbors <- function(data, k){
    nearest <- nn2(data, data, k, searchtype = "standard")
    return(nearest[[1]])
}

jaccard_coeff <- function(idx) {
    .Call('Rphenograph_jaccard_coeff', PACKAGE = 'Rphenograph', idx)
}
