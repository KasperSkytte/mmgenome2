#' Title
#'
#' @param scaffolds 
#' @param network 
#' @param original_data 
#' @param min_connections 
#' @param include_connections 
#'
#' @return
#' @export
#' 
#' @import igraph
#' @import dplyr
#'
#' @examples
mmextract_network <- function(scaffolds, 
                              network,
                              original_data,
                              min_connections = 2,
                              include_connections = "direct") {
  if(is.data.frame(scaffolds)) {
    scaffolds <- scaffolds[[1]]
  } else if(!is.vector(scaffolds) | !is.data.frame(scaffolds)) {
    stop("Scaffolds must be provided either as a vector or as a dataframe, where the first column contains the scaffold names.")
  }
  
  if (include_connections == "direct"){
    ns <- dplyr::filter(network, (network[[1]] %in% scaffolds | network[[2]] %in% scaffolds) & network$connections >= min_connections)
    out <- dplyr::filter(original_data, original_data[[1]] %in% ns$scaffold1 | original_data[[1]] %in% ns$scaffold2 | original_data[[1]] %in% scaffolds)
  }
  
  if (include_connections == "all"){
    snetwork <- dplyr::filter(network, connections >= min_connections)
    g <- igraph::graph.data.frame(snetwork, directed = F)
    g.clust <- igraph::clusters(g)
    clusters <- cbind.data.frame(igraph::V(g)$name,g.clust$membership)
    colnames(clusters) <- c("scaffold", "cluster")  
    ext.clusters <- dplyr::filter(clusters, clusters$scaffold %in% scaffolds)
    ext.scaffolds <- dplyr::filter(clusters, clusters$cluster %in% ext.clusters$cluster)
    out <- dplyr::filter(original_data, original_data[[1]] %in% scaffolds | original_data[[1]] %in% ext.scaffolds$scaffold)
  }
  return(out)  
}