#' Find connections of a subset of scaffolds in a larger set of scaffolds
#'
#' @param mm (\emph{required}) A dataframe loaded with \code{\link{mmload}} in which to find connections from \code{scaffolds}.
#' @param scaffolds (\emph{required}) The scaffolds from which to find connections in \code{mm}. Must be a vector of scaffold names or a dataframe with scaffold names in the first column.
#' @param network (\emph{required}) Paired-end or mate-pair connections between scaffolds in long format. The first and second columns must contain all connected scaffold pairs and the third column the number of connections. 
#' @param min_connections Filter all scaffold pairs with equal to or less than this number of connections before the extraction. (\emph{Default: } \code{2})
#' @param include_connections The connections to include. One of the following:
#' \describe{
#'   \item{\code{"direct"}}{Extract only scaffolds from \code{mm} that are directly connected to any of the scaffolds in \code{scaffolds}. (\emph{default})}
#'   \item{\code{"all"}}{Same as \code{"direct"}, except that any further connections to connected scaffolds are also included, continuing until no further connections are found.}
#' }
#' @export
#' 
#' @import igraph
#' @import dplyr
#'
mmexpand_network <- function(mm,
                             scaffolds, 
                             network,
                             min_connections = 2,
                             include_connections = "direct") {
  if(is.data.frame(scaffolds)) {
    scaffolds <- scaffolds[[1]]
  } else if(!is.vector(scaffolds) | !is.data.frame(scaffolds)) {
    stop("Scaffolds must be provided either as a vector or as a dataframe, where the first column contains the scaffold names.")
  }
  
  if (include_connections == "direct"){
    ns <- dplyr::filter(network, (network[[1]] %in% scaffolds | network[[2]] %in% scaffolds) & network[[3]] >= min_connections)
    out <- dplyr::filter(mm, mm[[1]] %in% ns$scaffold1 | mm[[1]] %in% ns$scaffold2 | mm[[1]] %in% scaffolds)
  }
  
  if (include_connections == "all"){
    snetwork <- dplyr::filter(network, connections >= min_connections)
    g <- igraph::graph.data.frame(snetwork, directed = F)
    g.clust <- igraph::clusters(g)
    clusters <- cbind.data.frame(igraph::V(g)$name,g.clust$membership)
    colnames(clusters) <- c("scaffold", "cluster")  
    ext.clusters <- dplyr::filter(clusters, clusters$scaffold %in% scaffolds)
    ext.scaffolds <- dplyr::filter(clusters, clusters$cluster %in% ext.clusters$cluster)
    out <- dplyr::filter(mm, mm[[1]] %in% scaffolds | mm[[1]] %in% ext.scaffolds$scaffold)
  }
  return(out)  
}