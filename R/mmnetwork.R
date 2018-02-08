#' Title
#'
#' @param mm 
#' @param network 
#' @param n_connections 
#' @param color_by 
#' @param locator 
#' @param selection 
#' @param highlight_labels 
#' @param highlight_color 
#' @param color_scale_log10 
#' @param links_scale 
#' @param scaffold_labels 
#' @param print_nolinks 
#' @param seed 
#'
#' @export
#' 
#' @import ggplot2
#' @import igraph
#' @import dplyr
#' @import sp
mmnetwork <- function(mm,
                      network, 
                      n_connections = 2,
                      color_by = NULL,
                      locator = FALSE,
                      selection = NULL,
                      highlight_labels = NULL,
                      highlight_color = "darkred",
                      color_scale_log10 = FALSE,
                      links_scale = 1,
                      scaffold_labels = FALSE,
                      print_nolinks = FALSE,
                      seed = .mm_Random.seed) {
  #Checks and error messages before anything else
  if(isTRUE(locator) & !is.null(selection))
    stop("Using the locator and highlighting a selection at the same time is not supported.")
  
  ## Subset the network
  snetwork <- dplyr::filter(network, network[[1]] %in% mm[[1]] & network[[2]] %in% mm[[1]] & network[["connections"]] >= n_connections)
  
  ## Convert to graph 
  g <- igraph::graph.data.frame(snetwork, directed = F)
  
  ## Calculate a layout
  if (!is.null(seed)) set.seed(seed)
  layout <- igraph::layout_with_fr(g)

  ## Extract layout coordinates
  gpoints <- merge(data.frame("scaffold" = igraph::V(g)$name, "x" = layout[,1], "y" = layout[,2]), mm, by = 1)
  
  ## Extract link coordinates 
  links <- merge(snetwork, gpoints[,1:3], by.x = "scaffold1", by.y = "scaffold")
  links <- merge(links, gpoints[,1:3], by.x = "scaffold2", by.y = "scaffold")
  colnames(links)[4:7] <- c("x", "y", "xend", "yend")
  
  p <- ggplot(data = gpoints, 
              aes_string(x = "x",
                         y = "y", 
                         size = "length")) +
    geom_segment(data = links,
                 aes_string(x = "x",
                            y = "y",
                            xend = "xend",
                            yend = "yend"),
                 color = "darkgrey", 
                 size = log10(links[["connections"]])*links_scale) +
    scale_size_area(name= "Scaffold length", max_size=20)
  
  if(!is.null(color_by)) {
    if(is.factor(mm[[color_by]]) | is.character(mm[[color_by]])) {
      p <- p +
        geom_point(alpha = 0.1, color = "black") +
        geom_point(data = subset(gpoints, gpoints[[color_by]] != "NA"), shape = 1) +
        guides(colour = guide_legend(override.aes = list(alpha = 1, size = 5, shape = 19)))
    } else if(is.numeric(mm[[color_by]])) {
      p <- p +
        geom_point(alpha = 0.7) +
        ifelse(isTRUE(color_scale_log10), 
               scale_colour_gradientn(colours = c("red", "green", "blue"), trans = "log10"), 
               scale_colour_gradientn(colours = c("red", "green", "blue")))
    }
  } else {
    p <- p +
      geom_point(alpha = 0.1, color = "black")
  }
  
  if (isTRUE(scaffold_labels)) {
    p <- p + 
      geom_text(label = gpoints[[1]], 
                size = 4,
                color = "black")
  }
  
  p <- p + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.key = element_blank())
  
  ### Highlight selected scaffolds
  if (!is.null(highlight_labels)) {
    if(is.data.frame(highlight_labels)) {
      highlight_labels <- highlight_labels[[1]]
    } 
    scaffolds <- as.character(highlight_labels)
    d <- dplyr::filter(gpoints, scaffold %in% scaffolds)
    p <- p +
      geom_text(data = d,
                color = highlight_color,
                size = 4,
                label = d[["scaffold"]])
  }
  
  if(isTRUE(print_nolinks)){
    nolinks <- mm[["scaffold"]][!(mm[["scaffold"]] %in% gpoints$scaffold)]
    cat("The following scaffolds have no links to other scaffolds:\n")
    cat(paste0(nolinks, collapse = ", "))
  }
  
  ##### Locator and selection #####
  if(isTRUE(locator)) {
    points <- mmlocator(p)
  }
  if(isTRUE(locator) | !is.null(selection)) {
    if(!isTRUE(locator) & !is.null(selection)) {
      points <- selection
    }
    p$selection <- points
    in_polygon <- sp::point.in.polygon(point.x = p$data[[colnames(points)[1]]],
                                       point.y = p$data[[colnames(points)[2]]],
                                       pol.x = points[[1]],
                                       pol.y = points[[2]],
                                       mode.checked = TRUE)
    scaffolds <- as.character(p$data$scaffold[which(in_polygon > 0)])
    p$data_in_selection <- mm[which(mm[[1]] %in% scaffolds),]
    p <- p + 
      geom_point(data = points,
                 aes_(x = points[[1]],
                      y = points[[2]]), 
                 color = "black",
                 size = 2,
                 inherit.aes = FALSE, 
                 na.rm = TRUE) +
      geom_polygon(data = points,
                   aes_(x = points[[1]],
                        y = points[[2]]),
                   fill = NA,
                   size = 0.5,
                   lty = 2,
                   color = "darkred",
                   inherit.aes = FALSE, 
                   na.rm = TRUE)
  }
  return(p)
}
