#' @title Visualise metagenomes in various ways
#'
#' @param mm (\emph{required}) A dataframe loaded with \code{\link{mmload}}.
#' @param x (\emph{required}) The variable from \code{mm} to plot on the first axis.
#' @param y (\emph{required}) The variable from \code{mm} to plot on the second axis.
#' @param min_length Remove scaffolds with a length at or below this value before plotting. (\emph{Default: } \code{0}) 
#' @param locator (\emph{Logical}) When \code{TRUE}, left-clicks in the plot are captured and the exact x/y-coordinates of the mouse clicks are returned. These coordinates can be used to highlight a selection of scaffolds in the plot, and also be used with \code{\link{mmextract}} to extract all scaffolds within the selection from the data. (\emph{Default: } \code{FALSE})
#' @param selection A 2-column dataframe with the x and y coordinates of points with which to draw a polygon onto the plot to highlight a selected region. A selection can be obtained by using the locator feature (by \code{locator = TRUE}). (\emph{Default: } \code{NULL})
#' @param x_scale Log10-scale (\code{"log10"}) or square-root scale \code{"sqrt"} the x axis. (\emph{Default: } \code{NULL})
#' @param x_limits Axis limits of the x axis. (\emph{Default: } \code{NULL})
#' @param y_scale Log10-scale (\code{"log10"}) or square-root scale \code{"sqrt"} the y axis. (\emph{Default: } \code{NULL})
#' @param y_limits Axis limits of the y axis. (\emph{Default: } \code{NULL})
#' @param color_by Color the scaffolds by a variable in \code{mm}. (\emph{Default: } \code{NULL})
#' @param alpha The transparancy of the scaffold points, where 0 is invisible and 1 is opaque. (\emph{Default: } \code{0.1})
#' @param scaffold_labels Add labels of a selection of scaffolds by providing either a character vector of scaffold names, or a dataframe with scaffold names in the first column. If set to \code{TRUE} then \emph{all} scaffolds will be labelled. (\emph{Default: } \code{FALSE})
#' @param fixed_size A fixed size for all scaffolds if set. If \code{NULL} then the scaffolds are scaled by length. (\emph{Default: } \code{NULL})
#' @param size_scale A factor to scale the sizes of the scaffolds plotted. Only applies when \code{fixed_size} is set to \code{NULL} and the scaffolds are scaled by length. (\emph{Default: } \code{1})
#' @param shared_genes (\emph{Logical}) If \code{TRUE}, lines will be drawn between scaffolds with any shared gene(s). (\emph{Default: } \code{FALSE})
#'
#' @export
#' 
#' @import ggplot2
#' @import tidyr
#' @import tibble
mmplot <- function(mm,
                   x,
                   y, 
                   min_length = 0,
                   color_by = NULL,
                   locator = FALSE,
                   selection = NULL,
                   x_scale = NULL,
                   x_limits = NULL,
                   y_scale = NULL,
                   y_limits = NULL,
                   alpha = 0.1,
                   scaffold_labels = FALSE,
                   fixed_size = NULL,
                   size_scale = 1,
                   shared_genes = FALSE) {
  #Checks and error messages before anything else
  if(isTRUE(locator) & !is.null(selection))
    stop("Using the locator and highlighting a selection at the same time is not supported.")
  
  if(!is.null(selection)) {
    selection <- as.data.frame(selection)
    if(ncol(selection) != 2)
      stop("A selection must be provided as a 2-column data frame or matrix containing the x- (first column) and y- (second column) coordinates of the points in the selection.")
  }
  
  if(!is.null(x_scale) | !is.null(y_scale)) {
    if(any(!c(x_scale, y_scale) %in% c("sqrt", "log10")) | length(x_scale) > 1 | length(y_scale) > 1)
      stop("Axis scales must be either \"log10\" or \"sqrt\"")
  }
  
  #filter based on minimum length
  mm <- subset(mm, length >= min_length)
  
  ##### base plot #####
  if(is.null(color_by)) {
    p <- ggplot(mm, aes_string(x = x, y = y, size = "length"))
  } else if(!is.null(color_by)) {
    p <- ggplot(mm, aes_string(x = x, y = y, size = "length", color = color_by))
  }
  #geom_point when fixed_size is set
  if(is.null(fixed_size)) {
    #if color_by is set and is numeric
    if(ifelse(!is.null(color_by), is.numeric(mm[[color_by]]), FALSE)) {
      p <- p + 
        geom_point(alpha = alpha, na.rm = TRUE)
    } else {
      p <- p + 
        geom_point(alpha = alpha, color = "black", na.rm = TRUE) #all other variables than numerics are "black"
    }
    p <- p + 
      scale_size_area(name = "Scaffold length", max_size = 20*size_scale) #resize according to the set size_scale
    #if color_by is set and is factor or character
    if(ifelse(!is.null(color_by), is.factor(mm[[color_by]]) | is.character(mm[[color_by]]), FALSE)) {
      p <- p +
        geom_point(data = subset(mm, mm[[color_by]] != "NA"), shape = 1, alpha = 0.7, na.rm = TRUE) + 
        guides(colour = guide_legend(override.aes = list(alpha = 1, size = 5, shape = 19)))
    }
  } else if(!is.null(fixed_size)) {
    #geom_point when fixed_size is NOT set
    p <- p + 
      geom_point(alpha = alpha, color = 'black', size = fixed_size, na.rm = TRUE) 
    #if color_by is set and is factor or character
    if(ifelse(!is.null(color_by), is.factor(mm[[color_by]]) | is.character(mm[[color_by]]), FALSE)) {
      p <- p +
        geom_point(data = subset(mm, mm[[color_by]] != "NA"), shape = 1, alpha = 0.7, size = fixed_size, na.rm = TRUE) + 
        guides(colour = guide_legend(override.aes = list(alpha = 1, size = 5, shape = 19)))
    }
  }
  
  #theme adjustments
  p <- p +
    theme(panel.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "grey95"),
          axis.line = element_line(color = "black"),
          axis.ticks = element_line(color = "black"),
          axis.text = element_text(color = "black"),
          legend.key = element_blank() 
          )
  
  #add scaffold names
  if(isTRUE(scaffold_labels)) {
    p <- p + geom_text(label = mm[[1]], size = 4, color = "black")
  } else if(is.vector(scaffold_labels) | is.data.frame(scaffold_labels)) {
    if(is.data.frame(scaffold_labels)) {
      scaffold_labels <- as.character(scaffold_labels[[1]])
    }
    labels_data <- subset(mm, mm[[1]] %in% as.character(scaffold_labels))
    p <- p + geom_text(data = labels_data, 
                       aes_(x = labels_data[[x]], 
                            y = labels_data[[y]],
                            label = labels_data[[1]]), 
                       size = 4, 
                       color = "black",
                       inherit.aes = FALSE)
  }
  
  #axis scales and limits
  if(!is.null(x_scale)) {
    if(x_scale == "log10") {
      p <- p + scale_x_log10(limits = x_limits)
    } else if(x_scale == "sqrt") {
      p <- p + scale_x_sqrt(limits = x_limits)
    }
  }
  if(!is.null(y_scale)) {
    if(y_scale == "log10") {
      p <- p + scale_y_log10(limits = y_limits)
    } else if(y_scale == "sqrt") {
      p <- p + scale_y_sqrt(limits = y_limits)
    }
  }
  
  #axis limits when no custom scale
  if(!is.null(x_limits) & is.null(x_scale)) {
    suppressMessages(p <- p + xlim(x_limits))
  }
  if(!is.null(y_limits) & is.null(y_scale)) {
    suppressMessages(p <- p + ylim(y_limits))
  }
  
  ##### network #####
  #extract connections between scaffolds
  #if (!is.null(network)) {
  #  snetwork <- subset(network, network$scaffold1 %in% data$scaffolds$scaffold & network$scaffold2 %in% data$scaffolds$scaffold & connections >= nconnections)
  #  
  #  links <- merge(snetwork, data$scaffolds[,c("scaffold",x,y)], by.x = "scaffold1", by.y = "scaffold") 
  #  links <- merge(links, data$scaffolds[,c("scaffold",x,y)], by.x = "scaffold2", by.y = "scaffold") 
  #  colnames(links)[4:7] <- c("x","y","xend","yend") 
  #  p <- p +  
  #    geom_segment(data = links, aes(x = x, y = y, xend = xend, yend = yend), color = "darkgrey", size = 1, alpha = 0.5) +
  #    geom_point(data = links, aes(x = x, y = y), size = 2, color = "darkgrey") +
  #    geom_point(data = links, aes(x = xend, y = yend), size = 2, color = "darkgrey")
  #}
  
  ##### Highlight selected scaffolds #####
  
  #asdqweqwe
  #qwqeqwe
  #qweqweqwe
  
  ##### Plot duplicates #####
  if (isTRUE(shared_genes)) {
    eg <- tidyr::separate_rows(mm[!is.na(mm[,"geneID"]),c("scaffold", "geneID")], "geneID")
    df <- eg[which(duplicated(eg[[2]]) | duplicated(eg[[2]], fromLast=TRUE)),] 
    
    #see https://stackoverflow.com/questions/48407650/how-do-i-get-all-pairs-of-values-in-a-variable-based-on-shared-values-in-a-diffe
    split <- split(df$scaffold, df$geneID)
    shared_genes <- split[which(sapply(split, length) > 1)] %>% 
      lapply(function(split) {
        sort(split) %>%
          combn(m = 2) %>%
          t()
      }) %>%
      do.call(what = rbind) %>%
      unique() %>%
      tibble::as_tibble()
    colnames(shared_genes) <- c("scaffold1", "scaffold2")
    
    segment_coords <- merge(shared_genes, mm[,c("scaffold", x, y)], by.x = "scaffold1", by.y = "scaffold") 
    segment_coords <- merge(segment_coords, mm[,c("scaffold", x, y)], by.x = "scaffold2", by.y = "scaffold") 
    colnames(segment_coords)[3:6] <- c("x","y","xend","yend") 
    
    p <- p +
      geom_segment(data = segment_coords, aes(x = x, y = y, xend = xend, yend = yend), color = "darkred", size = 1) +
      geom_point(data = segment_coords, aes(x = x, y = y), size = 2, color = "darkred", na.rm = TRUE) +
      geom_point(data = segment_coords, aes(x = xend, y = yend), size = 2, color = "darkred", na.rm = TRUE)
  }
  
  ##### Locator and selection #####
  if(isTRUE(locator)) {
    points <- mmlocator(p, x_scale, y_scale)
  }
  if(isTRUE(locator) | !is.null(selection)) {
    if(!isTRUE(locator) & !is.null(selection)) {
      points <- selection
    }
    p$selection <- points
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
