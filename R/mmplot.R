#' mmplot
#'
#' @param mm 
#' @param x 
#' @param y 
#' @param color_by 
#' @param alpha 
#' @param fixed_size 
#' @param size_scale 
#' @param min_length 
#'
#' @return
#' @export
#'
#' @examples
mmplot <- function(mm,
                   x,
                   y,
                   color_by = NULL,
                   alpha = 0.1,
                   fixed_size = NULL,
                   size_scale = 1, 
                   min_length = 0) {
  #Check if data has been loaded with mmload
  if(!any(class(mm) == "mm"))
    stop("The provided data is not of class \"mm\". Please use the mm_load() function to load data before using any mmgenome function.")
  
  #filter based on minimum length
  mm <- subset(mm, length >= min_length)
  
  ##### base plot #####
  if(is.null(color_by)) {
    p <- ggplot(mm, aes_string(x = x, y = y, size = "length"))
  } else if(!is.null(color_by)) {
    p <- ggplot(mm, aes_string(x = x, y = y, size = "length", color = color_by))
  }
  if(is.null(fixed_size)) {
    if(ifelse(!is.null(color_by), is.numeric(mm[[color_by]]), FALSE)) {
      p <- p + 
        geom_point(alpha = alpha)
    } else {
      p <- p + 
        geom_point(alpha = alpha, color = "black")
    }
    p <- p + 
      scale_size_area(name = "Scaffold length", max_size = 20*size_scale)
    if(ifelse(!is.null(color_by), is.factor(mm[[color_by]]) | is.character(mm[[color_by]]), FALSE)) {
      p <- p +
        geom_point(data = subset(mm, mm[[color_by]] != "NA"), shape = 1, alpha = 0.7) + 
        guides(colour = guide_legend(override.aes = list(alpha = 1, size = 5, shape = 19)))
    }
  } else if(!is.null(fixed_size)) {
    p <- p + 
      geom_point(alpha = alpha, color = 'black', size = fixed_size) 
    if(ifelse(!is.null(color_by), is.factor(mm[[color_by]]) | is.character(mm[[color_by]]), FALSE)) {
      p <- p +
        geom_point(data = subset(mm, mm[[color_by]] != "NA"), shape = 1, alpha = 0.7, size = fixed_size) + 
        guides(colour = guide_legend(override.aes = list(alpha = 1, size = 5, shape = 19)))
    }
  }
  
  p <- p +
    theme(panel.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "grey95"),
          axis.line = element_line(color = "black"),
          axis.ticks = element_line(color = "black"),
          axis.text = element_text(color = "black"),
          legend.key = element_blank()
          ) 
  return(p)
}