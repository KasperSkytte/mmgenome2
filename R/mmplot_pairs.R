#' Plot multiple combinations of variables in a pairs plot
#'
#' @param mm (\emph{required}) A dataframe loaded with \code{\link{mmload}}.
#' @param variables A vector of 3 or more variable names in \code{mm} to plot on each axis. If NULL, the default, then all coverage variables will be plotted. (\emph{Default: } \code{NULL})
#' @param textsize The text size of the axis titles.
#' @param ... Arguments passed on to \code{\link{mmplot}}, eg. \code{color_by}, \code{min_length}, axis scales and more, see help("mmplot").
#'
#' @export
#' 
#' @import cowplot
#' @import dplyr
#' @import ggplot2
#'
mmplot_pairs <- function(mm, 
                         variables = NULL, 
                         textsize = 5,
                         ...) {
  if(is.null(variables)) {
    variables <- names(dplyr::select(mm, dplyr::starts_with("cov_")))
  }
  if(length(variables) == 2)
    stop("Plotting only two variables is better done with mmplot")
  if(!is.character(variables))
    stop("The variables to plot must be provided as a character vector with 3 or more variable names")
  ## Make a blank plot
  emp <- data.frame(x = 0, y = 0)
  
  pblank <-  ggplot(emp, aes(x,y)) + 
    geom_blank() +
    theme_bw() +
    theme(plot.margin = margin(0,0,0,0),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_blank())
  
  ## Iterate through the different combinations of plots
  temp <- list()
  for (i in 1:length(variables)){
    for (j in 1:length(variables)){
      if (i < j){
        p <- mmplot(mm, 
                    x = variables[j], 
                    y = variables[i], 
                    ...) + 
          theme(plot.margin = margin(3,3,0,0, unit = "pt"), 
                legend.position = "none",
                axis.title.x = element_blank(), 
                axis.title.y = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.line = element_blank(),
                panel.border = element_rect(fill = NA, size = 0.5, linetype = 3, color = "gray75")
                )
      }    
      if (i == j){ 
        p <- pblank + geom_text(label = variables[i], size = textsize)
      }
      if(i > j){
        p <- pblank 
      }
      plotnr <- paste0("x_",variables[j],"y_",variables[i])
      temp[plotnr] <- list(p)
    }
  }
  ncol <- temp %>%
    length() %>%
    sqrt() %>%
    floor() %>%
    as.integer()
  for(i in 0:(ncol-2)) {
    plotID <- 2+(ncol+1)*i
    temp[[plotID]] <- temp[[plotID]] +
      theme(axis.text = element_text(),
            axis.ticks = element_line())
  }
  do.call(cowplot::plot_grid, c(temp, ncol = ncol, align = "hv"))
}