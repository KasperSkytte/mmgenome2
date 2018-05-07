#' @title Plot multiple combinations of variables in a pairs plot
#' 
#' @description Plots multiple variables from the given mm object in a grid plot with all pairs of variables using \code{\link{mmplot}}. 
#'
#' @param mm (\emph{required}) A dataframe loaded with \code{\link{mmload}}.
#' @param variables A vector of 3 or more variable names in \code{mm} to plot on each axis. If NULL, the default, then all coverage variables will be plotted as well as GC content. (\emph{Default: } \code{NULL})
#' @param textsize The text size of the axis titles. (\emph{Default: } \code{5})
#' @param axis_ticks Hide or show axis ticks on both axes. (\emph{Default: } \code{TRUE})
#' @param ... Arguments passed on to \code{\link{mmplot}}, eg. \code{color_by}, \code{min_length}, axis scales and more, see \code{help("mmplot")}.
#'
#' @export
#' 
#' @return A ggplot object. Note that mmgenome2 hides all warnings produced by ggplot objects.
#' 
#' @importFrom magrittr %>%
#' @importFrom cowplot plot_grid
#' @importFrom dplyr select starts_with
#' @import ggplot2
#' 
#' @examples
#' library(mmgenome2)
#' data(mmgenome2)
#' mmgenome2
#' mmplot_pairs(mmgenome2,
#'              variables = c("cov_C13.11.14", "cov_C13.11.25", "cov_C13.12.03", "cov_C14.01.09"),
#'              min_length = 10000,
#'              color_by = "taxonomy",
#'              x_scale = "log10",
#'              y_scale = "log10",
#'              x_limits = c(1, NA),
#'              y_limits = c(1, NA))
#' 
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
#' @author Soren M. Karst \email{smk@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}
mmplot_pairs <- function(mm, 
                         variables = NULL, 
                         textsize = 5,
                         axis_ticks = TRUE,
                         ...) {
  args <- list(...)
  if(any(names(args) == "locator")) {
    if(isTRUE(args[["locator"]]))
      stop("mmplot_pairs does not support the locator feature", call. = FALSE)
  }
  if(any(names(args) == "selection")) {
    if(!is.null(args[["selection"]]))
      stop("mmplot_pairs cannot highlight a \"selection\", use mmplot instead", call. = FALSE)
  }
  if(is.null(variables)) {
    variables <- names(dplyr::select(mm, dplyr::starts_with("cov_"), "gc"))
  }
  if(length(variables) == 2)
    stop("Plotting only two variables is better done with mmplot", call. = FALSE)
  if(!is.character(variables))
    stop("The variables to plot must be provided as a character vector with 3 or more variable names", call. = FALSE)
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
  if(isTRUE(axis_ticks)) {
    for(i in 0:(ncol-2)) {
      plotID <- 2+(ncol+1)*i
      temp[[plotID]] <- temp[[plotID]] +
        theme(axis.text = element_text(),
              axis.ticks = element_line())
    }
  }
  suppressWarnings(do.call(cowplot::plot_grid, c(temp, ncol = ncol, align = "hv")))
}