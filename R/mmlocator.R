#' @title Locate points in a ggplot2 plot
#' 
#' @description Internal function used within mmplot and mmplot_network
#'
#' @param plot ggplot2 object.
#' @param x_scale x axis scale, \code{"sqrt"} or \code{"log10"}.
#' @param y_scale y axis scale, \code{"sqrt"} or \code{"log10"}.
#'
#' @return A data frame with the x/y coordinates of the mousepositions clicked in the ggplot2 plot.
#' 
#' @import ggplot2
#' @import grid
#' @import clipr
#' 
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
#' @author Rasmus Hansen Kirkegaard \email{rhk@@bio.aau.dk}
#' @author Soren M. Karst \email{smk@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}
mmlocator <- function(plot, x_scale = NULL, y_scale = NULL) {
  #build ggplot object to be able to extract axis ranges, and show the plot
  suppressWarnings(ggobj <- ggplot2::ggplot_build(plot))
  suppressWarnings(print(ggobj$plot))
  #find the active panel
  panel <- unlist(grid::current.vpTree())
  panel <- unname(panel[grep("\\.name$", names(panel))])
  panel <- grep("panel", panel, fixed = TRUE, value = TRUE)
  if(length(panel) != 1){
    stop("no ggplot detected in current device")
  }
  grid::seekViewport(panel, recording = FALSE)
  #store the relative positions clicked (relative to the plot window 0-1) in a data frame and plot the positions as darkred dots. Repeat until user clicks "finish"
  points <- data.frame(x = as.numeric(), y = as.numeric())
  repeat {
    point <- grid::grid.locator(unit = "native")
    if(is.null(point))
      break
    grid::grid.points(point[1], point[2], pch = 16, gp = grid::gpar(cex = 0.5, col = "darkred"))
    points <- rbind(points, as.numeric(point))
  }
  
  ### Extract axis ranges. WARNING: keep an eye on changes in ggplot2 behind the scenes, accessing ranges has changed a few times already
  xrange <- ggobj$layout$panel_params[[1]][["x.range"]]
  yrange <- ggobj$layout$panel_params[[1]][["y.range"]]
  #rescale the relative positions by the axis ranges
  points[,1] <- (xrange[2]-xrange[1]) * points[,1] + xrange[1]
  points[,2] <- (yrange[2]-yrange[1]) * points[,2] + yrange[1]
  #adjust points if log10 or square-root scales have been applied in the plot
  if(!is.null(x_scale)) {
    if(x_scale == "log10") {
      points[,1] <- 10^points[,1]
    } else if(x_scale == "sqrt") {
      points[,1] <- points[,1]^2
    }
  }
  if(!is.null(y_scale)) {
    if(y_scale == "log10") {
      points[,2] <- 10^points[,2]
    } else if(y_scale == "sqrt") {
      points[,2] <- points[,2]^2
    }
  }
  #print the points as data frame code ready to copy/paste and store as a selection
  colnames(points) <- c(ggobj$plot$mapping$x, ggobj$plot$mapping$y)
  selection <- paste0("data.frame(", 
                      colnames(points[1]), 
                      " = ", 
                      paste0(round(points[1], 3)),
                      ", ",
                      colnames(points[2]),
                      " = ",
                      paste0(round(points[2], 3)),
                      ")"
  )
  message(paste0("Selection:\n", selection))
  userChoice <- readline(prompt = "Do you want to copy the selection to clipboard? (y/n or ENTER/ESC): ")
  if(tolower(userChoice) %in% c("y", "", "yes")) {
    clipr::write_clip(selection)
  }
  #and return the selection
  return(points)
}
