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
#' @importFrom shiny actionButton div fillPage icon observeEvent p plotOutput reactiveValues renderPlot runApp shinyApp stopApp
#' @importFrom clipr write_clip
#' 
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
#' @author Rasmus Hansen Kirkegaard \email{rhk@@bio.aau.dk}
#' @author Soren M. Karst \email{smk@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}
mmlocator <- function(plot, x_scale = NULL, y_scale = NULL) {
  axisnames <- c(plot[["mapping"]][["x"]], plot[["mapping"]][["y"]])
  app <- shinyApp(
    ui = fillPage(padding = c(5,5,50),
                  div(style="display: inline-block;",
                      actionButton(
                        inputId = "undo",
                        label = "Undo",
                        icon = icon("undo")
                      ),
                      actionButton(
                        inputId = "stop",
                        label = "Finish",
                        icon = icon("check")
                      )
                  ),
                  p(),
                  plotOutput("plot", 
                             click = "plot_click",
                             height = "100%",
                             width = "100%")
    ),
    server = function(input, output, session) {
      clickData <- reactiveValues(
        x = numeric(),
        y = numeric()
      )
      
      observeEvent(input$plot_click, {
        new_x <- input$plot_click$x
        new_y <- input$plot_click$y
        
        #adjust points if square-root scale have been applied in the plot
        if(!is.null(x_scale)) {
          if(x_scale == "sqrt")
            new_x <- new_x^2
        }
        if(!is.null(y_scale)) {
          if(y_scale == "sqrt")
            new_y <- new_y^2
        }
        
        #add the click coordinates to clickData
        clickData[["x"]] <- append(clickData[["x"]], new_x)
        clickData[["y"]] <- append(clickData[["y"]], new_y)
        
        #Update selection object in global environment
        df <- get(".current_selection", envir = globalenv())
        assign(".current_selection", 
               rbind(df,
                     data.frame(x = new_x,
                                y = new_y)),
               envir = globalenv())
        #message("point clicked:\n",
        #        "\t", colnames(df)[1], " = ", round(df[nrow(df), 1], 3),
        #        "\t", colnames(df)[2], " = ", round(df[nrow(df), 2], 3))
      })
      
      observeEvent(input$undo, {
        #remove the last point clicked
        lastPoint <- length(clickData[["x"]])
        clickData[["x"]] <- clickData[["x"]][-lastPoint]
        clickData[["y"]] <- clickData[["y"]][-lastPoint]
        #Update selection object in global environment
        df <- get(".current_selection", envir = globalenv())
        assign(".current_selection", df[-lastPoint,], envir = globalenv())
      })
      
      observeEvent(input$stop, {
        stopApp()
      })
      
      output$plot <- renderPlot({
        p <- plot + 
          geom_point(data = data.frame(x = clickData[["x"]],
                                       y = clickData[["y"]]),
                     aes_string(x = "x",
                                y = "y"),
                     color = "black",
                     inherit.aes = FALSE,
                     na.rm = TRUE) +
          geom_polygon(data = data.frame(x = clickData[["x"]],
                                         y = clickData[["y"]]),
                       aes_string(x = "x",
                                  y = "y"),
                       fill = NA,
                       size = 0.5,
                       lty = 2,
                       color = "darkred",
                       inherit.aes = FALSE, 
                       na.rm = TRUE)
        return(p)
      })
    }, 
    onStart = function() {
      df <- data.frame(x = numeric(), y = numeric())
      colnames(df) <- axisnames
      assign(".current_selection", df, envir = globalenv())
    }
  )
  suppressWarnings(
    runApp(app, 
           quiet = TRUE,
           launch.browser = rstudioapi::viewer))
  selection <- paste0("data.frame(", 
                      colnames(.current_selection[1]), 
                      " = ", 
                      paste0(round(.current_selection[1], 3)),
                      ", ",
                      colnames(.current_selection[2]),
                      " = ",
                      paste0(round(.current_selection[2], 3)),
                      ")"
  )
  message(paste0("Selection:\n", selection))
  userChoice <- readline(prompt = "Do you want to copy the selection to clipboard? (y/n or ENTER/ESC): ")
  if(tolower(userChoice) %in% c("y", "", "yes")) {
    clipr::write_clip(selection)
  }
  return(.current_selection)
}
