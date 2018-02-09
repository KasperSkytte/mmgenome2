#' Generates a ggplot2 style pairs plot
#'
#' Generates a ggplot2 style pairs plot.
#'
#' @usage mmplot_pairs(mm, variables)
#'
#' @param mm (required) The dataframe to plot.
#' @param variables (required) The vaiables to plot as a vector.
#' @param textsize Size of the legend text (deafult = 8).
#' @param ... Arguments passed on to mmplot()
#' 
#' @return A pairs plot object.
#' 
#' @export
#' @author Soren M. Karst \email{smk@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}
#' 
#' @examples
#' data(rocco)
#' 
#' mmplot_pairs(mm = d, variables = c("C13.12.03", "C14.01.09", "C13.11.25", "gc", "PC2"), color = "essential", log = c("C13.12.03", "C14.01.09", "C13.11.25"), minlength = 10000)
#' 
#' 

mmplot_pairs <- function(mm, variables, textsize = 8, log = NULL, ...){  
  
  ## Make a blank plot
  emp <- data.frame(x = 0, y = 0)
  
  pblank <-  ggplot(emp, aes(x,y)) + 
    geom_blank() +
    theme_bw() +z
    theme(panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(), 
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks=element_blank(),
          panel.border=element_blank())
  
  ## Iterate through the different combinations of plots
  
  temp <- list()
  
  for (i in 1:length(variables)){
    for (j in 1:length(variables)){
      if (i < j){
        p <- mmplot(mm, 
                    x = variables[j], 
                    y = variables[i], 
                    ...) + 
          theme(plot.margin=unit(c(0,0,0,0),"cm"), 
                legend.position = "none",
                panel.grid.minor=element_blank(), 
                panel.grid.major=element_blank(), 
                axis.title.x=element_blank(), 
                axis.title.y=element_blank(),
                axis.text.x=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks=element_blank())#,
                #panel.border=element_blank())
        if(!is.null(log)){
          if(variables[j] %in% log & !(variables[i] %in% log)){
            p <- p + scale_x_log10()
          }
          if(variables[i] %in% log & !(variables[j] %in% log)){
            p <- p + scale_y_log10()
          }
          if(variables[i] %in% log & variables[j] %in% log){
            p <- p + scale_x_log10() + scale_y_log10()
          }
        }
      }    
      if (i == j){ 
        p <- pblank + geom_text(label=variables[i], size = textsize)
      }
      if(i > j){
        p <- pblank 
      }
      plotnr <- paste("x",variables[j],"y",variables[i], sep = "")
      temp[plotnr] <- list(p)
    }
  }
  
  do.call(gridExtra::grid.arrange, c(temp, ncol=floor(sqrt(length(temp)))))
}