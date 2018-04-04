#' @title Plot all coverage profiles 
#' 
#' @description To look for patterns... WIP
#'
#' @param mm 
#' @param normalise 
#' @param alpha 
#' @param color_by 
#' @param color_vector 
#' @param color_scale_log10 
#' @param y_scale_log10 
#' @param plot_lines 
#'
#' @return A ggplot2 object.
#' 
#' @export
#' @import ggplot2
#' @importFrom dplyr left_join
#' @importFrom data.table melt
#' @importFrom magrittr %>% %<>%
#'
#' @examples
#' library(mmgenome2)
#' data(mmgenome2)
#' mmgenome2
#' selection <- data.frame(cov_C13.12.03 = c(7.676, 5.165, 6.386, 10.933), 
#'                         cov_C14.01.09 = c(24.852, 32.545, 53.062, 38.52))
#' mmgenome2_extraction <- mmextract(mmgenome2, 
#'                                   min_length = 3000,
#'                                   selection = selection,
#'                                   inverse = FALSE)
#' mmplot_cov_profiles(mmgenome2_extraction,
#'                     color_by = "taxonomy")
mmplot_cov_profiles <- function(mm,
                                normalise = TRUE,
                                alpha = 0.6,
                                color_by = NULL,
                                color_vector = c("blue", "green", "red"),
                                color_scale_log10 = FALSE,
                                y_scale_log10 = FALSE,
                                plot_lines = TRUE) {
  if(length(color_by) > 1)
    stop("color_by must be of length 1", call. = FALSE)
  cov_variables <- which(substring(tolower(colnames(mm)), 1, 4) == "cov_")
  if(isTRUE(normalise))
    mm[,cov_variables] %<>% lapply(function(x) {x/max(x)*100})
  colnames(mm)[1] <- "scaffold"
  gg <- data.table::melt(mm,
                         id.vars = 1,
                         measure.vars = cov_variables,
                         value.name = "Coverage") %>%
    dplyr::left_join(mm[,c(1, which(colnames(mm) %in% color_by))], by = "scaffold")
  if(!is.null(color_by)) {
    p <- ggplot(gg, 
                aes_string(x = "variable",
                           y = "Coverage",
                           color = color_by))
    if(is.numeric(mm[[color_by]])) {
      if(!isTRUE(color_scale_log10)) {
        p <- p + 
          scale_color_gradientn(colors = color_vector, breaks = c(20,40,60,80))
      } else if(isTRUE(color_scale_log10)) {
        p <- p + 
          scale_color_gradientn(colors = color_vector, trans = "log10", breaks = c(20,40,60,80))
      }
    }
  } else if(is.null(color_by)) {
    p <- ggplot(gg, 
                aes_string(x = "variable",
                           y = "Coverage"))
  }
  if(isTRUE(y_scale_log10))
    p <- p + scale_y_log10()
  if(isTRUE(plot_lines))
    p <- p + geom_line(aes_string(group = "scaffold"), alpha = alpha)
  if(isTRUE(normalise))
    p <- p + ylab("Coverage (normalised to 100)")
  
  p <- p + 
    geom_point(alpha = alpha)
  
  #theme adjustments
  p <- p +
    theme(axis.text.x = element_text(angle = 90),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          legend.key = element_blank()
    ) +
    xlab("")
  return(p)
}
