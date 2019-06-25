#' @title Custom ggplot object print to suppress warnings
#' @description Simply wraps \code{ggplot2:::print.ggplot()} in \code{suppressWarnings()}. The main reason for this is to hide warnings produced by log10-transformation in \code{mmplot}.
#'
#' @param ... All arguments are directly passed to \code{ggplot2:::print.ggplot()}
#' @import ggplot2
#' @export
print.ggplot <- function(...) {
  suppressWarnings(ggplot2:::print.ggplot(...))
}
