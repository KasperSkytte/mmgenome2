#' @title Compare coverage profiles
#'
#' @description Plots all coverage profiles in the data for comparison. This is useful for identifying scaffold patterns across multiple coverage profiles and spot potential contaminants.
#'
#' @param mm (\emph{required}) A dataframe loaded with \code{\link{mmload}}.
#' @param normalise (\emph{Logical}) Normalises each coverage column in \code{mm} to its highest value and multiplies by 100. (\emph{Default: } \code{TRUE})
#' @param alpha The transparancy of the scaffold points and lines, where 0 is invisible and 1 is opaque. (\emph{Default: } \code{0.6})
#' @param color_by Color the scaffolds by a variable in \code{mm}. (\emph{Default: } \code{NULL})
#' @param color_vector The colors from which to generate a color gradient when \code{color_by} is set and the variable is continuous. Any number of colors can be used. (\emph{Default: } \code{c("blue", "green", "red")})
#' @param color_scale_log10 (\emph{Logical}) Log10-scale the color gradient when \code{color_by} is set and the variable is continuous. (\emph{Default: } \code{FALSE})
#' @param y_scale_log10 (\emph{Logical}) Log10-scale the y axis. (\emph{Default: } \code{FALSE})
#' @param plot_lines (\emph{Logical}) Connect scaffolds with lines. (\emph{Default: } \code{TRUE})
#' @param interactive_plot (\emph{Logical}) Return an interactive \code{plotly} plot or not. (\emph{Default: } \code{FALSE})
#' @param mm_normalise A dataframe loaded with \code{\link{mmload}} to normalise by.
#' @param mm_normalise_index A coverage column to serve as an index for normalisation.
#'
#' @return A ggplot or plotly object. Note that mmgenome2 hides all warnings produced by ggplot objects.
#'
#' @export
#' @import ggplot2
#' @importFrom dplyr left_join
#' @importFrom data.table melt data.table
#' @importFrom magrittr %>% %<>%
#' @importFrom tidyr unite
#'
#' @examples
#' library(mmgenome2)
#' data(mmgenome2)
#' mmgenome2
#' selection <- data.frame(
#'   cov_C13.11.25 = c(7.2, 16.2, 25.2, 23.3, 10.1),
#'   cov_C14.01.09 = c(47, 77, 52.8, 29.5, 22.1)
#' )
#' mmgenome2_extraction <- mmextract(mmgenome2,
#'   min_length = 3000,
#'   selection = selection,
#'   inverse = FALSE
#' )
#' mmplot_cov_profiles(mmgenome2_extraction,
#'   color_by = "taxonomy",
#'   normalise = FALSE
#' )
mmplot_cov_profiles <- function(mm,
                                normalise = FALSE,
                                alpha = 0.6,
                                color_by = NULL,
                                color_vector = c("blue", "green", "red"),
                                color_scale_log10 = FALSE,
                                y_scale_log10 = FALSE,
                                plot_lines = TRUE,
                                interactive_plot = FALSE,
                                mm_normalise = NULL,
                                mm_normalise_index = NULL) {
  # can currently only color by 1 variable
  if (length(color_by) > 1) {
    stop("color_by must be of length 1", call. = FALSE)
  }

  # find which columns are coverage profiles
  cov_variables <- which(substring(tolower(colnames(mm)), 1, 4) == "cov_")

  # normalise to the highest value in each coverage profile
  if (isTRUE(normalise & is.null(mm_normalise))) {
    mm[, cov_variables] %<>% lapply(function(x) {
      x / max(x) * 100
    })
  }
  # normalise based on the amount of data for each sample
  if (isTRUE(normalise) & (!is.null(mm_normalise)) & (!is.null(mm_normalise_index))) {
    data_vol <- (mm_normalise$length * mm_normalise %>% select(starts_with("cov_"))) %>% summarise_all(funs(sum = sum))
    data_vol_normalised <- data_vol / as.numeric((data_vol %>% select(paste0(mm_normalise_index, "_sum"))))
    cov_variables <- which(substring(tolower(colnames(mm_normalise)), 1, 4) == "cov_")
    mm[, cov_variables] <- t(t(mm[, cov_variables]) / (as.numeric(data_vol_normalised)))
  }
  colnames(mm)[1] <- "scaffold"

  # make a data frame suited for ggplot2 and merge with the color_by variable if supplied
  gg <- data.table::melt(
    data.table::data.table(mm),
    id.vars = 1,
    measure.vars = cov_variables,
    value.name = "Coverage"
  ) %>%
    dplyr::left_join(mm[, c(1, 2, which(colnames(mm) %in% color_by)), drop = FALSE], by = "scaffold")

  # color_by
  if (!is.null(color_by)) {
    p <- ggplot(
      gg,
      aes_string(
        x = "variable",
        y = "Coverage",
        color = color_by
      )
    )

    # if numeric set custom colorscale and breaks if colored by GC
    if (is.numeric(mm[[color_by]])) {
      p <- p +
        scale_color_gradientn(
          colors = color_vector,
          trans = if (isTRUE(color_scale_log10)) "log10" else "identity",
          breaks = if (color_by == "gc") c(20, 40, 60, 80) else waiver()
        )
    }
  } else if (is.null(color_by)) {
    p <- ggplot(
      gg,
      aes_string(
        x = "variable",
        y = "Coverage"
      )
    )
  }

  # theme adjustments
  p <- p +
    theme(
      axis.text.x = element_text(angle = 90),
      axis.line = element_line(color = "black"),
      panel.background = element_blank(),
      panel.border = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey95"),
      legend.key = element_blank()
    ) +
    xlab("")

  # log10 scale y axis
  if (isTRUE(y_scale_log10)) {
    p <- p + scale_y_log10()
  }

  # plot lines
  if (isTRUE(plot_lines)) {
    p <- p + geom_line(aes_string(group = "scaffold"), alpha = alpha)
  }

  # adjust y axis label when normalise = TRUE
  if (isTRUE(normalise) & (is.null(mm_normalise)) & (is.null(mm_normalise_index))) {
    p <- p + ylab("Coverage (normalised to 100)")
  }
  # adjust
  if (isTRUE(normalise) & (!is.null(mm_normalise)) & (!is.null(mm_normalise_index))) {
    p <- p + ylab("Coverage\n(normalised by data volume)")
  }

  # add points to plot and generate plotly hover labels if interactive, return plot
  if (isTRUE(interactive_plot)) {
    checkReqPkgs(c("plotly", "purrr"))
    data_plotly <- gg %>%
      purrr::imap(~ paste(.y, .x, sep = ": ")) %>%
      as.data.frame() %>%
      tidyr::unite("test", sep = "<br>") %>%
      unlist(use.names = FALSE) %>%
      unique()
    p <- p +
      geom_point(
        alpha = alpha,
        aes(text = data_plotly)
      )
    plotly <- plotly::ggplotly(p, tooltip = "text") %>%
      plotly::layout(margin = list(
        l = 50,
        r = 0,
        b = 150,
        t = 20
      ))
    return(plotly)
  } else {
    p <- p +
      geom_point(alpha = alpha)
    return(p)
  }
}
