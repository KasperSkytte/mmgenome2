% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mmplot_cov_profiles.R
\name{mmplot_cov_profiles}
\alias{mmplot_cov_profiles}
\title{Compare coverage profiles}
\usage{
mmplot_cov_profiles(mm, normalise = TRUE, alpha = 0.6, color_by = NULL,
  color_vector = c("blue", "green", "red"), color_scale_log10 = FALSE,
  y_scale_log10 = FALSE, plot_lines = TRUE)
}
\arguments{
\item{mm}{(\emph{required}) A dataframe loaded with \code{\link{mmload}}.}

\item{normalise}{(\emph{Logical}) Normalises each coverage column in \code{mm} to its highest value and multiplies by 100. (\emph{Default: } \code{TRUE})}

\item{alpha}{The transparancy of the scaffold points and lines, where 0 is invisible and 1 is opaque. (\emph{Default: } \code{0.6})}

\item{color_by}{Color the scaffolds by a variable in \code{mm}. (\emph{Default: } \code{NULL})}

\item{color_vector}{The colors from which to generate a color gradient when \code{color_by} is set and the variable is continuous. Any number of colors can be used. (\emph{Default: } \code{c("blue", "green", "red")})}

\item{color_scale_log10}{(\emph{Logical}) Log10-scale the color gradient when \code{color_by} is set and the variable is continuous. (\emph{Default: } \code{FALSE})}

\item{y_scale_log10}{(\emph{Logical}) Log10-scale the y axis. (\emph{Default: } \code{FALSE})}

\item{plot_lines}{(\emph{Logical}) Connect scaffolds with lines. (\emph{Default: } \code{TRUE})}
}
\value{
A ggplot2 object.
}
\description{
Plots all coverage profiles in the data for comparison. This is useful for identifying scaffold patterns across multiple coverage profiles and spot potential contaminants.
}
\examples{
library(mmgenome2)
data(mmgenome2)
mmgenome2
selection <- data.frame(cov_C13.12.03 = c(7.676, 5.165, 6.386, 10.933), 
                        cov_C14.01.09 = c(24.852, 32.545, 53.062, 38.52))
mmgenome2_extraction <- mmextract(mmgenome2, 
                                  min_length = 3000,
                                  selection = selection,
                                  inverse = FALSE)
mmplot_cov_profiles(mmgenome2_extraction,
                    color_by = "taxonomy")
}