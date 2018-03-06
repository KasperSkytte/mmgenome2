#' @title Print basic statistics of metagenome data
#' 
#' @description Calculates some basic statistics like the N50 length, minimum-, mean-, and maximum scaffold lengths, mean GC content, information about essential genes, and more...
#'
#' @param mm (\emph{required}) A dataframe loaded with \code{\link{mmload}}.
#' @export
#'
#' @details 
#' The returned stats are calculated as follows:
#' \itemize{
#'    \item Scaffolds: The number of different scaffolds in the assembly.
#'    \item N50: The shortest sequence length at 50% of the assembly.
#'    \item Length.total: The total size of the assembly.
#'    \item Length.max: The size of the largest scaffold in the assembly.
#'    \item Length.mean: The average scaffold size in the assembly.
#'    \item Length.min: The size of the smallest scaffold in the assembly.
#'    \item weighted_GC_mean: The average GC content in the assembly, weighted by scaffold sizes.
#'    \item cov_*: The average coverage of each coverage variable in mm, weighted by scaffold sizes. (Only columns starting with "cov_" will be shown)
#'    \item Ess.genes.total: The total number of essential genes, if any have been loaded. 
#'    \item Ess.genes.unique: The number of unique essential genes, if any have been loaded. 
#' }
#' 
#' @return A dataframe with the calculated stats.
#'
#' @importFrom dplyr select starts_with
#' @importFrom tidyr separate_rows
#' 
#' @examples 
#' library(mmgenome2)
#' data(mmgenome2)
#' mmgenome2
#' mmstats(mmgenome2)
#' 
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
#' @author Soren M. Karst \email{smk@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}
mmstats <- function(mm) {
  mm$length <- as.numeric(mm$length)
  lengthTotal <- sum(mm$length)
  
  #N50
  lengths <- sort(as.integer(mm$length), decreasing = TRUE)
  N50 <- lengths[which(cumsum(lengths) >= lengthTotal/2)[1]]
  
  #essential genes
  if(any(colnames(mm) == "geneID")) {
    ess.genes <- tidyr::separate_rows(mm[!is.na(mm[,"geneID"]),c("scaffold", "geneID")], "geneID")[["geneID"]]
  } else
    ess.genes <- NULL
  
  stats <- c(Scaffolds = length(mm$scaffold),
             N50 = N50,
             Length.total = lengthTotal,
             Length.max = max(mm$length),
             Length.mean = round(mean(mm$length) ,2),
             Length.min = min(mm$length),
             weighted_GC_mean = round(sum(mm$gc*mm$length)/lengthTotal, 2),
             unlist(lapply(dplyr::select(mm, dplyr::starts_with("cov_")), function(x) {round(sum(x*mm$length)/lengthTotal, 2)})),
             Ess.genes.total = length(ess.genes),
             Ess.genes.unique = length(unique(ess.genes))
  )
  df <- data.frame(stats)
  colnames(df)[1] <- "General stats"
  df[,1] <- as.character(df[,1])
  return(df)
}