#' @title Print summary of basic statistics
#'
#' @param mm (\emph{required}) A dataframe loaded with \code{\link{mmload}}.
#' @export
#'
#' @import dplyr
#' @import tidyr
mmstats <- function(mm) {
  mm$length <- as.numeric(mm$length)
  lengthTotal <- sum(mm$length)
  
  #N50
  lengths <- sort(as.integer(mm$length), decreasing = TRUE)
  N50 <- lengths[which(cumsum(lengths) >= lengthTotal/2)[1]]
  
  #essential genes
  ess.genes <- tidyr::separate_rows(mm[!is.na(mm[,"geneID"]),c("scaffold", "geneID")], "geneID")[["geneID"]]
  
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