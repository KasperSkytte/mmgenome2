#' mmstats
#'
#' @param mm 
#' @param compare 
#'
#' @return
#' @export
#'
#' @examples
mmstats <- function(mm, compare) {
  attributes <- get(class(mm)[grepl(".mmID_", class(mm))], envir = .GlobalEnv)
  mm$length <- as.numeric(mm$length)
  lengthTotal <- sum(mm$length)
  
  #N50
  lengths <- sort(as.integer(mm$length), decreasing = TRUE)
  N50 <- lengths[which(cumsum(lengths) >= lengthTotal/2)[1]]
  
  #fix coverage names
  coverage <- unlist(lapply(mm[,attributes[["coverageCols"]]], function(x) {round(sum(x*mm$length)/lengthTotal, 2)}))
  names(coverage) <- paste0("cov_", names(coverage))
  
  stats <- c(Scaffolds = length(mm$scaffold),
             N50 = N50,
             Length.total = lengthTotal,
             Length.max = max(mm$length),
             Length.mean = round(mean(mm$length) ,2),
             Length.min = min(mm$length),
             weighted_GC_mean = round(sum(mm$gc_pct*mm$length)/lengthTotal, 2),
             coverage,
             Ess.genes.total = attributes[["total_Ess.genes"]],
             Ess.genes.unique = attributes[["unique_Ess.genes"]]
  )
  df <- data.frame(stats)
  colnames(df)[1] <- "General stats"
  df[,1] <- as.character(df[,1])
  return(df)
}