#' @title Print basic statistics of metagenome data
#' 
#' @description Calculates some basic statistics like the N50 length, minimum-, mean-, and maximum scaffold lengths, mean GC content, information about essential genes, and more...
#'
#' @param mm (\emph{required}) A dataframe loaded with \code{\link{mmload}}.
#' @param original_data If \code{mm} is a subset/extraction of another dataframe loaded with \code{\link{mmload}}, then provide here the original dataframe from which the extraction originates to compare the extraction to the original data, see examples. (\emph{Default: } \code{NULL})
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
#' @importFrom magrittr %<>% %>% 
#' 
#' @examples 
#' library(mmgenome2)
#' data(mmgenome2)
#' mmstats(mmgenome2)
#' 
#' #Compare an extraction with the original data from which the extraction originates:
#' selection <- data.frame(cov_C13.11.25 = c(7.2, 16.2, 25.2, 23.3, 10.1),
#'                         cov_C14.01.09 = c(47, 77, 52.8, 29.5, 22.1))
#' mmgenome2_extraction <- mmextract(mmgenome2, 
#'                                   selection = selection)
#'                                   
#' mmstats(mmgenome2_extraction, original_data = mmgenome2)
#' 
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
#' @author Soren M. Karst \email{smk@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}
mmstats <- function(mm, 
                    original_data = NULL) {
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
  
  #calculate weighted coverage
  calcCov <- function(x) {
    unlist(lapply(dplyr::select(x, dplyr::starts_with("cov_")), 
                  function(cov) {
                    nonEmptyIDs <- which(complete.cases(cov))
                    nonEmpty_mm <- x[nonEmptyIDs,]
                    round(sum(cov[nonEmptyIDs]*nonEmpty_mm$length)/sum(nonEmpty_mm$length), 2)
                  }))
  }
  
  #make a named vector with the stats
  stats <- c(Scaffolds = length(mm$scaffold),
             N50 = N50,
             Length.total = lengthTotal,
             Length.max = max(mm$length),
             Length.mean = round(mean(mm$length) ,2),
             Length.min = min(mm$length),
             weighted_GC_mean = round(sum(mm$gc*mm$length)/lengthTotal, 2),
             if(is.null(original_data)) {
               calcCov(mm)
             } else if(!is.null(original_data)) {
               #If original_data is supplied calculate the coverage profiles
               #as a fraction of the original data coverage profiles
               round(lengthTotal/sum(original_data$length)*(calcCov(mm)/calcCov(original_data))*100, 2) %>% 
                 `names<-`(paste0(names(.), " (% of original)"))
             },
             Ess.genes.total = length(ess.genes),
             Ess.genes.unique = length(unique(ess.genes))
  )
  #save the stats in a data frame for export
  stats %<>% 
    as.data.frame() %>%
    `colnames<-`("General stats")
  
  #prettyprint the stats data frame
  stats[,1] %>% 
    as.character() %>%
    prettyNum(big.mark = " ", 
              big.interval = 3, 
              drop0trailing = T,
              nsmall = 2) %>%
    trimws() %>%
    StrAlign(sep = ".") %>%
    as.data.frame() %>%
    `colnames<-`("General stats") %>%
    `rownames<-`(rownames(stats)) %>%
    print.data.frame(
      quote = FALSE,
      right = TRUE,
      row.names = TRUE,
      max = nrow(.))
  
  #return the stats data frame invisibly
  invisible(stats)
}
