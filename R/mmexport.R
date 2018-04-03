#' @title Export DNA sequences of selected scaffolds
#' 
#' @description Writes the DNA sequences of a given set of scaffolds from the assembly to a FASTA file. 
#'
#' @param scaffolds (\emph{required}) Either a vector of scaffold names or a dataframe loaded with \code{\link{mmload}} containing the scaffolds to extract from the \code{assembly}.
#' @param assembly (\emph{required}) The assembly to subset and export. Must be an object of class \code{"DNAStringSet"} as loaded with \code{\link{readDNAStringSet}}. By default any object named "assembly" in the global environment will be used.
#' @param file The file path and file name of the file to write the sequences to. (\emph{Default: } \code{"exported_assembly.fa"}) 
#'
#' @export
#' 
#' @importFrom Biostrings writeXStringSet
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
#' mmgenome2_extraction
#' \dontrun{
#' mmexport(mmgenome2_extraction,
#'          assembly = assembly,
#'          file = "bins/exported_assembly.fa")
#' }
#' 
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
#' @author Soren M. Karst \email{smk@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}
mmexport <- function(scaffolds, 
                     assembly = get("assembly", envir = .GlobalEnv), 
                     file = "exported_assembly.fa") {
  if(!any(class(assembly) == "DNAStringSet"))
    stop("The provided assembly is not of class \"DNAStringSet\". Use Biostrings::readDNAStringSet(path, format = \"fasta\") to load the assembly.", call. = FALSE)
  if(is.data.frame(scaffolds)) {
    scaffolds <- as.character(scaffolds[[1]])
  } else if(!any(is.vector(scaffolds), is.data.frame(scaffolds))) {
    stop("Scaffolds must be provided either as a vector, or as a dataframe, where the first column contains the scaffold names.", call. = FALSE)
  }
  Biostrings::writeXStringSet(assembly[scaffolds], filepath = file)
}
