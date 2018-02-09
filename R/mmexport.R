#' @title Export DNA sequences of selected scaffolds
#' 
#' @description Writes the DNA sequences of a set of scaffolds from the assembly to a FASTA file. 
#'
#' @param scaffolds (\emph{required}) Either a vector of scaffold names or a dataframe loaded with \code{\link{mmload}} containing the scaffolds to extract from the \code{assembly}.
#' @param assembly (\emph{required}) The assembly to subset and export. Must be an object of class \code{"DNAStringSet"} as loaded with \code{\link[Biostrings]{readDNAStringSet}}. 
#' @param file The file path and file name of the file to write the sequences to. (\emph{Default: } \code{"exported_assembly.fa"}) 
#'
#' @export
#' 
#' @import Biostrings
mmexport <- function(mm, assembly = assembly, file = "exported_assembly.fa") {
  if(!any(class(assembly) == "DNAStringSet"))
    stop("The provided assembly is not of class \"DNAStringSet\". Use Biostrings::readDNAStringSet(path, format = \"fasta\") to load the assembly.")
  if(is.data.frame(scaffolds)) {
    scaffolds <- scaffolds[[1]]
  } else if(!is.vector(scaffolds) | !is.data.frame(scaffolds)) {
    stop("The scaffolds to export from the assembly must be either provided as a vector or as a dataframe, where the first column contains scaffold names.")
  }
  Biostrings::writeXStringSet(assembly[scaffolds], filepath = file)
}
