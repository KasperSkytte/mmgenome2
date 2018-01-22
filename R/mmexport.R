#' mmexport
#'
#' @param mm 
#' @param assembly 
#' @param file 
#'
#' @return
#' @export
#'
#' @examples
mmexport <- function(mm, assembly, file = "exported_assembly.fa") {
  if(!any(class(mm) == "mm"))
    stop("The provided data is not a \"mm\" class data frame. Use mmload() to load various data into a \"mm\" data frame. (or class(df) <- append(class(df), \"mm\") if you know what you are doing)")
  if(!any(class(assembly) == "DNAStringSet"))
    stop("The provided assembly is not of class \"DNAStringSet\". Use Biostrings::readDNAStringSet(path, format = \"fasta\") to load the assembly.")
  Biostrings::writeXStringSet(assembly, filepath = file)
}
