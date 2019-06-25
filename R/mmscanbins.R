#' @title Retrieve scaffold names from all genome bins in a folder
#' @description Reads all FASTA files (genome bins) in a folder and returns all unique scaffold names. This is useful to get an overview of which scaffolds are already extracted from the metagenome by for example highlighting them in plots.
#'
#' @param binfolder (\emph{required}) Path to the folder to scan.
#' @param namesOnly Return only the scaffold names (\code{TRUE}) or a complete list with the DNA sequences of the scaffolds in each file (\code{FALSE}). (\emph{Default: } \code{TRUE})
#'
#' @return A character vector (if \code{namesOnly = TRUE}) or a list (if \code{namesOnly = FALSE}).
#' @export
#' @importFrom tools file_path_sans_ext file_ext
#' @importFrom Biostrings readDNAStringSet
#' @examples
#' \dontrun{
#' # Use mmextract and mmexport to extract scaffolds from the assembly and
#' # write their sequences to a FASTA file per bin. Then the mmscanbins() function
#' # can be used to retrieve the names of all scaffolds already extracted, and
#' # then for example highlight them in a plot:
#'
#' binned_scaffolds <- mmscanbins("path/to/folder", namesOnly = TRUE)
#' binned_scaffolds
#' mmplot(mmgenome2,
#'   min_length = 10000,
#'   x = "cov_C13.11.25",
#'   y = "cov_C14.01.09",
#'   color_by = "taxonomy",
#'   highlight_scaffolds = binned_scaffolds,
#'   highlight_color = "darkred",
#'   x_scale = "log10",
#'   y_scale = "log10"
#' )
#' }
mmscanbins <- function(binfolder, namesOnly = TRUE) {
  # check if binfolder is actually a path to a folder
  if (!is.character(binfolder)) {
    stop("That does not look like a path to a folder", call. = FALSE)
  }
  if (length(binfolder) > 1) {
    stop("binfolder must only contain 1 folder path", call. = FALSE)
  }
  if (!dir.exists(binfolder)) {
    stop("The folder \"", binfolder, "\" does not exist", call. = FALSE)
  }
  # get paths to all files in the folder
  filepaths <- list.files(
    path = binfolder,
    full.names = TRUE,
    all.files = FALSE,
    recursive = FALSE,
    ignore.case = TRUE
  )
  # subset to only FASTA files
  filepaths <- filepaths[which(tolower(tools::file_ext(filepaths)) %in% c("fasta", "fa"))]
  # if no FASTA files return error
  if (length(filepaths) > 0) {
    # read all FASTA files into a list
    bins <- list()
    filenames <- basename(filepaths)
    for (i in 1:length(filepaths)) {
      bins[[tools::file_path_sans_ext(filenames[i])]] <- Biostrings::readDNAStringSet(filepaths[[i]], format = "fasta")
    }
    # extract unique scaffold names and return them as a character vector, or return a list with the DNA
    if (isTRUE(namesOnly)) {
      scaffolds <- lapply(bins, names) %>%
        unlist() %>%
        unique() %>%
        as.character()
      return(scaffolds)
    } else if (!isTRUE(namesOnly)) {
      return(bins)
    }
  } else {
    stop("No FASTA files found in the folder \"", binfolder, "\"", call. = FALSE)
  }
}
