#' @title Merge link files
#' 
#' @description Searches a folder for all files with a filename ending with _link, and then all the files found are loaded and merged into one data frame. This is done by searching the files for all unique scaffold pairs and summing the number of connections per pair. The files must contain 3 columns called \code{"scaffold1"}, \code{"scaffold2"}, and \code{"connections"}.
#'
#' @param path (\emph{required}) Path to a folder.
#' @param verbose (\emph{Logical}) Whether to print status messages of the process. (\emph{Default: } \code{TRUE}) 
#'
#' @return A tibble
#' @export
#'
#' @examples
#' \dontrun{
#'   links <- mmmerge_linkfiles("path/to/folder")
#' }
mmmerge_linkfiles <- function(path, verbose = TRUE) {
  filepaths <- list.files(path = path, 
                          full.names = TRUE,
                          all.files = FALSE, 
                          recursive = FALSE,
                          ignore.case = TRUE)
  filepaths <- filepaths[grepl("*._link$", tools::file_path_sans_ext(filepaths))]
  if(length(filepaths) > 0) {
    filenames <- basename(filepaths)
    links <- list()
    for (i in 1:length(filenames)) {
      links[[stringr::str_remove(tools::file_path_sans_ext(filenames)[i], "_link$")]] <- data.table::fread(filepaths[i], data.table = FALSE)[1:3]
    }
  } else
    stop("No files with a filename ending with \"_link\" were found in the folder \"", path, "\"", call. = FALSE)
  if(!all(lapply(links, names) %>% unlist %>% unique == c("scaffold1", "scaffold2", "connections")))
    stop("All link files must have exactly 3 columns named \"scaffold1\", \"scaffold2\", and \"connections\", in that order.", call. = FALSE)
  if(isTRUE(verbose))
    message(paste0("Merging the following ", length(links), " link files found in the folder \"", path, "\":\n", paste0(filenames, collapse = "\n")))
  mergedlinks <- bind_rows(links) %>%
    group_by(scaffold1, scaffold2) %>%
    summarise_at("connections", sum) %>%
    ungroup() %>%
    arrange(desc(connections))
  if(isTRUE(verbose))
    message(paste0("\nTotal number of unique scaffold pairs: ", nrow(mergedlinks)))
  return(mergedlinks)
}
