#' @title Merge mm object with more data
#' 
#' @description Internal function used in mmload. \code{y} can be a named vector, dataframe, or a list containing any of these types of data to be merged with \code{x}.
#'
#' @param x mm object
#' @param y Object to merge with \code{x}
#' @param type Character string defining the type of data being merged
#' 
#' @importFrom dplyr left_join intersect filter
#' @importFrom tibble enframe
#'
#' @return A tibble.
#' 
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
mmmerge <- function(x, y, type) {
  #must be a data frame, named atomic vector, or a list of data frames and/or named vectors
  if(any(class(y) %in% c("list", "data.frame", "tbl", "tbl_df")) | is.atomic(y) | is.factor(y)) {
    #wrap non-lists in a list to work with the for-loop
    if(any(!class(y) %in% "list"))
      y <- list(y)
    #merge mm with each element in the provided list
    for(i in 1:length(y)) {
      string <- ifelse(length(y) == 1, paste0("'", type, "'"), paste0("'", type,"'", " element ", i))
      if(is.factor(y[[i]]))
        y[[i]] <- as.character(y[[i]])
      #first column must be the sequence names
      if(any(class(y[[i]]) %in% c("data.frame", "tbl", "tbl_df")) & length(y[[i]]) < 2)
        stop(paste0(string, " not accepted: Data frames must contain at least 2 columns where the first column contains the sequence names exactly matching those of the assembly."))
      
      #vectors must be named to be able to merge with mm
      if(is.atomic(y[[i]]) & is.null(names(y[[i]])))
        stop(paste0(string, " not accepted: The vector is not a named vector. The vector elements must be named by sequence names exactly matching those of the assembly."))
      
      #column names are preserved from data frames, but not from vectors. Use names of the provided list, or else create a dummy name
      if(is.atomic(y[[i]]) & !is.null(names(y[[i]]))) {
        y[[i]] <- tibble::enframe(y[[i]], name = "scaffold", value = ifelse((is.null(names(y)) | names(y)[[i]] == ""), paste0(type, i), names(y)[[i]]))
      }
      
      #merge x and y[[i]] by scaffold
      colnames(y[[i]])[1] <- "scaffold" #first columns must have same name
      y[[i]][1] <- lapply(y[[i]][1], as.character) #and must be character
      sharedScaffolds <- dplyr::intersect(x$scaffold, y[[i]][["scaffold"]]) #which scaffolds are shared between x and y[[i]]
      
      #print missing or excess scaffolds between x and y[[i]]
      if(!all(x$scaffold %in% y[[i]][["scaffold"]])) {
        missingScaffolds <- dplyr::filter(x, !scaffold %in% sharedScaffolds)[[1]]
        warning(paste0("Only ", length(sharedScaffolds), " of all ", length(x$scaffold), " scaffolds in the assembly match in ", string,  ". The following ", length(missingScaffolds), " scaffolds are missing:\n\"", paste(missingScaffolds, collapse = "\", \""), "\""))
      } else if(!all(y[[i]][["scaffold"]] %in% x$scaffold)) {
        excessScaffolds <- filter(y, !scaffold %in% sharedScaffolds)[[1]]
        warning(paste0(string, " contains more scaffolds than the assembly. The following ", length(excessScaffolds), " scaffolds have not been loaded:\n\"", paste(excessScaffolds, collapse = "\", \""), "\""))
      } else if (!any(x$scaffold %in% y[[i]][["scaffold"]]))
        #no match sucks
        stop("No scaffold names match between the assembly and ", string, ". ")
      x <- dplyr::left_join(x, 
                            y[[i]], 
                            by = "scaffold")
    }
  } else
    stop("Data must be provided as a data frame, named vector, or a list of multiple data frames and/or named vectors.")
  return(x)
}