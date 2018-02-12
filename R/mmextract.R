#' Extract all scaffolds within a selection polygon
#'
#' @param mm (\emph{required}) A dataframe loaded with \code{\link{mmload}}.
#' @param selection (\emph{required}) A 2-column dataframe with the x and y coordinates of a selection of points in an \code{mmplot}. The column names of the provided dataframe must match column names in \code{mm}. 
#' @param min_length Filter all scaffolds with a length at or below this value before the extraction. (\emph{Default: } \code{0})
#' @param inverse (\emph{Logical}) If \code{TRUE}, then the scaffolds within the \code{selection} are instead removed. (\emph{Default: } \code{FALSE})
#'
#' @export
#' 
#' @import sp
#' @import dplyr
mmextract <-  function(mm, 
                       selection,
                       min_length = 0,
                       inverse = FALSE) {
  if(!any(colnames(selection) %in% colnames(mm)))
    stop("Could not find any variable names in mm matching those in the selection.")
  
  #filter based on minimum length
  mm <- dplyr::filter(mm, length >= min_length)
  
  #return scaffolds only in the selection polygon
  in_polygon <- sp::point.in.polygon(point.x = mm[[colnames(selection)[1]]],
                                     point.y = mm[[colnames(selection)[2]]],
                                     pol.x = selection[[1]],
                                     pol.y = selection[[2]],
                                     mode.checked = TRUE)
  ifelse(isTRUE(inverse), 
         mm <- dplyr::filter(mm, in_polygon == 0), 
         mm <- dplyr::filter(mm, in_polygon > 0))
  return(mm)
}