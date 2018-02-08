#' Title
#'
#' @param mm 
#' @param selection 
#' @param min_length 
#' @param inverse 
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