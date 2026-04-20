#' @title Filter peaks in retention range
#'
#' @description
#' This function reduces the original data frame to retention times that are biological relevant (without the noise at the beginning and end).
#' This region needs to be defined by the user based on length of their elution time and chromatograms.
#'
#' @usage interestRegion(x, start, end)
#' @param x data frame
#' @param start numeric value from which RT should the peaks be picked
#' @param end numeric value to which RT should the peaks be picked
#'
#' @details
#' The format of the data frame should be: Peak names|abundances of samples|m.z|RT|...
#'
#'
#'
#' @returns A reduced data frame.
#'
#' @examples
#' neg_small <- interestRegion(neg, 1.4, 25)
#'
#'
#' @export

interestRegion <- function(x, start, end){
  x.new <- x[x$RT >= start & x$RT <= end,]
  return(x.new)
}
