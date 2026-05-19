#' @title Data format control
#' @description
#' Given the data table and number of samples, the functions checks if the table
#' is in good format for the rest of MetaboPeak functions
#'
#' @usage checkFormat(x, n)
#' @param x data frame
#' @param n number of samples
#'
#' @details
#' This function validates the structure of LC-MS peak tables.
#'
#' It checks whether:
#' \itemize{
#'   \item the input is a data.frame
#'   \item required columns ("m.z", "RT") are present
#'   \item the number of sample columns matches \code{n}
#'   \item sample abundance columns are numeric
#' }
#'
#' @returns Logical scalar with the value TRUE if everything is in order, or an
#' error specifying where is problem.
#'
#' @export
#'
#' @examples
#' \donttest{
#' data("neg", package = "MetaboPeak")
#' checkFormat(neg, 48)
#' }


checkFormat <- function(x, n){

  # must be data frame
  if (!is.data.frame(x)) {
    stop("Input must be a data.frame")
  }

  # required columns
  if (!"m.z" %in% colnames(x)) {
    stop("Column 'm.z' is missing")
  }

  if (!"RT" %in% colnames(x)) {
    stop("Column 'RT' is missing")
  }

  # enough columns?
  if (ncol(x) < n + 3) {
    stop("Number of sample columns does not match n")
  }

  # sample columns should be numeric
  sample_cols <- 2:(n+1)

  if (!all(sapply(x[, sample_cols], is.numeric))) {
    stop("Sample abundance columns must be numeric")
  }

  return(TRUE)
}
