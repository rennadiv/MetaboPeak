#' Filter peaks in a retention time range
#'
#' @description
#' This function reduces the original data frame to retention times (RT) that
#' are biologically relevant, removing background noise from the beginning
#' and the end of the run.
#'
#' @usage interestRegion(x, start, end)
#'
#' @param x A data frame containing chromatography data.
#' @param start A numeric value representing the starting RT window limit.
#' @param end A numeric value representing the ending RT window limit.
#'
#' @details
#' The input data frame `x` must contain a column explicitly named \code{RT}.
#' A typical format looks like:
#' \code{Peak names | abundances of samples | m.z | RT | ...}
#'
#' @returns A filtered data frame containing only rows within the specified RT range.
#'
#' @examples
#' # Assuming 'neg' is a data frame with an 'RT' column:
#' # neg_small <- interestRegion(neg, 1.4, 25)
#'
#' @export

interestRegion <- function(x, start, end) {
  # Safety check: Ensure start and end are numbers
  if (!is.numeric(start) || !is.numeric(end)) {
    stop("Both 'start' and 'end' must be numeric values.")
  }

  # Safety check: Ensure the RT column exists
  if (!"RT" %in% colnames(x)) {
    stop("The input data frame does not have a column named 'RT'.")
  }

  # Filter the data frame
  x.new <- x[x$RT >= start & x$RT <= end, ]

  return(x.new)
}
