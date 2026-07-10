#' Data format control
#'
#' @description
#' Given the data table and number of samples, the function checks if the table
#' is in the correct format for downstream MetaboPeak functions.
#'
#' @usage checkFormat(x, n)
#' @param x A data frame.
#' @param n An integer specifying the expected number of sample columns.
#'
#' @details
#' This function validates the structure of LC-MS peak tables.
#'
#' It checks whether:
#' \itemize{
#'   \item The input is a data.frame (or tibble).
#'   \item Required columns ("m.z", "RT") are present.
#'   \item The remaining columns match the expected sample count \code{n}.
#'   \item Sample abundance columns contain strictly numeric values.
#' }
#'
#' @returns Logical scalar with the value \code{TRUE} if everything is in order,
#' or throws an error specifying where the problem lies.
#'
#' @export
#'
#' @examples
#' \donttest{
#' data("neg", package = "MetaboPeak")
#' checkFormat(neg, 48)
#' }
checkFormat <- function(x, n) {

  # 1. Must be a data frame or tibble
  if (!is.data.frame(x)) {
    stop("Input must be a data.frame")
  }

  # 2. Check for required metadata columns
  if (!"m.z" %in% colnames(x)) {
    stop("Column 'm.z' is missing")
  }

  if (!"RT" %in% colnames(x)) {
    stop("Column 'RT' is missing")
  }

  # 3. Dynamic Column Identification
  # Instead of hardcoding 2:(n+1), we explicitly isolate sample columns
  all_cols <- colnames(x)
  metadata_cols <- c("m.z", "RT", "Peak names", "peak_id") # common labels to exclude

  # Sample columns are everything that isn't explicitly m.z or RT
  # (and optionally dropping an ID column if it exists)
  sample_cols <- all_cols[!(all_cols %in% c("m.z", "RT"))]

  # If there is a peak ID column at the start, drop it from the sample list
  if (length(sample_cols) > n && sample_cols[1] %in% c("Peak names", "peak_id", "ID")) {
    sample_cols <- sample_cols[-1]
  }

  # 4. Verify the number of sample columns matches 'n' exactly
  if (length(sample_cols) != n) {
    stop(paste0("Expected ", n, " sample columns, but found ", length(sample_cols), "."))
  }

  # 5. Check that all identified sample columns are strictly numeric
  # Loop works perfectly for both base R data.frames and tidyverse tibbles
  for (col in sample_cols) {
    if (!is.numeric(x[[col]])) {
      stop(paste0("Sample abundance column '", col, "' must be numeric"))
    }
  }

  return(TRUE)
}
