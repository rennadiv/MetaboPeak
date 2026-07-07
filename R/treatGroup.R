#' @title Unique code for combined treatment
#' @description
#' Based on the given abbreviation for the names of treatment, it makes an unique code for each treatment combination to which the sample is
#' exposed.
#'
#' @usage treatGroup(x, number.of.treatments, my_abbreviations)
#' @param x data frame with treatment info
#' @param number.of.treatments integer of treatment columns (types of treatment)
#' @param my_abbreviations character vector of the abbreviations you want to use in alphabetical order (optional)
#'
#' @details
#' The format of the data frame should be: sample name|treatment groups|...
#'
#' The vector of abbreviations has to be in alphabetical order (small cases first). If missing, then the treatment code will be formed
#' from the names of the treatments combined by underscore.
#'
#'
#' @returns x with new column with the treatment code
#'
#' @examples
#' abb <- c('B','S','A','E','D','W','n','N')
#' treatGroup(t_info, 4, abb)
#'
#' treatGroup(t_info, 4)
#'
#' @export
#'


treatGroup <- function(x, number.of.treatments, my_abbreviations = NULL) {
  # 1. Safely identify the treatment columns (Columns 2 to number.of.treatments + 1)
  target_cols <- 2:(number.of.treatments + 1)

  # 2. Extract just the columns we want to work with
  treatment_data <- x[, target_cols, drop = FALSE]

  if (is.null(my_abbreviations)) {
    # If no abbreviations, paste columns together separated by "_"
    x$Treatment <- interaction(treatment_data, sep = "_")

  } else {
    # If abbreviations exist, map them over the factor levels efficiently
    # Split abbreviations into a list matching the number of levels per column
    col_level_counts <- sapply(treatment_data, function(col) length(unique(col)))

    # Check if the user provided the correct total number of abbreviations
    if (sum(col_level_counts) != length(my_abbreviations)) {
      stop("The number of abbreviations does not match the total number of unique factor levels.")
    }

    # Split the abbreviation vector into a list for each column
    abbrev_groups <- split(my_abbreviations, rep(1:length(col_level_counts), col_level_counts))

    # Apply abbreviations to each column safely
    mapped_data <- as.data.frame(mapply(function(col, abbrev) {
      factor(col, labels = abbrev)
    }, treatment_data, abbrev_groups, SIMPLIFY = FALSE))

    # Combine the abbreviated columns together without any separators
    x$Treatment <- interaction(mapped_data, sep = "")
  }

  return(x)
}

