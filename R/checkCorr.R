#' @title Correlation of two peaks at the same retention time
#'
#' @description
#' Calculates the Spearman correlation coefficient between two specific peaks defined by their
#' mass-to-charge ratio (m/z) and target Retention Time (RT). The peaks can be queried from
#' a single dataset or across two different datasets (e.g., matching positive and negative ion modes).
#'
#' @usage checkCorr(x, y = NULL, masses, RT, n, tolerance = 0.05)
#'
#' @param x Data frame. Format must be: Peak identifiers in column 1, followed by \code{n} sample columns, "m.z", and "RT".
#' @param y Data frame (optional). A second data frame with compatible dimensions to \code{x}. If \code{NULL}, defaults to \code{y = x}.
#' @param masses Numeric vector of length 2 containing the target m/z values to compare \code{c(mass1, mass2)}.
#' @param RT Numeric. The target retention time around which to search for the peaks.
#' @param n Numeric. The exact number of sample abundance columns present in the dataset(s).
#' @param tolerance Numeric. Mass tolerance window (in Da) used to search for the closest m/z match. Default is \code{0.05}.
#'
#' @details
#' The function locates the closest peak matching \code{masses[1]} in dataset \code{x} and \code{masses[2]}
#' in dataset \code{y} within the specified mass tolerance. If multiple matching peaks are found, it isolates
#' the single feature closest to the target retention time (\code{RT}). Missing values are handled using
#' pairwise complete observations.
#'
#' @returns A numeric value indicating the Spearman correlation coefficient between the two target peaks.
#'
#' @examples
#' # Check correlation within the same dataset
#' # checkCorr(pos, masses = c(595.2, 449.08), RT = 11.8, n = 48)
#'
#' # Check cross-mode correlation between positive and negative datasets
#' # checkCorr(pos, neg, masses = c(595.2, 593.15), RT = 11.8, n = 48)
#'
#' @importFrom stats cor
#' @export
checkCorr <- function(x, y = NULL, masses, RT, n, tolerance = 0.05) {

  # Set default data mapping behavior if y is omitted
  if (is.null(y)) {
    y <- x
  }

  # Identify target abundance coordinates
  sample_cols <- 2:(n + 1)
  mass1 <- masses[1]
  mass2 <- masses[2]

  # ----------------------------------------------------------------------------
  # 1. PROCESS FIRST PEAK (Dataset X)
  # ----------------------------------------------------------------------------
  # Find all peaks within the mass tolerance window
  matches_x <- x[abs(x$m.z - mass1) <= tolerance, , drop = FALSE]

  if (nrow(matches_x) == 0) {
    stop(paste("The first mass (", mass1, ") was not found within a tolerance of ", tolerance, " Da.", sep = ""))
  }

  # From matches, isolate the one closest to the target Retention Time
  best_row_x <- matches_x[which.min(abs(matches_x$RT - RT)), , drop = FALSE]

  # Extract the abundance values as a clean numeric vector
  vec1 <- as.numeric(best_row_x[, sample_cols])

  # ----------------------------------------------------------------------------
  # 2. PROCESS SECOND PEAK (Dataset Y)
  # ----------------------------------------------------------------------------
  # Find all peaks within the mass tolerance window
  matches_y <- y[abs(y$m.z - mass2) <= tolerance, , drop = FALSE]

  if (nrow(matches_y) == 0) {
    stop(paste("The second mass (", mass2, ") was not found within a tolerance of ", tolerance, " Da.", sep = ""))
  }

  # Fix the original typo: Evaluate against matches_y$RT, not x!
  best_row_y <- matches_y[which.min(abs(matches_y$RT - RT)), , drop = FALSE]

  # Extract the abundance values as a clean numeric vector
  vec2 <- as.numeric(best_row_y[, sample_cols])

  # ----------------------------------------------------------------------------
  # 3. COMPUTE SPEARMAN CORRELATION
  # ----------------------------------------------------------------------------
  # Replace legacy 0 markers with NA to handle missing pairs elegantly
  vec1[vec1 == 0] <- NA
  vec2[vec2 == 0] <- NA

  # Compute the correlation coefficient
  out_corr <- stats::cor(
    x = vec1,
    y = vec2,
    method = 'spearman',
    use = 'pairwise.complete.obs'
  )

  return(out_corr)
}
