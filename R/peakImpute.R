#' Impute missing and low-intensity values
#'
#' @description
#' Replaces missing values and/or values below a specified threshold with a
#' fraction of the minimum detected abundance for each peak.
#'
#' The replacement value is calculated independently for each peak as:
#' \code{min(valid values) * fraction}
#' where valid values are those greater than \code{threshold}.
#'
#' @usage peakImpute(x, threshold, fraction, replace_low, replace_na)
#'
#' @param x Data frame. First column contains peak IDs, followed by abundance
#'   columns. Any additional annotation columns (e.g. m.z, RT) are preserved.
#' @param threshold Numeric. Values less than or equal to this threshold are
#'   considered undetected. Default = 1000.
#' @param fraction Numeric. Fraction of the minimum detected value used for
#'   imputation. Default = 0.1.
#' @param replace_low Logical. Replace values <= threshold. Default = TRUE.
#' @param replace_na Logical. Replace NA values. Default = TRUE.
#'
#' @returns A data frame with imputed abundance values, preserving original row and column layouts.
#'
#' @examples
#' # Example data setup assumed (neg dataset)
#' # imputed_data <- peakImpute(neg, threshold = 1000, fraction = 0.1)
#'
#' @export
peakImpute <- function(
    x,
    threshold = 1000,
    fraction = 0.1,
    replace_low = TRUE,
    replace_na = TRUE
) {

  # 1. Coerce to standard data.frame to eliminate modern Tidyverse/Tibble assignment bugs
  out <- as.data.frame(x)

  # 2. Identify numeric tracks safely
  abundance_cols <- sapply(out, is.numeric)

  # CRITICAL PROTECTION: The first column contains identities (Peak IDs)
  # and must NEVER be treated as an abundance matrix even if formatted as numbers.
  abundance_cols[1] <- FALSE

  # Exclude known metadata attributes explicitly
  abundance_cols[match(c("m.z", "RT"), names(out), nomatch = 0)] <- FALSE

  # 3. Process matrix operations via standard base R arrays
  abundance <- as.matrix(out[, abundance_cols, drop = FALSE])

  abundance <- t(apply(abundance, 1, function(v) {

    valid <- v[v > threshold & !is.na(v)]

    # If the feature is entirely missing across the cohort, return as is
    if (length(valid) == 0) {
      return(v)
    }

    # Calculate localized LOD fraction replacement value
    replacement <- min(valid) * fraction

    if (replace_low) {
      v[v <= threshold & !is.na(v)] <- replacement
    }

    if (replace_na) {
      v[is.na(v)] <- replacement
    }

    return(v)
  }))

  # 4. Safely push numeric matrices back into coordinates
  out[, abundance_cols] <- abundance

  return(out)
}
