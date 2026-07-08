#' @title Peak information and statistics
#'
#' @description
#' This function computes comprehensive descriptive statistics for mass spectrometry peaks
#' across all sample columns (including Retention Time, m/z, percentage of NAs, maximum
#' intensity, mean, standard deviation, and CV coefficient). Statistics can be returned
#' for all peaks or filtered down to a single target peak identifier.
#'
#' @usage peakInfo(x, n, ID)
#'
#' @param x Data frame. Format must be: Peak identifiers in column 1, followed by \code{n} sample columns, "m.z", and "RT".
#' @param n Numeric. The exact number of sample abundance columns present in the dataset.
#' @param ID Character or Numeric (optional). A specific peak identifier (from column 1) to retrieve data for.
#'
#' @details
#' The input abundance data values equal to 0 are internally transformed to \code{NA} to
#' guarantee accurate calculation of true missingness, means, and variances.
#'
#' @returns A clean data frame containing statistical summaries with maintained numeric data types.
#'
#' @examples
#' # Calculate statistics for all features assuming 48 samples
#' # stats_summary <- peakInfo(neg, n = 48)
#'
#' # Isolate data for specific peak ID "Peak_20"
#' # single_peak <- peakInfo(neg, n = 48, ID = "20")
#'
#' @importFrom matrixStats rowMeans2 rowSds rowMaxs
#' @export


peakInfo <- function(x, n, ID = NULL) {

  # --- 1. CORE ACCELERATION CHECKS ---
  if (!requireNamespace("matrixStats", quietly = TRUE)) {
    stop("Please install the 'matrixStats' package to run this optimized function.")
  }

  # --- 2. DATA SEGMENTATION & CLEANING ---
  # Extract the abundance columns cleanly (Columns 2 to n + 1)
  sample_cols <- 2:(n + 1)
  z_mat <- as.matrix(x[, sample_cols, drop = FALSE])

  # Treat 0 values as missing values (NA) for true analytical calculations
  z_mat[z_mat == 0] <- NA

  # Extract explicit structural metadata features
  peak_ids <- as.character(x[, 1])
  m.z      <- as.numeric(x[, 'm.z'])
  RT       <- as.numeric(x[, 'RT'])

  # --- 3. HIGH SPEED SUMMARY MATH CALCULATIONS ---
  # Math calculation of missingness fraction per row (Fixed the length denominator bug)
  pctNAs <- round(rowMeans(is.na(z_mat)), 2)

  # Fast, compiled C-level matrix stats loops
  Max  <- matrixStats::rowMaxs(z_mat, na.rm = TRUE)
  Mean <- round(matrixStats::rowMeans2(z_mat, na.rm = TRUE), 4)
  SD   <- round(matrixStats::rowSds(z_mat, na.rm = TRUE), 4)
  CV   <- round(SD / Mean, 3)

  # Clean up mathematical infinities or NaNs resulting from completely empty rows
  Max[is.infinite(Max)] <- NA

  # --- 4. ASSEMBLE CLEAN NON-DESTRUCTIVE DATA FRAME ---
  # Creating standard data frame columns preserves true data types (Numbers stay numbers)
  x_new <- data.frame(
    ID     = peak_ids,
    m.z    = m.z,
    RT     = RT,
    pctNAs = pctNAs,
    Max    = Max,
    Mean   = Mean,
    SD     = SD,
    CV     = CV,
    stringsAsFactors = FALSE
  )

  # Set row names to easily match matching records
  rownames(x_new) <- peak_ids

  # --- 5. TARGET SELECTION FILTER ---
  if (is.null(ID)) {
    return(x_new)
  } else {
    target_id <- as.character(ID)
    if (!target_id %in% rownames(x_new)) {
      stop(paste("Requested ID '", target_id, "' could not be found in the dataset.", sep = ""))
    }
    return(x_new[target_id, , drop = FALSE])
  }
}
