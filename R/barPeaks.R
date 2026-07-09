#' @title Barplots of the 3 lowest variable peaks based on mass (m/z)
#'
#' @description
#' Calculates the Coefficient of Variation (CV) for all peaks matching a target mass
#' within a specified tolerance window. It then displays side-by-side or stacked barplots
#' showing sample abundance profiles for the up to 3 peaks with the lowest CV coefficients
#' (i.e., the most stable, least variable features).
#'
#' @usage barPeaks(x, mass, n, tolerance = 0.05)
#'
#' @param x Data frame. Format must be: Peak identifiers in column 1, followed by \code{n} sample columns, "m.z", and "RT".
#' @param mass Numeric or Character. The target mass-to-charge ratio (m/z) value to search for.
#' @param n Numeric. The exact number of sample abundance columns present in the dataset.
#' @param tolerance Numeric. Mass tolerance window (in Da) around the target m/z. Default is \code{0.05}.
#'
#' @details
#' Zero abundance values are treated as missing values (\code{NA}) to prevent data skewing
#' during mean and variance calculations.
#'
#' @returns A named list containing:
#' \itemize{
#'   \item \code{peak_names}: Character vector of the top selected peak IDs.
#'   \item \code{mids}: List containing the bar midpoint coordinates from each barplot.
#'   \item \code{CV}: Numeric vector of the calculated CV coefficients for the top peaks.
#'   \item \code{data}: Data frame of the raw abundance values used for plotting.
#' }
#'
#' @examples
#' # barPeaks(neg, '355.10', 48)
#'
#' @importFrom graphics barplot par
#' @importFrom matrixStats rowMeans2 rowSds
#'
#' @export


barPeaks <- function(x, mass, n, tolerance = 0.05) {

  # ---- 1. CORE ACCELERATION & VALUE CHECKS ----
  if (!is.data.frame(x)) {
    stop("Input parameter 'x' must be a structured data.frame.")
  }
  if (!is.numeric(n) || length(n) != 1) {
    stop("Parameter 'n' must be a single numeric value representing total samples.")
  }
  if (!requireNamespace("matrixStats", quietly = TRUE)) {
    stop("Please install the 'matrixStats' package to run this optimized function.")
  }

  # ---- 2. DYNAMIC LAYOUT MANAGEMENT ----
  # Capture user's active graphical environment settings
  old_par <- graphics::par(no.readonly = TRUE)
  # Automatically reset layouts back to the original state upon function termination
  on.exit(graphics::par(old_par), add = TRUE)

  # ---- 3. TARGET FEATURE FILTERING ----
  target_mass <- as.numeric(mass)
  sample_cols <- 2:(n + 1)

  # Filter rows falling cleanly within the target mass tolerance envelope
  matched_rows <- x[abs(x$m.z - target_mass) <= tolerance, , drop = FALSE]

  if (nrow(matched_rows) == 0) {
    stop(paste("No peaks were found matching mass ", target_mass, " within a tolerance of ", tolerance, " Da.", sep = ""))
  }

  # Separate raw abundance dimensions
  abundance_mat <- as.matrix(matched_rows[, sample_cols, drop = FALSE])
  abundance_mat[abundance_mat == 0] <- NA

  # ---- 4. SPEED-OPTIMIZED METRIC CALCULATIONS ----
  Mean <- matrixStats::rowMeans2(abundance_mat, na.rm = TRUE)
  SD   <- matrixStats::rowSds(abundance_mat, na.rm = TRUE)
  CV   <- SD / Mean

  # Handle instances where a completely empty peak yields NaN metrics
  CV[is.nan(CV)] <- Inf

  # Sort ascending to discover the absolute lowest variance peaks
  ord <- order(CV)
  keep_count <- min(3, length(ord))
  top_idx <- ord[1:keep_count]

  # Isolate selected plotting elements
  final_rows <- matched_rows[top_idx, , drop = FALSE]
  final_mat  <- abundance_mat[top_idx, , drop = FALSE]
  final_cv   <- CV[top_idx]
  final_ids  <- as.character(final_rows[, 1])

  # ---- 5. HIGH-VISIBILITY PLOTTING INTERFACE ----
  # Set layout matrix based dynamically on total matching records found (max 3)
  graphics::par(mfrow = c(keep_count, 1), mar = c(3, 4, 3, 2) + 0.1)

  mids_list <- list()

  for (i in seq_len(keep_count)) {
    # Generate clean bar chart distributions
    mids <- graphics::barplot(
      height   = final_mat[i, ],
      main     = paste("ID: ", final_ids[i],
                       " | m/z =", round(final_rows$m.z[i], 4),
                       " | RT =", round(final_rows$RT[i], 2),
                       " | CV =", round(final_cv[i], 3)),
      col      = "#1f77b4",
      border   = "white",
      las      = 1,
      names.arg = colnames(final_mat)
    )

    mids_list[[i]] <- mids
  }

  # ---- 6. COMPACT PAYLOAD RETURN ----
  return(list(
    peak_names = final_ids,
    mids       = mids_list,
    CV         = final_cv,
    data       = as.data.frame(final_mat)
  ))
}
