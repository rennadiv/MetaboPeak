#' @title Find the sample name with the highest abundance for a given mass
#'
#' @description
#' Searches the dataset for a target mass-to-charge ratio (m/z) and identifies the specific
#' sample column that contains the absolute highest abundance value for that feature.
#'
#' @usage highestAbund(x, mass, n, tolerance = 0.02)
#'
#' @param x Data frame. Format must be: Peak identifiers in column 1, followed by \code{n} sample columns, "m.z", and "RT".
#' @param mass Numeric or Character. The target mass-to-charge ratio (m/z) value to search for.
#' @param n Numeric. The exact number of sample abundance columns present in the dataset.
#' @param tolerance Numeric. Mass tolerance window (in Da) around the target m/z. Default is \code{0.02}.
#'
#' @details
#' If multiple peaks are found within the specified mass tolerance window, the function
#' scans all of them and extracts the sample name corresponding to the absolute maximum
#' intensity across all matching rows, preventing unexpected function crashes.
#'
#' @returns A character string representing the name of the sample with the highest peak intensity.
#'
#' @examples
#' # Find the sample where m/z 593.128 is most abundant
#' # top_sample <- highestAbund(neg, mass = 593.128, n = 48)
#'
#' @export


highestAbund <- function(x, mass, n, tolerance = 0.02) {

  # ---- 1. CORE CHECKS ----
  if (!is.data.frame(x)) {
    stop("Input parameter 'x' must be a structured data.frame.")
  }
  if (!is.numeric(n) || length(n) != 1) {
    stop("Parameter 'n' must be a single numeric value representing total samples.")
  }

  # ---- 2. MASS TARGET FILTERING ----
  target_mass <- as.numeric(mass)
  sample_cols <- 2:(n + 1)

  # Find all peaks falling inside the tolerance envelope
  matched_rows <- x[abs(x$m.z - target_mass) <= tolerance, , drop = FALSE]

  if (nrow(matched_rows) == 0) {
    stop(paste("The target mass (", target_mass, ") was not found within a tolerance of ", tolerance, " Da.", sep = ""))
  }

  # ---- 3. EXTRACT THE ABSOLUTE MAXIMUM COORDINATE ----
  # Isolate only the sample abundance matrix space
  abundance_mat <- as.matrix(matched_rows[, sample_cols, drop = FALSE])

  # Replace any 0 values with NA so they don't interfere with math checks
  abundance_mat[abundance_mat == 0] <- NA

  # Find the exact 2D index (row, column) of the highest overall value
  max_coord <- which(abundance_mat == max(abundance_mat, na.rm = TRUE), arr.ind = TRUE)

  # Extract the column name using the column index from our coordinate
  # (Using the first match if there is an exact tie down to the decimal point)
  best_sample_index <- max_coord[1, "col"]
  top_sample_name <- colnames(abundance_mat)[best_sample_index]

  return(top_sample_name)
}
