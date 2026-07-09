#' @title Treatment groups with complete cases for a specific peak
#'
#' @description
#' Identifies which specific treatment groups have 100% complete data observations (zero missing
#' values or zeroes) for a target peak specified by its mass (m/z) and retention time (RT).
#'
#' @usage treatCases(x, y, mass, RT, n, mz_tolerance = 0.02, rt_tolerance = 0.1)
#'
#' @param x Data frame. Format must be: Peak identifiers in column 1, followed by \code{n} sample columns, "m.z", and "RT".
#' @param y Data frame. Treatment mapping table. Column 1 must match the sample column names in \code{x}. Must contain a column named \code{treatment}.
#' @param mass Numeric or Character. The target mass-to-charge ratio (m/z) value to search for.
#' @param RT Numeric. The target retention time (RT) value to search for.
#' @param n Numeric. The exact number of sample abundance columns present in the dataset.
#' @param mz_tolerance Numeric. Mass tolerance window (in Da) around the target m/z. Default is \code{0.02}.
#' @param rt_tolerance Numeric. Retention time window (in minutes) around the target RT. Default is \code{0.1}.
#'
#' @details
#' If multiple features reside within the overlapping mass and retention time tolerance envelopes,
#' the function isolates the single closest mathematical match to prevent output text crashes.
#'
#' @returns A character vector containing the names of the treatment groups where the peak is fully observed.
#'
#' @examples
#' # Find treatments where m/z 447.128 at RT 5.25 has complete cases
#' # groups <- treatCases(neg, t_info, mass = 447.128, RT = 5.25, n = 48)
#'
#' @export

treatCases <- function(x, y, mass, RT, n, mz_tolerance = 0.02, rt_tolerance = 0.1) {

  # ---- 1. CORE DATA CHECKS ----
  if (!is.data.frame(x)) {
    stop("Input parameter 'x' must be a structured data.frame.")
  }
  if (!is.data.frame(y)) {
    stop("Input parameter 'y' must be a structured data.frame.")
  }
  if (!"Treatment" %in% colnames(y)) {
    stop("The treatment mapping data frame 'y' must contain a column named 'treatment'.")
  }

  target_mass <- as.numeric(mass)
  target_rt   <- as.numeric(RT)
  sample_cols <- 2:(n + 1)

  # ---- 2. NUMERICAL FILTERING FOR METABOLITE FEATURE ----
  # Filter rows falling inside both the mass and retention time tolerance envelopes
  matched_rows <- x[abs(x$m.z - target_mass) <= mz_tolerance &
                      abs(x$RT - target_rt) <= rt_tolerance, , drop = FALSE]

  if (nrow(matched_rows) == 0) {
    stop(paste("No peak found matching m/z = ", target_mass, " and RT = ", target_rt, " within specified tolerances.", sep = ""))
  }

  # If multiple peaks match, safely pick the single closest one by combined distance
  if (nrow(matched_rows) > 1) {
    warning("Multiple peak entries found matching these coordinates. Selecting the single closest feature calculation.")
    mz_dist <- abs(matched_rows$m.z - target_mass)
    rt_dist <- abs(matched_rows$RT - target_rt)
    # Standard normalized Euclidean distance step
    total_dist <- (mz_dist / mz_tolerance) + (rt_dist / rt_tolerance)
    matched_rows <- matched_rows[which.min(total_dist), , drop = FALSE]
  }

  # ---- 3. CALCULATE GROUP MISSINGNESS PATTERNS ----
  # Extract abundances for the matched peak row
  peak_abundances <- as.numeric(matched_rows[1, sample_cols])

  # Treat legacy 0 values as missing values (NA)
  peak_abundances[peak_abundances == 0] <- NA

  # Pair sample labels to their respective abundance measurements explicitly
  sample_names <- colnames(x)[sample_cols]
  abundance_vector <- setNames(peak_abundances, sample_names)

  # Align the treatment information row tracking array to match our active abundance matrix sample sequence
  sample_order_in_y <- match(sample_names, as.character(y[, 1]))
  aligned_y <- y[sample_order_in_y, , drop = FALSE]

  treatments <- as.character(aligned_y$Treatment)

  # Group sample indices by their treatment label mappings
  group_indices <- split(seq_along(abundance_vector), treatments)

  # Check completeness across groups simultaneously using vectorized logical blocks
  complete_groups <- sapply(names(group_indices), function(g_name) {
    cols <- group_indices[[g_name]]
    group_values <- abundance_vector[cols]

    # A group is complete if it has exactly zero NA components
    all(!is.na(group_values))
  })

  # Filter down to return only the group strings that evaluated to TRUE
  elements <- names(complete_groups)[complete_groups]

  return(elements)
}
