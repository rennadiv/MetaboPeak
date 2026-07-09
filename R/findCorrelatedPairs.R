#' @title Find and rank highly correlated peak pairs within RT windows
#'
#' @description
#' Automatically scans a dataset (or compares two datasets) to find pairs of peaks that elute
#' within a shared Retention Time (RT) window and exhibit a strong Spearman correlation. This
#' is highly useful for identifying related metabolic features like adducts, isotopes, or fragments.
#'
#' @usage findCorrelatedPairs(x, y = NULL, n, rt_window = 0.05, cor_threshold = 0.85)
#'
#' @param x Data frame. Format must be: Peak identifiers in column 1, followed by \code{n} sample columns, "m.z", and "RT".
#' @param y Data frame (optional). A second data frame with compatible dimensions to \code{x} (e.g., negative mode). If \code{NULL}, defaults to \code{y = x}.
#' @param n Numeric. The exact number of sample abundance columns present in the dataset(s).
#' @param rt_window Numeric. The maximum allowable retention time difference (in minutes) between two peaks to be compared. Default is \code{0.05}.
#' @param cor_threshold Numeric. The minimum Spearman correlation coefficient required to report a pair. Default is \code{0.85}.
#'
#' @returns A data frame sorted by correlation strength, showing the IDs, m/z values, RT values, RT differences, and correlation scores of matching pairs.
#'
#' @importFrom stats cor
#' @export


findCorrelatedPairs <- function(x, y = NULL, n, rt_window = 0.05, cor_threshold = 0.85) {

  # 1. SETUP DATASETS
  same_dataset <- is.null(y)
  if (same_dataset) {
    y <- x
  }

  sample_cols <- 2:(n + 1)

  # Clean up zero values to NAs for accurate calculations
  mat_x <- as.matrix(x[, sample_cols, drop = FALSE])
  mat_x[mat_x == 0] <- NA

  mat_y <- as.matrix(y[, sample_cols, drop = FALSE])
  mat_y[mat_y == 0] <- NA

  # Extract metadata arrays
  ids_x <- as.character(x[, 1])
  mz_x  <- as.numeric(x$m.z)
  rt_x  <- as.numeric(x$RT)

  ids_y <- as.character(y[, 1])
  mz_y  <- as.numeric(y$m.z)
  rt_y  <- as.numeric(y$RT)

  # 2. FIND CANDIDATE PAIRS CLOSES IN TIME (Non-overlapping sliding approach)
  # Instead of calculating millions of correlations, we find which index rows are actually near each other in time
  pairs_list <- list()

  for (i in seq_along(rt_x)) {
    # Find indices in Y that are within the RT window
    matching_y_indices <- which(abs(rt_y - rt_x[i]) <= rt_window)

    # If comparing a dataset with itself, avoid duplicate pairs (A-B and B-A) and self-comparison (A-A)
    if (same_dataset) {
      matching_y_indices <- matching_y_indices[matching_y_indices > i]
    }

    if (length(matching_y_indices) > 0) {
      pairs_list[[length(pairs_list) + 1]] <- data.frame(
        idx_x = i,
        idx_y = matching_y_indices,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(pairs_list) == 0) {
    message("No peak pairs were found within the specified Retention Time window.")
    return(data.frame())
  }

  # Combine candidate index mappings
  all_candidates <- do.call(rbind, pairs_list)

  # 3. CALCULATE CORRELATIONS FOR THE CANDIDATES
  cor_scores <- numeric(nrow(all_candidates))

  # Calculate row-by-row Spearman correlations for valid candidates
  for (k in seq_len(nrow(all_candidates))) {
    row_x <- all_candidates$idx_x[k]
    row_y <- all_candidates$idx_y[k]

    cor_scores[k] <- stats::cor(
      x = mat_x[row_x, ],
      y = mat_y[row_y, ],
      method = "spearman",
      use = "pairwise.complete.obs"
    )
  }

  # 4. FILTER AND FORMAT THE RESULTS
  all_candidates$correlation <- cor_scores

  # Keep only pairs exceeding our threshold (ignoring missing math values)
  valid_pairs <- all_candidates[!is.na(all_candidates$correlation) & all_candidates$correlation >= cor_threshold, ]

  if (nrow(valid_pairs) == 0) {
    message("No pairs exceeded the specified correlation threshold.")
    return(data.frame())
  }

  # Construct final summary output frame using original metadata values
  results <- data.frame(
    Peak_X = ids_x[valid_pairs$idx_x],
    Peak_Y = ids_y[valid_pairs$idx_y],
    m.z_X  = mz_x[valid_pairs$idx_x],
    m.z_Y  = mz_y[valid_pairs$idx_y],
    RT_X   = rt_x[valid_pairs$idx_x],
    RT_Y   = rt_y[valid_pairs$idx_y],
    RT_Diff = abs(rt_x[valid_pairs$idx_x] - rt_y[valid_pairs$idx_y]),
    Correlation = valid_pairs$correlation,
    stringsAsFactors = FALSE
  )

  # Order results by strongest correlation first
  results <- results[order(-results$Correlation), ]
  rownames(results) <- NULL

  return(results)
}
