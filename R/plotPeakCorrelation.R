#' @title Plot intensity profiles of two correlated peaks
#'
#' @description
#' Generates a line plot comparing the intensity profiles of two peaks across all samples.
#' This is ideal for visually verifying co-elution, tracking adducts, or confirming cross-mode peak links.
#'
#' @usage plotPeakCorrelation(x, y = NULL, peak_id_x, peak_id_y, n)
#'
#' @param x Data frame. Format must be: Peak identifiers in column 1, followed by \code{n} sample columns, "m.z", and "RT".
#' @param y Data frame (optional). A second data frame if comparing across modes (e.g., negative mode). If \code{NULL}, defaults to \code{y = x}.
#' @param peak_id_x Character. The unique ID of the first peak (from column 1 of \code{x}).
#' @param peak_id_y Character. The unique ID of the second peak (from column 1 of \code{y}).
#' @param n Numeric. The exact number of sample abundance columns present in the dataset(s).
#'
#' @returns None. Generates a clean line plot in the active R graphics device.
#'
#' @importFrom stats cor
#' @importFrom graphics plot lines points legend
#' @export

plotPeakCorrelation <- function(x, y = NULL, peak_id_x, peak_id_y, n) {

  if (is.null(y)) {
    y <- x
  }

  sample_cols <- 2:(n + 1)

  # 1. Isolate the target rows
  row_x <- x[x[, 1] == peak_id_x, , drop = FALSE]
  row_y <- y[y[, 1] == peak_id_y, , drop = FALSE]

  if (nrow(row_x) == 0) stop(paste("Peak ID '", peak_id_x, "' not found in dataset x.", sep = ""))
  if (nrow(row_y) == 0) stop(paste("Peak ID '", peak_id_y, "' not found in dataset y.", sep = ""))

  # 2. Extract values and handle 0/NA
  vec_x <- as.numeric(row_x[, sample_cols])
  vec_y <- as.numeric(row_y[, sample_cols])

  vec_x[vec_x == 0] <- NA
  vec_y[vec_y == 0] <- NA

  # Calculate correlation for the plot subtitle
  current_cor <- stats::cor(vec_x, vec_y, method = "spearman", use = "pairwise.complete.obs")

  # 3. Setup clean plot boundaries using normalized scales so different intensities overlap nicely
  # We scale them 0 to 1 just for trend visualization
  scale_vector <- function(v) {
    if (all(is.na(v)) || max(v, na.rm = TRUE) == min(v, na.rm = TRUE)) return(v)
    (v - min(v, na.rm = TRUE)) / (max(v, na.rm = TRUE) - min(v, na.rm = TRUE))
  }

  scaled_x <- scale_vector(vec_x)
  scaled_y <- scale_vector(vec_y)

  # 4. Generate the Base R Plot
  graphics::plot(
    1:n, scaled_x,
    type = "b", col = "#1f77b4", pch = 16, lwd = 2,
    ylim = c(0, 1.1), xaxt = "n",
    xlab = "Samples", ylab = "Relative Intensity (Normalized)",
    main = paste("Profile Comparison:", peak_id_x, "vs", peak_id_y)
  )

  graphics::lines(1:n, scaled_y, type = "b", col = "#ff7f0e", pch = 17, lwd = 2)

  # Add sample index labels to the X-axis cleanly
  graphics::axis(1, at = 1:n, labels = 1:n)

  # Add metadata tracking subtitle details
  graphics::mtext(
    paste("Spearman r =", round(current_cor, 3),
          " | m/z X:", round(row_x$m.z, 3), " RT X:", round(row_x$RT, 2),
          " | m/z Y:", round(row_y$m.z, 3), " RT Y:", round(row_y$RT, 2)),
    side = 3, line = 0.25, cex = 0.8, col = "gray30"
  )

  # Add legend
  graphics::legend(
    "topright",
    legend = c(peak_id_x, peak_id_y),
    col = c("#1f77b4", "#ff7f0e"),
    lty = 1, lwd = 2, pch = c(16, 17),
    bty = "n"
  )
}
