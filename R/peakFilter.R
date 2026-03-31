#' @title Filter by NA, CV coefficient and/or correlated RT
#' @description
#' This function provides 3 types of filters. It needs the abundances table, treatment info with combined treatment code under the column
#' treatment (function treatGroup can be used to make this code) and some thresholds.
#'
#' 1) based on NA values: it calculates which peaks (rows) have more NAs than
#' threshold and if there are no complete cases for any treatment group, than it removes those rows.
#'
#' 2) based on CV coefficient: it calculates which peaks (rows) have bigger CV coefficient than
#' threshold and if there are no complete cases for any treatment group, than it removes those rows.
#'
#' 3) based on correlated RT: it rounds the RT and orders them, than based on thresholds groups similar retention times and calculates correlation
#' (spearman). From this correlated group it picks the one row with max m/z. (it might take some time)
#'
#' @usage peakFilter(x, y, fNA, fCV, fRT)
#' @param x data frame
#' @param y treatment info data table
#' @param fNA character vector
#' @param fCV character vector
#' @param fRT character vector
#'
#' @details
#' The format of the data frame should be: Peak names|abundances of samples|m.z|RT|...
#'
#' Treatment info must contain samples name, that corresponds with x in the 1. column and column treatment
#' (code for combination of treatments, use [treatGroup()])
#'
#' fNA and fCV can be either defined as F/FALSE or with 2 argument vector ('T', threshold). fRT can be either defined as F/FALSE, or
#' with 4 argument vector ('T', how many decimal places should RT be rounded, range value in which to look for similar RT, correlation
#' threshold - how well should the rows be correlated to take them as the same)
#'
#'
#' @returns A reduced data frame.
#'
#' @examples
#' peakFilter(neg, t_info_group, fNA = c('T',0.5), fCV = F, fRT = F)
#' peakFilter(neg, t_info_group, fNA = c('T',0.5), fCV = c('T', 0.4), fRT = F)
#' peakFilter(neg, t_info_group, fNA = c('T',0.5), fCV = c('T', 0.4), fRT = c('T',3,0.2,0.95)
#'
#'
#' @import dplyr
#' @importFrom igraph graph_from_data_frame
#' @importFrom igraph components
#' @importFrom stats cor
#'
#' @export
#'

peakFilter_fast <- function(
    x, y,
    fNA = c('T', 0.5),
    fCV = c('T', 0.8),
    fRT = c('T', 3, 0.02, 0.95)
) {

  # --- PREP ---
  x[x == 0] <- NA
  if (is.na(x[1,1])) x[1,1] <- 0

  rownames(x) <- x[,1]

  Z <- x[, y[,1], drop = FALSE]
  rownames(Z) <- x[,1]
  colnames(Z) <- paste(colnames(Z), y$treatment, sep = '_')

  y$treatment <- as.factor(y$treatment)

  # Convert to matrix for speed
  Z <- as.matrix(Z)

  # Precompute group column indices
  group_cols <- split(seq_along(y$treatment), y$treatment)

  # Helper: check if row has complete cases in ANY group
  has_complete_group <- function(mat) {
    keep_matrix <- sapply(group_cols, function(cols) {
      rowSums(!is.na(mat[, cols, drop = FALSE])) == length(cols)
    })
    rowSums(keep_matrix) > 0
  }

  # --- 1) NA FILTER ---
  if (fNA[1] == 'T') {
    threshold <- as.numeric(fNA[2])

    pctNAs <- rowMeans(is.na(Z))
    trash <- pctNAs > threshold

    keep_any_group <- has_complete_group(Z)

    keep <- !trash | keep_any_group
    Z <- Z[keep, , drop = FALSE]
  }

  # --- 2) CV FILTER ---
  if (fCV[1] == 'T') {
    threshold <- as.numeric(fCV[2])

    Mean <- rowMeans(Z, na.rm = TRUE)
    SD   <- apply(Z, 1, sd, na.rm = TRUE)
    CV   <- SD / Mean

    trash <- CV > threshold

    keep_any_group <- has_complete_group(Z)

    keep <- !trash | keep_any_group
    Z <- Z[keep, , drop = FALSE]
  }

  # --- 3) RT FILTER ---
  if (fRT[1] == 'T') {

    digits     <- as.numeric(fRT[2])
    rt_window  <- as.numeric(fRT[3])
    cor_thresh <- as.numeric(fRT[4])

    rt <- round(x[rownames(Z), "RT"], digits)
    mz <- x[rownames(Z), "m.z"]

    # Pre-sort once
    ord <- order(rt)
    rt <- rt[ord]
    Z  <- Z[ord, , drop = FALSE]
    mz <- mz[ord]

    keep <- rep(TRUE, nrow(Z))

    # Sliding window grouping (no repeated scanning)
    start <- 1
    n <- length(rt)

    while (start <= n) {
      end <- start

      # expand window
      while (end < n && (rt[end+1] - rt[start]) <= rt_window) {
        end <- end + 1
      }

      idx <- start:end

      if (length(idx) > 1) {
        sub <- Z[idx, , drop = FALSE]

        cmat <- suppressWarnings(
          cor(t(sub), method = "spearman", use = "pairwise.complete.obs")
        )

        cmat[lower.tri(cmat, diag = TRUE)] <- NA

        pairs <- which(cmat > cor_thresh, arr.ind = TRUE)

        if (nrow(pairs) > 0) {
          # simple grouping (avoid igraph)
          groups <- split(pairs[,1], pairs[,2])

          for (g in groups) {
            rows <- unique(c(g))

            # compute row max efficiently
            vals <- apply(Z[idx[rows], , drop = FALSE], 1, max, na.rm = TRUE)

            best_local <- rows[which.max(vals)]
            remove_local <- setdiff(rows, best_local)

            keep[idx[remove_local]] <- FALSE
          }
        }
      }

      start <- end + 1
    }

    Z <- Z[keep, , drop = FALSE]
  }

  # --- CLEAN COLUMN NAMES ---
  colnames(Z) <- sub("_.*$", "", colnames(Z))

  # --- RETURN FORMAT ---
  Z_out <- data.frame(ID = rownames(Z), Z, check.names = FALSE)
  Z_out$m.z <- x[Z_out$ID, "m.z"]
  Z_out$RT  <- x[Z_out$ID, "RT"]

  return(Z_out)
}
