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
#' @usage peakFilter(x, y, fNA, fCV, fRT, parallel, n_cors)
#' @param x data frame
#' @param y treatment info data table
#' @param fNA character vector
#' @param fCV character vector
#' @param fRT character vector
#' @parallel whether to check multiple groups at the same time on multiple CPU cores
#' @n_cors number of CPU cores your function should use when running in parallel
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
#' parallel can be either defined as F/FALSE or T/TRUE
#'
#' n_cores can be set as NULL and the computer will use number of available cores -1, or it can be set as any number that would work best for
#' used computer.
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
#' @import matrixStats
#' @importFrom igraph graph_from_data_frame
#' @importFrom igraph components
#' @importFrom stats cor
#'
#' @export
#'

peakFilter <- function(
    x, y,
    fNA = c('T', 0.5),
    fCV = c('T', 0.8),
    fRT = c('T', 3, 0.02, 0.95),
    parallel = FALSE,
    n_cores = NULL
) {

  # --- SETUP PARALLEL ---
  if (parallel) {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      stop("Please install 'future.apply' for parallel execution")
    }
    if (is.null(n_cores)) {
      n_cores <- future::availableCores() - 1
    }
    future::plan(future::multisession, workers = n_cores)
  }

  # --- LOAD FAST MATRIX OPS ---
  if (!requireNamespace("matrixStats", quietly = TRUE)) {
    stop("Please install 'matrixStats'")
  }

  # --- PREP ---
  x[x == 0] <- NA
  if (is.na(x[1,1])) x[1,1] <- 0

  rownames(x) <- x[,1]

  Z <- x[, y[,1], drop = FALSE]
  rownames(Z) <- x[,1]
  colnames(Z) <- paste(colnames(Z), y$treatment, sep = '_')

  y$treatment <- as.factor(y$treatment)
  Z <- as.matrix(Z)

  group_cols <- split(seq_along(y$treatment), y$treatment)

  has_complete_group <- function(mat) {
    keep_matrix <- sapply(group_cols, function(cols) {
      rowSums(!is.na(mat[, cols, drop = FALSE])) == length(cols)
    })
    rowSums(keep_matrix) > 0
  }

  # --- NA FILTER ---
  if (fNA[1] == 'T') {
    threshold <- as.numeric(fNA[2])

    pctNAs <- rowMeans(is.na(Z))
    trash <- pctNAs > threshold

    keep <- !trash | has_complete_group(Z)
    Z <- Z[keep, , drop = FALSE]
  }

  # --- CV FILTER ---
  if (fCV[1] == 'T') {
    threshold <- as.numeric(fCV[2])

    Mean <- matrixStats::rowMeans2(Z, na.rm = TRUE)
    SD   <- matrixStats::rowSds(Z, na.rm = TRUE)
    CV   <- SD / Mean

    trash <- CV > threshold

    keep <- !trash | has_complete_group(Z)
    Z <- Z[keep, , drop = FALSE]
  }

  # --- RT FILTER ---
  if (fRT[1] == 'T') {

    digits     <- as.numeric(fRT[2])
    rt_window  <- as.numeric(fRT[3])
    cor_thresh <- as.numeric(fRT[4])

    rt <- round(x[rownames(Z), "RT"], digits)
    mz <- x[rownames(Z), "m.z"]

    ord <- order(rt)
    rt <- rt[ord]
    Z  <- Z[ord, , drop = FALSE]
    mz <- mz[ord]

    # --- build RT windows ---
    groups <- list()
    start <- 1
    n <- length(rt)

    while (start <= n) {
      end <- start
      while (end < n && (rt[end+1] - rt[start]) <= rt_window) {
        end <- end + 1
      }
      groups[[length(groups) + 1]] <- start:end
      start <- end + 1
    }

    # --- FUNCTION TO PROCESS ONE GROUP ---
    process_group <- function(idx) {

      if (length(idx) <= 1) return(rep(TRUE, length(idx)))

      sub <- Z[idx, , drop = FALSE]

      cmat <- suppressWarnings(
        cor(t(sub), method = "spearman", use = "pairwise.complete.obs")
      )

      cmat[lower.tri(cmat, diag = TRUE)] <- NA
      pairs <- which(cmat > cor_thresh, arr.ind = TRUE)

      keep_local <- rep(TRUE, length(idx))

      if (nrow(pairs) > 0) {
        groups_local <- split(pairs[,1], pairs[,2])

        for (g in groups_local) {
          rows <- unique(g)

          vals <- matrixStats::rowMaxs(sub[rows, , drop = FALSE], na.rm = TRUE)
          best <- rows[which.max(vals)]

          keep_local[setdiff(rows, best)] <- FALSE
        }
      }

      keep_local
    }

    # --- APPLY (parallel or not) ---
    if (parallel) {
      keep_list <- future.apply::future_lapply(groups, process_group)
    } else {
      keep_list <- lapply(groups, process_group)
    }

    # --- COMBINE RESULTS ---
    keep <- rep(FALSE, nrow(Z))

    for (i in seq_along(groups)) {
      keep[groups[[i]]] <- keep_list[[i]]
    }

    Z <- Z[keep, , drop = FALSE]
  }

  # --- CLEAN COLUMN NAMES ---
  colnames(Z) <- sub("_.*$", "", colnames(Z))

  # --- OUTPUT ---
  Z_out <- data.frame(ID = rownames(Z), Z, check.names = FALSE)
  Z_out$m.z <- x[Z_out$ID, "m.z"]
  Z_out$RT  <- x[Z_out$ID, "RT"]

  return(Z_out)
}
