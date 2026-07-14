#' @title Filter by NA, CV coefficient and/or correlated RT
#'
#' @description
#' This function provides 3 types of filters for metabolomics peak data. It utilizes an
#' abundances table, sample-to-treatment mapping data, and customizable parameter thresholds.
#'
#' @usage peakFilter(x, y, fNA, fCV, fRT, protect_complete_groups, parallel, n_cores)
#'
#' @param x Data frame. Format must be: Peak names/IDs in column 1, followed by sample abundance columns, "m.z", and "RT".
#' @param y Treatment mapping data frame. Column 1 must match the sample column names in \code{x}. Must contain a column named \code{Treatment}.
#' @param fNA Vector. Either \code{FALSE} or a 2-element vector: \code{c('T', threshold_fraction)}. Default is \code{c('T', 0.5)}.
#' @param fCV Vector. Either \code{FALSE} or a 2-element vector: \code{c('T', threshold_cv)}. Default is \code{c('T', 0.8)}.
#' @param fRT Vector. Either \code{FALSE} or a 4-element vector: \code{c('T', rounding_digits, rt_window, correlation_threshold)}. Default is \code{c('T', 3, 0.02, 0.95)}.
#' @param protect_complete_groups Logical. If \code{TRUE}, peaks with 100\% complete measurements in at least one group are protected from being filtered out. Default is \code{TRUE}.
#' @param parallel Logical. Set to \code{TRUE} to check multiple retention time windows concurrently across multiple CPU cores. Default is \code{FALSE}.
#' @param n_cores Numeric. Number of CPU cores to utilize. If \code{NULL}, defaults to available cores minus 1.
#'
#' @returns A reduced data frame containing all original columns of \code{x}, filtered by the selected criteria.
#'
#' @import dplyr
#' @import matrixStats
#' @importFrom igraph graph_from_data_frame components
#' @importFrom stats cor
#'
#' @export
peakFilter <- function(
    x, y,
    fNA = c('T', 0.5),
    fCV = c('T', 0.8),
    fRT = c('T', 3, 0.02, 0.95),
    protect_complete_groups = TRUE,
    parallel = FALSE,
    n_cores = NULL
) {

  # --- REQUIRE CORE PACKAGES ---
  if (!requireNamespace("matrixStats", quietly = TRUE)) {
    stop("Please install 'matrixStats' to run this function.")
  }

  # --- PARALLEL SESSION MANAGEMENT ---
  if (parallel) {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      stop("Please install 'future.apply' for parallel execution capabilities.")
    }
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)

    if (is.null(n_cores)) {
      n_cores <- max(1, future::availableCores() - 1)
    }
    future::plan(future::multisession, workers = n_cores)
  }

  # --- INITIAL PROCESSING & CLEANING ---
  x_clean <- as.data.frame(x)
  y <- as.data.frame(y)

  sample_names <- as.character(y[, 1])
  abundance_matrix <- as.matrix(x_clean[, sample_names, drop = FALSE])

  # Ensure zeros are turned into explicit NA tokens
  abundance_matrix[abundance_matrix == 0] <- NA

  # Track initial state for summary calculations
  initial_peaks <- nrow(abundance_matrix)
  removed_na    <- 0
  removed_cv    <- 0
  removed_rt    <- 0

  # Ensure row identities match IDs
  rownames(abundance_matrix) <- as.character(x_clean[, 1])

  # Cache treatment groupings
  y$Treatment <- as.factor(y$Treatment)
  group_cols <- split(seq_along(y$Treatment), y$Treatment)

  # Helper function to detect rows with complete cases within any single treatment group
  has_complete_group <- function(mat) {
    if (!protect_complete_groups) {
      return(rep(FALSE, nrow(mat))) # Turn off protection if requested
    }
    keep_matrix <- sapply(group_cols, function(cols) {
      rowSums(!is.na(mat[, cols, drop = FALSE])) == length(cols)
    })
    if (is.vector(keep_matrix)) return(keep_matrix)
    rowSums(keep_matrix) > 0
  }

  # --- 1) NA FILTERING SCENARIO ---
  if (is.character(fNA) && fNA[1] == 'T') {
    # FIX: Explicitly target index [2] to avoid generating NAs from character values
    threshold <- as.numeric(fNA[2])
    pctNAs <- rowMeans(is.na(abundance_matrix))

    # Calculate group protections based on the masked NA matrix
    protected <- has_complete_group(abundance_matrix)
    keep_na <- (pctNAs <= threshold) | protected

    # Clean up logical NAs if any exist
    keep_na[is.na(keep_na)] <- FALSE

    removed_na <- sum(!keep_na)
    abundance_matrix <- abundance_matrix[keep_na, , drop = FALSE]
  }

  # --- 2) CV FILTERING SCENARIO ---
  if (is.character(fCV) && fCV[1] == 'T') {
    # FIX: Explicitly target index [2] to extract the threshold value safely
    threshold <- as.numeric(fCV[2])

    protected_by_completeness <- has_complete_group(abundance_matrix)

    group_cv_pass <- sapply(group_cols, function(cols) {
      sub_mat <- abundance_matrix[, cols, drop = FALSE]
      g_mean  <- matrixStats::rowMeans2(sub_mat, na.rm = TRUE)
      g_sd    <- matrixStats::rowSds(sub_mat, na.rm = TRUE)
      g_cv    <- g_sd / g_mean

      # Pass if group is clean OR completely empty (handled by NA sweep anyway)
      is.na(g_cv) | (g_cv <= threshold)
    })

    if (is.vector(group_cv_pass)) {
      keep_cv <- group_cv_pass | protected_by_completeness
    } else {
      keep_cv <- (rowSums(group_cv_pass, na.rm = TRUE) > 0) | protected_by_completeness
    }

    keep_cv[is.na(keep_cv)] <- FALSE
    removed_cv <- sum(!keep_cv)
    abundance_matrix <- abundance_matrix[keep_cv, , drop = FALSE]
  }

  # --- 3) RT FILTERING SCENARIO ---
  if (is.character(fRT) && fRT[1] == 'T') {
    if (!requireNamespace("igraph", quietly = TRUE)) {
      stop("Please install the 'igraph' library to run the RT structural network filter.")
    }

    # FIX: Extract indices [2], [3], and [4] safely from the input vector
    digits     <- as.numeric(fRT[2])
    rt_window  <- as.numeric(fRT[3])
    cor_thresh <- as.numeric(fRT[4])

    current_ids <- rownames(abundance_matrix)
    matched_indices <- match(current_ids, x_clean[, 1])

    rt <- round(x_clean[matched_indices, "RT"], digits)
    mz <- x_clean[matched_indices, "m.z"]

    ord         <- order(rt)
    rt          <- rt[ord]
    mz          <- mz[ord]
    current_ids <- current_ids[ord]

    groups <- list()
    start  <- 1
    n      <- length(rt)

    while (start <= n) {
      end <- start
      while (end < n && (rt[end + 1] - rt[start]) <= rt_window) {
        end <- end + 1
      }
      groups[[length(groups) + 1]] <- start:end
      start <- end + 1
    }

    process_group <- function(idx) {
      if (length(idx) <= 1) return(current_ids[idx])

      sub_ids <- current_ids[idx]
      sub_mat <- abundance_matrix[sub_ids, , drop = FALSE]
      sub_mz  <- mz[idx]
      names(sub_mz) <- sub_ids

      cmat <- suppressWarnings(
        cor(t(sub_mat), method = "spearman", use = "pairwise.complete.obs")
      )

      cmat[is.na(cmat)] <- 0
      rownames(cmat) <- sub_ids
      colnames(cmat) <- sub_ids

      cmat[lower.tri(cmat, diag = TRUE)] <- NA
      pairs <- which(cmat > cor_thresh, arr.ind = TRUE)

      if (nrow(pairs) == 0) return(sub_ids)

      edges <- data.frame(
        from = sub_ids[pairs[, 1]],
        to   = sub_ids[pairs[, 2]],
        stringsAsFactors = FALSE
      )

      g <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = sub_ids)
      comp <- igraph::components(g)

      kept_from_group <- sapply(split(names(comp$membership), comp$membership), function(cluster_elements) {
        cluster_mzs <- sub_mz[cluster_elements]
        names(cluster_mzs)[which.max(cluster_mzs)]
      })

      return(unname(kept_from_group))
    }

    if (parallel) {
      final_list <- future.apply::future_lapply(groups, process_group)
    } else {
      final_list <- lapply(groups, process_group)
    }

    final_kept_ids <- unique(unlist(final_list))

    removed_rt <- nrow(abundance_matrix) - length(final_kept_ids)
    abundance_matrix <- abundance_matrix[rownames(abundance_matrix) %in% final_kept_ids, , drop = FALSE]
  }

  # --- 4) PRINT DYNAMIC FILTER SUMMARY TO CONSOLE ---
  final_peaks <- nrow(abundance_matrix)

  cat("\n=== MetaboPeak Filter Summary ===\n")
  cat(paste0("Initial peak count:         ", initial_peaks, "\n"))

  if (!is.logical(fNA) || fNA != FALSE) {
    cat(paste0("Removed by NA filter:       ", removed_na, "\n"))
  }

  if (!is.logical(fCV) || fCV != FALSE) {
    cat(paste0("Removed by CV filter:       ", removed_cv, "\n"))
  }

  if (!is.logical(fRT) || fRT != FALSE) {
    cat(paste0("Removed by RT network loop: ", removed_rt, "\n"))
  }

  cat(paste0("Final peak count:           ", final_peaks, "\n"))
  cat("=================================\n\n")

  output_df <- x_clean[x_clean[, 1] %in% rownames(abundance_matrix), , drop = FALSE]
  return(output_df)
}
