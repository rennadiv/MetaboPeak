#' @title Filter by NA, CV coefficient and/or correlated RT
#'
#' @description
#' This function provides 3 types of filters for metabolomics peak data. It utilizes an
#' abundances table, sample-to-treatment mapping data, and customizable parameter thresholds.
#'
#' 1) Based on NA values: Removes peaks exceeding the missingness threshold, unless
#' the peak has complete measurements (no NAs) in at least one treatment group.
#'
#' 2) Based on CV coefficient: Removes peaks with a Coefficient of Variation exceeding
#' the threshold, unless the peak has complete measurements in at least one treatment group.
#'
#' 3) Based on correlated RT: Groups peaks with similar retention times (RT). Within each
#' group, it constructs a network graph of pairs exceeding a Spearman correlation threshold.
#' From each connected network component, it retains only the single peak with the maximum m/z value.
#'
#' @usage peakFilter(x, y, fNA, fCV, fRT, parallel, n_cores)
#'
#' @param x Data frame. Format must be: Peak names/IDs in column 1, followed by sample abundance columns, "m.z", and "RT".
#' @param y Treatment mapping data frame. Column 1 must match the sample column names in \code{x}. Must contain a column named \code{Treatment}.
#' @param fNA Vector. Either \code{FALSE} or a 2-element vector: \code{c('T', threshold_fraction)}. Default is \code{c('T', 0.5)}.
#' @param fCV Vector. Either \code{FALSE} or a 2-element vector: \code{c('T', threshold_cv)}. Default is \code{c('T', 0.8)}.
#' @param fRT Vector. Either \code{FALSE} or a 4-element vector: \code{c('T', rounding_digits, rt_window, correlation_threshold)}. Default is \code{c('T', 3, 0.02, 0.95)}.
#' @param parallel Logical. Set to \code{TRUE} to check multiple retention time windows concurrently across multiple CPU cores. Default is \code{FALSE}.
#' @param n_cores Numeric. Number of CPU cores to utilize. If \code{NULL}, defaults to available cores minus 1.
#'
#' @returns A reduced data frame containing all original columns of \code{x}, filtered by the selected criteria.
#'
#' @examples
#' # Example data setup assumed (neg dataset and t_info table)
#' # peakFilter(neg, t_info, fNA = c('T', 0.5), fCV = FALSE, fRT = FALSE)
#' # peakFilter(neg, t_info, fNA = c('T', 0.5), fCV = c('T', 0.4), fRT = FALSE)
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

    # Save user's current parallel configuration
    old_plan <- future::plan()
    # Explicitly revert to normal sequential R on exit to release background memory/CPU resources
    on.exit(future::plan(old_plan), add = TRUE)

    if (is.null(n_cores)) {
      n_cores <- max(1, future::availableCores() - 1)
    }
    future::plan(future::multisession, workers = n_cores)
  }

  # --- INITIAL PROCESSING & CLEANING ---
  # Work on a copy to prevent overriding global variables
  x_clean <- x

  # Ensure index values of 0 are treated as missing values (NA)
  sample_names <- as.character(y[, 1])
  abundance_matrix <- as.matrix(x_clean[, sample_names, drop = FALSE])
  abundance_matrix[abundance_matrix == 0] <- NA

  # Ensure row identities are cleanly mapped to character coordinates
  rownames(abundance_matrix) <- as.character(x_clean[, 1])

  # Cache treatment groupings
  y$Treatment <- as.factor(y$Treatment)
  group_cols <- split(seq_along(y$Treatment), y$Treatment)

  # Helper function to detect rows with complete cases within any single treatment group
  has_complete_group <- function(mat) {
    keep_matrix <- sapply(group_cols, function(cols) {
      rowSums(!is.na(mat[, cols, drop = FALSE])) == length(cols)
    })
    if (is.vector(keep_matrix)) return(keep_matrix) # handling single row outputs safely
    rowSums(keep_matrix) > 0
  }

  # --- 1) NA FILTERING SCENARIO ---
  if (is.character(fNA) && fNA[1] == 'T') {
    threshold <- as.numeric(fNA[2])
    pctNAs <- rowMeans(is.na(abundance_matrix))

    # Keep the row if it passes the threshold OR contains a completely full treatment category
    keep_na <- (pctNAs <= threshold) | has_complete_group(abundance_matrix)
    abundance_matrix <- abundance_matrix[keep_na, , drop = FALSE]
  }

  # --- 2) CV FILTERING SCENARIO ---
  if (is.character(fCV) && fCV[1] == 'T') {
    threshold <- as.numeric(fCV[2])

    # Calculate group-specific CVs to protect true biomarkers
    # without letting high-noise global peaks slip through
    group_cv_pass <- sapply(group_cols, function(cols) {
      sub_mat <- abundance_matrix[, cols, drop = FALSE]
      g_mean  <- matrixStats::rowMeans2(sub_mat, na.rm = TRUE)
      g_sd    <- matrixStats::rowSds(sub_mat, na.rm = TRUE)
      g_cv    <- g_sd / g_mean

      # A group passes if its internal variance is stable (under threshold)
      # OR if it's mostly empty (handled by NA filter anyway, so we don't punish it here)
      is.na(g_cv) | (g_cv <= threshold)
    })

    # Handle single-row edge case safety mapping
    if (is.vector(group_cv_pass)) {
      keep_cv <- any(group_cv_pass)
    } else {
      keep_cv <- rowSums(group_cv_pass) > 0
    }

    abundance_matrix <- abundance_matrix[keep_cv, , drop = FALSE]
  }

  # --- 3) RT FILTERING SCENARIO ---
  if (is.character(fRT) && fRT[1] == 'T') {
    if (!requireNamespace("igraph", quietly = TRUE)) {
      stop("Please install the 'igraph' library to run the RT structural network filter.")
    }

    digits     <- as.numeric(fRT[2])
    rt_window  <- as.numeric(fRT[3])
    cor_thresh <- as.numeric(fRT[4])

    # Subset matching data coordinates smoothly
    current_ids <- rownames(abundance_matrix)
    matched_indices <- match(current_ids, x_clean[, 1])

    rt <- round(x_clean[matched_indices, "RT"], digits)
    mz <- x_clean[matched_indices, "m.z"]

    # Order rows to track sequential RT steps smoothly
    ord  <- order(rt)
    rt   <- rt[ord]
    mz   <- mz[ord]
    abundance_matrix <- abundance_matrix[ord, , drop = FALSE]
    current_ids      <- current_ids[ord]

    # Generate sliding window grouping arrays
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

    # Isolated processing loop for finding network redundancies
    process_group <- function(idx) {
      if (length(idx) <= 1) return(current_ids[idx])

      sub_mat <- abundance_matrix[idx, , drop = FALSE]
      sub_ids <- current_ids[idx]
      sub_mz  <- mz[idx]

      # Compute pairwise Spearman metrics
      cmat <- suppressWarnings(
        cor(t(sub_mat), method = "spearman", use = "pairwise.complete.obs")
      )

      rownames(cmat) <- sub_ids
      colnames(cmat) <- sub_ids

      # Isolate strong correlation matches
      cmat[lower.tri(cmat, diag = TRUE)] <- NA
      pairs <- which(cmat > cor_thresh, arr.ind = TRUE)

      if (nrow(pairs) == 0) return(sub_ids)

      # Build dynamic network tracking maps using igraph
      edges <- data.frame(
        from = sub_ids[pairs[, 1]],
        to   = sub_ids[pairs[, 2]],
        stringsAsFactors = FALSE
      )

      g <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = sub_ids)
      comp <- igraph::components(g)

      # For each interconnected cluster, extract the node containing the max m/z value
      peak_map <- data.frame(id = sub_ids, cluster = comp$membership, mz = sub_mz, stringsAsFactors = FALSE)

      retained_ids <- peak_map %>%
        dplyr::group_by(cluster) %>%
        dplyr::slice_max(order_by = mz, n = 1, with_ties = FALSE) %>%
        dplyr::pull(id)

      return(retained_ids)
    }

    # Execution pathways
    if (parallel) {
      kept_ids_list <- future.apply::future_lapply(groups, process_group)
    } else {
      kept_ids_list <- lapply(groups, process_group)
    }

    final_kept_ids <- unlist(kept_ids_list)
    abundance_matrix <- abundance_matrix[rownames(abundance_matrix) %in% final_kept_ids, , drop = FALSE]
  }

  # --- CLEAN RETAINED ROWSETS AND RETURN ORIGINAL SCHEMA ---
  final_row_names <- rownames(abundance_matrix)
  output_data <- x[x[, 1] %in% final_row_names, , drop = FALSE]

  # Ensure the output row order precisely aligns with our filtered matrix state
  output_data <- output_data[match(final_row_names, output_data[, 1]), , drop = FALSE]

  return(output_data)
}
