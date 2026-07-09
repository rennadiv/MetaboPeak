#' @title Annotate peaks based on an in-house reference library
#'
#' @description
#' Matches experimental mass-to-charge ratios (m/z) and retention times (RT) against a reference
#' standard library using specified tolerance windows.
#'
#' @usage libraryHits(x, l, RT.range = 0.25, mz.range = 0.009, unique = TRUE)
#'
#' @param x Data frame. Experimental data. Format must be: Peak identifiers in column 1, followed by sample abundance columns, "m.z", and "RT".
#' @param l Data frame. Reference library. Format must be: Compound names/identifiers in column 1, "m.z", and "RT".
#' @param RT.range Numeric. Maximum allowable absolute retention time difference (in minutes) for a match. Default is \code{0.25}.
#' @param mz.range Numeric. Maximum allowable absolute mass-to-charge difference (in Da) for a match. Default is \code{0.009}.
#' @param unique Logical. If \code{TRUE}, when multiple experimental peaks match a single library compound, only the one closest in Retention Time is kept. Default is \code{TRUE}.
#'
#' @details
#' Elements that do not match any compound within the library parameters are cleanly categorized as \code{"unknown"}.
#'
#' @returns A named list containing 3 structured elements:
#' \itemize{
#'   \item \code{lib_small}: Data frame mapping successful hits, displaying time deviations, compound identities, library standards, and sample abundances.
#'   \item \code{lib_all}: Data frame mirroring the entire input matrix structure \code{x}, with an added \code{Compound} identity column.
#'   \item \code{unique_names}: Character vector listing all uniquely annotated compound assignments.
#' }
#'
#' @import dplyr
#' @export
libraryHits <- function(x, l, RT.range = 0.25, mz.range = 0.009, unique = TRUE) {

  # ---- 1. SANITIZE DATA STRUCTURING ----
  mat1 <- as.data.frame(l)
  mat2 <- as.data.frame(x)

  # Check if columns exist before using them
  if (!"m.z" %in% colnames(mat2) || !"RT" %in% colnames(mat2)) {
    stop("Experimental data frame 'x' must contain exact columns named 'm.z' and 'RT'.")
  }
  if (!"m.z" %in% colnames(mat1) || !"RT" %in% colnames(mat1)) {
    stop("Library data frame 'l' must contain exact columns named 'm.z' and 'RT'.")
  }

  # FORCE true numeric vectors to completely fix the binary operator math error
  mat2$m.z <- as.numeric(as.character(mat2$m.z))
  mat2$RT  <- as.numeric(as.character(mat2$RT))
  mat1$m.z <- as.numeric(as.character(mat1$m.z))
  mat1$RT  <- as.numeric(as.character(mat1$RT))

  # Standardize the first column name of the library for consistency
  colnames(mat1)[1] <- "Compound"
  mat1$Compound <- as.character(mat1$Compound)

  # ---- 2. REFERENCE PAIRWISE SEARCH LOOP ----
  res <- vector("list", nrow(mat1))

  for (i in seq_len(nrow(mat1))) {
    # Now guaranteed to be purely numeric math subtraction!
    mz.diff <- mat2$m.z - mat1$m.z[i]
    rt.diff <- mat2$RT - mat1$RT[i]

    # Isolate vector index hits matching tolerance conditions
    j <- which(abs(rt.diff) < RT.range & abs(mz.diff) < mz.range)

    # Handle multiple candidate hits if unique mode is enabled
    if (length(j) > 1 && unique) {
      best_hit <- j[which.min(abs(rt.diff[j]))]
      res[[i]] <- best_hit
    } else if (length(j) == 0) {
      res[[i]] <- integer(0)
    } else {
      res[[i]] <- j
    }
  }

  # ---- 3. BUILD THE 'LIB_SMALL' TARGET SUMMARY DATAFRAME ----
  lengths <- sapply(res, length)

  if (sum(lengths) == 0) {
    message("No reference library hits were discovered using the current tolerance thresholds.")

    # Rebuild a clean fallback lib_all structure if there are zero hits
    lib_all_empty <- mat2
    lib_all_empty <- cbind(
      Alignment.ID = mat2[, 1],
      Compound = "unknown",
      mat2[, -1, drop = FALSE],
      stringsAsFactors = FALSE
    )
    colnames(lib_all_empty)[1] <- colnames(mat2)[1]

    return(list(
      lib_small = data.frame(),
      lib_all = lib_all_empty,
      unique_names = character(0)
    ))
  }

  n <- rep.int(seq_len(length(res)), lengths)
  m <- unlist(res)

  # Isolate metadata column exclusions smoothly
  meta_cols <- which(colnames(mat2) %in% c("m.z", "RT"))

  lib_small <- cbind(
    mat1[n, , drop = FALSE],
    RTdiff   = mat2$RT[m] - mat1$RT[n],
    Peak_m.z = mat2$m.z[m],
    Peak_RT  = mat2$RT[m],
    mat2[m, -meta_cols, drop = FALSE],
    stringsAsFactors = FALSE
  )

  unique_names <- unique(lib_small$Compound)

  # ---- 4. FAST VECTORIZED DE-DUPLICATION AND JOINING FOR 'LIB_ALL' ----
  # Calculate peak appearance frequencies cleanly using dplyr syntax rules
  df_counts <- lib_small %>%
    dplyr::add_count(Peak_m.z, name = "Peak_frequency") %>%
    dplyr::mutate(Matching_Key_m.z = Peak_m.z) # Anchor clean linking column

  # Build a lightweight, distinct key-value mapping frame
  mapping_key <- df_counts %>%
    dplyr::group_by(Matching_Key_m.z) %>%
    dplyr::summarise(
      Compound_Assigned = dplyr::first(Compound),
      .groups = "drop"
    )

  # Join the mapping key directly back to the original experimental frame
  mat2_match <- mat2 %>%
    dplyr::left_join(mapping_key, by = c("m.z" = "Matching_Key_m.z")) %>%
    dplyr::mutate(Compound = ifelse(is.na(Compound_Assigned), "unknown", Compound_Assigned))

  # Re-order variables to match your exact output layout criteria
  lib_all <- mat2_match %>%
    dplyr::select(Alignment.ID = 1, Compound, dplyr::everything(), -Compound_Assigned)

  # ---- 5. RETURN COMBINED SUMMARY PAYLOAD ----
  output <- list(
    lib_small    = lib_small,
    lib_all      = lib_all,
    unique_names = unique_names
  )

  return(output)
}
