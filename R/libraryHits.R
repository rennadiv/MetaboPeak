#' @title Annotation peaks based on in house library
#'
#' @description
#' This function combines data frame with m/z and RT with given library.
#' Based on the ranges of m/z and RT it tries to match m/z,RT of the library with the detected.
#'
#' @usage libraryHits(x, l, RT.range, mz.range, unique)
#' @param x data frame
#' @param l library (data frame)
#' @param RT.range numeric value determining the range of RT
#' @param mz.range numeric value determining the range of mz
#' @param unique logical, if there should not be replicates with the same compound name
#'
#' @details
#' The format of the data frame (x) should be: Peak names|abundances of samples|m.z|RT|...
#'
#' The format of the library: Compound names|m.z|RT|...
#'
#' @returns
#' list of 3 objects (lib_small, lib_all and unique_name)
#'
#' * lib_small = data frame in the format: RT difference|Compound name|m.z|RT|abundances of just annotated samples
#'
#' * lib_all = data frame in the format: Compound name|x (all the samples, unknown = 'unknown')
#'
#' * unique_names = string vector of all the compound names, that have been annotated
#'
#' @examples
#' data_library <- libraryHits(neg, pkg_lib, RT.range = 0.25, mz.range = 0.009, unique = T)
#'
#' @import dplyr
#'
#' @export


libraryHits <- function(x, l, RT.range = 0.25, mz.range = 0.009, unique = T){
  mat1 <- l
  mat2 <- x

  lib_small <- NULL

  res <- vector("list",nrow(mat1))
  for(i in 1:nrow(mat1)) {
    mz.diff <- mat2$m.z - mat1$m.z[i]
    rt.diff <- mat2$RT - mat1$RT[i]
    j <- which(rt.diff > -RT.range & rt.diff < RT.range &
                 mz.diff > -mz.range & mz.diff < mz.range)
    if (length(j) > 1 & unique == T){
      p <- which.min(abs(rt.diff[j]))
      j <- j[p]
      res[[i]] <- j
    } else if (length(j) == 0) {
      res[[i]] <- NULL
    } else {
      res[[i]] <- j
    }
  }

  n <- rep.int(1:length(res), sapply(res, length))
  m <- unlist(res)
  o <- which(colnames(mat2) == 'm.z' | colnames(mat2) == 'RT')
  lib_small <- cbind(mat1[n,], RTdiff = mat2$RT[m] - mat1$RT[n], Peak_m.z = mat2$m.z[m], Peak_RT = mat2$RT[m], mat2[m,-o])
  colnames(lib_small)[1] <- 'Compound'
  unique_names <- unique(lib_small$Compound)

  df <- lib_small %>%
    add_count(lib_small$Peak_m.z, sort = F, name = 'Peak_frequency')
  k <- c()
  for (i in 1:nrow(x)){
    w <- which(x$m.z[i] == df$`neg_lib$Peak_m.z`)
    if (length(w) == 0){
      k_new <- 'unknown'
    } else if (length(w) > 0 & length(df$Peak_frequency[w]) == 1){
      k_new <- df$Compound[w]
    } else {
      k_new <- df$group[w][1]
    }
    k <- c(k,k_new)
  }
  lib_all <- cbind(Alignment.ID = x[,1],Compound = k,x[,-1])

  output <- list('lib_small' = lib_small, 'lib_all' <- lib_all, 'unique_names' <- unique_names)
  return(output)
}
