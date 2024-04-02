#' @title Annotation peaks based on in house library
#'
#' @description
#' This function combines data frame with m/z and RT with given library.
#' Based on the ranges of m/z and RT it tries to match m/z,RT of the library with the detected.
#'
#' @usage libraryHits(x, l, RT.range, mz.range, unique)
#' @param x data frame
#' @param l library (data frame)
#' @param RT.range numeric vector with 2 values
#' @param mz.range numeric vector with 2 values
#' @param unique logical, if there should not be replicates with the same compoud name
#'
#' @details
#' The format of the data frame (x) should be: Peak names|abundances of samples|m.z|RT|...
#'
#' The format of the library: Compound names|m.z|RT|...
#'
#' @returns data frame in the fromat: RT difference|Compound name|m.z|RT|abundances of samples
#'
#' @examples
#' data_new <- libraryHits(neg, pkg_lib, RT.range = c(-0.25, 0.25), mz.range = c(-0.009, 0.009), unique = T)
#'
#' @export


libraryHits <- function(x, l, RT.range = c(-0.25, 0.25), mz.range = c(-0.009, 0.009), unique = T){
  mat1 <- l
  mat2 <- x

  out <- NULL

  res <- vector("list",nrow(mat1))
  for(i in 1:nrow(mat1)) {
    mz.diff <- mat2$m.z - mat1$m.z[i]
    rt.diff <- mat2$RT - mat1$RT[i]
    j <- which(rt.diff > RT.range[1] & rt.diff < RT.range[2] &
                 mz.diff > mz.range[1] & mz.diff < mz.range[2])
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
  out <- cbind(mat1[n,], RTdiff = mat2$RT[m] - mat1$RT[n], Peak_m.z = mat2$m.z[m], Peak_RT = mat2$RT[m], mat2[m,-o])
  #out <- cbind(RTdiff=mat2$RT[m] - mat1$RT[n],mat1[n, -c(4:7)], mat2[m, -c(157,158)])
  return(out)
}
