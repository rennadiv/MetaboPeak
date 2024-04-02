#' @title Unique code for combined treatment
#' @description
#' Based on the given abbreviation for the names of treatment, it makes an unique code for each treatment combination to which the sample is
#' exposed.
#'
#' @usage treatGroup(x, treatment_names, my_abbreviations, number.of.treatments)
#' @param x data frame with treatment info
#' @param treatment_names character vector of the treatment names from the x
#' @param my_abbreviations character vector of the abbreviations you want to use
#' @param number.of.treatments integer of treatment columns (types of treatment)
#'
#' @details
#' * The format of the data frame should be: sample name|treatment groups|...
#'
#'
#' @returns x with new column with the treatment code
#'
#' @examples
#' tr <- c('S','B','AC','EC','D','W','N-','N+')
#' abb <- c('S','B','A','E','D','W','n','N')
#' treatGroup(t_info, tr, abb, 4)
#'
#' @export
#'


treatGroup <- function (x, treatment_names, my_abbreviations, number.of.treatments){
  new_el <- vector('list', number.of.treatments)
  g <- c()
  j = 0
  for (i in 1:number.of.treatments+1){
    x[,i] <- as.factor(x[,i])
    for (k in 1:(length(levels(x[,i]))-1)) {
      new_el[[i-1]] <- ifelse(x[,i] == tr[j+k], abb[j+k], abb[j+k+1])
      j <- j+2
    }
    g <- paste(g,new_el[[i-1]], sep = '')
  }
  x$treatment <- as.factor(g)
  return(x)

}




