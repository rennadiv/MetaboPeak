#' @title Unique code for combined treatment
#' @description
#' Based on the given abbreviation for the names of treatment, it makes an unique code for each treatment combination to which the sample is
#' exposed.
#'
#' @usage treatGroup(x, number.of.treatments, my_abbreviations)
#' @param x data frame with treatment info
#' @param number.of.treatments integer of treatment columns (types of treatment)
#' @param my_abbreviations character vector of the abbreviations you want to use in alphabetical order (optional)
#'
#' @details
#' The format of the data frame should be: sample name|treatment groups|...
#'
#' The vector of abbreviations has to be in alphabetical order (small cases first). If missing, then the treatment code will be formed
#' from the names of the treatments combined by underscore.
#'
#'
#' @returns x with new column with the treatment code
#'
#' @examples
#' abb <- c('B','S','A','E','D','W','n','N')
#' treatGroup(t_info, 4, abb)
#'
#' treatGroup(t_info, 4)
#'
#' @export
#'


treatGroup <- function (x, number.of.treatments, my_abbreviations){
  if (missing(my_abbreviations)){
    n <- 1:number.of.treatments+1
    x.new <- apply(x[,n],1,paste,collapse = "_")
    x$treatment <- x.new
    x$treatment <- as.factor(x$treatment)
    return(x)
  } else {
  new_el <- vector('list', number.of.treatments)
  g <- c()
  j = 0
  for (i in 1:number.of.treatments+1){
    x[,i] <- as.factor(x[,i])
    l <- length(levels(x[,i]))
    for (k in 1:length(levels(x[,i]))) {
      levels(x[,i])[k] <- my_abbreviations[j+k]
      new_el[[i-1]] <- as.character(x[,i])
    }
    j <- j+l
    g <- paste(g,new_el[[i-1]], sep = '')
  }
  x$treatment <- as.factor(g)
  return(x)
  }
}

