#' @title Correlation of two peaks in the same retention time.
#' @description
#' This function gives spearman correlation coefficient of two peaks, defined by their mass (m/z).
#' They can be from the same data frame, or two different data frames. The combining factor is their retention time (RT).
#'
#' @usage checkCorr(x, y = NULL, masses, RT, n)
#' @param x data frame
#' @param y NULL (default) or a data frame with compatible dimensions to x. The default is equivalent to y = x
#' @param masses vector of peak masses
#' @param RT retention time
#' @param n number of samples
#'
#' @details
#' The format of the data frame (x,y) should be: Peak names|abundances of samples|m.z|RT|...
#'
#' The masses can be either written as numeric (rounded) or as character (not rounded).
#'
#' The correlation is computed by spearman coefficiant and if there are missing values than the correlation between each pair of variables is computed using all complete pairs of observations on those variables.
#'
#' @returns Correlation of the two peaks.
#'
#' @examples
#' checkCorr(pos, masses = c(595.2, 449.08), RT = 11.8, n= 48)
#' checkCorr(pos, neg, masses = c(595.2, 593.15), RT = 11.8, n= 48)
#'
#' @importFrom stats cor
#' @export
#'

checkCorr <- function(x, y = NULL, masses, RT, n){
  n1 <- n+1
  data <- x
  if (!is.null(y)){
    data2 <- y
  } else {
    data2 <- x
  }
  mass1 <- masses[1]
  mass2 <- masses[2]
  decimalplaces <- function(x) {
    if ((x %% 1) != 0) {
      nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
    } else {
      return(0)
    }
  }
  f <- decimalplaces(as.numeric(mass1))
  g <- decimalplaces(as.numeric(mass2))

  df1 <- data[which(round(data$m.z,f) == mass1), c(2:n1,n1+1,n1+2)]
  df1 <- df1[which.min(df1$RT - RT), c(1:n)]

  df2 <- data2[which(round(data2$m.z,g) == mass2), c(2:n1,n1+1, n1+2)]
  df2 <- df2[which.min(df2$RT - RT), c(1:n)]

  if (nrow(df1) == 0) {
    print('First mass is not in the data table')
  } else if (nrow(df2) == 0) {
    print('Second mass is not in the data table')
  } else {
    c1 <- as.vector(t(df1))
    c2 <- as.vector(t(df2))
    if (sum(is.na(c1)) > 0 | sum(is.na(c2)) > 0){
      cor(c1,c2,method = 'spearman', use = 'pairwise.complete.obs')
    } else {
      cor(c1, c2, method = 'spearman')
    }
  }
}
