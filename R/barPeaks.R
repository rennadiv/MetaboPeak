#' @title Barplots of 3 lowest variable peaks based on mass (m/z).
#' @description
#' The function calculates CV coefficient (SD/Mean) of each peak with the given mass (m/z)
#' and plots abundances of all samples for the max 3 lowest CV coefficient.
#'
#' @usage barPeaks(x, mass, n)
#' @param x data frame
#' @param mass peak mass (m/z)
#' @param n number of samples
#'
#' @details
#' * The format of the data frame should be: Peak names|abundances of samples|m.z|RT|...
#'
#' * The mass should be put as a string.
#'
#' @returns Bar plots of max 3 peaks with similar m/z and their RT and CV coefficient.
#'
#' @examples
#' barPeaks(neg, '355.10', 48)
#'
#' @export


barPeaks <- function(x, mass, n){
  number.of.samples <- n
  data <- x
  n3 <- number.of.samples+3
  a <- data[substr(data$m.z, 1, nchar(mass)) == mass, c(2:n3)]
  aa <- a[,-c(number.of.samples+1, number.of.samples+2)]
  Mean <- round(rowMeans(aa, na.rm = T),4) # mean of all the rows on 4 decimal places
  SD <- round(apply(aa, 1, sd, na.rm=TRUE),4) # standard deviation for each row
  CV <- round(SD/Mean, 3) # CV coefficient for each row
  if (nrow(a) > 2){
    p.plot <- par(mfrow = c(3,1))
    a <- a[rank(CV) < 4,]
    aa <- aa[rank(CV) < 4,]
    for (i in c(1:3)) {
      barplot(as.matrix(aa[i,]), main = paste('m/z =', round(a[i, number.of.samples+1],3),
                                              ',RT =',round(a[i, number.of.samples+2],2),
                                              ',CV =', CV[rank(CV) < 4][i]), names.arg = rep('',n))
    }
  } else if (nrow(a) == 0) {
    print('This mass is not in the data frame.')
  } else {
    p.plot <- par(mfrow = c(nrow(a),1))
    for (i in c(1:nrow(a))) {
      barplot(as.matrix(aa[i,]), main = paste('m/z =', round(a[i, number.of.samples+1],3),
                                              ',RT =',round(a[i,number.of.samples+2],2),
                                              ',CV =', CV[rank(CV) < 4][i]), names.arg = rep('',n))
    }
  }
}
