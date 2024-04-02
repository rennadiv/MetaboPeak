#' @title Name of the sample in which m/z is highest.
#' @description
#' Given the mass (m/z) it looks for the highest value of abundance throughout the samples and returns name of the sample.
#'
#' @usage highestAbund(x, mass, n)
#' @param x data frame
#' @param mass peak mass (m/z)
#' @param n number of samples
#'
#' @details
#' The format of the data frame should be: Peak names|abundances of samples|m.z|RT|...
#'
#' The mass should be put as a string.
#'
#' @returns A string of the sample name.
#'
#' @examples
#' highestAbund(neg,'593.12805',48)
#' highestAbund(pos,'192.05482',48)
#'
#' @export

highestAbund <- function(x, mass, n){
  data <- x
  number.of.samples <- n
  n1 <- number.of.samples+1
  a <- data[substr(data$m.z, 1, nchar(mass)) == mass, c(2:n1)]
  if (nrow(a) > 1) {
    print('The m/z is not unique. Use more decimal places.')
  } else {
    colnames(a)[which(a == max(a, na.rm = T))]
  }
}
