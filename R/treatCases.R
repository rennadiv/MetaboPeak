#' @title Treatments of specific peak in complete cases.
#' @description
#' Based on m/z and RT this function will return treatments in which the complete cases are.
#'
#' @usage treatCases(x, y, mass, RT, n)
#' @param x data frame
#' @param y data frame - treatment info
#' @param mass peak mass (m/z)
#' @param RT retention time
#' @param n number of samples
#'
#' @details
#' The format of the data frame x should be: Peak names|abundances of samples|m.z|RT|...
#'
#' The data frame y needs column named treatment which contains the code for all treatments together (for this function treatGroup can be used).
#' Also the first column should corespond to the names of samples in data frame x.
#'
#' @returns A string of the treatments code.
#'
#' @examples
#' treatCases(neg,t_info_group,447.12851,5.252,48)
#' treatCases(neg,t_info_group,447.12872,12.386,48)
#' treatCases(neg,t_info_group,447.22214,13.226,48)
#'
#' @export

treatCases <- function(x, y, mass, RT, n){
  x[x == 0] <- NA
  if (is.na(x[1,1])){
    x[1,1] = 0
  }
  peak <- x[substr(x$m.z, 1, nchar(mass)) == mass & substr(x$RT, 1, nchar(RT)) == RT,2:(n+1)]
  if (nrow(peak) > 1){
    print('The combination of m/z and RT is not unique.')
  } else {
    elements <- c()
    z <- peak %>% dplyr::select(y[,1])
    colnames(z) <- paste(colnames(z), y$treatment, sep = '_')
    y$treatment <- as.factor(y$treatment)
    q <- levels(y$treatment)
    for (i in q) {
     if (complete.cases(z[,grepl(i, colnames(z))]) == T) new.elements <- i else next
     elements <- c(elements, new.elements)
     }
  }
  return(elements)
  }
