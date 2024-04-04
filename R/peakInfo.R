#' @title Peak information and statistics
#' @description
#' This function returns information about individual peaks within the whole data table (RT, m/z, % of NAs, Max intensity, Mean, SD and CV coefficient).
#' If ID is provided, than it returns the same information, but for just this one peak.
#'
#'
#' @usage peakInfo(x, n, ID)
#' @param x data frame
#' @param n number of samples
#' @param ID numeric (optional)
#'
#'
#' @details
#' * The format of the data frame x should be: Peak names|abundances of samples|m.z|RT|...
#'
#'
#' @returns A data frame with all the peak information for all or for one
#'
#' @examples
#' z <- peakInfo(neg,48)
#' z[1:5,]
#'
#' peakInfo(neg,48,20)
#'
#' @export
#'

peakInfo <- function(x, n, ID){
  x[x == 0] <- NA
  if (is.na(x[1,1])){
    x[1,1] = 0
  }
  rownames(x) <- x[,1]
  m.z <- x[,'m.z']
  RT <- x[,'RT']
  z <- x[,2:(n+1)]
  pctNAs <- round(rowSums(is.na(z))/length(z),2) # percentage of all NA in every row
  f1 <- function(x) (max(x, na.rm = T)) # function to not use NAs in getting maximum of a row
  Max <- apply(z, 1, f1) # taking maximum of a row
  Mean <- round(rowMeans(z, na.rm = T),4) # mean of all the rows on 4 decimal places
  SD <- round(apply(z, 1, sd, na.rm=TRUE),4) # standard deviation for each row
  CV <- round(SD/Mean, 3) # CV coefficient for each row
  x.new <- as.data.frame(cbind(ID = x[,1], m.z, RT, pctNAs, Max, Mean, SD, CV))
  if (missing(ID)){
    return(x.new)
  } else {
    peak <- x.new[which(rownames(x.new) == ID),]
    return(peak)
  }
}
