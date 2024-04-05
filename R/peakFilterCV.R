#' @title Filter by CV coefficient
#' @description
#' This function provides filter by CV coefficient and can be used with or without treatment information.
#'
#'
#' @usage peakFilterCV(x, cut, y)
#' @param x data frame
#' @param cut numeric value of the cutoff
#' @param y treatment info data table (optional)
#'
#' @details
#' The format of the data frame should be: Peak names|abundances of samples|m.z|RT|...
#'
#' Treatment info must contain samples name, that corresponds with x in the 1. column and column treatment
#' (code for combination of treatments, use [treatGroup()])
#'
#'
#' @returns A reduced data frame.
#'
#' @examples
#' peakFilterCV(neg, 0.8)
#' peakFilterCV(neg, 0.8, t_info_group)
#'
#' @import dplyr
#' @import igraph
#' @export
#'

peakFilterCV <- function(x, cut, y){

  x[x == 0] <- NA
  if (is.na(x[1,1])){
    x[1,1] = 0
  }


  ## Function which looks for empty vectors
  is.empty <- function(x) length(x)==0
  ## Function to not use NAs in getting maximum of a row
  f1 <- function(x) (max(x, na.rm = T))

  Mean <- round(rowMeans(z, na.rm = T),4) # mean of all the rows on 4 decimal places
  SD <- round(apply(z, 1, sd, na.rm=TRUE),4) # standard deviation for each row


  if (missing(y)){
    x$CV <- round(SD/Mean, 3) # CV coefficient for each row
    delete <- rownames(x)[which(x$CV > as.numeric(cut))]

    o <- which(colnames(x) == 'CV')
    x.CV <- x[! (rownames(x) %in% delete), -o]
    return(x.NA)
  } else {
    y$treatment <- as.factor(y$treatment)

    Z <- x %>% dplyr::select(y[,1])
    rownames(Z) <- x[,1]
    colnames(Z) <- paste(colnames(Z), y$treatment, sep = '_')

    z <- x
    z$CV <- round(SD/Mean, 3) # CV coefficient for each row
    trash.names <- rownames(z)[which(z$CV > as.numeric(cut))]
    Look.at <- vector('list',length(trash.names))
    names(Look.at) <- trash.names
    leave <- c()

    for (i in levels(y.$treatment)) {
      new.elements <- trash.names[complete.cases(z[trash.names,grepl(i, colnames(z))])] # going through all the treatments and looking if the peak is in all the replicants
      if (!is.empty(new.elements)){ # if there are some, then write it in the list
        for (j in new.elements) {
          Look.at[[j]] <- c(Look.at[[j]],i)
        }
      }
      leave <- c(leave, new.elements) # also write the peak number here (it's making duplicates - need to use unique)
    }

    Look.at_CV <- Look.at[unique(leave)]
    length(unique(leave))
    delete <- setdiff(trash.names, unique(leave)) # vector of clusters to delete

    o <- which(colnames(z) == 'CV')
    if (is.empty(delete)) {
      z.CV <- z[,-o]
    } else {
      z.CV <- z[! (rownames(z) %in% delete), -o]
    }
    colnames(z.CV) <- sapply(colnames(z.CV), function(x) strsplit(x,"_")[[1]][[1]], USE.NAMES = F)
    return(z.CV)
  }
}
