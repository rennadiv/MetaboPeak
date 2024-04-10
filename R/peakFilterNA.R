#' @title Filter by missing values
#' @description
#' This function provides filter by missing values and can be used with or without treatment information.
#'
#' @usage peakFilterNA(x, cut, y)
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
#' peakFilterNA(neg, 0.5)
#' peakFilterNA(neg, 0.5, t_info_group)
#'
#' @import dplyr
#' @export
#'

peakFilterNA <- function(x, cut, y){

  x[x == 0] <- NA
  if (is.na(x[1,1])){
    x[1,1] = 0
  }


  ## Function which looks for empty vectors
  is.empty <- function(x) length(x)==0
  ## Function to not use NAs in getting maximum of a row
  f1 <- NULL
  y. <- NULL

  f1 <- function(x) (max(x, na.rm = T))

  if (missing(y)){
    x$pctNAs <- round(rowSums(is.na(x))/length(x),2)
    delete <- rownames(x)[which(x$pctNAs > as.numeric(cut))]

    o <- which(colnames(x) == 'pctNAs')
    x.NA <- x[! (rownames(x) %in% delete), -o]
    return(x.NA)
  } else {
    y$treatment <- as.factor(y$treatment)

    Z <- x %>% dplyr::select(y[,1])
    rownames(Z) <- x[,1]
    colnames(Z) <- paste(colnames(Z), y$treatment, sep = '_')

    z <- x
    z$pctNAs <- round(rowSums(is.na(z))/length(z),2)
    trash.names <- rownames(z)[which(z$pctNAs > as.numeric(cut))]
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

    Look.at_NA <- Look.at[unique(leave)]
    length(unique(leave))
    delete <- setdiff(trash.names, unique(leave)) # vector of clusters to delete

    o <- which(colnames(z) == 'pctNAs')
    if (is.empty(delete)) {
      z.NA <- z[,-o]
    } else {
      z.NA <- z[! (rownames(z) %in% delete), -o]
    }
    colnames(z.NA) <- sapply(colnames(z.NA), function(x) strsplit(x,"_")[[1]][[1]], USE.NAMES = F)
    return(z.NA)
  }
}
