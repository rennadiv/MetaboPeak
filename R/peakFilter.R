#' @title Filter by NA, CV coefficient and/or correlated RT
#' @description
#' This function provides 3 types of filters. It needs the abundances table, treatment info with combined treatment code under the column
#' treatment (function treatGroup can be used to make this code) and some thresholds.
#'
#' 1) based on NA values: it calculates which peaks (rows) have more NAs than
#' threshold and if there are no complete cases for any treatment group, than it removes those rows.
#'
#' 2) based on CV coefficient: it calculates which peaks (rows) have bigger CV coefficient than
#' threshold and if there are no complete cases for any treatment group, than it removes those rows.
#'
#' 3) based on correlated RT: it rounds the RT and orders them, than based on thresholds groups similar retention times and calculates correlation
#' (spearman). From this correlated group it picks the one row with max m/z. (it might take some time)
#'
#' @usage peakFilter(x, y, fNA, fCV, fRT)
#' @param x data frame
#' @param y treatment info data table
#' @param fNA character vector
#' @param fCV character vector
#' @param fRT character vector
#'
#' @details
#' The format of the data frame should be: Peak names|abundances of samples|m.z|RT|...
#'
#' Treatment info must contain samples name, that corresponds with x in the 1. column and column treatment
#' (code for combination of treatments, use [treatGroup()])
#'
#' fNA and fCV can be either defined as F/FALSE or with 2 argument vector ('T', threshold). fRT can be either defined as F/FALSE, or
#' with 4 argument vector ('T', how many decimal places should RT be rounded, range value in which to look for similar RT, correlation
#' threshold - how well should the rows be correlated to take them as the same)
#'
#'
#' @returns A reduced data frame.
#'
#' @examples
#' peakFilter(neg, t_info_group, fNA = c('T',0.5), fCV = F, fRT = F)
#' peakFilter(neg, t_info_group, fNA = c('T',0.5), fCV = c('T', 0.4), fRT = F)
#' peakFilter(neg, t_info_group, fNA = c('T',0.5), fCV = c('T', 0.4), fRT = c('T',3,0.2,0.95)
#'
#'
#' @import dplyr
#' @importFrom igraph graph_from_data_frame
#' @importFrom igraph components
#' @importFrom stats cor
#'
#' @export
#'

peakFilter <- function (x, y, fNA = c('T',0.5), fCV = c('T', 0.8), fRT = c('T',3,0.02,0.95)){
  x[x == 0] <- NA
  if (is.na(x[1,1])){
    x[1,1] = 0
  }
  Z <- x %>% dplyr::select(y[,1])
  rownames(Z) <- x[,1]
  colnames(Z) <- paste(colnames(Z), y$treatment, sep = '_')

  y$treatment <- as.factor(y$treatment)

  ## Function which looks for empty vectors
  is.empty <- function(x) length(x)==0
  ## Function to not use NAs in getting maximum of a row
  f1 <- function(x) (max(x, na.rm = T))

  fce.NA <- function(z, y. = y, fNA. = fNA){
    z$pctNAs <- round(rowSums(is.na(z))/length(z),2)
    trash.names <- rownames(z)[which(z$pctNAs > as.numeric(fNA.[2]))]
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

  fce.CV <- function(z, y.=y, fCV.=fCV){
    Mean <- round(rowMeans(z, na.rm = T),4) # mean of all the rows on 4 decimal places
    SD <- round(apply(z, 1, sd, na.rm=TRUE),4) # standard deviation for each row
    z$CV <- round(SD/Mean, 3) # CV coefficient for each row

    trash.names <- rownames(z)[which(z$CV > as.numeric(fCV.[2]))]
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

  fce.RT <- function(z,x.=x, fRT.=fRT){
    x <- x.[rownames(z),]
    data <- data.frame('ID' = rownames(z), 'RT' = round(x$RT, as.numeric(fRT.[2])), 'mass' = x$m.z) # data table from the given one (df given by user)
    data_sorted <- data[order(data$RT),] # sorted by RT
    data_RTs <- data_sorted$RT # RTs in their own vector
    w <- length(unique(data_RTs)) # length of that vector
    my_list <- vector('list',w) # list where all the similar cluster would be stored
    names(my_list) <- sort(unique(data_RTs)) # naming the elements of the list by unique RT
    RTs <- c()
    keep <- c()
    trash <- c()
    z_filtered <- c()
    for (i in 1:w) {
      # for every RT, find all the RT with the values +-r (r = given by user)
      RTs <- data_sorted[data_RTs >= as.numeric(names(my_list[i]))-as.numeric(fRT.[3]) & data_RTs <= as.numeric(names(my_list[i]))+as.numeric(fRT[3]),]
      my_list[[i]] <- RTs # putting those RT into the list
      RTs <- c() # clearing for next loop

      if (length(my_list[[i]]$ID) == 1) next # if there are no similar RT, then jump to the next iteration
      # calculate correlation on the samples (rows taken are the ones with similar RT)
      var.corelation <- cor(as.matrix(t(z[my_list[[i]]$ID,])), method = 'spearman', use = "pairwise.complete.obs")
      # prevent duplicated pairs
      var.corelation <- var.corelation*lower.tri(var.corelation)
      # take those over the threshold (given by user)
      check.corelation <- which(var.corelation > as.numeric(fRT.[4]), arr.ind=TRUE)

      # get correlated groups:
      graph.cor <- graph_from_data_frame(check.corelation, directed = FALSE)
      groups.cor <- split(unique(as.vector(check.corelation)), components(graph.cor)$membership)
      cor_groups <- lapply(groups.cor,FUN=function(list.cor){rownames(var.corelation)[list.cor]})
      if (length(cor_groups) == 0) next # if there are no correlated peaks above threshold, than next iteration
      for (j in 1:length(cor_groups)) { # looking through all the groups in cor_groups
        z$Max <- apply(z, 1, f1)
        table <- z[match(cor_groups[[j]], rownames(z)),] # table with the correlated peaks and more columns
        keep_new <- rownames(table)[which.max(table$Max)] # the one peak to keep based on Max value (taking the highest from all)
        trash_new <- rownames(table)[-which.max(table$Max)] # the rest of the peaks
        keep <- c(keep,keep_new)
        keep <- unique(keep)
        trash <- c(trash, trash_new)
        trash <- unique(trash)
      }
    }
    o <- which(colnames(z) == 'Max')
    z.RT <- z[! rownames(z) %in% trash, -o] # new data frame, without the removed rows (peaks)
    colnames(z.RT) <- sapply(colnames(z.RT), function(x) strsplit(x,"_")[[1]][[1]], USE.NAMES = F)
    return(z.RT)
  }

  if (fNA[1] == 'T' & fCV[1] == 'T' & fRT[1] == 'T'){
    z.NA <- fce.NA(Z)
    z.CV <- fce.CV(z.NA)
    Z <- fce.RT(z.CV)
  } else if (fNA[1] == 'T' & fCV[1] == 'T' & fRT[1] != 'T'){
    z.NA <- fce.NA(Z)
    Z <- fce.CV(z.NA)
  } else if (fNA[1] == 'T' & fCV[1] != 'T' & fRT[1] == 'T'){
    z.NA <- fce.NA(Z)
    Z <- fce.RT(z.NA)
  } else if (fNA[1] == 'T' & fCV[1] != 'T' & fRT[1] != 'T'){
    Z <- fce.NA(Z)
  } else if (fNA[1] != 'T' & fCV[1] == 'T' & fRT[1] == 'T'){
    z.CV <- fce.CV(Z)
    Z <- fce.RT(z.CV)
  } else if (fNA[1] != 'T' & fCV[1] == 'T' & fRT[1] != 'T'){
    Z <- fce.CV(Z)
  } else if (fNA[1] != 'T' & fCV[1] != 'T' & fRT[1] == 'T'){
    Z <- fce.RT(Z)
  } else {
    Z <- x %>% dplyr::select(y[,1])
    rownames(Z) <- x[,1]
    colnames(Z) <- paste(colnames(Z), y$treatment, sep = '_')
  }
  Z <- cbind(rownames(Z),Z)
  Z$m.z <- x[rownames(Z),'m.z']
  Z$RT <- x[rownames(Z),'RT']
  return(Z)
}

