#' @title Filter by similar RT and their correlation
#' @description
#' This function provides filter by RT and their correlation.
#'
#' @usage peakFilterNA(x, dec, RTrange, cut)
#' @param x data frame
#' @param dec numeric value, on how many decimal places the RT should be rounded
#' @param RTrange numeric value, which defines the range (+-) of similar RT
#' @param cut numeric value of the correlation cutoff
#'
#' @details
#' The format of the data frame should be: Peak names|abundances of samples|m.z|RT|...
#'
#'
#'
#' @returns A reduced data frame.
#'
#' @examples
#' peakFilterRT(neg,3,0.2,0.95)
#'
#'
#' @import dplyr
#' @import igraph
#' @export
#'

peakFilterRT <- function(x, dec, RTrange, cut){
  data <- data.frame('ID' = rownames(x), 'RT' = round(x$RT, as.numeric(dec)), 'mass' = x$m.z) # data table from the given one (df given by user)
  data_sorted <- data[order(data$RT),] # sorted by RT
  data_RTs <- data_sorted$RT # RTs in their own vector
  w <- length(unique(data_RTs)) # length of that vector
  my_list <- vector('list',w) # list where all the similar cluster would be stored
  names(my_list) <- sort(unique(data_RTs)) # naming the elements of the list by unique RT
  RTs <- c()
  keep <- c()
  trash <- c()
  z_filtered <- c()
  z <- x
  for (i in 1:w) {
    # for every RT, find all the RT with the values +-r (r = given by user)
    RTs <- data_sorted[data_RTs >= as.numeric(names(my_list[i]))-as.numeric(RTrange) & data_RTs <= as.numeric(names(my_list[i]))+as.numeric(RTrange),]
    my_list[[i]] <- RTs # putting those RT into the list
    RTs <- c() # clearing for next loop

    if (length(my_list[[i]]$ID) == 1) next # if there are no similar RT, then jump to the next iteration
    # calculate correlation on the samples (rows taken are the ones with similar RT)
    var.corelation <- cor(as.matrix(t(z[my_list[[i]]$ID,])), method = 'spearman', use = "pairwise.complete.obs")
    # prevent duplicated pairs
    var.corelation <- var.corelation*lower.tri(var.corelation)
    # take those over the threshold (given by user)
    check.corelation <- which(var.corelation > as.numeric(cut), arr.ind=TRUE)

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
