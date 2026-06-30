## HEADER --------------------
##
## Script name: prepare_neg
##
## Purpose of script: Import the original metabolomics datasets, clean and standardize them,
##                    and create the example datasets distributed with the package.
##
## Author: Renata Divinová
##
## Date Created: 2026-06-30
##
## Copyright (c) Renata Divinová, 2026
## Email: divinova.r@czechglobe.cz
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------

rm(list=ls()) # removing all from global environment

# Load packages -----------------------------------------------------------


# Load data ---------------------------------------------------------------

neg <- read.csv('/home/renca/PRACE/RD_packages/MetaboPeak/code/data-raw/raw/negMSMS.csv', sep = '\t', dec = '.')


# Polishing ---------------------------------------------------------------

# renaming colnames

neg <- neg[,c(1,4:51,3,2)]
y1 <- sapply(colnames(neg[2:49]), function(x) strsplit(x,"_")[[1]][[10]], USE.NAMES = F)
y2 <- sapply(colnames(neg[2:49]), function(x) strsplit(x,"_")[[1]][[11]], USE.NAMES = F)
y <- paste(y1,y2, sep = '')
colnames(neg)[2:49] <- y

colnames(neg)[1] <- "peak_id"


# Saving data set ---------------------------------------------------------

usethis::use_data(
  neg,
  overwrite = TRUE
)
