## HEADER --------------------
##
## Script name: prepare_pos
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

pos <- read.csv('/home/renca/PRACE/RD_packages/MetaboPeak/code/data-raw/raw/posMSMS.csv', sep = '\t', dec = '.')


# Polishing ---------------------------------------------------------------

# renaming colnames

pos <- pos[,c(1,4:51,3,2)]
x1 <- sapply(colnames(pos[2:49]), function(x) strsplit(x,"_")[[1]][[10]], USE.NAMES = F)
x2 <- sapply(colnames(pos[2:49]), function(x) strsplit(x,"_")[[1]][[11]], USE.NAMES = F)
x <- paste(x1,x2, sep = '')
colnames(pos)[2:49] <- x

colnames(pos)[1] <- "peak_id"

# Saving data set ---------------------------------------------------------

usethis::use_data(
  pos,
  overwrite = TRUE
)
