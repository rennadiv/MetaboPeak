## HEADER --------------------
##
## Script name: prepare_treatment
##
## Purpose of script: Import the treatment dataset and add an unique code for the treatment combination.
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

t_info <- read.csv('/home/renca/PRACE/RD_packages/MetaboPeak/code/data-raw/raw/treatment_info.csv', sep = ',', dec = '.')


# Polishing ---------------------------------------------------------------

abb <- c('B','S','A','E','D','W','n','N')
t_info <- treatGroup(t_info, 4, abb)

# Saving data set ---------------------------------------------------------

usethis::use_data(
  t_info,
  overwrite = TRUE
)

