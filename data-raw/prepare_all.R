## HEADER --------------------
##
## Script name: prepare_all
##
## Purpose of script: Prepare all the datasets, including MSMS negative and positive mode and treatment dataset.
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

source("data-raw/prepare_neg.R")
source("data-raw/prepare_pos.R")
source("data-raw/prepare_treatment.R")
