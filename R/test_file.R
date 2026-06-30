## HEADER --------------------
##
## Script name: Script for testing
##
## Purpose of script: To check all the functions and use it on the package datasets
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

# SET WORKING DIRECTORY -----------------------------
cat("SETTING WORKING DIRECTORY...\n\n", sep = "")
wd <- "/home/renca/PRACE/RD_packages/MetaboPeak/code/"
setwd(wd)
cat("WORKING DIRECTORY HAS BEEN SET TO: ", wd, sep = "")

# INSTALL PACKAGES & LOAD LIBRARIES -----------------
cat("INSTALLING PACKAGES & LOADING LIBRARIES... \n\n", sep = "")
packages <- c("tidyverse", "stringr", "readxl") # list of packages to load
n_packages <- length(packages) # count how many packages are required

new.pkg <- packages[!(packages %in% installed.packages())] # determine which packages aren't installed

# install missing packages
if(length(new.pkg)){
  install.packages(new.pkg)
}

# load all requried libraries
for(n in 1:n_packages){
  cat("Loading Library #", n, " of ", n_packages, "... Currently Loading: ", packages[n], "\n", sep = "")
  lib_load <- paste("library(\"",packages[n],"\")", sep = "") # create string of text for loading each library
  eval(parse(text = lib_load)) # evaluate the string to load the library
}
