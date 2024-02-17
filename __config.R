# -------------------------------------------------------------
# CONFIGURATION
# -------------------------------------------------------------

# working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set current dir as WD

# -------------------------------------------------------------
# 1. Packages Imports
# -------------------------------------------------------------
library("dada2"); packageVersion("dada2")
library("phyloseq");packageVersion("phyloseq")

library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("stringr")
library("vegan")
library("gridExtra")
library("tidyverse")

#source("fonctions.R")


# -------------------------------------------------------------
# 2. Config
# -------------------------------------------------------------

# folders 
data_folder  = file.path(getwd(),  gene)
output_results =  file.path(data_folder, "results_v1") # desired output folder name
rawdata =  file.path(data_folder, "rawdata") # folder containing fastq files


# Folders creation
dir.create(file.path(data_folder), showWarnings = FALSE)
dir.create(file.path(output_results), showWarnings = FALSE)
dir.create(file.path(rawdata), showWarnings = FALSE)


# load tab files
taxa_file = file.path(data_folder, "taxa.rds")
seqtab_nochim_file = file.path(data_folder, "seqtab_nochim.rds")
dbs = list("16S" = "~/silva_nr99_v138.1_wSpecies_train_set.fa.gz")

# list of samples to exclude
exclusion_list = c("H201","H202","H203", "Undetermined")