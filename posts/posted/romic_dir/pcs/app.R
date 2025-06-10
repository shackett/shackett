library(shiny)
library(romic)

tidy_omics <- readRDS("gasch2K.Rds")
romic::app_pcs(tidy_omics)
