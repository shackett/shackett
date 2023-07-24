library(dplyr)
library(shiny)

tidy_omics <- readRDS("heatshock_gasch2K.Rds")
romic::app_heatmap(tidy_omics)