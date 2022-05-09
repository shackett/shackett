library(dplyr)
library(impulse)
library(shiny)

tidy_omics <- readRDS(domics_paths$local_path[domics_paths$type == app_type])
romic::app_heatmap(tidy_omics)