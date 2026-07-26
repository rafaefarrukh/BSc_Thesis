###################################################################################################
# GLOBAL
###################################################################################################

# LIBRARIES #####################################

# shiny dashboard
library(shiny)
library(shinyWidgets)
library(bslib)

# data manipulation
library(tidyverse)

# visualization
library(plotly)
library(DT)
library(r3dmol)
library(visNetwork)

# IMPORT DATA ###################################

# data objects
annotated_genome <- read_rds("data/annotated_genome.rds") # data_prep/1_genome.R
gem_network <- read_rds("data/gem_network.rds") # data_prep/3_gem.R

# all step 7 csv results
load("data/step_7_results.RData") # data_prep/csv_to_rdata.R

# ui files
lapply(list.files(path = "ui", pattern = "\\.R$", full.names = TRUE), source)

# server files
lapply(list.files(path = "server", pattern = "\\.R$", full.names = TRUE), source)
###################################################################################################