# Required libraries
library(rprojroot)
library(readr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggthemes)
theme_set(theme_clean()+theme(legend.position = "top"))
library(StatOrdPattHxC)

data("LinfLsup")

# Define the base path of the project using rprojroot
base_path <- find_root(has_file("README.md"))
plot_path <- file.path(base_path, "Plots")
data_path <- file.path(base_path, "Data", "csv")

# create Plots path if not exists
if(!dir.exists(plot_path)) {
  dir.create(plot_path)
}

# create Data/csv if not exists. Please copy all the csv files here in this directory
if(!dir.exists(data_path)) {
  dir.create(data_path, recursive = TRUE)
}

print("done")