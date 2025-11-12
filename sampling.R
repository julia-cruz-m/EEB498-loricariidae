# Load libraries
library(boot)
library(tidyverse)

# Read in CSV
pca = read.csv('bootstrap_PCA.csv')

# N genera (species richness)
N <- 3

# Function to sample unique genera then sample a species from each genus
random_species <- function(data, indices) {
  # Sample N unique genera
  selected_genera <- sample(unique(data$genus), N, replace = FALSE)
  # Filter data to only those genera
  sampled <- data %>%
    filter(genus %in% selected_genera) %>%
    group_by(genus) %>%
    # Sample one species per genus
    slice_sample(n = 1) %>%
    ungroup()
  return(sampled)
}

# 10 replicates with 3 genera/species each
N_rep <- 10 # N replicates
boot_list <- replicate(
  N_rep, 
  random_species(pca), 
  simplify = FALSE)
