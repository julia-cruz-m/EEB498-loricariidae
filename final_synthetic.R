# Final synthetic data

library(boot)
library(tidyverse)
library(purrr)
library(geometry)

pca = read.csv("bootstrap_PCA.csv")

# Make an empty dataframe that will be populated with all the synthetic data
NND_all <- tibble()

# For loop to apply each function for each species richness (N) from 2-20
for (N in 2:20) {
  # Sampling ----
  random_species <- function(data, N) {
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
  
  # 1000 replicates with N genera/species each
  N_rep <- 1000
  boot_list <- replicate(
    N_rep, 
    random_species(pca, N), 
    simplify = FALSE)
  
  # NND & sdNND ----
  nnd <- function(df) {
    # Select coordinate columns (e.g. Comp1–Comp6 -> 6D)
    coords <- df %>%
      select(starts_with("Comp")) %>%
      select(1:6)
    # Compute pairwise Euclidean distance as a matrix
    dist_matrix <- as.matrix(dist(coords))
    # Find nearest neighbour distance for each species (smallest non-zero)
    nnd <- apply(dist_matrix, 1, function(x) min(x[x > 0]))
    # Return average NND & sdNND for that tibble
    tibble(
      avg_nnd = mean(nnd),
      sd_nnd = sd(nnd)
    )
  }
  # Apply function to each replicate (each tibble)
  nnd_df <- map_dfr(boot_list,
                    nnd, .id = "replicate"
  )
  
  # Centroid ----
  n_centroid <- map_dfr(seq_along(boot_list), function(i) {
    centroid <- boot_list[[i]]
    centroid %>%
      summarise(
        replicate = i,
        comp1_cen = mean(Comp1, na.rm = TRUE),
        comp2_cen = mean(Comp2, na.rm = TRUE),
        comp3_cen = mean(Comp3, na.rm = TRUE),
        comp4_cen = mean(Comp4, na.rm = TRUE),
        comp5_cen = mean(Comp5, na.rm = TRUE),
        comp6_cen = mean(Comp6, na.rm = TRUE)
      )
  })
  
  # Average CD ----
  
  # Add centroid as a column to boot_list
  boot_list <- map2(
    boot_list,
    split(n_centroid, n_centroid$replicate),
    ~ mutate(.x,
             comp1_cen = .y$comp1_cen,
             comp2_cen = .y$comp2_cen,
             comp3_cen = .y$comp3_cen,
             comp4_cen = .y$comp4_cen,
             comp5_cen = .y$comp5_cen,
             comp6_cen = .y$comp6_cen)
  )
  
  avg_cd <- function(df) {
    # Select coordinate columns (e.g. Comp1–Comp6 -> 6D)
    coords <- df %>% 
      select(starts_with("Comp")) %>% 
      select(1:6)
    centroids <- df %>% 
      slice(1) %>% 
      select(starts_with("comp")) %>% 
      select(ends_with("_cen"))
    # Calculate Euclidean distance from each species to centroid
    dists <- sqrt(rowSums((as.matrix(coords) - as.numeric(centroids))^2))
    # Average distance across species
    mean(dists)
  }
  
  # Apply function to each tibble
  n_cd <- tibble(
    replicate = seq_along(boot_list),
    cd = map_dbl(boot_list, avg_cd)
  )
  
  # Hypervolume ----
  # Create empty dataframe, which will be populated later
  hv_df <- tibble(
    replicate = seq_along(boot_list),
    hv3d = NA_real_,
    hv6d = NA_real_
  )
  
  for (i in seq_along(boot_list)) {
    df_i <- boot_list[[i]]
    
    # 3D hypervolume (N > 3)
    if (N >= 4 && nrow(df_i) >= 4) {
      dimensions_3 <- df_i %>% 
        select(Comp1, Comp2, Comp3) %>% 
        as.matrix()
      
      # Calculation
      hv3 <- tryCatch({
        convhulln(dimensions_3, "FA")$vol
      }, error = function(e) NA_real_) # Return NA if there is an error
      
      # Add into dataframe
      hv_df$hv3d[i] <- hv3
    }
    
    # 6D hypervolume (N > 6)
    if (N >= 7 && nrow(df_i) >= 7) {
      dimensions_6 <- df_i %>%
        select(Comp1, Comp2, Comp3, Comp4, Comp5, Comp6) %>%
        as.matrix()
      
      hv6 <- tryCatch({
        convhulln(dimensions_6, "FA")$vol
      }, error = function(e) NA_real_)
      
      # Add into dataframe
      hv_df$hv6d[i] <- hv6
    }
  }
  # Cleaning ----
  # Adding species richness as a column
  nnd_df <- nnd_df %>%
    mutate(sp_rich = N) %>% 
    mutate(replicate = as.integer(replicate))
  
  # Adding centroid, CD, and hypervolume to original NND dataframe
  final_df <- nnd_df %>%
    left_join(n_centroid, by = "replicate") %>%
    left_join(n_cd, by = "replicate") %>% 
    left_join(hv_df, by = "replicate")
  
  NND_all <- bind_rows(NND_all, final_df)
}

final_all <- NND_all %>% 
  select(sp_rich, replicate, avg_nnd, sd_nnd, hv3d, hv6d, cd, starts_with("comp"))

write.csv(final_all, "final_all.csv")
