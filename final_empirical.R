# Final empirical

library(boot)
library(tidyverse)
library(purrr)
library(geometry)

pca = read.csv("bootstrap_PCA.csv")
assemblages = read.csv("assemblages.csv") %>%  # Table with genus and assemblage columns
  select(genus, assemblage) %>% 
  filter(assemblage != "")

# Add community sorting to pca dataframe
pca <- assemblages %>% 
  left_join(pca, by = "genus", relationship = "many-to-many")

assemblage_list <- pca %>%
  # Convert dataframe into a list of lists by assemblage
  group_split(assemblage) %>% 
  set_names(unique(pca$assemblage)) %>% 
  # Sample one individual per genus if there are multiple
  map(~ .x %>%  
        group_by(genus) %>% 
        slice_sample(n = 1) %>% 
        ungroup()
  )

# Empty dataframe that will be populated with all empirical data
empirical_all <- tibble()

# Metric Calculations ----

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
nnd_df <- map_dfr(assemblage_list,
                  nnd, .id = "assemblage"
)

# Centroid ----
n_centroid <- map_dfr(names(assemblage_list), function(nm) {
  centroid <- assemblage_list[[nm]]
  centroid %>%
    summarise(
      assemblage = nm,
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
assemblage_list <- map2(
  assemblage_list,
  split(n_centroid, n_centroid$assemblage),
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
  assemblage = names(assemblage_list),
  cd = map_dbl(assemblage_list, avg_cd)
)

# Hypervolume ----
# Create empty dataframe, which will be populated later
hv_df <- tibble(
  assemblage = names(assemblage_list),
  hv3d = NA_real_,
  hv6d = NA_real_
)

for (i in seq_along(assemblage_list)) {
  df_i <- assemblage_list[[i]]
  
  # 3D hypervolume (N > 3)
  if (nrow(df_i) >= 4) {
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
  if (nrow(df_i) >= 7) {
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
  mutate(
    sp_rich = sapply(assemblage_list, nrow))

# Adding centroid, CD, and hypervolume to original NND dataframe
final_df <- nnd_df %>%
  left_join(n_centroid, by = "assemblage") %>%
  left_join(n_cd, by = "assemblage") %>% 
  left_join(hv_df, by = "assemblage")

empirical_all <- bind_rows(empirical_all, final_df) %>% 
  select(sp_rich, assemblage, avg_nnd, sd_nnd, hv3d, hv6d, cd, starts_with("comp"))

write.csv(empirical_all, "empirical_all.csv")
