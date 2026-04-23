# Morphospace Metrics for GBIF Data as Empirical Assemblages

library(boot)
library(tidyverse)
library(purrr)
library(geometry)

assemblages = read.csv("gbif_empirical.csv") %>% 
  select(-sp_rich)
pca = read.csv("bootstrap_PCA.csv")

# Add community sorting to pca dataframe
pca_assemblages <- assemblages %>%
  left_join(
    pca,
    by = c("genus", "species"),
    relationship = "many-to-many"
    ) %>% 
  filter(!is.na(genus))

# The coding chunk below is for occurrences where we don't have PCA data for the species but for the genus

# Setting seed (optional)
set.seed(123)

# Define PCA component columns
comp_cols <- paste0("Comp", 1:6)

# Fill missing species
pca_assemblages_filled <- pca_assemblages %>%
  rowwise() %>%
  mutate(
    # If species is NA, replace species + all Comp columns
    across(
      c(species, all_of(comp_cols)),
      ~ {
        if (is.na(species)) {
          
          # Get all species from this genus in reference PCA
          candidates <- pca %>% 
            filter(genus == cur_data()$genus)
          
          if (nrow(candidates) > 0) {
            sampled <- candidates %>% slice_sample(n = 1)
            
            if (cur_column() == "species") {
              return(sampled$species)
            } else {
              return(sampled[[cur_column()]])
            }
          }
        }
        .
      }
    )
  ) %>%
  ungroup()

# Fill remaining NA columns
pca_assemblages_comp <- pca_assemblages_filled %>%
  mutate(row_id = row_number())

for (i in seq_len(nrow(pca_assemblages_comp))) {
  
  # Skip rows without NA in Comp columns
  if (!any(is.na(pca_assemblages_comp[i, comp_cols]))) next
  
  g <- pca_assemblages_comp$genus[i]
  
  candidates <- pca %>%
    filter(
      genus == g,
      if_all(all_of(comp_cols), ~ !is.na(.))
    )
  
  if (nrow(candidates) == 0) next
  
  sampled <- candidates %>% slice_sample(n = 1)
  
  for (comp in comp_cols) {
    if (is.na(pca_assemblages_comp[[comp]][i])) {
      pca_assemblages_comp[[comp]][i] <- sampled[[comp]]
    }
  }
}

pca_assemblages_comp <- pca_assemblages_comp %>%
  select(-row_id)

# Ensuring unique species in each cell
pca_clean <- pca_assemblages_comp %>%
  filter(!is.na(Comp1)) %>% 
  select(-Comp7,-Comp8,-Comp9,-Comp10,-Comp11,-Comp12,-Comp13,-Comp14,-Comp15,-Comp16,-Comp17,-Comp18,-Comp19,-Comp20,-Comp21,-Comp22,-Comp23,-Comp24,-Comp25,-Comp26,-Comp27,-Comp28,-Comp29,-Comp30) %>% 
  group_by(cell_id, genus, species) %>%
  slice(1) %>%
  ungroup()

# Calculate species richness
richness <- pca_clean %>%
  group_by(cell_id) %>%
  mutate(sp_rich = n(), .groups = "drop") %>% 
  filter(sp_rich > 1)

assemblage_list <- richness %>%
  group_split() %>% 
  rlang::set_names(map_chr(., ~ as.character(unique(.x$cell_id))))

# Keeping a dataframe of metadata to re-attach to the final table later
meta <- pca_assemblages %>% 
  distinct(assemblage = as.character(cell_id), cell_long = cen_long, cell_lat = cen_lat, country, locality)

# Empty dataframe that will be populated with all empirical data
empirical_all <- tibble()

# Metric Calculations ----

# NND & sdNND ----

nnd <- function(df) {
  coords <- df %>%
    select(starts_with("Comp")) %>%
    select(1:6)
  
  n_species <- nrow(coords)
  
  # Compute pairwise distances
  dist_matrix <- as.matrix(dist(coords))
  
  # Remove self distances
  diag(dist_matrix) <- NA
  
  # If all distances are NA or zero → return NA
  if (all(is.na(dist_matrix))) {
    return(tibble(
      avg_nnd = NA_real_,
      sd_nnd  = NA_real_
    ))
  }
  
  # Compute NND for each species
  nnd_vals <- apply(dist_matrix, 1, function(x) {
    min_val <- min(x, na.rm = TRUE)
    # If min_val is infinite or NaN, set NA
    if (!is.finite(min_val)) return(NA_real_)
    min_val
  })
  
  tibble(
    avg_nnd = mean(nnd_vals, na.rm = TRUE),
    sd_nnd  = sd(nnd_vals, na.rm = TRUE)
  )
}

# Apply function to each replicate (each tibble)
nnd_df <- map_dfr(assemblage_list, nnd, .id = "assemblage")

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
  left_join(hv_df, by = "assemblage") %>% 
  left_join(meta, by = "assemblage")

empirical_species <- pca_clean %>% 
  mutate(assemblage = as.character(cell_id)) %>% 
  select(-country,-locality) %>% 
  left_join(final_df, by = "assemblage", relationship = "many-to-many")

empirical_all <- bind_rows(empirical_species, final_df) %>% 
  select(sp_rich, genus, species, assemblage, cell_long, cell_lat, country, locality, avg_nnd, sd_nnd, hv3d, hv6d, cd, starts_with("Comp")) %>% 
  filter(!is.na(genus))

write.csv(empirical_all, "empirical_gbif.csv")
