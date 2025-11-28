# R Script to clean data to get ready for bootstrapping

library(tidyverse)

PCA = read.csv("Loricariidae_PCA_species.csv") %>% 
  select(id = X, starts_with("Comp")) # Changing a column name

genus_species = read.csv("lujan_ct.csv") %>% 
  filter(Family == "Loricariidae") %>% # Only loricariids
  # Changing column names and selecting only necessary columns
  select(id = Catalog.Number.for.Common.Use, binomial = Updated.Binomial.for.Common.Use) %>% 
  # Replace spaces with underscores to match PCA data
  mutate(id = gsub(" ", "_", trimws(id))) %>%
  # Filtering this dataset to only include ids that are present in the PCA dataframe
  filter(id %in% PCA$id)

# Join dataframes together
merged_df <- genus_species %>%
  inner_join(PCA, by = "id") %>% 
  # Separate binomial names into genus and species columns, respectively
  separate(binomial, into = c("genus", "species"), sep = " ", extra = "merge", fill = "right")

write_csv(merged_df, "bootstrap_PCA.csv")
