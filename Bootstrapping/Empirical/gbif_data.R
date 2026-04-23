# GBIF Occurrence Data Collection

library(rgbif)
library(sf)
library(terra)
library(tidyverse)
library(rnaturalearth)
library(dplyr)

file.edit("~/.Renviron")
# Paste the account information below (lines 14-16) directly into the .Renviron file that opens after the line above (should only have to do once)

# Fill in according to GBIF account info
GBIF_USER=username
GBIF_PWD=password
GBIF_EMAIL=email

# Verify login into GBIF (should print the inputted user & email)
Sys.getenv("GBIF_USER")
Sys.getenv("GBIF_EMAIL")

# Get GBIF taxon key to filter for Loricariidae only
tax_key <- name_backbone(
  name = "Loricariidae",
  rank = "family"
)$usageKey

print(tax_key)

# Occurrence data download (with coordinates only)
occ_dl <- occ_download(
 pred("taxonKey", tax_key),
 pred("hasCoordinate", TRUE),
 format = "SIMPLE_CSV"
)

# Check if download is done before proceeding with next line
occ_download_meta(occ_dl)

occ_file <- occ_download_get(occ_dl, overwrite = TRUE)
occ <- occ_download_import(occ_file)

# Convert occurrences to an object in sf
occ_sf <- occ %>%
  # Filter for occurrences w/ coordinates only
  filter(!is.na(decimalLatitude), !is.na(decimalLongitude)) %>% 
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

# Clip map to South America only (time-saving)
sa <- ne_countries(
  continent = "South America",
  returnclass = "sf"
)

occ_sf <- st_join(occ_sf, sa, join = st_within)
occ_sf <- occ_sf[!is.na(occ_sf$name), ]

# Project occurrences onto a metric coordinate system (the UTM most appropriate for the region)
occ_proj <- st_transform(occ_sf, crs = 6933)

# Grid of 1km quadrats
coords <- st_coordinates(occ_proj)

# 1km raster

r <- rast(
  xmin = min(coords[,1]),
  xmax = max(coords[,1]),
  ymin = min(coords[,2]),
  ymax = max(coords[,2]),
  resolution = 1000,
  crs = "EPSG:6933"
)

cell_id <- cellFromXY(r, coords)

# Assign each occurrence to a grid cell
occ_df <- occ_sf %>%
  st_drop_geometry() %>%
  mutate(cell_id = cell_id) %>%
  filter(!is.na(cell_id))

# Making assemblages
# Keep unique species per cell
assemblages_unique <- occ_df %>%
  select(cell_id, binomial = species, country = countryCode, locality) %>%
  filter(!is.na(binomial), binomial != "") %>%
  distinct(cell_id, binomial, .keep_all = TRUE)

# Calculate richness per cell
richness <- assemblages_unique %>%
  group_by(cell_id) %>%
  summarise(sp_rich = n(), .groups = "drop")

# Join richness back
empirical_assemblages <- assemblages_unique %>%
  left_join(richness, by = "cell_id") %>%
  mutate(binomial = as.character(binomial)) %>%
  separate(
    col = binomial,
    into = c("genus", "species"),
    sep = " ",
    extra = "merge",
    fill = "right"
  )

# Getting centroid coordinates for each cell to identify each by latitude and longitude
centroids <- xyFromCell(r, unique(empirical_assemblages$cell_id))

centroids_df <- tibble(
  cell_id = unique(empirical_assemblages$cell_id),
  centroid_x = centroids[,1],
  centroid_y = centroids[,2]
) %>%
  st_as_sf(coords = c("centroid_x", "centroid_y"), crs = 6933) %>%
  st_transform(4326) %>%
  mutate(
    cen_long = st_coordinates(.)[,1],
    cen_lat = st_coordinates(.)[,2]
  ) %>%
  st_drop_geometry()

# Preparing & writing final .csv output
empirical <- empirical_assemblages %>% 
  left_join(centroids_df, by = "cell_id") %>% 
  select(
    sp_rich,
    cell_id,
    cen_long,
    cen_lat,
    genus,
    species,
    country,
    locality
  ) %>% 
  filter(sp_rich > "1") %>% 
  arrange(sp_rich, cell_id)

write_csv(empirical, "gbif_empirical.csv")
