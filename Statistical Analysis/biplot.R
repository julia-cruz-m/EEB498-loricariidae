# Biplot creation

library(plyr)
library(dplyr)
library(tidyverse)
library(ggfortify)
library(cowplot)

assemblages = read_csv("pca_clean.csv") %>% 
  select(-.groups, assemblage = cell_id)
metrics = read_csv("empirical_gbif.csv") %>% 
  select(assemblage, avg_nnd, sd_nnd, hv3d, hv6d, cd)

joined = left_join(assemblages, metrics, by = "assemblage") %>% 
  select(-id)

# Convex hull function
# find_hull <- function(df) df[chull(df$Comp1, df$Comp2), ]

# Add filter template (highest species richness)
filtered <- assemblages %>% 
  filter(assemblage == "8991625") %>% 
  mutate(group = "filter")

hulls_filtered <- filtered %>%
  group_by(group) %>%
  slice(chull(Comp1, Comp2)) %>%
  ungroup()

ggplot(data = filtered, aes(x = Comp1, y = Comp2, colour = group, fill = group)) +
  geom_point() + 
  geom_polygon(data = hulls_filtered, alpha = 0.5) +
  labs(x = "PC1", y = "PC2") +
  theme_cowplot() +
  # Colour change
  scale_colour_manual(values = "lightblue3") +
  scale_fill_manual(values = "lightblue")

# Filter: lowest HV3D
filtered1 <- assemblages %>% 
  filter(assemblage == "19099035") %>% 
  mutate(group = "filter")

hulls_filtered1 <- filtered1 %>%
  group_by(group) %>%
  slice(chull(Comp1, Comp2)) %>%
  ungroup()

ggplot(data = filtered1, aes(x = Comp1, y = Comp2, colour = group, fill = group)) +
  geom_point() + 
  geom_polygon(data = hulls_filtered1, alpha = 0.5) +
  labs(x = "PC1", y = "PC2") +
  theme_cowplot() +
  # Colour change
  scale_colour_manual(values = "lightblue3") +
  scale_fill_manual(values = "lightblue")

# Filter: highest HV3D
filtered2 <- assemblages %>% 
  filter(assemblage == "4449715") %>% 
  mutate(group = "filter")

hulls_filtered2 <- filtered2 %>%
  group_by(group) %>%
  slice(chull(Comp1, Comp2)) %>%
  ungroup()

ggplot(data = filtered2, aes(x = Comp1, y = Comp2, colour = group, fill = group)) +
  geom_point() + 
  geom_polygon(data = hulls_filtered2, alpha = 0.5) +
  labs(x = "PC1", y = "PC2") +
  theme_cowplot() +
  # Colour change
  scale_colour_manual(values = "lightblue3") +
  scale_fill_manual(values = "lightblue")

# Highest HV6D
filtered1 <- assemblages %>% 
  filter(assemblage == "8991625") %>% 
  mutate(group = "filter")

hulls_filtered1 <- filtered2 %>%
  group_by(group) %>%
  slice(chull(Comp1, Comp2)) %>%
  ungroup()

ggplot(data = filtered1, aes(x = Comp1, y = Comp2, colour = group, fill = group)) +
  geom_point() + 
  geom_polygon(data = hulls_filtered1, alpha = 0.5) +
  labs(x = "PC1", y = "PC2") +
  theme_cowplot() +
  # Colour change
  scale_colour_manual(values = "lightblue3") +
  scale_fill_manual(values = "lightblue")

# Lowest HV6D
filtered2 <- assemblages %>% 
  filter(assemblage == "12627820") %>% 
  mutate(group = "filter")

hulls_filtered2 <- filtered2 %>%
  group_by(group) %>%
  slice(chull(Comp1, Comp2)) %>%
  ungroup()

ggplot(data = filtered2, aes(x = Comp1, y = Comp2, colour = group, fill = group)) +
  geom_point() + 
  geom_polygon(data = hulls_filtered2, alpha = 0.5) +
  labs(x = "PC1", y = "PC2") +
  theme_cowplot() +
  # Colour change
  scale_colour_manual(values = "lightblue3") +
  scale_fill_manual(values = "lightblue") +
