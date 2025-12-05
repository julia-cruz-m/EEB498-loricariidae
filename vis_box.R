# Box Plots
library(tidyverse)
library(cowplot)

synthetic = read.csv("final_all.csv") %>% 
  mutate(type = "Synthetic")
empirical = read.csv("empirical_all.csv") %>% 
  mutate(type = "Empirical")

final = bind_rows(synthetic, empirical) %>%
  filter(sp_rich < 14) # Because currently don't have empirical assemblage richness > 13

# NND
ggplot(final, aes(x = type, y = avg_nnd, fill = type)) + 
  geom_boxplot() +
  geom_signif(comparisons = list(c("Empirical", "Synthetic")), 
              map_signif_level=TRUE) +
  labs(x = "Assemblage Type", y = "Average NND") +
  scale_fill_manual(values = c("lightblue2", "gray")) +
  theme_cowplot() +
  theme(legend.position = "none")

# sdNND
ggplot(final, aes(x = type, y = sd_nnd, fill = type)) + 
  geom_boxplot() +
  geom_signif(comparisons = list(c("Empirical", "Synthetic")), 
              map_signif_level=TRUE) +
  labs(x = "Assemblage Type", y = "Evenness (sdNND)") +
  scale_fill_manual(values = c("lightblue2", "lightgray")) +
  theme_cowplot() +
  theme(legend.position = "none")

# CD
ggplot(final, aes(x = type, y = cd, fill = type)) + 
  geom_boxplot() +
  geom_signif(comparisons = list(c("Empirical", "Synthetic")), 
              map_signif_level=TRUE) +
  labs(x = "Assemblage Type", y = "Average CD") +
  scale_fill_manual(values = c("lightblue2", "lightgray")) +
  theme_cowplot() +
  theme(legend.position = "none")

# Adjusting so that species richness does not confound w/ results
# Currently empirical only has a max richness of 13 whereas synthetic has 20
rich_adj <- final %>% 
  filter(sp_rich < 14)

# HV3D
ggplot(rich_adj, aes(x = type, y = hv3d, fill = type)) + 
  geom_boxplot() +
  geom_signif(comparisons = list(c("Empirical", "Synthetic")), 
              map_signif_level=TRUE) +
  labs(x = "Assemblage Type", y = "3D Hypervolume") +
  scale_fill_manual(values = c("lightblue2", "lightgray")) +
  theme_cowplot() +
  theme(legend.position = "none")

# HV6D
ggplot(rich_adj, aes(x = type, y = hv6d, fill = type)) + 
  geom_boxplot() +
  geom_signif(comparisons = list(c("Empirical", "Synthetic")), 
              map_signif_level=TRUE) +
  labs(x = "Assemblage Type", y = "6D Hypervolume") +
  scale_fill_manual(values = c("lightblue2", "gray")) +
  theme_cowplot() +
  theme(legend.position = "none")
