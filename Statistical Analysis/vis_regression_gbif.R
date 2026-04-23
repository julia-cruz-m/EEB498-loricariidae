# Linear Regression Visualization & Tests

library(tidyverse)
library(cowplot)
library(scales)
library(broom)
library(ggeffects)

synthetic = read.csv("final_all_1.csv") %>% 
  # Adding assemblage type column (Random/Empirical)
  mutate(type = "Random")

empirical = read.csv("empirical_gbif.csv") %>% 
  # Adding assemblage type column
  mutate(type = "Empirical", sd_nnd = ifelse(sd_nnd == 0, NA, sd_nnd)) %>% 
  # 1 row per assemblage
  distinct(assemblage, .keep_all = TRUE)

# Code below counts how many Empirical assemblages there are
(n_distinct(empirical$genus))

# Combining synthetic + empirical assemblage data into one dataframe
final = bind_rows(synthetic, empirical) %>% 
  filter(sp_rich > 1, sp_rich <= 25)

# Assemblage type column as a factor
final$type <- factor(final$type)

# NND *** ----

# Making sure nnd is not 0 or infinite
nnd_nonzero = final %>% 
  filter(avg_nnd > 0) %>% 
  filter(is.finite(avg_nnd))

ggplot(data = nnd_nonzero, aes(x = log(sp_rich), y = avg_nnd, color = type)) +
  geom_smooth(
    method = "glm",
    method.args = list(family = Gamma(link = "log")),
    se = TRUE) +
  labs(y = "Average NND", 
       x = "log(Species Richness)", 
       colour = "Assemblage Type",
       title = "C") +
  theme_cowplot() +
  scale_colour_manual(breaks = c("Empirical", "Random"), values = c("blue", "gray50")) +
  coord_cartesian(
    xlim = range(log(nnd_nonzero$sp_rich)),
    clip = "off"
  ) +
  theme(plot.margin = margin(5.5, 40, 5.5, 5.5)) +
  
  # Intercept difference caption
  annotate("text", x = Inf, y = Inf,
           label = "Intercept~difference:~italic(p) < 0.001", hjust = 1, vjust = 1, size = 3.5, parse = TRUE) +
  # Slope difference caption
  annotate("text", x = Inf, y = 0.26,
           label = "Slope~difference:~italic(p) < 0.001",hjust = 1, vjust = 1, size = 3.5, parse = TRUE) +
  # Random assemblage label
  annotate("text", x = 3, y = 0.12,
           label = "Random", colour = "gray50") +
  # Empirical assemblage label
  annotate("text", x = 2.6, y = 0.06,
           label = "Empirical", colour = "blue") +
  theme(legend.position = "none")


# Summary of slope difference b/w assemblage types
nnd <- glm(avg_nnd ~ log(sp_rich) * type,
                 data = nnd_nonzero,
                 family = Gamma(link = "log"))

summary(nnd) # ***
summary(nnd)$coefficients # For table

# Predicted values
nnd_pred <- ggpredict(nnd, terms = c("sp_rich [all]", "type"))

nnd_effect <- as.data.frame(nnd_pred) %>%
  dplyr::select(x, predicted, group) %>%
  tidyr::pivot_wider(names_from = group, values_from = predicted) %>%
  dplyr::mutate(percent_diff = (Empirical / Random - 1) * 100)

# Effect size summary
nnd_summary <- nnd_effect %>%
  summarise(
    min_diff = min(percent_diff),
    max_diff = max(percent_diff),
    mean_diff = mean(percent_diff)
  )

nnd_summary

# Min & max % differences
nnd_effect[c(1, nrow(nnd_effect)), ]

# Max divergence
nnd_max <- nnd_effect %>%
  mutate(abs_diff = abs(percent_diff)) %>%
  slice_max(abs_diff, n = 1)

nnd_max

# Min divergence
nnd_min <- nnd_effect %>% 
  mutate(abs_diff = abs(percent_diff)) %>% 
  slice_min(abs_diff, n = 1)

nnd_min

# % change in slope
nnd_slope <- nnd_effect %>%
  summarise(
    emp_change = (last(Empirical) - first(Empirical)) / first(Empirical) * 100,
    rand_change = (last(Random) - first(Random)) / first(Random) * 100
  )

nnd_slope

# sdNND *** ----
sdnnd_nonzero = final %>% 
  filter(sd_nnd > 0) %>% 
  filter(sp_rich >= 3)

ggplot(data = sdnnd_nonzero, aes(x = log(sp_rich), y = sd_nnd, color = type)) +
  geom_smooth(method = glm, 
              method.args = list(family = Gamma(link = "log")), 
              se = T) +
  labs(y = "sdNND", 
       x = "log(Species Richness)", 
       colour = "Assemblage Type",
       title = "D") +
  theme_cowplot() +
  scale_colour_manual(breaks = c("Random", "Empirical"), values = c("gray50", "blue")) +
  coord_cartesian(
    xlim = range(log(sdnnd_nonzero$sp_rich)),
    clip = "off"
  ) +
  theme(plot.margin = margin(5.5, 40, 5.5, 5.5)) +
  
  # Intercept difference caption
  annotate("text", x = Inf, y = Inf,
           label = "Intercept~difference:~italic(p) < 0.001", hjust = 1, vjust = 1, size = 3.5, parse = TRUE) +
  # Slope difference caption
  annotate("text", x = Inf, y = 0.095,
           label = "Slope~difference:~italic(p) < 0.001",hjust = 1, vjust = 1, size = 3.5, parse = TRUE) +
  # Random label
  annotate("text", x = 3, y = 0.047,
           label = "Random", colour = "gray50") +
  # Empirical label
  annotate("text", x = 3, y = 0.081,
           label = "Empirical", colour = "blue") +
  theme(legend.position = "none")

sdnnd <- glm(sd_nnd ~ log(sp_rich) * type, data = sdnnd_nonzero, 
             family = Gamma(link = "log"))
summary(sdnnd) # ***
summary(sdnnd)$coefficients

# Predicted values
sdnnd_pred <- ggpredict(sdnnd, terms = c("sp_rich [all]", "type"))

sdnnd_effect <- as.data.frame(sdnnd_pred) %>%
  dplyr::select(x, predicted, group) %>%
  tidyr::pivot_wider(names_from = group, values_from = predicted) %>%
  dplyr::mutate(percent_diff = (Empirical / Random - 1) * 100)

# Effect size summary
sdnnd_summary <- sdnnd_effect %>%
  summarise(
    min_diff = min(percent_diff),
    max_diff = max(percent_diff),
    mean_diff = mean(percent_diff)
  )

sdnnd_summary

# Min & max % differences
sdnnd_effect[c(1, nrow(sdnnd_effect)), ]

# Max divergence
sdnnd_max <- sdnnd_effect %>%
  mutate(abs_diff = abs(percent_diff)) %>%
  slice_max(abs_diff, n = 1)

sdnnd_max

# Slope differences
sdnnd_slope <- sdnnd_effect %>%
  summarise(
    emp_change = (last(Empirical) - first(Empirical)) / first(Empirical) * 100,
    rand_change = (last(Random) - first(Random)) / first(Random) * 100
  )

sdnnd_slope

# CD *** ----
ggplot(data = final, aes(x = log(sp_rich), y = cd, color = type)) +
  geom_smooth(
    method = "glm",
    method.args = list(family = Gamma(link = "log")),
    se = TRUE) +
  labs(y = "Average CD", 
       x = "log(Species Richness)", 
       colour = "Assemblage Type",
       title = "E") +
  theme_cowplot() +
  scale_colour_manual(breaks = c("Empirical", "Random"), values = c("blue", "gray50"))  +
  coord_cartesian(
    xlim = range(log(final$sp_rich)),
    clip = "off"
  ) +
  theme(plot.margin = margin(5.5, 40, 5.5, 5.5)) +
  annotate("text", x = Inf, y = Inf,
           label = "Intercept~difference:~italic(p) < 0.001", hjust = 1, vjust = 1, size = 3.5, parse = TRUE) +
  annotate("text", x = Inf, y = 0.258,
           label = "Slope~difference:~italic(p)==0.417",hjust = 1, vjust = 1, size = 3.5, parse = TRUE) +
  annotate("text", x = 3, y = 0.218,
           label = "Random", colour = "gray50") +
  annotate("text", x = 3, y = 0.2467,
           label = "Empirical", colour = "blue") +
  theme(legend.position = "none")

cd <- glm(cd ~ log(sp_rich) * type,
              data = final,
              family = Gamma(link = "log"))

summary(cd) # non-sig p = 0.417 # ***
summary(cd)$coefficients

# Predicted values
cd_pred <- ggpredict(cd, terms = c("sp_rich [all]", "type"))

cd_effect <- as.data.frame(cd_pred) %>%
  dplyr::select(x, predicted, group) %>%
  tidyr::pivot_wider(names_from = group, values_from = predicted) %>%
  dplyr::mutate(percent_diff = (Empirical / Random - 1) * 100)

# Effect size summary
cd_summary <- cd_effect %>%
  summarise(
    min_diff = min(percent_diff),
    max_diff = max(percent_diff),
    mean_diff = mean(percent_diff)
  )

cd_summary

# Min & max % differences
cd_effect[c(1, nrow(cd_effect)), ]

# Max divergence
cd_max <- cd_effect %>%
  mutate(abs_diff = abs(percent_diff)) %>%
  slice_max(abs_diff, n = 1)

cd_max

# Slope differences
cd_slope <- cd_effect %>%
  summarise(
    emp_change = (last(Empirical) - first(Empirical)) / first(Empirical) * 100,
    rand_change = (last(Random) - first(Random)) / first(Random) * 100
  )

cd_slope

# HV3D *** ----
ggplot(data = final %>% 
         filter(sp_rich > 3), aes(x = log(sp_rich), y = hv3d, color = type)) +
  geom_smooth(
    data = final %>% filter(log(sp_rich) <= 2.9 & (log(sp_rich) >= 1.3)),
    method = "glm",
    method.args = list(family = Gamma(link = "log")),
    se = T) +
  labs(y = "3D Hypervolume", 
       x = "log(Species Richness)", 
       colour = "Assemblage Type",
       title = "A") +
  theme_cowplot() +
  scale_colour_manual(breaks = c("Random", "Empirical"), values = c("gray50", "blue")) +
  coord_cartesian(
    xlim = c(1.4, 2.85),
    clip = "off") +
  theme(plot.margin = margin(5.5, 40, 5.5, 5.5)) +
  annotate("text", x = Inf, y = Inf,
           label = "Intercept~difference:~italic(p) < 0.001", hjust = 1, vjust = 1, size = 3.5, parse = TRUE) +
  annotate("text", x = Inf, y = 0.018,
           label = "Slope~difference:~italic(p) < 0.001",hjust = 1, vjust = 1, size = 3.5, parse = TRUE) +
  annotate("text", x = 2.7, y = 0.016,
           label = "Random", colour = "gray50") +
  annotate("text", x = 2.8, y = 0.007,
           label = "Empirical", colour = "blue") +
  theme(legend.position = "none")

hv3d <- glm(hv3d ~ log(sp_rich) * type,
                data = final,
                family = Gamma(link = "log")) # ***

summary(hv3d) # ***
summary(hv3d)$coefficients

# Predicted values
hv3d_pred <- ggpredict(hv3d, terms = c("sp_rich [4:25]", "type"))

hv3d_effect <- as.data.frame(hv3d_pred) %>%
  dplyr::select(x, predicted, group) %>%
  tidyr::pivot_wider(names_from = group, values_from = predicted) %>%
  dplyr::mutate(percent_diff = (Empirical / Random - 1) * 100)

# Effect size summary
hv3d_summary <- hv3d_effect %>%
  summarise(
    min_diff = min(percent_diff),
    max_diff = max(percent_diff),
    mean_diff = mean(percent_diff)
  )

hv3d_summary

# Min & max % differences
hv3d_effect[c(1, nrow(hv3d_effect)), ]

# Max divergence
hv3d_max <- hv3d_effect %>%
  mutate(abs_diff = abs(percent_diff)) %>%
  slice_max(abs_diff, n = 1)

hv3d_max

# Slope differences
hv3d_slope <- hv3d_effect %>%
  summarise(
    emp_change = (last(Empirical) - first(Empirical)) / first(Empirical) * 100,
    rand_change = (last(Random) - first(Random)) / first(Random) * 100
  )

hv3d_slope

# HV6D *** (code to change linetype is here) ----
ggplot(data = final %>% filter(sp_rich > 6), 
       aes(x = log(sp_rich), y = hv6d, color = type)) +
  geom_smooth(
    data = final %>% filter((log(sp_rich) <= 2.9 & log(sp_rich) >= 1.9)),
    method = "glm",
    method.args = list(family = Gamma(link = "log")),
    se = TRUE
  ) +
  labs(y = "6D Hypervolume", 
       x = "log(Species Richness)", 
       colour = "Assemblage Type",
       title = "B") +
  theme_cowplot() +
  scale_colour_manual(breaks = c("Random", "Empirical"), values = c("gray50", "blue")) +
  coord_cartesian(xlim = c(1.95, 2.89), clip = "off") +  # limit x
  theme(plot.margin = margin(5.5, 40, 5.5, 5.5)) +   # give space on right
  annotate("text", x = Inf, y = Inf,
           label = "Intercept~difference:~italic(p) < 0.001", hjust = 1, vjust = 1, size = 3.5, parse = TRUE) +
  annotate("text", x = Inf, y = 2.4e-06,
           label = "Slope~difference:~italic(p) < 0.001", hjust = 1, vjust = 1, size = 3.5, parse = TRUE) +
  annotate("text", x = 2.75, y = 1.9e-06,
           label = "Random", colour = "gray50") +
  annotate("text", x = 2.875, y = 5e-07,
           label = "Empirical", colour = "blue") +
  theme(legend.position = "none")

hv6d <- glm(hv6d ~ log(sp_rich) * type,
            data = final,
            family = Gamma(link = "log")) # ***

summary(hv6d) # ***
summary(hv6d)$coefficients

# Predicted values
hv6d_pred <- ggpredict(hv6d, terms = c("sp_rich [7:25]", "type"))

hv6d_effect <- as.data.frame(hv6d_pred) %>%
  dplyr::select(x, predicted, group) %>%
  tidyr::pivot_wider(names_from = group, values_from = predicted) %>%
  dplyr::mutate(percent_diff = (Empirical / Random - 1) * 100)

# Effect size summary
hv6d_summary <- hv6d_effect %>%
  summarise(
    min_diff = min(percent_diff),
    max_diff = max(percent_diff),
    mean_diff = mean(percent_diff)
  )

hv6d_summary

# Min & max % differences
hv6d_effect[c(1, nrow(hv6d_effect)), ]

# Max divergence
hv6d_max <- hv6d_effect %>%
  mutate(abs_diff = abs(percent_diff)) %>%
  slice_max(abs_diff, n = 1)

hv6d_max

# Slope differences
hv6d_slope <- hv6d_effect %>%
  summarise(
    emp_change = (last(Empirical) - first(Empirical)) / first(Empirical) * 100,
    rand_change = (last(Random) - first(Random)) / first(Random) * 100
  )

hv6d_slope
