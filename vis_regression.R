# Linear Regression Visualization & Tests

library(tidyverse)
library(cowplot)

synthetic = read.csv("final_all.csv") %>% 
  mutate(type = "Synthetic")
empirical = read.csv("empirical_all.csv") %>% 
  mutate(type = "Empirical")

final = bind_rows(synthetic, empirical)

# FINAL

# NND
ggplot(data = final, aes(x = sp_rich, y = avg_nnd, color = type)) +
  geom_point(data = final %>% filter(type == "Empirical")) +
  geom_smooth(method = glm, method.args = list(family = 'poisson'), se = F) +
  labs(y = "Average NND", x = "Species Richness", colour = "Assemblage Type",
       caption = "p < 0.05") +
  theme_cowplot() +
  scale_colour_manual(breaks = c("Synthetic", "Empirical"),
                      values = c("gray", "blue")) +
  theme(plot.caption = element_text(face = "italic", hjust = 1))

# Summary of slope difference b/w assemblage types
nnd <- lm(avg_nnd ~ sp_rich * type, data = final)
summary(nnd)

# Assumes slopes are the same; tests whether lines are above/below each other
nnd_int <- lm(avg_nnd ~ sp_rich + type, data = final)
summary(nnd_int)

# sdNND
ggplot(data = final, aes(x = sp_rich, y = sd_nnd, color = type)) +
  geom_point(data = final %>% filter(type == "Empirical")) +
  geom_smooth(method = glm, method.args = list(family = 'poisson'), se = F) +
  labs(y = "Evenness (sdNND)", x = "Species Richness", colour = "Assemblage Type",
       caption = "p = 0.384") +
  theme_cowplot() +
  scale_colour_manual(breaks = c("Synthetic", "Empirical"),
                      values = c("gray", "blue"))  +
  theme(plot.caption = element_text(face = "italic", hjust = 1))

sdnnd <- lm(sd_nnd ~ sp_rich * type, data = final)
summary(sdnnd)

sdnnd_int <- lm(sd_nnd ~ sp_rich + type, data = final)
summary(sdnnd_int)

# CD
ggplot(data = final, aes(x = sp_rich, y = cd, color = type)) +
  geom_point(data = final %>% filter(type == "Empirical")) +
  geom_smooth(method = glm, method.args = list(family = 'poisson'), se = F) +
  labs(y = "Average CD", x = "Species Richness", colour = "Assemblage Type",
       caption = "p = 0.217") +
  theme_cowplot() +
  scale_colour_manual(breaks = c("Synthetic", "Empirical"),
                      values = c("gray", "blue")) +
  theme(plot.caption = element_text(face = "italic", hjust = 1))

cd <- lm(cd ~ sp_rich * type, data = final)
summary(cd)

cd_int <- lm(cd ~ sp_rich + type, data = final)
summary(cd_int)

# HV3D
ggplot(data = final, aes(x = sp_rich, y = hv3d, color = type)) +
  geom_point(data = final %>% filter(type == "Empirical")) +
  geom_smooth(method = glm, method.args = list(family = 'poisson'), se = F) +
  labs(y = "3D Hypervolume", x = "Species Richness", colour = "Assemblage Type",
       caption = "p = 0.058") +
  theme_cowplot() +
  scale_colour_manual(breaks = c("Synthetic", "Empirical"),
                      values = c("gray", "blue")) +
  theme(plot.caption = element_text(face = "italic", hjust = 1))

hv3d <- lm(hv3d ~ sp_rich * type, data = final)
summary(hv3d)

hv3d_int <- lm(hv3d ~ sp_rich + type, data = final)
summary(hv3d_int)

# HV6D
ggplot(data = final, aes(x = sp_rich, y = hv6d, color = type)) +
  geom_point(data = final %>% filter(type == "Empirical")) +
  geom_smooth(method = glm, method.args = list(family = 'poisson'), se = F) +
  labs(y = "6D Hypervolume", x = "Species Richness", colour = "Assemblage Type",
       caption = "p = 0.973") +
  theme_cowplot() +
  scale_colour_manual(breaks = c("Synthetic", "Empirical"),
                      values = c("gray", "blue")) +
  theme(plot.caption = element_text(face = "italic", hjust = 1))

hv6d <- lm(hv6d ~ sp_rich * type, data = final)
summary(hv6d)

hv6d_int <- lm(hv6d ~ sp_rich + type, data = final)
summary(hv6d_int)

# Code to move legend to top right instead of right
#  + theme(legend.position = c(0.95, 0.95),
#    legend.justification = c("right", "top"))

# SYNTHETIC VISUALIZATION

# NND
ggplot(data = final, aes(x = sp_rich, y = avg_nnd, color = type)) +
#  geom_point(colour = "lightblue") +
  geom_smooth(method = glm, method.args = list(family = 'poisson'), se = F) +
  xlab("Species Richness") +
  ylab("Average NND") +
  theme_cowplot()

nnd <- glm(avg_nnd ~ sp_rich, 
             data = synthetic, 
             family = poisson)

summary(nnd)

# sdNND
hist(synthetic$sd_nnd)

# CD
ggplot(data = synthetic, aes(x = sp_rich, y = cd)) +
  geom_point(colour = "lightblue") +
  geom_smooth(method = glm, method.args = list(family = 'poisson'), colour = "navy", se = F) +
  xlab("Species Richness") +
  ylab("Average CD") +
  theme_cowplot()

# Hypervolume
# 3D
ggplot(data = synthetic, aes(x = sp_rich, y = hv3d)) +
  geom_point(colour = "lightblue") +
  geom_smooth(method = glm, method.args = list(family = 'poisson'), colour = "navy", se = F) +
  xlab("Species Richness") +
  ylab("3D Hypervolume") +
  theme_cowplot()

# 6D
ggplot(data = synthetic, aes(x = sp_rich, y = hv6d)) +
  geom_point(colour = "lightblue") +
  geom_smooth(method = glm, method.args = list(family = 'poisson'), colour = "navy", se = F) +
  xlab("Species Richness") +
  ylab("6D Hypervolume") +
  theme_cowplot()