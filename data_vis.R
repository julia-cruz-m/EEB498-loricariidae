# Visualization

library(tidyverse)
library(cowplot)

synthetic = read.csv("final_all.csv")

# NND
ggplot(data = synthetic, aes(x = sp_rich, y = avg_nnd)) +
  geom_point(colour = "lightblue") +
  geom_smooth(method = glm, method.args = list(family = 'poisson'), colour = "navy", se = F) +
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