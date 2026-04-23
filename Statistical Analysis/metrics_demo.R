# Metrics Demonstration Plot

library(plyr)
library(dplyr)
library(tidyverse)
library(ggfortify)
library(cowplot)
library(ggimage)
library(ggarrow)

species = read_csv("bootstrap_PCA.csv")

# Filter for target species only
filtered <- species %>% 
  filter((genus == "Hemipsilichthys" & species == "nimius") | 
           (genus == "Curculionichthys" & species == "pirba") |
           (genus == "Gymnotocinclus" & species == "anosteos") |
           (genus == "Pareiorhaphis" & species == "garapiply") |
           (genus == "Ancistrus" & species == "ranunculus") |
           (genus == "Avalithoxus" & species == "jantjae") |
           (genus == "Cordylancistrus" & species == "nsp Maranon 3") |
           #(genus == "Exastilithoxus" & species == "sp Cuao") |
           (genus == "Hypancistrus" & species == "debilittera") |
           #(genus == "Hypostomus" & species == "taphorni") |
           #(genus == "Lasiancistrus" & species == "schomburgkii") |
           #(genus == "Leptoancistrus" & species == "nsp Atrato") |
           #(genus == "Lithoxancistrus" & species == "orinoco") |
           (genus == "Panaqolus" & species == "albomacula") |
           (genus == "Panaque" & species == "sp. L191") |
           #(genus == "Peckoltia" & species == "relictum") |
           (genus == "Pseudocanthicus" & species == "fordii") |
           (genus == "Pseudolithoxus" & species == "tigris") |
           #(genus == "Pseudoqolus" & species == "koko") |
           #(genus == "Scobinancistrus" & species == "aureatus") |
           (genus == "Brochiloricaria" & species == "chauliodon") |
           (genus == "Crossoloricaria" & species == "variegata") |
           (genus == "Metaloricaria" & species == "paucidens") |
           (genus == "Pseudohemiodon" & species == "unillanos")) %>% # |
           #(genus == "Pseudancistrus" & species == "genisetiger")) %>% 
  mutate(Legend = "Hypervolume")

hulls_filtered <- filtered %>%
  group_by(Legend) %>%
  slice(chull(Comp1, Comp2)) %>%
  ungroup()

# Adding centroid
centroid <- data.frame(
  PC1 = mean(filtered$Comp1),
  PC2 = mean(filtered$Comp2)
)

# Adding arrows
library(FNN)

# Nearest neighbours
nn <- get.knn(filtered[,c("Comp1","Comp2")], k = 1)

filtered$nn <- nn$nn.index[,1]

nn_arrows <- data.frame(
  x = filtered$Comp1,
  y = filtered$Comp2,
  xend = filtered$Comp1[filtered$nn],
  yend = filtered$Comp2[filtered$nn]
)

nn_arrows$mutual <- mapply(function(i,j){
  filtered$nn[j] == i
}, 1:nrow(filtered), filtered$nn)

# Centroid arrows
centroid_arrows <- data.frame(
  x = filtered$Comp1,
  y = filtered$Comp2,
  xend = centroid$PC1,
  yend = centroid$PC2
)

ggplot(data = filtered, aes(x = Comp1, y = Comp2, colour = Legend, fill = Legend)) +
  geom_polygon(data = hulls_filtered, aes(colour = "Hypervolume", fill = "Hypervolume"), 
               alpha = 0.4,
               colour = NA) +
  labs(x = "PC1", y = "PC2") +
  theme_cowplot() +
  geom_image(
    data = tibble(Comp1 = 0.16, Comp2 = 0.16),
    aes(x = Comp1, y = Comp2, image = "l191.png"),
    size = 0.5,
    inherit.aes = FALSE
  ) +
  geom_image(
    data = tibble(Comp1 = -0.44, Comp2 = -0.1),
    aes(x = Comp1, y = Comp2, image = "maranon.png"),
    size = 0.5,
    inherit.aes = FALSE
  ) +
  geom_image(
    data = tibble(Comp1 = -0.15, Comp2 = 0.15),
    aes(x = Comp1, y = Comp2, image = "anosteos.png"),
    size = 0.5,
    inherit.aes = FALSE
  ) +
  geom_image(
    data = tibble(Comp1 = 0.45, Comp2 = -0.05),
    aes(x = Comp1, y = Comp2, image = "jantjae.png"),
    size = 0.5,
    inherit.aes = FALSE
  ) +
  geom_image(
    data = tibble(Comp1 = 0.45, Comp2 = -0.22),
    aes(x = Comp1, y = Comp2, image = "variegata.png"),
    size = 0.5,
    inherit.aes = FALSE
  ) +
  coord_cartesian(ylim = c(-0.25,0.3), xlim = c(-0.55, 0.55)) +
  geom_point(shape = "square", size = 3) +
  geom_arrow_segment(data = subset(nn_arrows, !mutual),
               aes(x=x, y=y, xend=xend, yend=yend, colour = "Nearest Neighbour Distance"),
               length_head = unit(2, "mm"),
               resect_head = unit(2, "mm"),
               resect_fins = unit(2, "mm"),
               inherit.aes = F) +
  geom_arrow_segment(
    data = subset(nn_arrows, mutual),
    aes(x = x, y = y, xend = xend, yend = yend, colour = "Nearest Neighbour Distance"),
    length_head = unit(2, "mm"),
    length_fins = unit(2,"mm"),
    resect_head = unit(0.75, "mm"),
    resect_fins = unit(2, "mm"),
    inherit.aes = FALSE
  ) +
  geom_arrow_segment(
    data = centroid_arrows,
    aes(x = x, y = y, xend = xend, yend = yend, colour = "Centroid Distance"),
    length_head = unit(2,"mm"),
    resect_head = unit(2, "mm"),
    resect_fins = unit(2, "mm"),
    alpha = 0.6,
    inherit.aes = FALSE
  ) +
  geom_point(data = centroid, aes(PC1, PC2),
           colour = "magenta", size = 3, shape = 16, stroke = 2, inherit.aes = F) +
  scale_colour_manual(
    name = "Legend",
    values = c(
      "Centroid Distance" = "magenta",
      "Nearest Neighbour Distance" = "navyblue",
      "Hypervolume" = "skyblue2"
    )
  ) +
  scale_fill_manual(
    name = "Legend",
    values = c(
      "Hypervolume" = "skyblue2",
      "Centroid Distance" = "magenta",
      "Nearest Neighbour Distance" = "navyblue"
    )
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(
        linetype = c(1, 0, 1),
        shape = c(NA, 22, NA),
        fill = c(NA, "skyblue2", NA),
        size = c(0.8, 6, 0.8),
        alpha = c(1, 0.4, 1)
      )
    ),
    fill = "none"
  ) +
  theme(legend.title = element_blank(),
        legend.position = c(0.04, 1), #
        legend.justification = c("left", "top"))

# Plain plot with just points
ggplot(data = filtered, aes(x = Comp1, y = Comp2)) +
  labs(x = "PC1", y = "PC2") +
  geom_point(colour = "skyblue2", size = 3, shape = "square") +
  coord_cartesian(ylim = c(-0.25,0.19), xlim = c(-0.5, 0.55)) +
  geom_image(
    data = tibble(Comp1 = 0.16, Comp2 = 0.16),
    aes(x = Comp1, y = Comp2, image = "l191.png"),
    size = 0.5,
    inherit.aes = FALSE) +
  geom_image(
    data = tibble(Comp1 = -0.44, Comp2 = -0.1),
    aes(x = Comp1, y = Comp2, image = "maranon.png"),
    size = 0.5,
    inherit.aes = FALSE) +
  geom_image(
    data = tibble(Comp1 = -0.15, Comp2 = 0.15),
    aes(x = Comp1, y = Comp2, image = "anosteos.png"),
    size = 0.5,
    inherit.aes = FALSE) +
  geom_image(
    data = tibble(Comp1 = 0.45, Comp2 = -0.05),
    aes(x = Comp1, y = Comp2, image = "jantjae.png"),
    size = 0.5,
    inherit.aes = FALSE) +
  geom_image(
    data = tibble(Comp1 = 0.45, Comp2 = -0.22),
    aes(x = Comp1, y = Comp2, image = "variegata.png"),
    size = 0.5,
    inherit.aes = FALSE) +
  theme_cowplot()

# METRIC BY METRIC VISUAL

# Hypervolume
ggplot(data = filtered, aes(x = Comp1, y = Comp2, colour = Legend, fill = Legend)) +
  geom_polygon(data = hulls_filtered, aes(colour = "Hypervolume", fill = "Hypervolume"), 
               alpha = 0.4,
               colour = NA) +
  labs(x = "PC1", y = "PC2") +
  theme_cowplot() +
  coord_cartesian(ylim = c(-0.25,0.19), xlim = c(-0.5, 0.55)) +
  geom_point(shape = "square", size = 3) +
  scale_color_manual(values = "skyblue2") +
  scale_fill_manual(values = "skyblue2") +
  theme(legend.position = "none") +
  geom_image(
    data = tibble(Comp1 = 0.16, Comp2 = 0.16),
    aes(x = Comp1, y = Comp2, image = "l191.png"),
    size = 0.5,
    inherit.aes = FALSE
  ) +
  geom_image(
    data = tibble(Comp1 = -0.44, Comp2 = -0.1),
    aes(x = Comp1, y = Comp2, image = "maranon.png"),
    size = 0.5,
    inherit.aes = FALSE
  ) +
  geom_image(
    data = tibble(Comp1 = -0.15, Comp2 = 0.15),
    aes(x = Comp1, y = Comp2, image = "anosteos.png"),
    size = 0.5,
    inherit.aes = FALSE
  ) +
  geom_image(
    data = tibble(Comp1 = 0.45, Comp2 = -0.05),
    aes(x = Comp1, y = Comp2, image = "jantjae.png"),
    size = 0.5,
    inherit.aes = FALSE
  ) +
  geom_image(
    data = tibble(Comp1 = 0.45, Comp2 = -0.22),
    aes(x = Comp1, y = Comp2, image = "variegata.png"),
    size = 0.5,
    inherit.aes = FALSE
  )

# NND
ggplot(data = filtered, aes(x = Comp1, y = Comp2, colour = Legend)) +
  labs(x = "PC1", y = "PC2") +
  theme_cowplot() +
  coord_cartesian(ylim = c(-0.25,0.19), xlim = c(-0.5, 0.55)) +
  geom_point(shape = "square", size = 3) +
  scale_color_manual(values = c("skyblue2", "navyblue")) +
  geom_arrow_segment(data = subset(nn_arrows, !mutual),
                     aes(x=x, y=y, xend=xend, yend=yend, colour = "Nearest Neighbour Distance"),
                     length_head = unit(2, "mm"),
                     resect_head = unit(2, "mm"),
                     resect_fins = unit(2, "mm"),
                     inherit.aes = F) +
  geom_image(
    data = tibble(Comp1 = 0.16, Comp2 = 0.16),
    aes(x = Comp1, y = Comp2, image = "l191.png"),
    size = 0.5,
    inherit.aes = FALSE
  ) +
  geom_image(
    data = tibble(Comp1 = -0.44, Comp2 = -0.1),
    aes(x = Comp1, y = Comp2, image = "maranon.png"),
    size = 0.5,
    inherit.aes = FALSE
  ) +
  geom_image(
    data = tibble(Comp1 = -0.15, Comp2 = 0.15),
    aes(x = Comp1, y = Comp2, image = "anosteos.png"),
    size = 0.5,
    inherit.aes = FALSE
  ) +
  geom_image(
    data = tibble(Comp1 = 0.45, Comp2 = -0.05),
    aes(x = Comp1, y = Comp2, image = "jantjae.png"),
    size = 0.5,
    inherit.aes = FALSE
  ) +
  geom_image(
    data = tibble(Comp1 = 0.45, Comp2 = -0.22),
    aes(x = Comp1, y = Comp2, image = "variegata.png"),
    size = 0.5,
    inherit.aes = FALSE
  ) +
  geom_arrow_segment(
    data = subset(nn_arrows, mutual),
    aes(x = x, y = y, xend = xend, yend = yend, colour = "Nearest Neighbour Distance"),
    length_head = unit(2, "mm"),
    length_fins = unit(2,"mm"),
    resect_head = unit(0.75, "mm"),
    resect_fins = unit(2, "mm"),
    inherit.aes = FALSE
  ) +
  theme(legend.position = "none")

# CD
ggplot(data = filtered, aes(x = Comp1, y = Comp2, colour = Legend)) +
  labs(x = "PC1", y = "PC2") +
  theme_cowplot() +
  coord_cartesian(ylim = c(-0.25,0.19), xlim = c(-0.5, 0.55)) +
  geom_point(shape = "square", size = 3) +
  geom_image(
    data = tibble(Comp1 = 0.16, Comp2 = 0.16),
    aes(x = Comp1, y = Comp2, image = "l191.png"),
    size = 0.5,
    inherit.aes = FALSE
  ) +
  geom_image(
    data = tibble(Comp1 = -0.44, Comp2 = -0.1),
    aes(x = Comp1, y = Comp2, image = "maranon.png"),
    size = 0.5,
    inherit.aes = FALSE
  ) +
  geom_image(
    data = tibble(Comp1 = -0.15, Comp2 = 0.15),
    aes(x = Comp1, y = Comp2, image = "anosteos.png"),
    size = 0.5,
    inherit.aes = FALSE
  ) +
  geom_image(
    data = tibble(Comp1 = 0.45, Comp2 = -0.05),
    aes(x = Comp1, y = Comp2, image = "jantjae.png"),
    size = 0.5,
    inherit.aes = FALSE
  ) +
  geom_image(
    data = tibble(Comp1 = 0.45, Comp2 = -0.22),
    aes(x = Comp1, y = Comp2, image = "variegata.png"),
    size = 0.5,
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = c("magenta", "skyblue2")) +
  geom_arrow_segment(
    data = centroid_arrows,
    aes(x = x, y = y, xend = xend, yend = yend, colour = "Centroid Distance"),
    length_head = unit(2,"mm"),
    resect_head = unit(2, "mm"),
    resect_fins = unit(2, "mm"),
    alpha = 0.6,
    inherit.aes = FALSE
  ) +
  geom_point(data = centroid, aes(PC1, PC2),
             colour = "magenta", size = 3, shape = 16, stroke = 2, inherit.aes = F) +
  theme(legend.position = "none")

nn_arrows <- nn_arrows %>%
  mutate(nnd = as.numeric(nn$nn.dist))

# sdNND
# First, categorize the nnd values
nn_arrows$NND_category <- cut(nn_arrows$nnd, 
                              breaks = quantile(nn_arrows$nnd, probs = 0:4 / 4), 
                              include.lowest = TRUE, 
                              labels = c("Low", "Medium-Low", "Medium-High", "High"))

# Now, you can use scale_colour_manual to manually assign colors based on these categories
ggplot(data = filtered, aes(x = Comp1, y = Comp2, colour = Legend)) +
  labs(x = "PC1", y = "PC2") +
  theme_cowplot() +
  coord_cartesian(ylim = c(-0.25,0.19), xlim = c(-0.5, 0.55)) +
  geom_point(shape = "square", size = 3, colour = "skyblue2") +  # For point colors, can be modified if needed
  geom_arrow_segment(data = subset(nn_arrows, !mutual),
                     aes(x = x, y = y, xend = xend, yend = yend, colour = NND_category),  # Use categories for color
                     length_head = unit(2, "mm"),
                     resect_head = unit(2, "mm"),
                     resect_fins = unit(2, "mm"),
                     inherit.aes = F) +
  geom_arrow_segment(
    data = subset(nn_arrows, mutual),
    aes(x = x, y = y, xend = xend, yend = yend, colour = NND_category),  # Use categories for color
    length_head = unit(2, "mm"),
    length_fins = unit(2, "mm"),
    resect_head = unit(0.75, "mm"),
    resect_fins = unit(2, "mm"),
    inherit.aes = FALSE
  ) +
  scale_colour_manual(values = c("Low" = "skyblue", 
                                 "Medium-Low" = "lightblue", 
                                 "Medium-High" = "blue", 
                                 "High" = "navyblue")) +  # Assign colors to the categories
  geom_image(
    data = tibble(Comp1 = 0.16, Comp2 = 0.16),
    aes(x = Comp1, y = Comp2, image = "l191.png"),
    size = 0.5,
    inherit.aes = FALSE
  ) +
  geom_image(
    data = tibble(Comp1 = -0.44, Comp2 = -0.1),
    aes(x = Comp1, y = Comp2, image = "maranon.png"),
    size = 0.5,
    inherit.aes = FALSE
  ) +
  geom_image(
    data = tibble(Comp1 = -0.15, Comp2 = 0.15),
    aes(x = Comp1, y = Comp2, image = "anosteos.png"),
    size = 0.5,
    inherit.aes = FALSE
  ) +
  geom_image(
    data = tibble(Comp1 = 0.45, Comp2 = -0.05),
    aes(x = Comp1, y = Comp2, image = "jantjae.png"),
    size = 0.5,
    inherit.aes = FALSE
  ) +
  geom_image(
    data = tibble(Comp1 = 0.45, Comp2 = -0.22),
    aes(x = Comp1, y = Comp2, image = "variegata.png"),
    size = 0.5,
    inherit.aes = FALSE
  ) +
  theme(legend.position = "none")
