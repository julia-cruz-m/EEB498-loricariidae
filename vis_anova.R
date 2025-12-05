# ANOVA Test + EMMEAN Visualization

library(tidyverse)
library(cowplot)
library(ggpubr)
library(rstatix)
library(broom)
library(emmeans)

synthetic = read.csv("final_all.csv") %>% 
  mutate(type = "Synthetic")
empirical = read.csv("empirical_all.csv") %>% 
  mutate(type = "Empirical")

final = bind_rows(synthetic, empirical)

# NND ----
# Test for slope homogeneity
homogen <- final %>% 
  anova_test(avg_nnd ~ type * sp_rich)

# Summary of slope difference b/w assemblage types
nnd <- lm(avg_nnd ~ sp_rich * type, data = final)
summary(nnd)

nnd_metrics <- augment(nnd) %>%
  select(-.hat, -.sigma, -.fitted) # Remove details
head(nnd_metrics, 3)

# THIS SECTION DOES NOT WORK BC SAMPLE SIZE IS TOO HIGH
# Assess normality of residuals using shapiro wilk test
# shapiro_test(nnd_metrics$.resid)

# Levine's test: homoegeneity of residual variances
lev <- nnd_metrics %>% 
  levene_test(.resid ~ type)

# Identify outliers
out <- nnd_metrics %>% 
  filter(abs(.std.resid) > 3) %>%
  as.data.frame()

# Remove outliers
no_out <- nnd_metrics %>% 
  filter(abs(.std.resid) < 3)

# ANCOVA
anc <- final %>% 
  anova_test(avg_nnd ~ sp_rich + type)
get_anova_table(anc) # Sig diff in avg_nnd b/w assemblage types

# Estimated adjusted means
means <- final %>% 
  emmeans_test(
    avg_nnd ~ type, covariate = sp_rich,
    p.adjust.method = "bonferroni"
  )
means

get_emmeans(means)

# Visualization
means <- means %>% 
  add_xy_position(x = "type", fun = "mean_se")
ggline(get_emmeans(means), x = "type", y = "emmean") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) + 
  stat_pvalue_manual(means, hide.ns = TRUE, tip.length = FALSE) +
  labs(subtitle = get_test_label(anc, detailed = TRUE),
       caption = get_pwc_label(means),
       y = "Estimated Marginal Mean",
       x = "Assemblage Type") +
  theme_cowplot()

# sdNND ----
# Test for slope homogeneity
sdhomogen <- final %>% 
  anova_test(sd_nnd ~ type * sp_rich)

# Summary of slope difference b/w assemblage types
sdnnd <- lm(sd_nnd ~ sp_rich * type, data = final)
summary(nnd)

sdnnd_metrics <- augment(sdnnd) %>%
  select(-.hat, -.sigma, -.fitted) # Remove details
head(sdnnd_metrics, 3)

# THIS SECTION DOES NOT WORK BC SAMPLE SIZE IS TOO HIGH
# Assess normality of residuals using shapiro wilk test
# shapiro_test(nnd_metrics$.resid)

# Levine's test: homoegeneity of residual variances
sdlev <- sdnnd_metrics %>% 
  levene_test(.resid ~ type)

# Identify outliers
sdout <- sdnnd_metrics %>% 
  filter(abs(.std.resid) > 3) %>%
  as.data.frame()

# Remove outliers
sdno_out <- sdnnd_metrics %>% 
  filter(abs(.std.resid) < 3)

# ANCOVA
sdanc <- final %>% 
  anova_test(sd_nnd ~ sp_rich + type)
get_anova_table(sdanc)

# Estimated adjusted means
sdmeans <- final %>% 
  emmeans_test(
    sd_nnd ~ type, covariate = sp_rich,
    p.adjust.method = "bonferroni"
  )
sdmeans

get_emmeans(sdmeans)

# Visualization
sdmeans <- sdmeans %>% 
  add_xy_position(x = "type", fun = "mean_se")
ggline(get_emmeans(sdmeans), x = "type", y = "emmean") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) + 
  stat_pvalue_manual(sdmeans, hide.ns = TRUE, tip.length = FALSE) +
  labs(subtitle = get_test_label(sdanc, detailed = TRUE),
       caption = get_pwc_label(sdmeans),
       y = "Estimated Marginal Mean",
       x = "Assemblage Type") +
  theme_cowplot()
# CD ----
# Test for slope homogeneity
cd_homogen <- final %>% 
  anova_test(cd ~ type * sp_rich)

# Summary of slope difference b/w assemblage types
cd <- lm(cd ~ sp_rich * type, data = final)
summary(cd)

cd_metrics <- augment(cd) %>%
  select(-.hat, -.sigma, -.fitted) # Remove details
head(cd_metrics, 3)

# THIS SECTION DOES NOT WORK BC SAMPLE SIZE IS TOO HIGH
# Assess normality of residuals using shapiro wilk test
# shapiro_test(nnd_metrics$.resid)

# Levine's test: homoegeneity of residual variances
cd_lev <- cd_metrics %>% 
  levene_test(.resid ~ type)

# Identify outliers
cd_out <- cd_metrics %>% 
  filter(abs(.std.resid) > 3) %>%
  as.data.frame()

# Remove outliers
cd_no_out <- cd_metrics %>% 
  filter(abs(.std.resid) < 3)

# ANCOVA
cd_anc <- final %>% 
  anova_test(cd ~ sp_rich + type)
get_anova_table(cd_anc)

# Estimated adjusted means
cd_means <- final %>% 
  emmeans_test(
    cd ~ type, covariate = sp_rich,
    p.adjust.method = "bonferroni"
  )
cd_means

get_emmeans(cd_means)

# Visualization
cd_means <- cd_means %>% 
  add_xy_position(x = "type", fun = "mean_se")
ggline(get_emmeans(cd_means), x = "type", y = "emmean") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) + 
  stat_pvalue_manual(cd_means, hide.ns = TRUE, tip.length = FALSE) +
  labs(subtitle = get_test_label(cd_anc, detailed = TRUE),
       caption = get_pwc_label(cd_means),
       y = "Estimated Marginal Mean",
       x = "Assemblage Type") +
  theme_cowplot()

# HV3D ----
# Test for slope homogeneity
hv3_homogen <- final %>% 
  anova_test(hv3d ~ type * sp_rich)

# Summary of slope difference b/w assemblage types
hv3 <- lm(hv3d ~ sp_rich * type, data = final)
summary(hv3)

hv3_metrics <- augment(hv3) %>%
  select(-.hat, -.sigma, -.fitted) # Remove details
head(hv3_metrics, 3)

# THIS SECTION DOES NOT WORK BC SAMPLE SIZE IS TOO HIGH
# Assess normality of residuals using shapiro wilk test
# shapiro_test(nnd_metrics$.resid)

# Levine's test: homoegeneity of residual variances
hv3_lev <- hv3_metrics %>% 
  levene_test(.resid ~ type)

# Identify outliers
hv3_out <- hv3_metrics %>% 
  filter(abs(.std.resid) > 3) %>%
  as.data.frame()

# Remove outliers
hv3_no_out <- hv3_metrics %>% 
  filter(abs(.std.resid) < 3)

# ANCOVA
hv3_anc <- final %>% 
  anova_test(hv3d ~ sp_rich + type)
get_anova_table(hv3_anc)

# Estimated adjusted means
hv3_means <- final %>% 
  emmeans_test(
    hv3d ~ type, covariate = sp_rich,
    p.adjust.method = "bonferroni"
  )
hv3_means

get_emmeans(hv3_means)

# Visualization
hv3_means <- hv3_means %>% 
  add_xy_position(x = "type", fun = "mean_se")
ggline(get_emmeans(hv3_means), x = "type", y = "emmean") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) + 
  stat_pvalue_manual(hv3_means, hide.ns = TRUE, tip.length = FALSE) +
  labs(subtitle = get_test_label(hv3_anc, detailed = TRUE),
       caption = get_pwc_label(hv3_means),
       y = "Estimated Marginal Mean",
       x = "Assemblage Type") +
  theme_cowplot()

# HV6D ----
# Test for slope homogeneity
hv6_homogen <- final %>% 
  anova_test(hv6d ~ type * sp_rich)

# Summary of slope difference b/w assemblage types
hv6 <- lm(hv6d ~ sp_rich * type, data = final)
summary(hv6)

hv6_metrics <- augment(hv6) %>%
  select(-.hat, -.sigma, -.fitted) # Remove details
head(hv6_metrics, 3)

# THIS SECTION DOES NOT WORK BC SAMPLE SIZE IS TOO HIGH
# Assess normality of residuals using shapiro wilk test
# shapiro_test(nnd_metrics$.resid)

# Levine's test: homoegeneity of residual variances
hv6_lev <- hv6_metrics %>% 
  levene_test(.resid ~ type)

# Identify outliers
hv6_out <- hv6_metrics %>% 
  filter(abs(.std.resid) > 3) %>%
  as.data.frame()

# Remove outliers
hv6_no_out <- hv6_metrics %>% 
  filter(abs(.std.resid) < 3)

# ANCOVA - THIS SECTION DOES NOT WORK
# I BELIEVE IT IS BECAUSE HV IS TOO CLOSE TO 0 (RESCALE?)
hv6_anc <- final %>% 
  anova_test(hv6d ~ sp_rich + type)
get_anova_table(hv6_anc)

# Estimated adjusted means
hv6_means <- final %>% 
  emmeans_test(
    hv6d ~ type, covariate = sp_rich,
    p.adjust.method = "bonferroni"
  )
hv6_means

get_emmeans(hv6_means)

# Visualization
hv6_means <- hv6_means %>% 
  add_xy_position(x = "type", fun = "mean_se")
ggline(get_emmeans(hv6_means), x = "type", y = "emmean") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) + 
  stat_pvalue_manual(hv6_means, hide.ns = TRUE, tip.length = FALSE) +
  labs(subtitle = get_test_label(hv6_anc, detailed = TRUE),
       caption = get_pwc_label(hv6_means),
       y = "Estimated Marginal Mean",
       x = "Assemblage Type"
  )

