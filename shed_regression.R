
# Regressions and code for associated figures for shedding duration and peak 
# viral load (proxied by minimum Ct) with both pre-season HAI titer on a 
# categorical basis (Fig. 3 in main text) and on a continuous basis (Supporting 
# Information Fig. 8).


library(tidyverse)
library(dplyr)
library(lme4)
library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(lmtest)
library(nlme)
library(broom.mixed)
library(lmerTest)
library(mice)

source('shed_model_run.R')
shed_data <- read.csv("mock_dataset.csv")

# Adding in a clean date based on the first date of collection
shed_data$npsdatecol <- as.Date(shed_data$npsdatecol, format = "%m/%d/%Y")
shed_data$date <- as.integer(shed_data$npsdatecol - as.Date("2017-01-16"))

# Adding in age categories
shed_data$age_cat <- NA 
shed_data$age_cat <- ifelse(shed_data$age <= 4, "<5", shed_data$age_cat)
shed_data$age_cat <- ifelse(shed_data$age >= 5 & shed_data$age <= 11, "5-11", shed_data$age_cat)
shed_data$age_cat <- ifelse(shed_data$age >= 12 & shed_data$age <= 18, "12-18", shed_data$age_cat)
shed_data$age_cat <- ifelse(shed_data$age >= 19 & shed_data$age <= 40, "19-40", shed_data$age_cat)
shed_data$age_cat <- ifelse(shed_data$age >= 41, ">40", shed_data$age_cat)
shed_data$age_cat <- relevel(as.factor(shed_data$age_cat), ref = ">40")

# One observation per infection in shed_data 
shed_data <- shed_data %>% 
  filter(npsh3 == 1)
shed_data <- distinct(shed_data, indid, .keep_all = TRUE)

# Create a total duration and transform the Ct variable from deltaCt (x) to peak Ct
params_df <- params_df %>%
  mutate(tot_dur_est = wp + wc) %>%
  mutate(peakct = 37 - x)

# Combine data sets 
regress_data <- params_df %>%
  left_join(shed_data, by = "indid")

# Categorical titer variable
regress_data$hai_cat <- NA
regress_data$hai_cat <- ifelse(regress_data$flua_h3n2 <40 , "Titer<40", regress_data$hai_cat)
regress_data$hai_cat <- ifelse(regress_data$flua_h3n2 >=40 , "Titer>=40", regress_data$hai_cat)

# DURATION with categorical titer ################################################
# Function: sample posterior draws and fit regression
set.seed(123)
fit_once <- function(df) {
  df_sample <- df %>%
    group_by(indid_inf) %>% # Used indid_inf to distinguish when an individual 
    #had multiple infections within one season, this is not present in the mock data set
    slice_sample(n = 1) %>%
    ungroup()
  
  model <- lm(tot_dur_est ~ hai_cat + age_cat + sex + 
                  hiv_status, data = df_sample)
}
    # in paper lmer was used, with (1|year) as random effect to adjust for A(H3N2) 
    # being observed in multiple years

# Run model across multiple resampled datasets
list_mods <- vector("list", 100)
for (i in 1:100) {
  list_mods[[i]] <- fit_once(regress_data)
}

# Pool the results (Rubin’s rules)
pooled_mods <- pool(list_mods)
# Extract pooled estimates
dur_cat_est <- summary(pooled_mods) %>%
  mutate(lower = estimate - 1.96 * std.error) %>%
  mutate(upper = estimate + 1.96 * std.error) 

# Add intercept to the estimates, so the results are centered around the intercept 
intercept <- dur_cat_est$estimate[which(dur_cat_est$term == '(Intercept)')]
dur_cat_est <- dur_cat_est %>%
  mutate(estimate_int = estimate + intercept, lower = lower + intercept, 
         upper = upper + intercept)

# Clean labels for the plot
dur_cat_est <- dur_cat_est %>%
  filter(term != '(Intercept)') %>%
  mutate(estlabel = round(estimate_int, digits = 2))
dur_cat_est$estlabel <- ifelse(dur_cat_est$p.value < 0.001, paste(dur_cat_est$estlabel, "***"), 
                               dur_cat_est$estlabel)
dur_cat_est$estlabel <- ifelse(dur_cat_est$p.value < 0.01 & dur_cat_est$p.value >= 0.001, 
                           paste(dur_cat_est$estlabel, "**"), dur_cat_est$estlabel)
dur_cat_est$estlabel <- ifelse(dur_cat_est$p.value < 0.05 & dur_cat_est$p.value >= 0.01, 
                           paste(dur_cat_est$estlabel, "*"), dur_cat_est$estlabel)

# Create plot 
dur_cat_est$term <- factor(dur_cat_est$term, 
                       levels = c("hai_catTiter>=40", "age_cat<5", "age_cat5-11", 
                                  "age_cat12-18", "age_cat19-40", "sexMale",                                
                                  "hiv_statusPositive"))
ref <- data.frame(term = c('hai_catTiter<40', 'age_cat>40', 'sexfemale', 'hivneg'), 
                  estlabel = '', lower = 0, upper = 0, estimate = intercept)
order <- c("hivneg","hiv_statusPositive","sexfemale","sexMale","age_cat>40",
           "age_cat19-40","age_cat12-18",  "age_cat5-11", "age_cat<5", 
           "hai_catTiter<40", "hai_catTiter>=40")
dur_cat_plot <- ggplot(dur_cat_est, aes(x = term, y = estimate_int, ymin = lower, ymax = upper, 
                                label = estlabel)) +  
  annotate('rect', xmin = 0, xmax = 2.5, ymin = 0, ymax = 16, fill = "#0d0887", alpha = 0.5) +
  annotate('rect', xmin = 2.5, xmax = 4.5, ymin = 0, ymax = 16, fill = "#9c179e", alpha = 0.5) +
  annotate('rect', xmin = 4.5, xmax = 9.5, ymin = 0, ymax = 16, fill = "#ed7953", alpha = 0.5) +
  annotate('rect', xmin = 9.5, xmax = 11.5, ymin = 0, ymax = 16, fill = "#f0f921", alpha = 0.5) +
  geom_text(size = 3, nudge_x = .4, color = "black") +
  geom_linerange(linewidth = 0.5) + geom_point(size = 1.5) +
  geom_point(data = ref, aes (x = term, y = estimate), color = 'black', fill = 'white', pch = 21) +
  coord_flip(ylim = c(0, 16)) + 
  scale_y_continuous(name = "Shedding duration (days)") + 
  scale_x_discrete(name = "", limits = order, labels = c("HIV Negative",
                                                         "HIV Positive", "Sex - Female", 
                                                         "Sex - Male", "Age >40",
                                                         "Age 19-40", "Age 12-18",
                                                         "Age 5-11", "Age <5", "Pre-Season HAI Titer < 40",
                                                         "Pre-Season HAI Titer >=40 ")) +
  labs(color = "") + 
  theme_classic() +
  geom_hline(yintercept = intercept, color = 'grey', alpha = 0.8)
print(dur_cat_plot)


# PEAK VIRAL LOAD with categorical titer #########################################
set.seed(123)
fit_once <- function(df) {
  df_sample <- df %>%
    group_by(indid_inf) %>% 
    slice_sample(n = 1) %>%
    ungroup()
  
  model <- lm(peakct ~ hai_cat + age_cat + sex + 
                hiv_status, data = df_sample)
}
list_mods <- vector("list", 100)
for (i in 1:100) {
  list_mods[[i]] <- fit_once(regress_data)
}
pooled_mods <- pool(list_mods)
ct_cat_est <- summary(pooled_mods) %>%
  mutate(lower = estimate - 1.96 * std.error) %>%
  mutate(upper = estimate + 1.96 * std.error) 
intercept <- ct_cat_est$estimate[which(ct_cat_est$term == '(Intercept)')]
ct_cat_est <- ct_cat_est %>%
  mutate(estimate_int = estimate + intercept, lower = lower + intercept, 
         upper = upper + intercept)
ct_cat_est <- ct_cat_est %>%
  filter(term != '(Intercept)') %>%
  mutate(estlabel = round(estimate_int, digits = 2))
ct_cat_est$estlabel <- ifelse(ct_cat_est$p.value < 0.001, paste(ct_cat_est$estlabel, "***"), 
                               ct_cat_est$estlabel)
ct_cat_est$estlabel <- ifelse(ct_cat_est$p.value < 0.01 & ct_cat_est$p.value >= 0.001, 
                               paste(ct_cat_est$estlabel, "**"), ct_cat_est$estlabel)
ct_cat_est$estlabel <- ifelse(ct_cat_est$p.value < 0.05 & ct_cat_est$p.value >= 0.01, 
                               paste(ct_cat_est$estlabel, "*"), ct_cat_est$estlabel)
ct_cat_est$term <- factor(ct_cat_est$term, 
                           levels = c("hai_catTiter>=40", "age_cat<5", "age_cat5-11", 
                                      "age_cat12-18", "age_cat19-40", "sexMale",                                
                                      "hiv_statusPositive"))
ref <- data.frame(term = c('hai_catTiter<40', 'age_cat>40', 'sexfemale', 'hivneg'), 
                  estlabel = '', lower = 0, upper = 0, estimate = intercept)
order <- c("hivneg","hiv_statusPositive","sexfemale","sexMale","age_cat>40",
           "age_cat19-40","age_cat12-18",  "age_cat5-11", "age_cat<5", 
           "hai_catTiter<40", "hai_catTiter>=40")
ct_cat_plot <- ggplot(ct_cat_est, aes(x = term, y = estimate_int, ymin = lower, ymax = upper, 
                                        label = estlabel)) +  
  annotate('rect', xmin = 0, xmax = 2.5, ymin = 5, ymax = 37, fill = "#0d0887", alpha = 0.5) +
  annotate('rect', xmin = 2.5, xmax = 4.5, ymin = 5, ymax = 37, fill = "#9c179e", alpha = 0.5) +
  annotate('rect', xmin = 4.5, xmax = 9.5, ymin = 5, ymax = 37, fill = "#ed7953", alpha = 0.5) +
  annotate('rect', xmin = 9.5, xmax = 11.5, ymin = 5, ymax = 37, fill = "#f0f921", alpha = 0.5) +
  geom_text(size = 3, nudge_x = .4, color = "black") +
  geom_linerange(linewidth = 0.5) + geom_point(size = 1.5) +
  geom_point(data = ref, aes (x = term, y = estimate), color = 'black', fill = 'white', pch = 21) +
  coord_flip(ylim = c(5, 37)) + 
  scale_y_continuous(name = "Minimum Ct Of Episode") + 
  scale_x_discrete(name = "", limits = order, labels = c("HIV Negative",
                                                         "HIV Positive", "Sex - Female", 
                                                         "Sex - Male", "Age >40",
                                                         "Age 19-40", "Age 12-18",
                                                         "Age 5-11", "Age <5", "Pre-Season HAI Titer < 40",
                                                         "Pre-Season HAI Titer >=40 ")) +
  labs(color = "", caption = "Lower Minimum Ct is a Proxy For Higher Peak Viral Shedding") + 
  theme_classic() +
  geom_hline(yintercept = intercept, color = 'grey', alpha = 0.8)
print(ct_cat_plot)

# DURATION w/ continuous titer ################################################
set.seed(123)
fit_once <- function(df) {
  df_sample <- df %>%
    group_by(indid_inf) %>%
    slice_sample(n = 1) %>%
    ungroup()
  
  model <- lm(tot_dur_est ~ flua_h3n2 + age_cat + sex + 
                hiv_status, data = df_sample)
}
list_mods <- vector("list", 100)
for (i in 1:100) {
  list_mods[[i]] <- fit_once(regress_data)
}
pooled_mods <- pool(list_mods)
dur_con_est <- summary(pooled_mods) %>%
  mutate(lower = estimate - 1.96 * std.error) %>%
  mutate(upper = estimate + 1.96 * std.error) 

intercept <- dur_con_est$estimate[which(dur_con_est$term == '(Intercept)')]
dur_con_est <- dur_con_est %>%
  mutate(estimate_int = estimate + intercept, lower = lower + intercept, 
         upper = upper + intercept)
dur_con_est <- dur_con_est %>%
  filter(term != '(Intercept)') %>%
  mutate(estlabel = round(estimate_int, digits = 2))
dur_con_est$estlabel <- ifelse(dur_con_est$p.value < 0.001, paste(dur_con_est$estlabel, "***"), 
                               dur_con_est$estlabel)
dur_con_est$estlabel <- ifelse(dur_con_est$p.value < 0.01 & dur_con_est$p.value >= 0.001, 
                               paste(dur_con_est$estlabel, "**"), dur_con_est$estlabel)
dur_con_est$estlabel <- ifelse(dur_con_est$p.value < 0.05 & dur_con_est$p.value >= 0.01, 
                               paste(dur_con_est$estlabel, "*"), dur_con_est$estlabel)
dur_con_est$term <- factor(dur_con_est$term, 
                           levels = c("flua_h3n2", "age_cat<5", "age_cat5-11", 
                                      "age_cat12-18", "age_cat19-40", "sexMale",                                
                                      "hiv_statusPositive"))
ref <- data.frame(term = c('age_cat>40', 'sexfemale', 'hivneg'), 
                  estlabel = '', lower = 0, upper = 0, estimate = intercept)
order <- c("hivneg","hiv_statusPositive","sexfemale","sexMale","age_cat>40",
           "age_cat19-40","age_cat12-18",  "age_cat5-11", "age_cat<5", 
           "flua_h3n2")
dur_con_plot <- ggplot(dur_con_est, aes(x = term, y = estimate_int, ymin = lower, ymax = upper, 
                                        label = estlabel)) +  
  annotate('rect', xmin = 0, xmax = 2.5, ymin = 0, ymax = 16, fill = "#0d0887", alpha = 0.5) +
  annotate('rect', xmin = 2.5, xmax = 4.5, ymin = 0, ymax = 16, fill = "#9c179e", alpha = 0.5) +
  annotate('rect', xmin = 4.5, xmax = 9.5, ymin = 0, ymax = 16, fill = "#ed7953", alpha = 0.5) +
  annotate('rect', xmin = 9.5, xmax = 10.5, ymin = 0, ymax = 16, fill = "#f0f921", alpha = 0.5) +
  geom_text(size = 3, nudge_x = .4, color = "black") +
  geom_linerange(linewidth = 0.5) + geom_point(size = 1.5) +
  geom_point(data = ref, aes (x = term, y = estimate), color = 'black', fill = 'white', pch = 21) +
  coord_flip(ylim = c(0, 16)) + 
  scale_y_continuous(name = "Shedding duration (days)") + 
  scale_x_discrete(name = "", limits = order, labels = c("HIV Negative",
                                                         "HIV Positive", "Sex - Female", 
                                                         "Sex - Male", "Age >40",
                                                         "Age 19-40", "Age 12-18",
                                                         "Age 5-11", "Age <5",
                                                         "Pre-Season HAI Titer - 
                                              Four-Fold Increase")) +
  labs(color = "") + 
  theme_classic() +
  geom_hline(yintercept = intercept, color = 'grey', alpha = 0.8)
print(dur_con_plot)


# PEAK VIRAL LOAD w/ continuous titer #########################################
set.seed(123)
fit_once <- function(df) {
  df_sample <- df %>%
    group_by(indid_inf) %>% 
    slice_sample(n = 1) %>%
    ungroup()
  
  model <- lm(peakct ~ flua_h3n2 + age_cat + sex + 
                hiv_status, data = df_sample)
}
list_mods <- vector("list", 100)
for (i in 1:100) {
  list_mods[[i]] <- fit_once(regress_data)
}
pooled_mods <- pool(list_mods)
ct_con_est <- summary(pooled_mods) %>%
  mutate(lower = estimate - 1.96 * std.error) %>%
  mutate(upper = estimate + 1.96 * std.error) 
intercept <- ct_con_est$estimate[which(ct_con_est$term == '(Intercept)')]
ct_con_est <- ct_con_est %>%
  mutate(estimate_int = estimate + intercept, lower = lower + intercept, 
         upper = upper + intercept)
ct_con_est <- ct_con_est %>%
  filter(term != '(Intercept)') %>%
  mutate(estlabel = round(estimate_int, digits = 2))
ct_con_est$estlabel <- ifelse(ct_con_est$p.value < 0.001, paste(ct_con_est$estlabel, "***"), 
                              ct_con_est$estlabel)
ct_con_est$estlabel <- ifelse(ct_con_est$p.value < 0.01 & ct_con_est$p.value >= 0.001, 
                              paste(ct_con_est$estlabel, "**"), ct_con_est$estlabel)
ct_con_est$estlabel <- ifelse(ct_con_est$p.value < 0.05 & ct_con_est$p.value >= 0.01, 
                              paste(ct_con_est$estlabel, "*"), ct_con_est$estlabel)

# Create plot (Figure 3)
ct_con_est$term <- factor(ct_con_est$term, 
                          levels = c("flua_h3n2", "age_cat<5", "age_cat5-11", 
                                     "age_cat12-18", "age_cat19-40", "sexMale",                                
                                     "hiv_statusPositive"))
ref <- data.frame(term = c('age_cat>40', 'sexfemale', 'hivneg'), 
                  estlabel = '', lower = 0, upper = 0, estimate = intercept)
order <- c("hivneg","hiv_statusPositive","sexfemale","sexMale","age_cat>40",
           "age_cat19-40","age_cat12-18",  "age_cat5-11", "age_cat<5", 
           "flua_h3n2")
ct_con_plot <- ggplot(ct_con_est, aes(x = term, y = estimate_int, ymin = lower, ymax = upper, 
                                      label = estlabel)) +  
  annotate('rect', xmin = 0, xmax = 2.5, ymin = 5, ymax = 37, fill = "#0d0887", alpha = 0.5) +
  annotate('rect', xmin = 2.5, xmax = 4.5, ymin = 5, ymax = 37, fill = "#9c179e", alpha = 0.5) +
  annotate('rect', xmin = 4.5, xmax = 9.5, ymin = 5, ymax = 37, fill = "#ed7953", alpha = 0.5) +
  annotate('rect', xmin = 9.5, xmax = 10.5, ymin = 5, ymax = 37, fill = "#f0f921", alpha = 0.5) +
  geom_text(size = 3, nudge_x = .4, color = "black") +
  geom_linerange(linewidth = 0.5) + geom_point(size = 1.5) +
  geom_point(data = ref, aes (x = term, y = estimate), color = 'black', fill = 'white', pch = 21) +
  coord_flip(ylim = c(5, 37)) + 
  scale_y_continuous(name = "Minimum Ct Of Episode") + 
  scale_x_discrete(name = "", limits = order, labels = c("HIV Negative",
                                                         "HIV Positive", "Sex - Female", 
                                                         "Sex - Male", "Age >40",
                                                         "Age 19-40", "Age 12-18",
                                                         "Age 5-11", "Age <5", "Pre-Season HAI Titer < 40",
                                                         "Pre-Season HAI Titer >=40 ")) +
  labs(color = "", caption = "Lower Minimum Ct is a Proxy For Higher Peak Viral Shedding") + 
  theme_classic() +
  geom_hline(yintercept = intercept, color = 'grey', alpha = 0.8)
print(ct_con_plot)
