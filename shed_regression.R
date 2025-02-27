
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

source('shed_model_run.R')
shed_data$age_cat <- relevel(as.factor(shed_data$age_cat), ref = ">40")

# DURATION w/ categorical titer ################################################
shed_data$hai_cat <- NA
shed_data$hai_cat <- ifelse(shed_data$flua_h3n2 <40 , "Titer<40", shed_data$hai_cat)
shed_data$hai_cat <- ifelse(shed_data$flua_h3n2 >=40 , "Titer>=40", shed_data$hai_cat)

dur_cat <- lm(tot_duration ~ hai_cat + age_cat + sex + hiv_status, 
               data = shed_data) 
    # in paper lmer was used, with (1|year) as random effect to adjust for A(H3N2) 
    # being observed in multiple years
dur_cat_est <- tidy(dur_cat, conf.int = TRUE)
dur_cat_est_int <- dur_cat_est$estimate[dur_cat_est$term == "(Intercept)"]
dur_cat_est <- dur_cat_est %>%
  mutate(estimate_int = estimate + dur_cat_est_int) %>% 
  mutate(conf.low = conf.low + dur_cat_est_int) %>%
  mutate(conf.high = conf.high + dur_cat_est_int) %>%
  filter(term != '(Intercept)') %>%
  mutate(estlabel = round(estimate_int, digits = 2))

dur_cat_est$estlabel <- ifelse(dur_cat_est$p.value < 0.001, paste(dur_cat_est$estlabel, "***"), 
                               dur_cat_est$estlabel)
dur_cat_est$estlabel <- ifelse(dur_cat_est$p.value < 0.01 & dur_cat_est$p.value >= 0.001, 
                           paste(dur_cat_est$estlabel, "**"), dur_cat_est$estlabel)
dur_cat_est$estlabel <- ifelse(dur_cat_est$p.value < 0.05 & dur_cat_est$p.value >= 0.01, 
                           paste(dur_cat_est$estlabel, "*"), dur_cat_est$estlabel)
dur_cat_est$term <- factor(dur_cat_est$term, 
                       levels = c("hai_catTiter>=40", "age_cat<5", "age_cat5-11", 
                                  "age_cat12-18", "age_cat19-40", "sexMale",                                
                                  "hiv_statusPositive"))
ref <- data.frame(term = c('hai_catTiter<40', 'age_cat>40', 'sexfemale', 'hivneg'), 
                  std.error = 0, statistic = 0, df = 0, p.value = 0,
                  conf.low = 0, conf.high = 0, estlabel = '', estimate = dur_cat_est_int)
order <- c("hivneg","hiv_statusPositive","sexfemale","sexMale","age>40",
           "age_cat19-40","age_cat12-18",  "age_cat5-11", "age_cat<5", 
           "hai_catTiter<40", "hai_catTiter>=40")
dur_cat_plot <- ggplot(dur_cat_est, aes(x = term, y = estimate_int, ymin = conf.low, ymax = conf.high, 
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
  geom_hline(yintercept = dur_cat_est_int, color = 'grey', alpha = 0.8)
print(dur_cat_plot)


# PEAK VIRAL LOAD w/ categorical titer #########################################
ct_cat <- lm(meanct ~ hai_cat + age_cat + sex + hiv_status, 
             data = shed_data) 
ct_cat_est <- tidy(ct_cat, conf.int = TRUE)
ct_cat_est_int <- ct_cat_est$estimate[ct_cat_est$term == "(Intercept)"]
ct_cat_est <- ct_cat_est %>%
  mutate(estimate_int = estimate + ct_cat_est_int) %>% 
  mutate(conf.low = conf.low + ct_cat_est_int) %>%
  mutate(conf.high = conf.high + ct_cat_est_int) %>%
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
                  std.error = 0, statistic = 0, df = 0, p.value = 0,
                  conf.low = 0, conf.high = 0, estlabel = '', estimate = ct_cat_est_int)
order <- c("hivneg","hiv_statusPositive","sexfemale","sexMale","age>40",
           "age_cat19-40","age_cat12-18",  "age_cat5-11", "age_cat<5", 
           "hai_catTiter<40", "hai_catTiter>=40")
ct_cat_plot <- ggplot(ct_cat_est, aes(x = term, y = estimate_int, ymin = conf.low, ymax = conf.high, 
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
  labs(color = "") + 
  theme_classic() +
  labs(color = "", caption = "Lower Minimum Ct is a Proxy For Higher Peak Viral Shedding") +
  geom_hline(yintercept = ct_cat_est_int, color = 'grey', alpha = 0.8)
print(ct_cat_plot)


# DURATION w/ categorical titer ################################################

# Dividing all log(HAI) values by the average so that average can be used as a 
# reference group 
AVEH3N2 <- mean(log(shed_data$flua_h3n2, base = 4))
shed_data$hai_con <- log(shed_data$flua_h3n2, base = 4) / AVEH3N2

dur_con <- lm(tot_duration ~ hai_con + age_cat + sex + hiv_status, 
              data = shed_data) 
dur_con_est <- tidy(dur_con, conf.int = TRUE)
dur_con_est_int <- dur_con_est$estimate[dur_con_est$term == "(Intercept)"]
dur_con_est <- dur_con_est %>%
  mutate(estimate_int = estimate + dur_con_est_int) %>% 
  mutate(conf.low = conf.low + dur_con_est_int) %>%
  mutate(conf.high = conf.high + dur_con_est_int) %>%
  filter(term != '(Intercept)') %>%
  mutate(estlabel = round(estimate_int, digits = 2))

dur_con_est$estlabel <- ifelse(dur_con_est$p.value < 0.001, paste(dur_con_est$estlabel, "***"), 
                               dur_con_est$estlabel)
dur_con_est$estlabel <- ifelse(dur_con_est$p.value < 0.01 & dur_con_est$p.value >= 0.001, 
                               paste(dur_con_est$estlabel, "**"), dur_con_est$estlabel)
dur_con_est$estlabel <- ifelse(dur_con_est$p.value < 0.05 & dur_con_est$p.value >= 0.01, 
                               paste(dur_con_est$estlabel, "*"), dur_con_est$estlabel)
dur_con_est$term <- factor(dur_con_est$term, 
                           levels = c("hai_con", "age_cat<5", "age_cat5-11", 
                                      "age_cat12-18", "age_cat19-40", "sexMale",                                
                                      "hiv_statusPositive"))
ref <- data.frame(term = c('hai_ave', 'age_cat>40', 'sexfemale', 'hivneg'), 
                  std.error = 0, statistic = 0, df = 0, p.value = 0,
                  conf.low = 0, conf.high = 0, estlabel = '', estimate = dur_con_est_int)
order <- c("hivneg","hiv_statusPositive","sexfemale","sexMale","age>40",
           "age_cat19-40","age_cat12-18",  "age_cat5-11", "age_cat<5", 
           "hai_ave", "hai_con")
dur_con_plot <- ggplot(dur_con_est, aes(x = term, y = estimate_int, ymin = conf.low, ymax = conf.high, 
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
                                                         "Age 5-11", "Age <5", "Pre-Season HAI Titer - 
                                                         Average",
                                                         "Pre-Season HAI Titer - 
                                              Four-Fold Increase")) +
  labs(color = "") + 
  theme_classic() +
  geom_hline(yintercept = dur_con_est_int, color = 'grey', alpha = 0.8)
print(dur_con_plot)


# PEAK VIRAL LOAD w/ categorical titer #########################################
ct_con <- lm(meanct ~ hai_con + age_cat + sex + hiv_status, 
              data = shed_data) 
ct_con_est <- tidy(ct_con, conf.int = TRUE)
ct_con_est_int <- ct_con_est$estimate[ct_con_est$term == "(Intercept)"]
ct_con_est <- ct_con_est %>%
  mutate(estimate_int = estimate + ct_con_est_int) %>% 
  mutate(conf.low = conf.low + ct_con_est_int) %>%
  mutate(conf.high = conf.high + ct_con_est_int) %>%
  filter(term != '(Intercept)') %>%
  mutate(estlabel = round(estimate_int, digits = 2))

ct_con_est$estlabel <- ifelse(ct_con_est$p.value < 0.001, paste(ct_con_est$estlabel, "***"), 
                              ct_con_est$estlabel)
ct_con_est$estlabel <- ifelse(ct_con_est$p.value < 0.01 & ct_con_est$p.value >= 0.001, 
                              paste(ct_con_est$estlabel, "**"), ct_con_est$estlabel)
ct_con_est$estlabel <- ifelse(ct_con_est$p.value < 0.05 & ct_con_est$p.value >= 0.01, 
                              paste(ct_con_est$estlabel, "*"), ct_con_est$estlabel)
ct_con_est$term <- factor(ct_con_est$term, 
                          levels = c("hai_con", "age_cat<5", "age_cat5-11", 
                                     "age_cat12-18", "age_cat19-40", "sexMale",                                
                                     "hiv_statusPositive"))
ref <- data.frame(term = c('hai_ave', 'age_cat>40', 'sexfemale', 'hivneg'), 
                  std.error = 0, statistic = 0, df = 0, p.value = 0,
                  conf.low = 0, conf.high = 0, estlabel = '', estimate = ct_con_est_int)
order <- c("hivneg","hiv_statusPositive","sexfemale","sexMale","age>40",
           "age_cat19-40","age_cat12-18",  "age_cat5-11", "age_cat<5", 
           "hai_ave", "hai_con")
ct_con_plot <- ggplot(ct_con_est, aes(x = term, y = estimate_int, ymin = conf.low, ymax = conf.high, 
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
                                                         "Age 5-11", "Age <5", "Pre-Season HAI Titer - 
                                                         Average",
                                                         "Pre-Season HAI Titer - 
                                              Four-Fold Increase")) +
  labs(color = "") + 
  theme_classic() +
  labs(color = "", caption = "Lower Minimum Ct is a Proxy For Higher Peak Viral Shedding") +
  geom_hline(yintercept = ct_cat_est_int, color = 'grey', alpha = 0.8)
print(ct_con_plot)