
# Regressions and code for associated figures for the household transmission model 
# with both pre-season HAI titer on a categorical basis (Fig 5. in main text) and 
# on a continuous basis (Supporting Information Fig. 9).

library(tidyverse)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(lme4)
library(sjPlot)
library(jtools)
library(merTools)
library(zoo)

source('shed_model_run.R')
shed_data$age_cat <- relevel(as.factor(shed_data$age_cat), ref = ">40")

# Need to formulate data in a matrix format for the household transmission model; 
# matrix has a observation row for every single day in the season of study for each participant
indid_l <- c(unique(shed_data$indid))
profile <- data.frame()
for (i in indid_l) {
  sub <- shed_data %>% filter(indid == i)
  profile <- rbind(profile, sub[1, ])
}
rm(i, indid_l, sub)

matrix <- data.frame(indid = rep(c(profile$indid), each = 286), date = rep(0:285, times = 500),
                     hh = rep(c(profile$hh_id), each = 286), site = rep(c(profile$site), each = 286), 
                     preseasonhai = rep(c(profile$flua_h3n2), each = 286),
                     age = rep(c(profile$age), each = 286), age_cat = rep(c(profile$age_cat), each = 286),
                     sex = rep(c(profile$sex), each = 286), hiv = rep(c(profile$hiv_status), each = 286), 
                     hh_cat = rep(c(profile$hh_size), each = 286))

# adding start date and viral load information for the primary infections 
inf_dat <- data.frame()
indid_inf <- unique(shed_data$indid_inf[which(!is.na(shed_data$indid_inf))])
for (n in indid_inf) {
  sub <- shed_data %>% filter(indid_inf == n)
  inf_dat <- rbind(inf_dat, sub[1, ])
}
rm(sub, indid_inf, n)

matrix$h3_present <- NA
matrix$prolif_load <- NA
matrix$clear_load <- NA
new <- data.frame()
indid_l <- c(inf_dat$indid_inf)
for (i in indid_l) {
  ind <- inf_dat$indid[which(inf_dat$indid_inf == i)]
  sub <- matrix %>% filter(indid == ind)
  s <- inf_dat$round_start[which(inf_dat$indid_inf == i)]
  sub$h3_present[1:(s-1)] <- 0
  # one day buffer
  sub$h3_present[s] <- 1   
  sub$h3_present[(s+1):286] <- NA
  
  sub$prolif_load <- ifelse(sub$date <= inf_dat$start_point[which(inf_dat$indid_inf == i)], 0, sub$prolif_load)
  sub$clear_load <- ifelse(sub$date <= inf_dat$start_point[which(inf_dat$indid_inf == i)], 0, sub$clear_load)
  sub$prolif_load <- ifelse(sub$date > inf_dat$start_point[which(inf_dat$indid_inf == i)] &
                              sub$date <= inf_dat$center_point[which(inf_dat$indid_inf == i)],
                            (sub$date - inf_dat$start_point[which(inf_dat$indid_inf == i)]) * 
                              inf_dat$prolif_slope[which(inf_dat$indid_inf == i)],
                            sub$prolif_load)
  sub$clear_load <- ifelse(sub$date > inf_dat$start_point[which(inf_dat$indid_inf == i)] &
                             sub$date <= inf_dat$center_point[which(inf_dat$indid_inf == i)],0, sub$clear_load)
  sub$clear_load <- ifelse(sub$date > inf_dat$center_point[which(inf_dat$indid_inf == i)] &
                             sub$date <= inf_dat$end_point[which(inf_dat$indid_inf == i)],
                           inf_dat$meanct[which(inf_dat$indid_inf == i)] + 
                             ((sub$date - inf_dat$center_point[which(inf_dat$indid_inf == i)]) * 
                                inf_dat$clear_slope[which(inf_dat$indid_inf == i)]),
                           sub$clear_load)
  sub$prolif_load <- ifelse(sub$date > inf_dat$center_point[which(inf_dat$indid_inf == i)] &
                              sub$date <= inf_dat$end_point[which(inf_dat$indid_inf == i)], 0, sub$prolif_load)
  sub$prolif_load <- ifelse(sub$date >= inf_dat$end_point[which(inf_dat$indid_inf == i)], 0, sub$prolif_load)
  sub$clear_load <- ifelse(sub$date >= inf_dat$end_point[which(inf_dat$indid_inf == i)], 0, sub$clear_load)
  new <- rbind(new, sub)
}

# no infections recorded
indid_one <- unique(c(inf_dat$indid))
noinf <- setdiff(unique(c(matrix$indid)), indid_one)
matrix$h3_present <- ifelse(matrix$indid %in% noinf, 0, matrix$h3_present)
matrix$prolif_load <- ifelse(matrix$indid %in% noinf, 0, matrix$prolif_load)
matrix$clear_load <- ifelse(matrix$indid %in% noinf, 0, matrix$clear_load)
noinf_df <- matrix %>% filter(indid %in% noinf)

matrix <- rbind(new, noinf_df) %>% 
  arrange(indid)
rm(i, s, sub, ind, indid_l, inf_dat, new, indid_one, noinf, noinf_df)

# household load
matrix$hh_prof_load <- NA
matrix$hh_clear_load <- NA
hh_l <- c(unique(matrix$hh))
for (h in hh_l) {
  for (d in 0:286){
    matrix$hh_prof_load[which(matrix$hh == h & matrix$date == d)] <- 
      sum(matrix$prolif_load[which(matrix$hh == h & matrix$date == d)])
    matrix$hh_clear_load[which(matrix$hh == h & matrix$date == d)] <- 
      sum(matrix$clear_load[which(matrix$hh == h & matrix$date == d)])
  }
}
rm(h, d, hh_l)
# remove individuals own household load 
matrix <- matrix %>%
  mutate(hh_prof_load = hh_prof_load - prolif_load) %>%
  mutate(hh_clear_load = hh_clear_load - clear_load) %>% 
  mutate(ind_lod_tot = prolif_load + clear_load)

# add proxy for community surveillance - number of infections present at the time
# split by study site and by year
matrix$community <- 0
matrix$community <- ifelse(matrix$ind_lod_tot > 0, matrix$community - 1, matrix$community)
for (s in c("Agincourt", "Klerksdorp")) {
  for (d in 0:286) { 
    matrix$community[which(matrix$date == d & matrix$site == s)] <- 
      (matrix$community[which(matrix$date == d & matrix$site == s)] + 
         length(which(matrix$ind_lod_tot[which(matrix$date == d & matrix$site == s)] > 0))) / 
      length(matrix$indid[which(matrix$date == d & matrix$site == s)])
  }
}
rm(d, s)
# moving window smoothing - 4 day window to span between samples then transform
matrix$comsmooth <- rollapply(matrix$community, width = 4, function(...) {round(mean(...), digits = 3)}, partial = TRUE)
matrix$comsmooth <- matrix$comsmooth * 100

# exclude observations after the infection point
matrix <- matrix %>% filter(!is.na(h3_present))

################################################################################
# Household transmission model regression and associated figure with titer on 
# a categorical scale (main text).

matrix$age_cat <- relevel(as.factor(matrix$age_cat), ref = ">40")
matrix <- matrix %>% mutate(hh_load = hh_prof_load + hh_clear_load)
matrix <- matrix %>% mutate(hh_loadscale = hh_load / 10)
matrix$hai_cat <- NA 
matrix$hai_cat <- ifelse(matrix$preseasonhai <40, "titer<40", matrix$hai_cat)
matrix$hai_cat <- ifelse(matrix$preseasonhai >=40, "titer>=40", matrix$hai_cat)
matrix$hh_cat <- ifelse(matrix$hh_cat == "5-Mar", "3-5", matrix$hh_cat)
matrix$hh_cat <- ifelse(matrix$hh_cat == "10-Jun", "6-10", matrix$hh_cat)

trans_mod_cat <- glm(h3_present ~ comsmooth + hh_loadscale + hai_cat + age_cat + 
                   sex + hiv + hh_cat,
                 data = matrix, family = poisson())
trans_cat_est <- data.frame(summary(trans_mod_cat)[["coefficients"]])%>%
  tibble::rownames_to_column("term") %>%
  mutate(LL = exp(Estimate - (1.96 * Std..Error))) %>%
  mutate(UL = exp(Estimate + (1.96 * Std..Error))) %>%
  mutate(Estimate = exp(Estimate)) %>%
  rename("pvalue" = 'Pr...z..') %>%
  dplyr::select(term, Estimate, LL, UL, pvalue) %>%
  mutate(estlabel = round(Estimate, digits = 3))

trans_cat_est$estlabel <- ifelse(trans_cat_est$pvalue < 0.001, paste(trans_cat_est$estlabel, "***"), 
                              trans_cat_est$estlabel)
trans_cat_est$estlabel <- ifelse(trans_cat_est$pvalue < 0.01 & trans_cat_est$pvalue >= 0.001, 
                              paste(trans_cat_est$estlabel, "**"), trans_cat_est$estlabel)
trans_cat_est$estlabel <- ifelse(trans_cat_est$pvalue < 0.05 & trans_cat_est$pvalue >= 0.01, 
                              paste(trans_cat_est$estlabel, "*"), trans_cat_est$estlabel)
trans_cat_est$term <- factor(trans_cat_est$term, 
                          levels = c("comsmooth", "hh_loadscale", "hai_cattiter>=40", 
                                     "age_cat<5", "age_cat5-11", 
                                     "age_cat12-18", "age_cat19-40", "sexMale", 
                                     "hivPositive", "hh_cat3-5", "hh_cat6-10"))

ref <- data.frame(term = c('hai_cattiter<40', 'age>40', 'sexfemale', 'hivneg', 'hh_cat>10'),
                  estimate = 1, std.error = 0, statistic = 0, df = 0, p.value = 0,
                  LL = 0, UL = 0, estlabel = '')

order <- c("hh_cat>10", "hh_cat6-10", "hh_cat3-5",
           "hivneg", "hivPositive", "sexfemale","sexMale","age>40",
           "age_cat19-40","age_cat12-18",  "age_cat5-11", "age_cat<5", 
           'hai_cattiter<40', "hai_cattiter>=40", "hh_loadscale", "comsmooth")

transcat_plot <- ggplot(trans_cat_est, aes(x = term, y = Estimate, ymin = LL, ymax = UL, 
                                        label = estlabel)) +
  scale_x_discrete(limits = order, labels = c("Household Size >10",
                                              "Household Size 6-10", 
                                              "Household Size 3-5", "HIV Negative", 
                                              "HIV Positive", "Sex - Female", 
                                              "Sex - Male", "Age >40",
                                              "Age 19-40", "Age 12-18",
                                              "Age 5-11", "Age <5", "Pre-Season HAI Titer <40", 
                                              "Pre-Season HAI Titer >=40",
                                              "Household FOI 
                                              (Step by 10 Ct)", "Community FOI")) +
  scale_y_continuous(trans = 'log2') + #plot hazard ratio on log2 scale
  coord_flip(ylim = c(0.04, 14)) + 
  annotate('rect', xmin = 0, xmax = 3.5, ymin = 0, ymax = 14, fill = "#0d0887", alpha = 0.5) +
  annotate('rect', xmin = 3.5, xmax = 5.5, ymin = 0, ymax = 14, fill = "#5c01a6", alpha = 0.5) +
  annotate('rect', xmin = 5.5, xmax = 7.5, ymin = 0, ymax = 14, fill = "#9c179e", alpha = 0.5) +
  annotate('rect', xmin = 7.5, xmax = 12.5, ymin = 0, ymax = 14,fill = "#cc4778", alpha = 0.5) +
  annotate('rect', xmin = 12.5, xmax = 14.5, ymin = 0, ymax = 14, fill = "#ed7953", alpha = 0.5) +
  annotate('rect', xmin = 14.5, xmax = 15.5, ymin = 0, ymax = 14, fill = "#fdb42f", alpha = 0.5) +
  annotate('rect', xmin = 15.5, xmax = 16.5, ymin = 0, ymax = 14, fill = "#f0f921", alpha = 0.5) +
  geom_hline(yintercept = 1, color = 'grey', alpha = 0.8) +
  geom_linerange(linewidth = 0.3) + geom_point(size = 1.5) +
  geom_text(size = 3, nudge_x = .4, color = "black") +
  geom_point(data = ref, aes (x = term, y = estimate), color = 'black', fill = 'white', pch = 21) +
  labs(y = "Hazard Ratio", x = "", color = "") +
  theme_classic() 
print(transcat_plot)


################################################################################
# Household transmission model regression and associated figure with titer on 
# a continuous scale (supporting information).
AVEH3N2 <- mean(log(profile$flua_h3n2, base = 4))
matrix$hai_con <- log(matrix$preseasonhai, base = 4) / AVEH3N2

trans_mod_con <- glm(h3_present ~ comsmooth + hh_loadscale + hai_con + age_cat + 
                       sex + hiv + hh_cat,
                     data = matrix, family = poisson())
trans_con_est <- data.frame(summary(trans_mod_cat)[["coefficients"]])%>%
  tibble::rownames_to_column("term") %>%
  mutate(LL = exp(Estimate - (1.96 * Std..Error))) %>%
  mutate(UL = exp(Estimate + (1.96 * Std..Error))) %>%
  mutate(Estimate = exp(Estimate)) %>%
  rename("pvalue" = 'Pr...z..') %>%
  dplyr::select(term, Estimate, LL, UL, pvalue) %>%
  mutate(estlabel = round(Estimate, digits = 3))

trans_con_est$estlabel <- ifelse(trans_con_est$pvalue < 0.001, paste(trans_con_est$estlabel, "***"), 
                                 trans_con_est$estlabel)
trans_con_est$estlabel <- ifelse(trans_con_est$pvalue < 0.01 & trans_con_est$pvalue >= 0.001, 
                                 paste(trans_con_est$estlabel, "**"), trans_con_est$estlabel)
trans_con_est$estlabel <- ifelse(trans_con_est$pvalue < 0.05 & trans_con_est$pvalue >= 0.01, 
                                 paste(trans_con_est$estlabel, "*"), trans_con_est$estlabel)
trans_con_est$term <- factor(trans_con_est$term, 
                             levels = c("comsmooth", "hh_loadscale", "hai_con", 
                                        "age_cat<5", "age_cat5-11", 
                                        "age_cat12-18", "age_cat19-40", "sexMale", 
                                        "hivPositive", "hh_cat3-5", "hh_cat6-10"))

ref <- data.frame(term = c('hai_ave', 'age>40', 'sexfemale', 'hivneg', 'hh_cat>10'),
                  estimate = 1, std.error = 0, statistic = 0, df = 0, p.value = 0,
                  LL = 0, UL = 0, estlabel = '')

order <- c("hh_cat>10", "hh_cat6-10", "hh_cat3-5",
           "hivneg", "hivPositive", "sexfemale","sexMale","age>40",
           "age_cat19-40","age_cat12-18",  "age_cat5-11", "age_cat<5", 
           'hai_ave', "hai_con", "hh_loadscale", "comsmooth")

transcon_plot <- ggplot(trans_con_est, aes(x = term, y = Estimate, ymin = LL, ymax = UL, 
                                           label = estlabel)) +
  scale_x_discrete(limits = order, labels = c("Household Size >10",
                                              "Household Size 6-10", 
                                              "Household Size 3-5", "HIV Negative", 
                                              "HIV Positive", "Sex - Female", 
                                              "Sex - Male", "Age >40",
                                              "Age 19-40", "Age 12-18",
                                              "Age 5-11", "Age <5", "Pre-Season HAI Titer - 
                                              Average", 
                                              "Pre-Season HAI Titer - 
                                              Four-Fold Increase",
                                              "Household FOI 
                                              (Step by 10 Ct)", "Community FOI")) +
  scale_y_continuous(trans = 'log2') + #plot hazard ratio on log2 scale
  coord_flip(ylim = c(0.04, 14)) + 
  annotate('rect', xmin = 0, xmax = 3.5, ymin = 0, ymax = 14, fill = "#0d0887", alpha = 0.5) +
  annotate('rect', xmin = 3.5, xmax = 5.5, ymin = 0, ymax = 14, fill = "#5c01a6", alpha = 0.5) +
  annotate('rect', xmin = 5.5, xmax = 7.5, ymin = 0, ymax = 14, fill = "#9c179e", alpha = 0.5) +
  annotate('rect', xmin = 7.5, xmax = 12.5, ymin = 0, ymax = 14,fill = "#cc4778", alpha = 0.5) +
  annotate('rect', xmin = 12.5, xmax = 14.5, ymin = 0, ymax = 14, fill = "#ed7953", alpha = 0.5) +
  annotate('rect', xmin = 14.5, xmax = 15.5, ymin = 0, ymax = 14, fill = "#fdb42f", alpha = 0.5) +
  annotate('rect', xmin = 15.5, xmax = 16.5, ymin = 0, ymax = 14, fill = "#f0f921", alpha = 0.5) +
  geom_hline(yintercept = 1, color = 'grey', alpha = 0.8) +
  geom_linerange(linewidth = 0.3) + geom_point(size = 1.5) +
  geom_text(size = 3, nudge_x = .4, color = "black") +
  geom_point(data = ref, aes (x = term, y = estimate), color = 'black', fill = 'white', pch = 21) +
  labs(y = "Hazard Ratio", x = "", color = "") +
  theme_classic() 
print(transcon_plot)



