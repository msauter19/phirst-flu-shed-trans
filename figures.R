
# Code to produce the figures included in the main text.

library(ggplot2)
library(cowplot)
library(dplyr)
library(lares)

################################################################################
# Fig. 1 Distribution of pre-season serum influenza antibody titers by age group

fig1_dat <- read.csv("Fig1_titer_by_age.csv")

h1_violin <- ggplot(fig1_dat, aes(x = age_cat, 
                                 y = log10(flua_h1n1pdm_gm))) +
  geom_violin(fill = "#9bc7e0", alpha = 0.8, color = "white", linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.4) +
  scale_y_continuous(breaks = c(log10(5), log10(80), log10(1280)),
    minor_breaks = c(1, log10(20), log10(40), log10(160), log10(320), log10(640)),
    labels = c("5", "80","1280"), limits = c(0.5, 3.2)) +
  labs(x = "Age at Enrollment (years)", y = "Pre-season HAI titer",
    title = "Pre-season A(H1N1)pdm09 HAI") +
  scale_x_discrete(limits = c("<5", "5-11", "12-18", "19-40", ">40")) +
  theme_light()

h3_violin <- ggplot(fig1_dat, aes(x = age_cat, 
                                  y = log10(flua_h3n2_gm))) +
  geom_violin(fill = "#28607f", alpha = 0.8, color = "white", linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.4) +
  scale_y_continuous(breaks = c(log10(5), log10(80), log10(1280)),
                     minor_breaks = c(1, log10(20), log10(40), log10(160), log10(320), log10(640)),
                     labels = c("5", "80","1280"), limits = c(0.5, 3.2)) +
  labs(x = "Age at Enrollment (years)", y = "Pre-season HAI titer",
       title = "Pre-season A(H3N2) HAI") +
  scale_x_discrete(limits = c("<5", "5-11", "12-18", "19-40", ">40")) +
  theme_light()

vic_violin <- ggplot(fig1_dat, aes(x = age_cat, 
                                  y = log10(flub_victoria_gm))) +
  geom_violin(fill = "#a5da76", alpha = 0.8, color = "white", linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.4) +
  scale_y_continuous(breaks = c(log10(5), log10(80), log10(1280)),
                     minor_breaks = c(1, log10(20), log10(40), log10(160), log10(320), log10(640)),
                     labels = c("5", "80","1280"), limits = c(0.5, 3.2)) +
  labs(x = "Age at Enrollment (years)", y = "Pre-season HAI titer",
       title = "Pre-season B/Victoria HAI") +
  scale_x_discrete(limits = c("<5", "5-11", "12-18", "19-40", ">40")) +
  theme_light()

yam_violin <- ggplot(fig1_dat, aes(x = age_cat, 
                                   y = log10(flub_yamagata_gm))) +
  geom_violin(fill = "#3c611b", alpha = 0.8, color = "white", linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.4) +
  scale_y_continuous(breaks = c(log10(5), log10(80), log10(1280)),
                     minor_breaks = c(1, log10(20), log10(40), log10(160), log10(320), log10(640)),
                     labels = c("5", "80","1280"), limits = c(0.5, 3.2)) +
  labs(x = "Age at Enrollment (years)", y = "Pre-season HAI titer",
       title = "Pre-season B/Yamagata HAI") +
  scale_x_discrete(limits = c("<5", "5-11", "12-18", "19-40", ">40")) +
  theme_light()

fig1_plot <- plot_grid(h1_violin, h3_violin, vic_violin, yam_violin, ncol = 2)

fig1_plot

################################################################################
# Fig. 2. Influenza viral shedding kinetics estimates for each subtype/lineage. 

fig2_dat <- read.csv("Fig2_shedding_kinetic_mean_estimates.csv")

h1_ind <- fig2_dat %>% filter(TYPE == "AH1")
h1_pmean <- -mean(h1_ind$meanwp)
h1_ctmean <- mean(h1_ind$meanct)
h1_cmean <- mean(h1_ind$meanwc)
h3_ind <- fig2_dat %>% filter(TYPE == "AH3")
h3_pmean <- -mean(h3_ind$meanwp)
h3_ctmean <- mean(h3_ind$meanct)
h3_cmean <- mean(h3_ind$meanwc)
vic_ind <- fig2_dat %>% filter(TYPE == "BVic")
vic_pmean <- -mean(vic_ind$meanwp)
vic_ctmean <- mean(vic_ind$meanct)
vic_cmean <- mean(vic_ind$meanwc)
yam_ind <- fig2_dat %>% filter(TYPE == "BYam")
yam_pmean <- -mean(yam_ind$meanwp)
yam_ctmean <- mean(yam_ind$meanct)
yam_cmean <- mean(yam_ind$meanwc)

means_df <- data.frame(x = c(h1_pmean, 0, h1_cmean, h3_pmean, 0, h3_cmean,
                             vic_pmean, 0, vic_cmean, yam_pmean, 0, yam_cmean), 
                       y = c(37, h1_ctmean, 37, 37, h3_ctmean, 37,
                             37, vic_ctmean, 37, 37, yam_ctmean, 37),
                       TYPE = c("AH1", "AH1", "AH1", 
                                "AH3", "AH3", "AH3",
                                "BVic", "BVic", "BVic",
                                "BYam", "BYam", "BYam"))
ci_df <- data.frame(x = c(-ci_var(h1_ind, meanwp)$upper_ci, 
                          0, ci_var(h1_ind, meanwc)$upper_ci, 
                          ci_var(h1_ind, meanwc)$lower_ci, 0, 
                          -ci_var(h1_ind, meanwp)$lower_ci,
                          -ci_var(h3_ind, meanwp)$upper_ci, 
                          0, ci_var(h3_ind, meanwc)$upper_ci, 
                          ci_var(h3_ind, meanwc)$lower_ci, 0, 
                          -ci_var(h3_ind, meanwp)$lower_ci,
                          -ci_var(vic_ind, meanwp)$upper_ci, 
                          0, ci_var(vic_ind, meanwc)$upper_ci, 
                          ci_var(vic_ind, meanwc)$lower_ci, 0, 
                          -ci_var(vic_ind, meanwp)$lower_ci,
                          -ci_var(yam_ind, meanwp)$upper_ci, 
                          0, ci_var(yam_ind, meanwc)$upper_ci, 
                          ci_var(yam_ind, meanwc)$lower_ci, 0, 
                          -ci_var(yam_ind, meanwp)$lower_ci), 
                    y = c(37, ci_var(h1_ind, meanct)$lower_ci, 37, 37,
                          ci_var(h1_ind, meanct)$upper_ci, 37,
                          37, ci_var(h3_ind, meanct)$lower_ci, 37, 37,
                          ci_var(h3_ind, meanct)$upper_ci, 37,
                          37, ci_var(vic_ind, meanct)$lower_ci, 37, 37,
                          ci_var(vic_ind, meanct)$upper_ci, 37,
                          37, ci_var(yam_ind, meanct)$lower_ci, 37, 37,
                          ci_var(yam_ind, meanct)$upper_ci, 37),
                    TYPE = c(rep('AH1', 6), rep('AH3', 6), 
                             rep('BVic', 6), rep('BYam', 6)))
iqr_df <- data.frame(x = c(-quantile(h1_ind$meanwp, 0.75, na.rm = T), 
                           0, quantile(h1_ind$meanwc, 0.75, na.rm = T), 
                           quantile(h1_ind$meanwc, 0.25, na.rm = T), 0, 
                           -quantile(h1_ind$meanwp, 0.25, na.rm = T),
                           -quantile(h3_ind$meanwp, 0.75, na.rm = T), 
                           0, quantile(h3_ind$meanwc, 0.75, na.rm = T), 
                           quantile(h3_ind$meanwc, 0.25, na.rm = T), 0, 
                           -quantile(h3_ind$meanwp, 0.25, na.rm = T),
                           -quantile(vic_ind$meanwp, 0.75, na.rm = T), 
                           0, quantile(vic_ind$meanwc, 0.75, na.rm = T), 
                           quantile(vic_ind$meanwc, 0.25, na.rm = T), 0, 
                           -quantile(vic_ind$meanwp, 0.25, na.rm = T),
                           -quantile(yam_ind$meanwp, 0.75, na.rm = T), 
                           0, quantile(yam_ind$meanwc, 0.75, na.rm = T), 
                           quantile(yam_ind$meanwc, 0.25, na.rm = T), 0, 
                           -quantile(yam_ind$meanwp, 0.25, na.rm = T)), 
                     y = c(37, quantile(h1_ind$meanct, 0.25, na.rm = T), 37, 37,
                           quantile(h1_ind$meanct, 0.75, na.rm = T), 37,
                           37, quantile(h3_ind$meanct, 0.25, na.rm = T), 37, 37,
                           quantile(h3_ind$meanct, 0.75, na.rm = T), 37,
                           37, quantile(vic_ind$meanct, 0.25, na.rm = T), 37, 37,
                           quantile(vic_ind$meanct, 0.75, na.rm = T), 37,
                           37, quantile(yam_ind$meanct, 0.25, na.rm = T), 37, 37,
                           quantile(yam_ind$meanct, 0.75, na.rm = T), 37),
                     TYPE = c(rep('AH1', 6), rep('AH3', 6), 
                              rep('BVic', 6), rep('BYam', 6)))

kin_plot <- ggplot(data = ci_df, aes(x = x, y = y, fill = TYPE)) + 
  geom_polygon() +
  geom_polygon(data = iqr_df, aes(x = x, y = y, fill = TYPE), alpha = 0.5) +
  geom_line(data = means_df, aes(x = x, y = y), linewidth = 0.7, color = "white") +
  scale_y_reverse() +
  coord_cartesian(ylim = c(37,15), xlim = c(-4,8), clip = "on") +
  geom_vline(aes(xintercept = 0), color = "black", alpha = 0.6) + 
  scale_fill_manual(values = c("#9bc7e0", "#28607f","#a5da76","#3c611b")) +
  labs(x = "Time since minimum Ct (days)", y = "Cycle threshold (Ct)") +
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4, 6, 8)) +
  facet_wrap(~TYPE, ncol = 2, labeller = 
               as_labeller(c('AH1' = "A(H1N1)pdm09", 'AH3' = "A(H3N2)",
                             'BVic' = "B/Victoria", 'BYam' = "B/Yamagata"))) +
  theme_bw(base_size = 12) + 
  theme(legend.position="none", strip.background = element_blank())

pbox <- ggplot(data = fig2_dat, aes(x = TYPE, y = meanwp, fill = TYPE)) + 
  geom_boxplot() + 
  scale_fill_manual(values = c(AH1 = "#9bc7e0", AH3 = "#28607f",
                               BVic = "#a5da76", BYam = "#3c611b")) +
  scale_x_discrete(limits = c("AH1", "AH3", "BVic", "BYam")) + 
  theme_classic(base_size = 12) + 
  theme(axis.ticks = element_blank(), axis.title.x = element_blank()) +
  labs(title = "Proliferation Duration", y = "Duration (days)") +
  guides(fill="none") +
  scale_y_continuous(trans = "log2", breaks = c(1, 2, 4, 8, 16, 32), 
                     limits = c(0.5, 32))

cbox <- ggplot(data = fig2_dat, aes(x = TYPE, y = meanwc, fill = TYPE)) + 
  geom_boxplot() + 
  scale_fill_manual(values = c(AH1 = "#9bc7e0", AH3 = "#28607f",
                               BVic = "#a5da76", BYam = "#3c611b")) +
  scale_x_discrete(limits = c("AH1", "AH3", "BVic", "BYam")) + 
  theme_classic(base_size = 12) + 
  theme(axis.ticks = element_blank(), axis.title.y = element_blank(),
        axis.title.x = element_blank()) + 
  labs(title = "Clearance Duration")+
  guides(fill="none") +
  scale_y_continuous(trans = "log2", breaks = c(1, 2, 4, 8, 16, 32),
                     limits = c(0.5, 32))

fig2_dat <- fig2_dat %>% mutate(tot_duration = meanwp + meanwc)

tbox <- ggplot(data = fig2_dat, aes(x = TYPE, y = tot_duration, fill = TYPE)) + 
  geom_boxplot() +
  scale_fill_manual(values = c(AH1 = "#9bc7e0", AH3 = "#28607f",
                               BVic = "#a5da76", BYam = "#3c611b")) +
  scale_x_discrete(limits = c("AH1", "AH3", "BVic", "BYam")) + 
  theme_classic(base_size = 12) + 
  theme(axis.ticks = element_blank(), axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  labs(title = "Total Duration")+
  guides(fill="none") +
  scale_y_continuous(trans = "log2", breaks = c(1, 2, 4, 8, 16, 32),
                     limits = c(0.5, 32))

ctbox <- ggplot(data = fig2_dat, aes(x = TYPE, y = meanct, fill = TYPE)) + 
  geom_boxplot() + 
  scale_fill_manual(values = c(AH1 = "#9bc7e0", AH3 = "#28607f",
                               BVic = "#a5da76", BYam = "#3c611b")) +
  scale_x_discrete(limits = c("AH1", "AH3", "BVic", "BYam")) + 
  theme_classic(base_size = 12) + 
  theme(axis.ticks = element_blank(), axis.title.x = element_blank()) +
  scale_y_reverse() +
  labs(title = "Minimum Ct", y = "Minimum Ct")+
  guides(fill="none")

text_size <- theme(axis.title = element_text(size = 11), axis.text  = element_text(size = 10),
  plot.title = element_text(size = 13, face = "plain"), 
  strip.text = element_text(size = 13, face = "plain"))
kin_plot <- kin_plot + text_size
pbox <- pbox + text_size
cbox <- cbox + text_size
tbox <- tbox + text_size
ctbox <- ctbox + text_size

box_plot <- plot_grid(pbox, cbox, tbox, ctbox, nrow = 1, labels = c("B", "C", "D", "E"), align = "h")
fig2_plot <- plot_grid(kin_plot, box_plot, labels = c("A", ""), ncol = 1, rel_heights = c(1.5,1))

fig2_plot

################################################################################
# Fig. 3. Predictors of influenza virus shedding characteristics by 
# subtype/lineage based on multivariable linear regression.

# Data provided are model-estimates. The model input includes individual-level 
# data that cannot be shared. The regression that produces these estimates can be
# simulated with pseudo data in shed_regression.R; this file also includes
# the code for data manipulation of the model outputs for figure generation.

fig3_dur_est <- read.csv('Fig3_shed_duration_model_estimates.csv')
fig3_dur_int <- read.csv('Fig3_shed_duration_model_intercepts.csv')

fig3_minct_est <- read.csv('Fig3_shed_min_ct_model_estimates.csv')
fig3_minct_int <- read.csv('Fig3_shed_min_ct_model_intercepts.csv')

order <- c("hivneg","cleanhivUnknown", "cleanhivPositive CD4 >=200",
           "cleanhivPositive CD4 <200","sexfemale","sexMale","age>40",
           "age_cat19-40","age_cat12-18",  "age_cat5-11", "age_cat<5", 
           "cutoff140Titer<40", "cutoff140Titer>=40")
ref <- data.frame(effect = rep('ref', 16), 
                   group = c(rep('A(H1N1)pdm09', 4), rep('A(H3N2)', 4),
                             rep('B/Yamagata', 4), rep('B/Victoria', 4)),
                   term = c(rep(c('cutoff140Titer<40', 'age>40', 'sexfemale', 
                                  'hivneg'), 4)), estlabel = '', lower = 0, upper = 0)

ref_dur <- ref %>% left_join(fig3_dur_int, by = "group")
fig3_dur_est$term <- factor(fig3_dur_est$term,
                            levels = c("cutoff140Titer>=40", "age_cat<5", 
                                       "age_cat5-11", "age_cat12-18", "age_cat19-40", 
                                       "sexMale","cleanhivPositive CD4 <200", 
                                       "cleanhivPositive CD4 >=200", "cleanhivUnknown"))

shed_dur_plot <- ggplot(fig3_dur_est, aes(x = term, y = estimate_int, ymin = lower, ymax = upper, 
                                      label = estlabel)) +  
  annotate('rect', xmin = 0, xmax = 4.5, ymin = 0, ymax = 16, fill = "#0d0887", alpha = 0.5) +
  annotate('rect', xmin = 4.5, xmax = 6.5, ymin = 0, ymax = 16, fill = "#9c179e", alpha = 0.5) +
  annotate('rect', xmin = 6.5, xmax = 11.5, ymin = 0, ymax = 16, fill = "#ed7953", alpha = 0.5) +
  annotate('rect', xmin = 11.5, xmax = 13.5, ymin = 0, ymax = 16, fill = "#f0f921", alpha = 0.5) +
  geom_text(size = 3, nudge_x = .4, color = "black") +
  geom_linerange(linewidth = 0.5) + geom_point(size = 1.5) +
  geom_point(data = ref_dur, aes (x = term, y = estimate), color = 'black', fill = 'white', pch = 21) +
  coord_flip(ylim = c(0, 16)) + 
  scale_y_continuous(name = "Shedding duration (days)") + 
  scale_x_discrete(name = "", limits = order, labels = c("HIV Negative",
                                                         "HIV Status Unknown", "PLWH CD4 >=200", 
                                                         "PLWH CD4 <200", "Sex - Female", 
                                                         "Sex - Male", "Age >40",
                                                         "Age 19-40", "Age 12-18",
                                                         "Age 5-11", "Age <5", "Pre-Season HAI Titer <40",
                                                         "Pre-Season HAI Titer >=40")) +
  facet_wrap(~group, nrow = 1) + 
  geom_hline(data = fig3_dur_int, aes(yintercept = estimate), color = 'grey', alpha = 0.8) +
  theme_classic() +
  theme(legend.position="none",
        strip.background = element_blank(), panel.spacing.x = unit(1.5, "lines"))

ref_ct <- ref %>% left_join(fig3_minct_int, by = "group")
fig3_minct_est$term <- factor(fig3_minct_est$term, 
                          levels = c("cutoff140Titer>=40", "age_cat<5", "age_cat5-11", 
                                     "age_cat12-18", "age_cat19-40", "sexMale",                                
                                     "cleanhivPositive CD4 <200", "cleanhivPositive CD4 >=200",   
                                     "cleanhivUnknown"))

shed_minct_plot <- ggplot(fig3_minct_est, aes(x = term, y = estimate_int, ymin = lower, ymax = upper, 
                                     label = estlabel)) +  
  annotate('rect', xmin = 0, xmax = 4.5, ymin = 5, ymax = 37, fill = "#0d0887", alpha = 0.5) +
  annotate('rect', xmin = 4.5, xmax = 6.5, ymin = 5, ymax = 37, fill = "#9c179e", alpha = 0.5) +
  annotate('rect', xmin = 6.5, xmax = 11.5, ymin = 5, ymax = 37, fill = "#ed7953", alpha = 0.5) +
  annotate('rect', xmin = 11.5, xmax = 13.5, ymin = 5, ymax = 37, fill = "#f0f921", alpha = 0.5) +
  geom_text(size = 3, nudge_x = .4, color = "black") +
  geom_linerange(linewidth = 0.5) + geom_point(size = 1.5) +
  geom_point(data = ref_ct, aes (x = term, y = estimate), color = 'black', fill = 'white', pch = 21) +
  coord_flip(ylim = c(5, 37)) + 
  scale_y_continuous(name = "Minimum Ct Of Episode") + 
  scale_x_discrete(name = "", limits = order, labels = c("HIV Negative",
                                                         "HIV Status Unknown", "PLWH CD4 >=200", 
                                                         "PLWH CD4 <200", "Sex - Female", 
                                                         "Sex - Male", "Age >40",
                                                         "Age 19-40", "Age 12-18",
                                                         "Age 5-11", "Age <5", "Pre-Season HAI Titer <40",
                                                         "Pre-Season HAI Titer >=40")) +
  labs(caption = "Lower Minimum Ct is a Proxy For Higher Peak Viral Shedding") +
  facet_wrap(~group, nrow = 1) + 
  geom_hline(data = fig3_minct_int, aes(yintercept = estimate), color = 'grey', alpha = 0.8) +
  theme_classic() +
  theme(legend.position="none",
        strip.background = element_blank(), panel.spacing.x = unit(1.5, "lines"))

fig3_plot <- plot_grid(shed_dur_plot, shed_minct_plot, ncol = 1)

fig3_plot

################################################################################
# Fig. 4. Visualization of key variables of the influenza household transmission 
# model. 

# Code provided does not include reference to the calculation of community 
# prevalence or total household viral load; these data were produced as part of 
# household transmission model which included individual-level data. This model, 
# and relevant calculations, can be simulated with pseudo data in 
# transmission_model.R.

fig4_dat <- read.csv('Fig4_transmission_model_visualization.csv')

fig4_dat$hh_load <- ifelse(fig4_dat$hh_load == 0, NA, fig4_dat$hh_load)

ag_house <- fig4_dat %>% filter(site == "Agincourt")
ag_house <- ag_house %>% arrange(round_start, age)
ag_house <- ag_house %>%
  mutate(indid = factor(indid, unique(indid))) %>%
  mutate(yeardate = as.Date(yeardate, format='%m/%d/%Y'))

agcom <- ggplot(ag_house)+
  geom_raster(aes(x = as.Date(yeardate), y = site, fill = comsmooth), hjust =0, vjust =0) +
  theme_classic() + 
  labs(x = NULL, y = NULL, fill = NULL, title = "Agincourt Site (Rural)") +
  scale_fill_gradientn(colours = c("#cfefd0", "#5ec962"),limits=c(0, 7), breaks = c(0, 6), na.value = "transparent") +
  scale_x_date(expand=c(0,0), limits = c(as.Date("2018-07-01"), as.Date("2018-11-01")), 
               date_breaks = "1 month") +
  scale_y_discrete(expand=c(0,0)) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.ontop = TRUE,
        panel.background = element_rect(fill = NA, color = 'black'),
        plot.margin = unit(c(1,1,1,1), "pt"),
        legend.key.height = unit(0.1, 'cm'),
        legend.key.width = unit(0.8, 'cm')) +
  geom_vline(xintercept = as.Date("2018-08-01"), color = "black", linetype = 2)+
  geom_vline(xintercept = as.Date("2018-09-01"), color = "black", linetype = 2)+
  geom_vline(xintercept = as.Date("2018-10-01"), color = "black", linetype = 2)

aghh <- ggplot(ag_house)+
  geom_raster(aes(x = as.Date(yeardate), y = site, fill = hh_load), hjust =0, vjust =0) +
  theme_classic() +
  labs(x = NULL, y = NULL, fill = NULL) +
  scale_fill_gradientn(colours = c("#d2daec", "#3b528b"),limits=c(0, 55), breaks = c(0, 25, 50), na.value = "transparent") +
  scale_x_date(expand=c(0,0), limits = c(as.Date("2018-07-01"), as.Date("2018-11-01")), 
               date_breaks = "1 month") +
  scale_y_discrete(expand=c(0,0)) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(), 
        panel.grid = element_blank(),
        panel.ontop = TRUE,
        panel.background = element_rect(fill = NA, color = 'black'),
        plot.margin = unit(c(1,1,1,1), "pt"),
        legend.key.height = unit(0.2, 'cm')) +
  geom_vline(xintercept = as.Date("2018-08-01"), color = "black", linetype = 2)+
  geom_vline(xintercept = as.Date("2018-09-01"), color = "black", linetype = 2)+
  geom_vline(xintercept = as.Date("2018-10-01"), color = "black", linetype = 2)

agind <- ggplot(ag_house) +
  geom_raster(aes(x = as.Date(yeardate), y = indid, fill = ind_lod_tot), hjust =0, vjust =0) +
  theme_classic() +
  labs(x = "Time of Year", y = NULL, fill = NULL) +
  scale_fill_gradientn(colours = c("#f2c0fe", "#440154"),limits=c(1, 37), na.value = "transparent") +
  scale_x_date(expand=c(0,0), limits = c(as.Date("2018-07-01"), as.Date("2018-11-01")), 
               date_breaks = "1 month", date_labels = "%b %d") +
  scale_y_discrete(expand=c(0,0), limits = rev) + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle=50, hjust = 0.9),
        panel.background = element_rect(fill = NA, color = 'black'),
        plot.margin = unit(c(1,1,1,1), "pt"),
        legend.key.height = unit(0.4, 'cm')) +
  geom_hline(yintercept = 0 + 0:5, colour = "black") +
  geom_vline(xintercept = as.Date("2018-08-01"), color = "black", linetype = 2)+
  geom_vline(xintercept = as.Date("2018-09-01"), color = "black", linetype = 2)+
  geom_vline(xintercept = as.Date("2018-10-01"), color = "black", linetype = 2)

kl_house <- fig4_dat %>% filter(site == "Klerksdorp")
kl_house <- kl_house %>% arrange(round_start, age)
kl_house <- kl_house %>%
  mutate(indid = factor(indid, unique(indid))) %>%
  mutate(yeardate = as.Date(yeardate, format='%m/%d/%Y'))

klcom <- ggplot(kl_house)+
  geom_raster(aes(x = as.Date(yeardate), y = site, fill = comsmooth), hjust =0, vjust =0) +
  theme_classic() + 
  labs(x = NULL, y = NULL, fill = NULL, title = "Klerksdorp Site (Urban)") +
  scale_fill_gradientn(colours = c("#cfefd0", "#5ec962"),limits=c(0, 7), breaks = c(0, 6), na.value = "transparent") +
  scale_x_date(expand=c(0,0), limits = c(as.Date("2018-07-01"), as.Date("2018-11-01")), 
               date_breaks = "1 month") +
  scale_y_discrete(expand=c(0,0)) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.ontop = TRUE,
        panel.background = element_rect(fill = NA, color = 'black'),
        legend.position = "none",
        plot.margin = unit(c(1, 1, 1, 1), "pt")) +
  geom_vline(xintercept = as.Date("2018-08-01"), color = "black", linetype = 2)+
  geom_vline(xintercept = as.Date("2018-09-01"), color = "black", linetype = 2)+
  geom_vline(xintercept = as.Date("2018-10-01"), color = "black", linetype = 2)

klhh <- ggplot(kl_house)+
  geom_raster(aes(x = as.Date(yeardate), y = site, fill = hh_load), hjust =0, vjust =0) +
  theme_classic() +
  labs(x = NULL, y = NULL, fill = NULL) +
  scale_fill_gradientn(colours = c("#d2daec", "#3b528b"),limits=c(0, 55), breaks = c(0, 25, 50), na.value = "transparent") +
  scale_x_date(expand=c(0,0), limits = c(as.Date("2018-07-01"), as.Date("2018-11-01")), 
               date_breaks = "1 month") +
  scale_y_discrete(expand=c(0,0)) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(), 
        panel.grid = element_blank(),
        panel.ontop = TRUE,
        panel.background = element_rect(fill = NA, color = 'black'),
        legend.position = "none",
        plot.margin = unit(c(1, 1, 1, 1), "pt")) +
  geom_vline(xintercept = as.Date("2018-08-01"), color = "black", linetype = 2)+
  geom_vline(xintercept = as.Date("2018-09-01"), color = "black", linetype = 2)+
  geom_vline(xintercept = as.Date("2018-10-01"), color = "black", linetype = 2)

klind <- ggplot(kl_house) +
  geom_raster(aes(x = as.Date(yeardate), y = indid, fill = ind_lod_tot), hjust =0, vjust =0) +
  theme_classic() +
  labs(x = "Time of Year", y = NULL, fill = NULL) +
  scale_fill_gradientn(colours = c("#f2c0fe", "#440154"),limits=c(1, 37), na.value = "transparent") +
  scale_x_date(expand=c(0,0), limits = c(as.Date("2018-07-01"), as.Date("2018-11-01")), 
               date_breaks = "1 month", date_labels = "%b %d") +
  scale_y_discrete(expand=c(0,0), limits = rev) + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle=50, hjust = 0.9),
        panel.background = element_rect(fill = NA, color = 'black'),
        legend.position = "none",
        plot.margin = unit(c(1, 1, 1, 1), "pt")) +
  geom_hline(yintercept = 0 + 0:5, colour = "black") +
  geom_vline(xintercept = as.Date("2018-08-01"), color = "black", linetype = 2)+
  geom_vline(xintercept = as.Date("2018-09-01"), color = "black", linetype = 2)+
  geom_vline(xintercept = as.Date("2018-10-01"), color = "black", linetype = 2)

################################################################################
# Fig. 5. Factors influencing the risk of influenza virus infection acquisition 
# by subtype/lineage based on multivariable logistic regression. 

# Data provided are model-estimates. The model input includes individual-level 
# data that cannot be shared. The regression that produces these estimates can be
# simulated with pseudo data in transmission_model.R; this file also includes
# the code for data manipulation of the model outputs for figure generation.

fig5_dat <- read.csv('Fig5_household_transmission_model_estimates.csv')

fig5_dat <- fig5_dat %>% filter(term != '(Intercept)')

fig5_dat$term <- factor(fig5_dat$term, 
                          levels = c("comsmooth", "hh_load", "titercutofftiter>=40", 
                                     "age_cat<5", "age_cat5-11", 
                                     "age_cat12-18", "age_cat19-40", "sexMale", "cleanhivPositive CD4 <200",
                                     "cleanhivPositive CD4 >=200", "cleanhivUnknown", 
                                     "hh_s"))

ref <- data.frame(effect = rep('ref', 16), 
                  group = c(rep('A(H1N1)pdm09', 4), rep('A(H3N2)', 4),
                            rep('B/Yamagata', 4), rep('B/Victoria', 4)),
                  term = c(rep(c('titercutofftiter<40', 'age>40', 'sexfemale', 'hivneg'), 4)),
                  estimate_log = 1, std.error = 0, statistic = 0, df = 0, p.value = 0,
                  lower = 0, upper = 0, estlabel = '')

order <- c("hh_s",
           "hivneg","cleanhivUnknown", "cleanhivPositive CD4 >=200",
           "cleanhivPositive CD4 <200","sexfemale","sexMale","age>40",
           "age_cat19-40","age_cat12-18",  "age_cat5-11", "age_cat<5", 
           'titercutofftiter<40', "titercutofftiter>=40", "hh_load", "comsmooth")

fig5_plot <- ggplot(fig5_dat, aes(x = term, y = estimate_log, ymin = lower, ymax = upper, 
                                        label = estlabel)) +
  scale_x_discrete(limits = order, labels = c("Household Size", "HIV Negative",
                                              "HIV Status Unknown", "PLWH CD4 >=200", 
                                              "PLWH CD4 <200", "Sex - Female", 
                                              "Sex - Male", "Age >40",
                                              "Age 19-40", "Age 12-18",
                                              "Age 5-11", "Age <5", "Pre-Season HAI Titer <40", 
                                              "Pre-Season HAI Titer >=40",
                                              "Household FOI", "Community FOI")) +
  scale_y_continuous(trans = 'log2') + 
  coord_flip(ylim = c(0.04, 16)) + 
  annotate('rect', xmin = 0, xmax = 1.5, ymin = 0, ymax = 16, fill = "#0d0887", alpha = 0.5) +
  annotate('rect', xmin = 1.5, xmax = 5.5, ymin = 0, ymax = 16, fill = "#5c01a6", alpha = 0.5) +
  annotate('rect', xmin = 5.5, xmax = 7.5, ymin = 0, ymax = 16, fill = "#9c179e", alpha = 0.5) +
  annotate('rect', xmin = 7.5, xmax = 12.5, ymin = 0, ymax = 16,fill = "#cc4778", alpha = 0.5) +
  annotate('rect', xmin = 12.5, xmax = 14.5, ymin = 0, ymax = 16, fill = "#ed7953", alpha = 0.5) +
  annotate('rect', xmin = 14.5, xmax = 15.5, ymin = 0, ymax = 16, fill = "#fdb42f", alpha = 0.5) +
  annotate('rect', xmin = 15.5, xmax = 16.5, ymin = 0, ymax = 16, fill = "#f0f921", alpha = 0.5) +
  geom_hline(yintercept = 1, color = 'grey', alpha = 0.8) +
  geom_linerange(linewidth = 0.3) + geom_point(size = 1.5) +
  geom_text(size = 3, nudge_x = .4, color = "black") +
  geom_point(data = ref, aes (x = term, y = estimate_log), color = 'black', fill = 'white', pch = 21) +
  labs(y = "Hazard Ratio", x = "") +
  facet_wrap(~group, nrow = 1) + 
  theme_classic() +
  theme(legend.position="none",
        strip.background = element_blank(), panel.spacing.x = unit(1.5, "lines"))

fig5_plot






