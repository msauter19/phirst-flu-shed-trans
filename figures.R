
# Code to produce the figures included in the main text. All figures were 
# replicated for each respective subtype/lineage observed in the paper. 
# For Fig. 3 see shed_regression.R and for Fig. 5 see transmission_model.R

source('shed_model_run.R')
library(lares)
library(gridExtra)
library(grid)


# Fig. 1 Distribution of the pre-season serum influenza antibody log base 
# 10 scale by the age of the participant at enrollment. 

titersdist_plot <- ggplot(data = shed_data, aes(x = age, y = log(flua_h3n2, base = 10))) +
  geom_point(color = "white", fill = "#28607f", alpha = 0.4, pch = 21) + 
  geom_smooth(color = "white", fill = "#28607f", alpha = 1,
              method="lm", formula=  y ~ splines::bs(x, knots = c(5, 12, 18, 40), degree = 3)) +
  scale_y_continuous(breaks = c(1,2,3), labels = c("1" = "10", "2" = "100", "3" = "1000")) +
  scale_x_continuous(breaks = c(0, 25, 50, 75), labels = c("0" = "0", "25" = "25", "50" = "50", "75" = "75+" )) +
  coord_cartesian(xlim = c(0,75), ylim = c(0.5,3.5)) +
  theme_light() +
  labs(x = "Age at Enrollment (years)", y = "HAI titer", title = "Pre-season A(H3N2) HAI")
print(titersdist_plot)

################################################################################
# Fig. 2 Influenza viral shedding trajectory estimates.
individual <- distinct(shed_data, indid_inf, .keep_all = TRUE) %>%
  filter(!(is.na(meanwc)))
raw <- shed_data %>% 
  filter(!(is.na(meanwc))) %>%
  mutate(DATE = date - center_point) 
pmean <- -mean(individual$meanwp)
ctmean <- mean(individual$meanct)
cmean <- mean(individual$meanwc)
means_df <- data.frame(x = c(pmean, 0, cmean), 
                       y = c(37, ctmean, 37))
ci_df <- data.frame(x = c(-ci_var(individual, meanwp)$upper_ci, 
                          0, ci_var(individual, meanwc)$upper_ci, 
                          ci_var(individual, meanwc)$lower_ci, 0, 
                          -ci_var(individual, meanwp)$lower_ci), 
                    y = c(37, ci_var(individual, meanct)$lower_ci, 37, 37,
                          ci_var(individual, meanct)$upper_ci, 37))
iqr_df <- data.frame(x = c(-quantile(individual$meanwp, 0.75, na.rm = T), 
                           0, quantile(individual$meanwc, 0.75, na.rm = T), 
                           quantile(individual$meanwc, 0.25, na.rm = T), 0, 
                           -quantile(individual$meanwp, 0.25, na.rm = T)), 
                     y = c(37, quantile(individual$meanct, 0.25, na.rm = T), 37, 37,
                           quantile(individual$meanct, 0.75, na.rm = T), 37))

kin_plot <- ggplot(data = raw, aes(x = DATE, y = npsh3ct)) + 
  geom_point(color = "white", fill = "black", alpha = 0.4, pch = 21) +
  geom_polygon(data = ci_df, aes(x = x, y = y), fill = "#28607f") +
  geom_polygon(data = iqr_df, aes(x = x, y = y), alpha = 0.5, fill = "#28607f") +
  geom_line(data = means_df, aes(x = x, y = y), linewidth = 0.7, color = "white") +
  scale_y_reverse() +
  coord_cartesian(ylim = c(37,15), xlim = c(-5,16), clip = "on") +
  geom_vline(aes(xintercept = 0), color = "black", alpha = 0.6) +
  labs(x = "Time since minimum Ct (days)", y = "Cycle threshold (Ct)") + 
  theme_bw() +
  theme(legend.position="none", strip.background = element_blank())
print(kin_plot)

pbox <- ggplot(data = individual, aes(y = log(meanwp, base = 2))) + 
  geom_boxplot(fill = "#28607f") +
  theme_classic() + 
  theme(axis.ticks = element_blank()) +
  labs(title = "Proliferation Duration", y = "Duration (days)", 
       x = NULL) +
  guides(fill="none") +
  scale_y_continuous(limits = c(0, 5), breaks = c(0, 1, 2, 3, 4, 5), 
                     labels=c("1", "2", "4", "8", "16", "32"))
print(pbox)

cbox <- ggplot(data = individual, aes(y = log(meanwc, base = 2))) + 
  geom_boxplot(fill = "#28607f") + 
  theme_classic() + 
  theme(axis.ticks = element_blank()) + 
  labs(title = "Clearance Duration", y = element_blank(), 
       x = NULL)+
  guides(fill="none") +
  scale_y_continuous(limits = c(0, 5), breaks = c(0, 1, 2, 3, 4, 5), 
                     labels=c("1", "2", "4", "8", "16", "32"))
print(cbox)

tbox <- ggplot(data = individual, aes(y = log(tot_duration, base = 2))) + 
  geom_boxplot(fill = "#28607f") +
  theme_classic() + 
  theme(axis.ticks = element_blank()) +
  labs(title = "Total Duration", y = element_blank(), 
       x = NULL)+
  guides(fill="none") +
  scale_y_continuous(limits = c(0, 5), breaks = c(0, 1, 2, 3, 4, 5), 
                     labels=c("1", "2", "4", "8", "16", "32"))
print(tbox)

ctbox <- ggplot(data = individual, aes(y = meanct)) + 
  geom_boxplot(fill = "#28607f") + 
  theme_classic() + 
  theme(axis.ticks = element_blank()) +
  scale_y_reverse() +
  labs(title = "Minimum Ct", y = "Minimum Ct", 
       x = NULL)+
  guides(fill="none")
print(ctbox)

################################################################################
# Fig. 4 Schematic representation of the impact of two key components of the 
# household influenza virus transmission model: community prevalence and 
# household viral load on the individual infection risk. 

# Need data in matrix form, similar to transmission model
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
indid_one <- unique(c(inf_dat$indid))
noinf <- setdiff(unique(c(matrix$indid)), indid_one)
matrix$h3_present <- ifelse(matrix$indid %in% noinf, 0, matrix$h3_present)
matrix$prolif_load <- ifelse(matrix$indid %in% noinf, 0, matrix$prolif_load)
matrix$clear_load <- ifelse(matrix$indid %in% noinf, 0, matrix$clear_load)
noinf_df <- matrix %>% filter(indid %in% noinf)
matrix <- rbind(new, noinf_df) %>% 
  arrange(indid)
rm(i, s, sub, ind, indid_l, inf_dat, new, indid_one, noinf, noinf_df)
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
matrix <- matrix %>%
  mutate(hh_prof_load = hh_prof_load - prolif_load) %>%
  mutate(hh_clear_load = hh_clear_load - clear_load) %>% 
  mutate(ind_lod_tot = prolif_load + clear_load)
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
matrix$comsmooth <- rollapply(matrix$community, width = 4, function(...) {round(mean(...), digits = 3)}, partial = TRUE)
matrix$comsmooth <- matrix$comsmooth * 100

# now combine matrix with indicators of the start from the shedding estimates
starts <- shed_data %>% 
  filter(!is.na(round_start) & infcluster == 1) %>%
  dplyr::select(indid, round_start)

figure_data <- matrix %>% 
  left_join(starts, by = "indid", relationship = "many-to-many") %>% 
  mutate(hh_load = hh_prof_load + hh_clear_load) 
figure_data$ind_lod_tot <- ifelse(figure_data$ind_lod_tot == 0, NA, figure_data$ind_lod_tot) 
figure_data$round_start <- ifelse(is.na(figure_data$round_start), 500, figure_data$round_start)
figure_data$yeardate <- as.Date(as.Date("1/16/2017",format='%m/%d/%Y') + figure_data$date, format = '%m/%d/%Y')
figure_data$comsmooth <- ifelse(figure_data$community == 0, NA, figure_data$comsmooth)


# agincourt 
agfam <- figure_data %>% filter(hh == "A113")
agfam <- agfam %>% arrange(round_start, age)
agfam <- agfam %>%
  mutate(indid = factor(indid, unique(indid)))
agfam$hh_load <- ifelse(agfam$hh_load == 0, NA, agfam$hh_load)

agcom <- ggplot(agfam)+
  geom_raster(aes(x = as.Date(yeardate), y = hh, fill = comsmooth), hjust =0, vjust =0) +
  theme_classic() +  
  labs(x = element_blank(), y = element_blank(), fill = element_blank()) +
  scale_fill_gradientn(colours = c("#cfefd0", "#5ec962"),limits=c(0, 4), breaks =c(0,2), na.value = "transparent") +
  scale_x_date(expand=c(0,0), limits = c(as.Date("2017-01-17"), as.Date("2017-11-01")), 
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
  geom_vline(xintercept = as.numeric(as.Date("2017-03-01")), color = "black", linetype = 2)+
  geom_vline(xintercept = as.numeric(as.Date("2017-06-01")), color = "black", linetype = 2)+
  geom_vline(xintercept = as.numeric(as.Date("2017-09-01")), color = "black", linetype = 2)

aghh <- ggplot(agfam)+
  geom_raster(aes(x = as.Date(yeardate), y = hh, fill = hh_load), hjust =0, vjust =0) +
  theme_classic() +
  labs(x = element_blank(), y = element_blank(), fill = element_blank()) +
  scale_fill_gradientn(colours = c("#d2daec", "#3b528b"),limits=c(0, 55), breaks = c(0, 25, 50), na.value = "transparent") +
  scale_x_date(expand=c(0,0), limits = c(as.Date("2017-01-17"), as.Date("2017-11-01")), 
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
  geom_vline(xintercept = as.numeric(as.Date("2017-03-01")), color = "black", linetype = 2)+
  geom_vline(xintercept = as.numeric(as.Date("2017-06-01")), color = "black", linetype = 2)+
  geom_vline(xintercept = as.numeric(as.Date("2017-09-01")), color = "black", linetype = 2)

agind <- ggplot(agfam) +
  geom_raster(aes(x = as.Date(yeardate), y = indid, fill = ind_lod_tot), hjust =0, vjust =0) +
  theme_classic() +
  labs(x = "Time of Year", y = element_blank(), fill = element_blank()) +
  scale_fill_gradientn(colours = c("#f2c0fe", "#440154"),limits=c(1, 37), na.value = "transparent") +
  scale_x_date(expand=c(0,0), limits = c(as.Date("2017-01-17"), as.Date("2017-11-01")), 
               date_breaks = "1 month", date_labels = "%b %d") +
  scale_y_discrete(expand=c(0,0), limits = rev) + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle=50, hjust = 0.9),
        panel.background = element_rect(fill = NA, color = 'black'),
        plot.margin = unit(c(1,1,1,1), "pt"),
        legend.key.height = unit(0.4, 'cm')) +
  geom_hline(yintercept = 0 + 0:20, colour = "black") +
  geom_vline(xintercept = as.numeric(as.Date("2017-03-01")), color = "black", linetype = 2)+
  geom_vline(xintercept = as.numeric(as.Date("2017-06-01")), color = "black", linetype = 2)+
  geom_vline(xintercept = as.numeric(as.Date("2017-09-01")), color = "black", linetype = 2)

blank <- ggplot() + theme_minimal()
agname <- textGrob('Agincourt Site (Rural)', gp = gpar(fontsize = 12.5), hjust=0.7)
comname <- textGrob('Community 
Prevalence
(%)', gp = gpar(fontsize = 10), hjust=0.4, vjust = 0.3)
hhname <- textGrob('Total Household 
Viral Load
(Ct)', gp = gpar(fontsize = 10), hjust = 0.4)
indname <- textGrob('Individual 
Viral Load
In Shared 
Household
(Ct)', gp = gpar(fontsize = 10), vjust = 0, hjust=0.4)
charac1 <- textGrob('Individual 
HAI Titer   Age', gp = gpar(fontsize = 9), hjust = 0.5, vjust = 0.9)

# above code replicated for klerksdorp household in Fig 4 and arranged together

grid.arrange(
  arrangeGrob(blank, blank, charac1, blank, nrow=4, heights = c(0.5,1,1,4)),
  arrangeGrob(agname, agcom, aghh, agind, nrow = 4, heights = c(0.5,1,1,4)),
  arrangeGrob(blank, comname, hhname, indname, nrow = 4, heights = c(0.5,1,1,4)),
  ncol = 3, widths = c(1.3,4,1.8))






