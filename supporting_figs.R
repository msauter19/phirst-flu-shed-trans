
# Code to produce the figures included in the supporting information. 
# For Fig. 8 see shed_regression.R and for Fig. 9 see transmission_model.R 

source('shed_model_run.R')

# Supporting Information Fig. 2 A sample of the estimated viral shedding kinetic 
# fits from the shedding model. 

plot_ct_fit <- function(params_df, lod, indiv_data, ntraces){
  params_df %>% 
    sample_n(ntraces) %>% 
    ggplot() + 
    # Plot traces:
    geom_segment(aes(x=-Inf, y=lod, xend=tp-wp, yend=lod), alpha=0.3, col="#73D055FF") + 
    geom_segment(aes(x=tp-wp, y=lod, xend=tp, yend=lod-x), alpha=0.3, col="#73D055FF") + 
    geom_segment(aes(x=tp, y=lod-x, xend=tp+wc, yend=lod), alpha=0.3, col="#73D055FF") + 
    geom_segment(aes(x=tp+wc, y=lod, xend=Inf, yend=lod), alpha=0.3, col="#73D055FF") + 
    geom_segment(aes(x=-Inf, y=lod, xend=tpmean-wpmean, yend=lod), linewidth = 1, col="#2D708EFF") + 
    geom_segment(aes(x=tpmean-wpmean, y=lod, xend=tpmean, yend=lod-xmean), linewidth = 1, col="#2D708EFF") + 
    geom_segment(aes(x=tpmean, y=lod-xmean, xend=tpmean+wcmean, yend=lod), linewidth = 1, col="#2D708EFF") + 
    geom_segment(aes(x=tpmean+wcmean, y=lod, xend=Inf, yend=lod), linewidth = 1, col="#2D708EFF") + 
    # Plot data:
    geom_point(data=indiv_data, aes(x=date, y=npsh3ct), size=3) + 
    theme_classic() + 
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + 
    labs(x="Time (days)", y="Ct") + 
    scale_y_reverse() + 
    facet_wrap(~id_clean)
}

params_df <- `paramsdf_5-11` 
params_df$tpmean <- NA
params_df$xmean <- NA
params_df$wpmean <- NA 
params_df$wcmean <- NA

for(i in 1:length(unique(`paramsdf_5-11`$id_clean))) {
  params_df$tpmean[which(params_df$id_clean == i)] <- mean(params_df$tp[which(params_df$id_clean == i)])
  params_df$xmean[which(params_df$id_clean == i)] <- mean(params_df$x[which(params_df$id_clean == i)])
  params_df$wpmean[which(params_df$id_clean == i)] <- mean(params_df$wp[which(params_df$id_clean == i)])
  params_df$wcmean[which(params_df$id_clean == i)] <- mean(params_df$wc[which(params_df$id_clean == i)])
}

indid_l <- na.omit(c(unique(mod_dat$indid_inf[which(mod_dat$age_cat == "5-11")])))
subset <- mod_dat %>%
  filter(indid_inf %in% indid_l)
id_list <- sort(unique(subset$indid_inf))
id_df <- data.frame(indid_inf = id_list, id_clean = 1:length(id_list), stringsAsFactors=FALSE)
subset <- left_join(subset, id_df, by="indid_inf")

params_df$id_clean <- as.numeric(params_df$id_clean)
fitsplot <- plot_ct_fit(params_df, 37, subset, ntraces = 500)
print(fitsplot)


# Supporting Information Fig. 3 - A sample of model trace plots for estimation 
# of the peak timing of viral shedding (tp). Here for the example age group of 5-11.

tpsampling <- traceplot(`ct_fit_5-11`, pars = c("tp"), inc_warmup = FALSE)
print(tpsampling)

# Supporting Information Fig. 4-6 Mosaic plot of nasopharyngeal swabs collected. 
# In paper, code replicated for each year of study observation, and altered to 
# include all subtypes/lineages, as well as to separate between negative swabs 
# and missing swabs. The figures in the paper also were edited to demonstrate the 
# specific time period in which the pre-season serum samples were drawn. 

shed_data$yeardate <- as.Date(as.Date("1/16/2017",format='%m/%d/%Y') + 
                                ((shed_data$funum - 1) * 3.5), format = '%m/%d/%Y')
shed_data$TYPE <- NA
shed_data$TYPE <- ifelse(shed_data$npsh3 == 1, "AH3", shed_data$TYPE)
shed_data$TYPE <- ifelse(is.na(shed_data$npsh3), "Negative", shed_data$TYPE)

mosaic <- ggplot(shed_data, aes(x = yeardate, y = indid, fill = TYPE, alpha = npsh3ct)) +
  geom_tile() +
  theme_minimal() +
  labs(x = "Collection Date", y = "Individual", 
       title = "Mock Cohort") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_rect(fill = NA, color = NA),
        panel.grid.major = element_line(color = NA),
        panel.grid.minor = element_line(color = NA)) +
  scale_x_date(expand=c(0,0), date_breaks = "1 month", date_labels = "%b %d") +
  guides(alpha="none") +
  scale_alpha(range = c(1, .2),limits=c(11, 37)) + 
  scale_fill_manual(name = element_blank(), labels = c("A(H3N2)", "PCR Negative"),
                    values=c(AH3 = "#28607f", Negative = "#F0F0F0"))
print(mosaic)