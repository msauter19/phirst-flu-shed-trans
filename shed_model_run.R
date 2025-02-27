
# Taking the mock data set, clustering infections, and then inputting into 
# model to estimate the kinetics of shedding (model written in STAN). 

library(tidyverse)
library(dplyr)
library(rstan)
library(purrr)
library(shinystan)
library(bayesplot)
rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())

# Mock data set
shed_data <- read.csv("mock_dataset.csv")
shed_data$npsh3ct <- as.double(shed_data$npsh3ct)

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

# Clustering algorithm to pair together samples in an infection, with new infections 
# defined by the presence of a new subtype/lineage (not relevant for this mock data set) 
# or two weeks of negative samples since the previous positive samples. 
clus <- shed_data %>%
  filter(npsh3 == 1)
indid_l <- unique(c(clus$indid))
clusters <- c()

for (n in indid_l) {
  #subset data for that individual
  sub <- clus %>%
    filter(indid == n) 
  if (nrow(sub) == 1) {
    clusters <- c(clusters, 1)
  }
  if (nrow(sub) > 1) {
    sub <- sub %>%
      arrange(funum) %>%
      dplyr::select(date)
    sub_dist <- dist(sub, method = "euclidean")
    hc_sub <- hclust(sub_dist, method = "single")
    
    clusters <- c(clusters, cutree(hc_sub, h = 14))
  }
}

clus$infcluster <- clusters
clus <- clus %>% 
  # making new variable to combine indid and the cluster/infection number 
  mutate(indid_inf = paste(indid, "_I", infcluster, sep = ""))
neg <- shed_data %>%
  filter(is.na(npsh3)) %>%
  mutate(infcluster = NA) %>% 
  mutate(indid_inf = NA)
shed_data <- clus %>%
  rbind(neg) %>%
  arrange(indid, funum)

rm(clus, hc_sub, neg, sub, clusters, indid_l, n, sub_dist)

# Now make a data set to input to the model that only includes individuals that 
# had a positive samples and then creates a two week bubble of observations 
# around the positive samples of an infection for the model. 
infect <- c(unique(shed_data$indid[which(shed_data$infcluster == 1)]))
model_dat <- shed_data %>%
  filter(indid %in% infect)

# Center the data around a pseudo peak time, with two weeks of extra negative samples on each side of the infection. 
mod_dat <- data.frame()
model_dat$center_date <- NA
model_dat$npsh3ct <- as.double(model_dat$npsh3ct)
for(i in infect) {
  sub <- model_dat %>%
    filter(indid == i)
  tp <- c(sub$date[which(sub$npsh3ct == min(sub$npsh3ct, na.rm = TRUE))])
  sub$center_date <- sub$npsdatecol[which(sub$date == tp[1])]
  sub$date = sub$date - tp[1]
  sub <- sub %>%
    filter(date[which(date == min(date[which(!(is.na(npsh3ct)))]))] - date < 14 &
             date - date[which(date == max(date[which(!(is.na(npsh3ct)))]))] < 14)
  mod_dat <- rbind(mod_dat, sub)
}
rm(i, sub, tp, infect, model_dat)

# Make a new data frame that will only contain the first event of each of the infections
one_dat <- data.frame()
indid_inf <- unique(mod_dat$indid_inf[which(!is.na(mod_dat$indid_inf))])
for (n in indid_inf) {
  sub <- mod_dat %>%
    filter(indid_inf == n)
  one_dat <- rbind(one_dat, sub[1, ])
}
rm(indid_inf, n, sub)

# Relevant functions for data parsing from shedding model
make_indiv_params_df <- function(extracted_params, parnames, n_indiv){
  out <- reduce(lapply(parnames, function(x) parseparam(extracted_params,x,n_indiv)), cbind) %>% 
    as_tibble %>% 
    mutate(iteration=1:n()) %>% 
    pivot_longer(-iteration) %>% 
    separate(name, c("param","id"), sep="_") %>%
    pivot_wider(id_cols = c("id","iteration"), names_from="param", values_from="value")
}
parseparam <- function(extracted_params, parname, n_indiv){
  as_tibble(setNames(as.data.frame(extracted_params[[parname]]), makenames(parname,n_indiv)))
}
makenames <- function(parname, n_indiv){
  unlist(lapply(1:n_indiv, function(x) paste0(parname, "_", x)))
}
make_shared_params_df <- function(extracted_params, parnames){
  out <- reduce(lapply(parnames, function(x) 
    as_tibble(setNames(as.data.frame(extracted_params[[x]]),x))
  ), cbind) %>%
    as_tibble() %>%
    mutate(iteration=1:n())
}


################################################################################
# Running model to estimate shedding kinetics. Run separately for each 
# subtype/lineage and age group (only age group for this mock data set). 

#type_l <- c("AH1", "AH3", "BVic", "BYam")
age_l <- c("<5", "5-11", "12-18", "19-40", ">40")

for(age in age_l) {
    
    indid_l <- na.omit(c(unique(mod_dat$indid_inf[which(mod_dat$age_cat == age)])))
    subset <- mod_dat %>%
      filter(indid_inf %in% indid_l)
    person_dat <- one_dat %>% 
      filter(indid_inf %in% indid_l)
    
    # must change clean_id to consecutive list of individuals
    id_list <- sort(unique(subset$indid_inf))
    id_df <- data.frame(indid_inf = id_list, id_clean = 1:length(id_list), stringsAsFactors=FALSE)
    subset <- left_join(subset, id_df, by="indid_inf")
    person_dat <- left_join(person_dat, id_df, by = "indid_inf") %>%
      dplyr::select(id_clean, everything())
    
    rm(id_df, id_list, indid_l)
    
    # parameters differed for type A and type B as indicated below #############
    #if (type == "AH3" | type == "AH1") {
      # for influenza a
    pars <- list(N = nrow(subset), 
                          n_id = length(unique(subset$id_clean)),
                          id = subset$id_clean,
                          t = subset$date, 
                          ct = subset$npsh3ct,
                          lod = as.numeric(37),
                          tpsd_p = as.numeric(2),
                          dctpmean_p = as.numeric(37/2),
                          dctpsd_p = as.numeric(37/6),
                          wpmax = as.numeric(15),
                          wpmean_p = as.numeric(2),
                          wpsd_p = as.numeric(1), 
                          wcmax = as.numeric(30),
                          wcmean_p = as.numeric(4),
                          wcsd_p = as.numeric(3),
                          sigma = as.numeric(1)
      )
    #}
    #if (type == "BVic" | type == "BYam") {
      # for influenza b
      #pars <- list(N = nrow(subset), 
      #             n_id = length(unique(subset$id_clean)),
      #             id = subset$id_clean,
      #             t = subset$date, 
      #             ct = subset$CT,
      #             lod = as.numeric(37),
      #             tpsd_p = as.numeric(2),
      #             dctpmean_p = as.numeric(37/2),
      #             dctpsd_p = as.numeric(37/6),
      #             wpmax = as.numeric(15),
      #             wpmean_p = as.numeric(2),
      #             wpsd_p = as.numeric(1), 
      #             wcmax = as.numeric(40),
      #             wcmean_p = as.numeric(5),
      #             wcsd_p = as.numeric(3),
      #             sigma = as.numeric(1)
      #)
    #}
    
    # running model ################################################################
    
    c_mod <- stanc("flu_ct_model.stan")
    ct_model <- stan_model(stanc_ret = c_mod)
    
    fit_startq <- Sys.time()
    ct_fit <- sampling(ct_model, 
                       data = list(
                         N = as.list(pars)$N, 
                         n_id = as.list(pars)$n_id,
                         id = as.list(pars)$id,
                         t = as.list(pars)$t, 
                         ct = as.list(pars)$ct, 
                         lod = as.list(pars)$lod,
                         tpsd_p = as.list(pars)$tpsd_p,
                         dctpmean_p = as.list(pars)$dctpmean_p,
                         dctpsd_p = as.list(pars)$dctpsd_p,
                         wpmax = as.list(pars)$wpmax,
                         wpmean_p = as.list(pars)$wpmean_p,
                         wpsd_p = as.list(pars)$wpsd_p,
                         wcmax = as.list(pars)$wcmax,
                         wcmean_p = as.list(pars)$wcmean_p,
                         wcsd_p = as.list(pars)$wcsd_p,
                         sigma = as.list(pars)$sigma), 
                       iter=1000, chains=4) #, control = list(adapt_delta = 0.9))
    fit_endq <- Sys.time()
    #print(paste0("Fit time: ",difftime(fit_endq, fit_startq, units="min")," mins"))
    rm(fit_endq, fit_startq)
    
    # saving data ##################################################################
    c <- get_sampler_params(ct_fit)
    summary(do.call(rbind, c), digits = 2)
    
    #print(ct_fit, pars = c("dctpmean", "dctpsd", "wpmean", "wpsd", "wcmean", "wcsd"))
    #print(ct_fit, pars = c("tp"))
    #print(ct_fit, pars = c("x"))
    #print(ct_fit, pars = c("wp"))
    #print(ct_fit, pars = c("wc"))
    #print(ct_fit, pars = c("mu"))
    
    params <- rstan::extract(ct_fit)
    indiv_params_df <- make_indiv_params_df(params, c("tp","x","wp","wc"), pars$n_id) %>%
      rename(id_clean = id)
    shared_params_df <- make_shared_params_df(params, c("dctpmean","wpmean","wcmean","dctpsd","wpsd","wcsd"))
    params_df <- indiv_params_df %>% 
      left_join(shared_params_df, by="iteration")
    
    person_dat$id_clean <- as.character(person_dat$id_clean)
    params_df <- params_df %>%
      left_join(person_dat, by="id_clean") 
    params_df$center_date = as.Date(params_df$center_date, format = "%m/%d/%Y")
    params_df$tp_trans = params_df$center_date + round(params_df$tp, digits = 0)
    
    
    ct_name <- paste("ct_fit", age, sep = "_")
    assign(ct_name, ct_fit)
    df_name <- paste("paramsdf", age, sep = "_")
    assign(df_name, params_df)
    
    rm(ct_fit, subset, paramsdf, indiv_params_df, shared_params_df, params, 
       params_df, person_dat, pars, c, c_mod, ct_model, ct_name, df_name)
}

################################################################################
# Now adding the new shedding estimates into data set for use of regressions and figures

shedding <- `paramsdf_<5` %>%
  rbind(`paramsdf_5-11`) %>%
  rbind(`paramsdf_12-18`) %>%
  rbind(`paramsdf_19-40`) %>%
  rbind(`paramsdf_>40`) %>%
  select(iteration, indid, indid_inf, tp, x, wp, wc)

indid_l <- c(unique(shedding$indid_inf))
means <- data.frame()
for(i in indid_l) {
  subset <- shedding %>% 
    filter(indid_inf == i) %>%
    mutate(meantp = mean(tp), 
           meanct = mean(37 - x), 
           meanwp = mean(wp),
           meanwc = mean(wc)) %>%
    filter(iteration == 1)
  means <- rbind(means, subset)
}
rm(subset, i, indid_l) 

id_list <- c(unique(shed_data$indid_inf))
id_list <- id_list[!is.na(id_list)]
date_dat <- data.frame()
for(i in id_list) {
  sub <- shed_data %>%
    filter(indid_inf == i) %>% 
    dplyr::select(indid, indid_inf, date, npsdatecol, npsh3ct)
  tp <- c(sub$date[which(sub$npsh3ct == max(sub$npsh3ct, na.rm = TRUE))])
  sub$center_date <- sub$date[which(sub$date == tp[1])]
  date_dat <- rbind(date_dat, sub[1,])
}
rm(id_list, sub, i, tp)

date_dat <- date_dat %>% 
  dplyr::select(indid_inf, center_date)
means <- means %>%
  left_join(date_dat, by = "indid_inf")

means <- means %>% 
  mutate(center_point = center_date + meantp) %>%
  mutate(start_point = center_point - meanwp) %>%
  mutate(end_point = center_point + meanwc) %>% 
  mutate(round_start = round(start_point)) %>% 
  mutate(tot_duration = meanwp + meanwc) %>%
  arrange(indid)

# calculating proliferation/clearance slope to use later to get day by day estimate of viral load
means <- means %>%
  mutate(prolif_slope = meanct / meanwp) %>%
  mutate(clear_slope = -meanct / meanwc) %>% 
  dplyr::select(indid_inf, center_date, center_point, start_point, end_point, round_start, 
                meanct, prolif_slope, clear_slope, meanwp, meanwc, tot_duration)

# combine with master data
shed_data <- shed_data %>% 
  left_join(means, by = "indid_inf")

rm(means, shedding, indid_l, date_dat, make_indiv_params_df, make_shared_params_df,
   makenames, parseparam, age, age_l)