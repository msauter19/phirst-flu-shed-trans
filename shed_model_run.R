
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

# Add in subtype/lineage 
shed_data$TYPE <- NA
shed_data$TYPE <- ifelse(!is.na(shed_data$npsh3ct), "AH3", shed_data$npsh3ct)

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

################################################################################
# Running model to estimate shedding kinetics. Run separately for each 
# subtype/lineage. 

type <- "AH3" # BYam AH1 BVic
mod_dat$CT <- mod_dat$npsh3ct
    
indid_l <- c(unique(mod_dat$indid_inf[which(mod_dat$TYPE == type)]))
subset <- mod_dat %>%
  filter(indid_inf %in% indid_l)
person_dat <- one_dat %>% 
  filter(indid_inf %in% indid_l)

# must change clean_id to consecutive list of individuals
id_list <- sort(unique(mod_dat$indid_inf))
id_df <- data.frame(indid_inf = id_list, id_clean = 1:length(id_list), stringsAsFactors=FALSE)
subset <- left_join(subset, id_df, by="indid_inf")
person_dat <- left_join(person_dat, id_df, by = "indid_inf") %>%
  dplyr::select(id_clean, everything())
    
rm(id_df, id_list, indid_l)
    
# parameters differed for type A and type B as indicated below #############
if (type == "AH3" | type == "AH1") {
  # for influenza a
  pars <- list(N = nrow(subset), 
               n_id = length(unique(subset$id_clean)),
               id = subset$id_clean,
               t = subset$date, 
               ct = subset$CT,
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
               sigma_prior_scale = as.numeric(5)
  )
}

if (type == "BVic" | type == "BYam") {
  # for influenza b
  pars <- list(N = nrow(subset), 
               n_id = length(unique(subset$id_clean)),
               id = subset$id_clean,
               t = subset$date, 
               ct = subset$CT,
               lod = as.numeric(37),
               tpsd_p = as.numeric(2),
               dctpmean_p = as.numeric(37/2),
               dctpsd_p = as.numeric(37/6),
               wpmax = as.numeric(15),
               wpmean_p = as.numeric(2),
               wpsd_p = as.numeric(1), 
               wcmax = as.numeric(40),
               wcmean_p = as.numeric(5),
               wcsd_p = as.numeric(3),
               sigma_prior_scale = as.numeric(5)
  )
} 
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
                     sigma_prior_scale = as.list(pars)$sigma_prior_scale), 
                   iter=1000, chains=4) #, control = list(adapt_delta = 0.9))
fit_endq <- Sys.time()
#print(paste0("Fit time: ",difftime(fit_endq, fit_startq, units="min")," mins"))
rm(fit_endq, fit_startq)
    
# saving data ##################################################################
c <- get_sampler_params(ct_fit)
params <- rstan::extract(ct_fit)
#summary(do.call(rbind, c), digits = 2)
#print(ct_fit, pars = c("dctpmean", "dctpsd", "wpmean", "wpsd", "wcmean", "wcsd"))
#print(ct_fit, pars = c("tp"))
#print(ct_fit, pars = c("x"))
#print(ct_fit, pars = c("wp"))
#print(ct_fit, pars = c("wc"))
#print(ct_fit, pars = c("mu"))
 
make_indiv_params_df <- function(extracted_params, parnames, n_indiv){
  out <- reduce(lapply(parnames, function(x) parseparam(extracted_params,x,n_indiv)), cbind) %>% 
    as_tibble %>% 
    mutate(iteration=1:n()) %>% 
    pivot_longer(-iteration) %>% 
    separate(name, c("param","id"), sep="_") %>%
    pivot_wider(id_cols = c("id","iteration"), names_from="param", values_from="value")
}
params_df <- make_indiv_params_df(params, c("tp","x","wp","wc"), pars$n_id) %>%
  rename(id_clean = id)
person_dat$id_clean <- as.character(person_dat$id_clean)
person_dat <- person_dat %>% select(id_clean, indid, indid_inf)
params_df <- params_df %>%
  left_join(person_dat, by="id_clean") 

add <- one_dat %>% select(indid_inf, center_date)
params_df <- params_df %>% 
  left_join(add, by = "indid_inf")

params_df$center_date = as.Date(params_df$center_date, format = "%m/%d/%Y")
params_df$tp_trans = params_df$center_date + round(params_df$tp, digits = 0)

rm(add, person_dat, shed_data, type, subset, params, mod_dat, one_dat, pars, 
   c, c_mod, ct_model, ct_fit)
