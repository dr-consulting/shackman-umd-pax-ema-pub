####################################################################################################
# Study 2 Modeling Script: Negative Mood Models

# Description: 
#   The analyses below involve a series of increasingly complex Bayesian multilevel regression 
#   models. The analyses addressed three broad research questions designed to provide a better 
#   understanding of the association between dispositional negativity and momentary 
#   negative affect:
#     1. What is a reasonable estimate of the tonic or "unique" association between dispositional 
#     negativity and momentary negative affect? 
#     2. What is a reasonable estimate of the association between dispositional negativity and 
#     momentary negative affect that can be attributed to differences in overall emotional context? 
#     3. What is a reasonable estimate of the association between dispositional negativity and 
#     momentary negative affect that can be attributed to reactivity to recent emotionally salient 
#     events?

# Modeling Notes: 
#   * Exploratory analyses revealed that negative mood ratings were approximately symmetrical in 
#   their distribution, thus a weakly informative normal prior was chosen for negative mood scores
#   * Missingness was addressed via multiple imputation (see imputation script elsewhere)
#   * To generate a more informative posterior predictive distribution in the missingness models, 
#   we included summary scores of the EMA
#   * NEG is mainly an aggregate of negative mood and combines high energy (anxious) and low 
#   energy (depressed) affective states. A multilevel factor analysis supported the combination of 
#   these two dimensions of negative momentary mood in a single composite. 
#   * Event ratings were individually mean-centered to maintain separation of between- and within-
#   subjects sources of variation in negative mood
####################################################################################################

#---------------------------------------------------------------------------------------------------
library(brms)
library(rstan)
library(rstanarm)
library(bayesplot)
library(tidyverse)
library(glue)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
#---------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------
# Location of repo stored locally
wd<-'~/github/ATNL/shackman-umd-pax-ema-pub'
data.folder <- '{wd}/Data' %>% glue()
study2.model <- '{wd}/Study_2_model_summaries' %>% glue()

# Will save very large posterior files from analyses (not recommended for git repo)
# For anyone attempting to reproduce these analyses be sure to identify a storage location with sufficient memory
posterior_save_dir <- "/media/dr-owner/HDD1"
study2.out <- '{posterior_save_dir}/EMA_S2_Bayesian_Posteriors/gamma' %>% glue()

# Also generally not recommended to store image files on GH... 
study2.graphics <- '{study2.out}/diagnostic_plots' %>% glue()
#--------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------
# Loading data from study 2 - PAX
load(paste0(data.folder, '/study2_data.RData'))
source('{wd}/utils.R' %>%  glue())

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
prior_config <- c(
    set_prior(
        'lognormal(2, .5)', 
        class = 'shape'
    ), 
    set_prior(
        'lognormal(-1.5, .5)', 
        class = 'sd'
    ), 
    set_prior(
        'lognormal(-.75, .55)', 
        class = "Intercept"
    )
)

# Add slopes to the priors
prior_config <- c(
    prior_config, 
    set_prior(
        'normal(0, 2)', 
        class = 'b'
    )
)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# A repeat of the model above, but with positive events 
S2_NEG_PosEvnt_DN_form <- bf(
    NEG ~ 1 + c.PosEvnt + c.DN + (1 + c.PosEvnt|ID)
)+Gamma(link = "log")

# Running model with priors (see above)
S2_NEG_PosEvnt_DN <- brm_multiple(S2_NEG_PosEvnt_DN_form,
                                  data = dat.study2_list, 
                                  chains = 3,
                                  cores = 3,
                                  iter = 15000,
                                  warmup = 10000, 
                                  control = list(adapt_delta = .99, 
                                                 max_treedepth = 15), 
                                  seed = 19610412, 
                                  save_all_pars = TRUE, 
                                  save_model = "S2_NEG_PosEvnt_DN", 
                                  open_progress = FALSE,
                                  refresh = 0, 
                                  prior = prior_config)

create_diagnostic_plots(S2_NEG_PosEvnt_DN, dat.study2_list, n_samples = 100, dir_path = study2.graphics)

sink(paste0(study2.model, '/S2_NEG_PosEvnt_DN.txt'))
print(summary(S2_NEG_PosEvnt_DN), digits = 5)
sink()

# Given size of files saved off posteriors to external hard drive
save(list = c('S2_NEG_PosEvnt_DN', 'S2_NEG_PosEvnt_DN_form'), 
     file = paste0(study2.out, "/S2_NEG_PosEvnt_DN.RData"))
remove(S2_NEG_PosEvnt_DN)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Swapping proportion of negative events for DN in this "two-level" model. Used to estimate unique
# and shared variance accounted for by proportion of negative events and DN in momentary mood.
S2_NEG_NegEvnt_prop.NegEvnt_form <- bf(
    NEG ~ 1 + c.NegEvnt + prop.NegEvnt + (1 + c.NegEvnt|ID)
)+Gamma(link = "log")

# Running model with priors (see above)
S2_NEG_NegEvnt_prop.NegEvnt <- brm_multiple(S2_NEG_NegEvnt_prop.NegEvnt_form,
                                            data = dat.study2_list, 
                                            chains = 3,
                                            cores = 3,
                                            iter = 15000,
                                            warmup = 10000, 
                                            control = list(adapt_delta = .99, 
                                                           max_treedepth = 15), 
                                            seed = 19610412, 
                                            save_all_pars = TRUE, 
                                            save_model = "S2_NEG_NegEvnt_prop.NegEvnt", 
                                            open_progress = FALSE,
                                            refresh = 0, 
                                            prior = prior_config)

create_diagnostic_plots(S2_NEG_NegEvnt_prop.NegEvnt, dat.study2_list, n_samples = 100, dir_path = study2.graphics)

sink(paste0(study2.model, '/S2_NEG_NegEvnt_prop.NegEvnt.txt'))
print(summary(S2_NEG_NegEvnt_prop.NegEvnt), digits = 5)
sink()

# Given size of files saved off posteriors to external hard drive
save(list = c('S2_NEG_NegEvnt_prop.NegEvnt', 'S2_NEG_NegEvnt_prop.NegEvnt_form'), 
     file = paste0(study2.out, "/S2_NEG_NegEvnt_prop.NegEvnt.RData"))
remove(S2_NEG_NegEvnt_prop.NegEvnt)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Final model - with reactivity effect included (in theory accounting for level-1 variance)
S2_NEG_NegEvnt_x_DN_prop.NegEvnt_form <- bf(
    NEG ~ 1 + c.NegEvnt * c.DN + prop.NegEvnt + (1 + c.NegEvnt|ID)
)+Gamma(link = "log")

# Running model with priors (see above)
S2_NEG_NegEvnt_x_DN_prop.NegEvnt <- brm_multiple(S2_NEG_NegEvnt_x_DN_prop.NegEvnt_form,
                                                 data = dat.study2_list, 
                                                 chains = 3,
                                                 cores = 3,
                                                 iter = 15000,
                                                 warmup = 10000, 
                                                 control = list(adapt_delta = .99, 
                                                                max_treedepth = 15), 
                                                 seed = 19610412, 
                                                 save_all_pars = TRUE, 
                                                 save_model = "S2_NEG_NegEvnt_x_DN_prop.NegEvnt", 
                                                 open_progress = FALSE,
                                                 refresh = 0, 
                                                 prior = prior_config)

create_diagnostic_plots(S2_NEG_NegEvnt_x_DN_prop.NegEvnt, dat.study2_list, n_samples = 100, dir_path = study2.graphics)

sink(paste0(study2.model, '/S2_NEG_NegEvnt_x_DN_prop.NegEvnt.txt'))
print(summary(S2_NEG_NegEvnt_x_DN_prop.NegEvnt), digits = 5)
sink()

# Given size of files saved off posteriors to external hard drive
save(list = c('S2_NEG_NegEvnt_x_DN_prop.NegEvnt', 'S2_NEG_NegEvnt_x_DN_prop.NegEvnt_form'), 
     file = paste0(study2.out, "/S2_NEG_NegEvnt_x_DN_prop.NegEvnt.RData"))
remove(S2_NEG_NegEvnt_DN_prop.NegEvnt)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Final model - with reactivity effect included (in theory accounting for level-1 variance)
S2_NEG_PosEvnt_x_DN_prop.PosEvnt_form <- bf(
    NEG ~ 1 + c.PosEvnt * c.DN + prop.PosEvnt + (1 + c.PosEvnt|ID)
)+Gamma(link = "log")

# Running model with priors (see above)
S2_NEG_PosEvnt_x_DN_prop.PosEvnt <- brm_multiple(S2_NEG_PosEvnt_x_DN_prop.PosEvnt_form,
                                                 data = dat.study2_list, 
                                                 chains = 3,
                                                 cores = 3,
                                                 iter = 15000,
                                                 warmup = 10000, 
                                                 control = list(adapt_delta = .99, 
                                                                max_treedepth = 15), 
                                                 seed = 19610412, 
                                                 save_all_pars = TRUE, 
                                                 save_model = "S2_NEG_PosEvnt_x_DN_prop.PosEvnt", 
                                                 open_progress = FALSE,
                                                 refresh = 0, 
                                                 prior = prior_config)

create_diagnostic_plots(S2_NEG_PosEvnt_x_DN_prop.PosEvnt, dat.study2_list, n_samples = 100, dir_path = study2.graphics)

sink(paste0(study2.model, '/S2_NEG_PosEvnt_x_DN_prop.PosEvnt.txt'))
print(summary(S2_NEG_PosEvnt_x_DN_prop.PosEvnt), digits = 5)
sink()

# Given size of files saved off posteriors to external hard drive
save(list = c('S2_NEG_PosEvnt_x_DN_prop.PosEvnt', 'S2_NEG_PosEvnt_x_DN_prop.PosEvnt_form'), 
     file = paste0(study2.out, "/S2_NEG_PosEvnt_x_DN_prop.PosEvnt.RData"))
remove(S2_NEG_PosEvnt_x_DN_prop.PosEvnt)
gc()
