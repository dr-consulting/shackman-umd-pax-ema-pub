###################################################################################################
# Study 1 Modeling Script

# Negative and Positive event MLMs
# Description: 
#   The analyses below examine a bivariate association between individual dispositional negativity 
#   scores and ratings of recent negative and positive events. Participants provided event ratings 
#   via repeated ecological momentary assessments completed on their smartphones. Dispositional 
#   negativity scores were created by aggregating the same set of items tapping trait anxiety and 
#   negative affect completed at two separate timepoints and in two different contexts. 

# Modeling Notes: 
#   * In this analysis, negative and positive event ratings have been transformed from a 5-point 
#   to a dichotomous scale. The value 1 represented "no" an event did not occur. All other values
#   were set to "yes" an event did occur. 
#   * Missingness was addressed using 10 imputed data sets. Script for imputation generation is
#   available in the repo
###################################################################################################

#--------------------------------------------------------------------------------------------------
# Loading relevant packages
library(brms)
library(rstan)
library(bayesplot)
library(tidyverse)
library(glue)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
#--------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------
# Location of repo stored locally
wd <- '~/github/ATNL/shackman-umd-pax-ema-pub'
data.folder <- '{wd}/Data' %>% glue()
study1.model <- '{wd}/Study_1_model_summaries' %>% glue()

# Will save very large posterior files from analyses (not recommended for git repo)
# For anyone attempting to reproduce these analyses be sure to identify a storage location with sufficient memory
posterior_save_dir <- "/media/dr-owner/HDD1"
study1.out <- '{posterior_save_dir}/EMA_S1_Bayesian_Posteriors' %>% glue()

# Also generally not recommended to store image files on GH... 
study1.graphics <- '{study1.out}/diagnostic_plots' %>% glue()
#--------------------------------------------------------------------------------------------------
# source local utility functions
source('{wd}/utils.R' %>%  glue())
#--------------------------------------------------------------------------------------------------
load('{data.folder}/study1_data_probationary.RData' %>% glue())

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-
# Simple Event Models:

# Substantive Questions: 
#   1. Is dispositional negativity associated with occurrence of recent negative events? 
#   2. Is dispositional negativity associated with occurrence of recent positive events? 
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-
#Negative Event Model: 
NegEvnt_DN <- bf(
  NegEvnt_dich ~ 1 + c.DN + (1|ID) 
) + bernoulli()

#Running model with priors (see above)
S1_NegEvnt_DN_dich <- brm_multiple(
  NegEvnt_DN,
  data = dat.study1_list, 
  chains = 3,
  iter = 15000,
  warmup = 10000, 
  control = list(adapt_delta = .99, 
                 max_treedepth = 15), 
  seed = 7201969, 
  save_all_pars = TRUE, 
  save_model = "S1_NegEvnt_DN_dich", 
  open_progress = FALSE,
  refresh = 0
)

save(list = c("S1_NegEvnt_DN_dich", "dat.study1_list", "dat.study1_model"), 
     file = '{study1.out}/S1_NegEvnt_DN_dich.RData' %>% glue()) 

sink('{study1.model}/S1_NegEvnt_DN_dich.txt' %>% glue())
print(summary(S1_NegEvnt_DN_dich), digits = 5)
sink()

create_diagnostic_plots(S1_NegEvnt_DN_dich, dat.study1_list, n_samples = 100, dir_path = study1.graphics)

# cleanup 
remove('S1_NegEvnt_DN_dich')
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-
#Positive Event Model: 
PosEvnt_DN <- bf(
  PosEvnt_dich ~ 1 + c.DN + (1|ID) 
) + bernoulli()

#Running model with priors (see above)
S1_PosEvnt_DN_dich <- brm_multiple(
  PosEvnt_DN,
  data = dat.study1_list, 
  chains = 3,
  iter = 15000,
  warmup = 10000, 
  control = list(adapt_delta = .99, 
                 max_treedepth = 15), 
  seed = 7201969, 
  save_all_pars = TRUE, 
  save_model = "S1_PosEvnt_DN_dich", 
  open_progress = FALSE,
  refresh = 0
)

save(list = c("S1_PosEvnt_DN_dich", "dat.study1_list", "dat.study1_model"), 
     file = '{study1.out}/S1_PosEvnt_DN_dich.RData' %>% glue()) 

sink('{study1.model}/S1_PosEvnt_DN_dich.txt' %>%  glue())
print(summary(S1_PosEvnt_DN_dich), digits = 5)
sink()

create_diagnostic_plots(S1_PosEvnt_DN_dich, dat.study1_list, n_samples = 100, dir_path = study1.graphics)

# cleanup 
remove('S1_PosEvnt_DN_dich')
gc()
