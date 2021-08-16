####################################################################################################
# Study 2 Modeling Script: Positive Mood Models

# Description: 
#   The analyses below involve a series of increasingly complex Bayesian multilevel regression 
#   models. The analyses addressed three broad research questions designed to provide a better 
#   understanding of the association between dispositional negativity and momentary 
#   positive affect:
#     1. What is a reasonable estimate of the tonic or "unique" association between dispositional 
#     negativity and momentary positive affect? 
#     2. What is a reasonable estimate of the association between dispositional negativity and 
#     momentary positive affect that can be attributed to differences in overall emotional context? 
#     3. What is a reasonable estimate of the association between dispositional negativity and 
#     momentary positive affect that can be attributed to reactivity to recent emotionally salient 
#     events?

# Modeling Notes: 
#   * Exploratory analyses revealed that positive mood ratings were approximately symmetrical in 
#   their distribution, thus a weakly informative normal prior was chosen for positive mood scores
#   * Missingness was addressed via multiple imputation (see imputation script elsewhere)
#   * To generate a more informative posterior predictive distribution in the missingness models, 
#   we included summary scores of the EMA
#   * POS is mainly an aggregate of positive mood and combines high energy (cheerful) and low 
#   energy (calm) affective states. A multilevel factor analysis supported the combination of these
#   two dimensions of positive momentary mood in a single composite. 
#   * Event ratings were individually mean-centered to maintain separation of between- and within-
#   subjects sources of variation in positive mood
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
study2.out <- '{posterior_save_dir}/EMA_S2_Bayesian_Posteriors/gaussian' %>% glue()

# Also generally not recommended to store image files on GH... 
study2.graphics <- '{study2.out}/diagnostic_plots' %>% glue()
#--------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------
# Loading data from study 2 - PAX
load(paste0(data.folder, '/study2_data.RData'))
source('{wd}/utils.R' %>%  glue())

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-
# Loading data from study 2 - PAX
load(paste0(data.folder, '/study2_data.RData'))

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Null Model - Intercept Only
# Initial null intercept model - will be necessary to generate final variance estimates
prior_config <- c(
  set_prior(
    'normal(0, 2)', 
    class = "Intercept"
  )
)

S2_POS_ucm_form <- bf(
  POS ~ 1 + (1|ID)
) + gaussian()

# Running model with priors (see above)
S2_POS_ucm <- brm_multiple(
  S2_POS_ucm_form,
  data = dat.study2_list, 
  chains = 3,
  cores = 3,
  iter = 15000,
  warmup = 10000, 
  control = list(adapt_delta = .99, 
                 max_treedepth = 15), 
  seed = 7201969, 
  save_all_pars = TRUE, 
  save_model = "S2_POS_ucm", 
  open_progress = FALSE,
  refresh = 0, 
  prior = prior_config
)

save(list = c("S2_POS_ucm", "S2_POS_ucm_form", "dat.study2_list", "dat.study2_model"), 
     file = '{study2.out}/S2_POS_ucm.RData' %>% glue()) 

sink(paste0(study2.model, '/S2_POS_ucm.txt'))
print(summary(S2_POS_ucm), digits = 5)
sink()

create_diagnostic_plots(S2_POS_ucm, dat.study2_list, n_samples = 100, dir_path = study2.graphics)

remove(list=c("S2_POS_ucm", "S2_POS_ucm_form"))
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-
# Model with a level-1 random effect for positive event rating (individually mean centered)
prior_config <- c(
  prior_config, 
  set_prior(
    'normal(0, 3)', 
    class = 'b'
  )
)

S2_POS_NegEvnt_form <- bf(
  POS ~ 1 + c.NegEvnt + (1 + c.NegEvnt|ID)
) + gaussian()

# Running model with priors (see above)
S2_POS_NegEvnt <- brm_multiple(
  S2_POS_NegEvnt_form,
  data = dat.study2_list, 
  chains = 3,
  cores = 3,
  iter = 15000,
  warmup = 10000, 
  control = list(adapt_delta = .99, 
                 max_treedepth = 15), 
  seed = 7201969, 
  save_all_pars = TRUE, 
  save_model = "S2_POS_NegEvnt", 
  open_progress = FALSE,
  refresh = 0, 
  prior = prior_config
)

save(list = c("S2_POS_NegEvnt", "S2_POS_NegEvnt_form", "dat.study2_list", "dat.study2_model"), 
     file = '{study2.out}/S2_POS_NegEvnt.RData' %>% glue()) 

sink(paste0(study2.model, '/S2_POS_NegEvnt.txt'))
print(summary(S2_POS_NegEvnt), digits = 5)
sink()

create_diagnostic_plots(S2_POS_NegEvnt, dat.study2_list, n_samples = 100, dir_path = study2.graphics)

remove(list=c("S2_POS_NegEvnt", "S2_POS_NegEvnt_form"))
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-
# Model with a level-1 random effect for positive event rating (individually mean-centered)
S2_POS_PosEvnt_form <- bf(
  POS ~ 1 + c.PosEvnt + (1 + c.PosEvnt|ID)
) + gaussian()

# Running model with priors (see above)
S2_POS_PosEvnt <- brm_multiple(
  S2_POS_PosEvnt_form,
  data = dat.study2_list, 
  chains = 3,
  cores = 3,
  iter = 15000,
  warmup = 10000, 
  control = list(adapt_delta = .99, 
                 max_treedepth = 15), 
  seed = 7201969, 
  save_all_pars = TRUE, 
  save_model = "S2_POS_PosEvnt", 
  open_progress = FALSE,
  refresh = 0, 
  prior = prior_config
)

save(list = c("S2_POS_PosEvnt", "S2_POS_PosEvnt_form", "dat.study2_list", "dat.study2_model"), 
     file = '{study2.out}/S2_POS_PosEvnt.RData' %>% glue()) 

sink(paste0(study2.model, '/S2_POS_PosEvnt.txt'))
print(summary(S2_POS_PosEvnt), digits = 5)
sink()

create_diagnostic_plots(S2_POS_PosEvnt, dat.study2_list, n_samples = 100, dir_path = study2.graphics)

remove(list=c("S2_POS_PosEvnt", "S2_POS_PosEvnt_form"))
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Model with level 1 random effect for recent positive event ratings and DN 
#   Note that due to the individually mean-centered level-1 effect, this model should capture an 
#   estimate of DN's total between-subjects "effect" on momentary mood. 
S2_POS_NegEvnt_DN_form <- bf(
  POS ~ 1 + c.NegEvnt + c.DN + (1 + c.NegEvnt|ID)
) + gaussian()

# Running model with priors (see above)
S2_POS_NegEvnt_DN <- brm_multiple(
  S2_POS_NegEvnt_DN_form,
  data = dat.study2_list, 
  chains = 3,
  cores = 3,
  iter = 15000,
  warmup = 10000, 
  control = list(adapt_delta = .99, 
                 max_treedepth = 15), 
  seed = 7201969, 
  save_all_pars = TRUE, 
  save_model = "S2_POS_NegEvnt_DN", 
  open_progress = FALSE,
  refresh = 0, 
  prior = prior_config
)

save(list = c("S2_POS_NegEvnt_DN", "S2_POS_NegEvnt_DN_form", "dat.study2_list", "dat.study2_model"), 
     file = '{study2.out}/S2_POS_NegEvnt_DN.RData' %>% glue()) 

sink(paste0(study2.model, '/S2_POS_NegEvnt_DN.txt'))
print(summary(S2_POS_NegEvnt_DN), digits = 5)
sink()

create_diagnostic_plots(S2_POS_NegEvnt_DN, dat.study2_list, n_samples = 100, dir_path = study2.graphics)

remove(list=c("S2_POS_NegEvnt_DN", "S2_POS_NegEvnt_DN_form"))
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Model with level 1 random effect for recent positive event ratings and DN 
#   Note that due to the individually mean-centered level-1 effect, this model should capture an 
#   estimate of DN's total between-subjects "effect" on momentary mood. Relatedly, accounting for 
#   slight differences related to the Bayesian approach, this effect should be effectively the same
#   as the one observed in the previous model. 
S2_POS_PosEvnt_DN_form <- bf(
  POS ~ 1 + c.PosEvnt + c.DN + (1 + c.PosEvnt|ID)
) + gaussian()

# Running model with priors (see above)
S2_POS_PosEvnt_DN <- brm_multiple(
  S2_POS_PosEvnt_DN_form,
  data = dat.study2_list, 
  chains = 3,
  cores = 3,
  iter = 15000,
  warmup = 10000, 
  control = list(adapt_delta = .99, 
                 max_treedepth = 15), 
  seed = 7201969, 
  save_all_pars = TRUE, 
  save_model = "S2_POS_PosEvnt_DN", 
  open_progress = FALSE,
  refresh = 0, 
  prior = prior_config
)

save(list = c("S2_POS_PosEvnt_DN", "S2_POS_PosEvnt_DN_form", "dat.study2_list", "dat.study2_model"), 
     file = '{study2.out}/S2_POS_PosEvnt_DN.RData' %>% glue()) 

sink(paste0(study2.model, '/S2_POS_PosEvnt_DN.txt'))
print(summary(S2_POS_PosEvnt_DN), digits = 5)
sink()

create_diagnostic_plots(S2_POS_PosEvnt_DN, dat.study2_list, n_samples = 100, dir_path = study2.graphics)

remove(list=c("S2_POS_PosEvnt_DN", "S2_POS_PosEvnt_DN_form"))
gc()


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# This model includes *average* positive event ratings as a measure of an individual's overall 
# positive event context. This is our proxy for overall exposure to more intense positive events. 
S2_POS_NegEvnt_Exp_form <- bf(
  POS ~ 1 + c.NegEvnt + prop.NegEvnt + (1 + c.NegEvnt|ID)
) + gaussian()

# Running model with priors (see above)
S2_POS_NegEvnt_Exp <- brm_multiple(
  S2_POS_NegEvnt_Exp_form,
  data = dat.study2_list, 
  chains = 3,
  cores = 3,
  iter = 15000,
  warmup = 10000, 
  control = list(adapt_delta = .99, 
                 max_treedepth = 15), 
  seed = 7201969, 
  save_all_pars = TRUE, 
  save_model = "S2_POS_NegEvnt_Exp", 
  open_progress = FALSE,
  refresh = 0, 
  prior = prior_config
)

save(list = c("S2_POS_NegEvnt_Exp", "S2_POS_NegEvnt_Exp_form", "dat.study2_list", "dat.study2_model"), 
     file = '{study2.out}/S2_POS_NegEvnt_Exp.RData' %>% glue()) 

sink(paste0(study2.model, '/S2_POS_NegEvnt_Exp.txt'))
print(summary(S2_POS_NegEvnt_Exp), digits = 5)
sink()

create_diagnostic_plots(S2_POS_NegEvnt_Exp, dat.study2_list, n_samples = 100, dir_path = study2.graphics)

remove(list=c("S2_POS_NegEvnt_Exp", "S2_POS_NegEvnt_Exp_form"))
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Model is the same as the one above, with the exception that the focus is on inclusion of between-
# subjects differences in average positive mood rating. The idea is similar though in that this is 
# our proxy for the overall contextual effect of reporting more intense positive events on average.
S2_POS_PosEvnt_Exp_form <- bf(
  POS ~ 1 + c.PosEvnt + prop.PosEvnt + (1 + c.PosEvnt|ID)
) + gaussian()

# Running model with priors (see above)
S2_POS_PosEvnt_Exp <- brm_multiple(
  S2_POS_PosEvnt_Exp_form,
  data = dat.study2_list, 
  chains = 3,
  cores = 3,
  iter = 15000,
  warmup = 10000, 
  control = list(adapt_delta = .99, 
                 max_treedepth = 15), 
  seed = 7201969, 
  save_all_pars = TRUE, 
  save_model = "S2_POS_PosEvnt_Exp", 
  open_progress = FALSE,
  refresh = 0, 
  prior = prior_config
)

save(list = c("S2_POS_PosEvnt_Exp", "S2_POS_PosEvnt_Exp_form", "dat.study2_list", "dat.study2_model"), 
     file = '{study2.out}/S2_POS_PosEvnt_Exp.RData' %>% glue()) 

sink(paste0(study2.model, '/S2_POS_PosEvnt_Exp.txt'))
print(summary(S2_POS_PosEvnt_Exp), digits = 5)
sink()

create_diagnostic_plots(S2_POS_PosEvnt_Exp, dat.study2_list, n_samples = 100, dir_path = study2.graphics)

remove(list=c("S2_POS_PosEvnt_Exp", "S2_POS_PosEvnt_Exp_form"))
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# This model in conjunction with the models including just DN or just average positive event ratings
# will be used to isolate the amount of "unique" variance attributable to DN, and to overall 
# positive contexts
S2_POS_NegEvnt_DN_Exp_form <- bf(
  POS ~ 1 + c.NegEvnt + prop.NegEvnt + c.DN + (1 + c.NegEvnt|ID)
) + gaussian()

# Running model with priors (see above)
S2_POS_NegEvnt_DN_Exp <- brm_multiple(
  S2_POS_NegEvnt_DN_Exp_form,
  data = dat.study2_list, 
  chains = 3,
  cores = 3,
  iter = 15000,
  warmup = 10000, 
  control = list(adapt_delta = .99, 
                 max_treedepth = 15), 
  seed = 7201969, 
  save_all_pars = TRUE, 
  save_model = "S2_POS_NegEvnt_DN_Exp", 
  open_progress = FALSE,
  refresh = 0, 
  prior = prior_config
)

save(list = c("S2_POS_NegEvnt_DN_Exp", "S2_POS_NegEvnt_DN_Exp_form", "dat.study2_list", "dat.study2_model"), 
     file = '{study2.out}/S2_POS_NegEvnt_DN_Exp.RData' %>% glue()) 

sink(paste0(study2.model, '/S2_POS_NegEvnt_DN_Exp.txt'))
print(summary(S2_POS_NegEvnt_DN_Exp), digits = 5)
sink()

create_diagnostic_plots(S2_POS_NegEvnt_DN_Exp, dat.study2_list, n_samples = 100, dir_path = study2.graphics)

remove(list=c("S2_POS_NegEvnt_DN_Exp", "S2_POS_NegEvnt_DN_Exp_form"))
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# This model in conjunction with the models including just DN or just average positive event ratings
# will be used to isolate the amount of "unique" variance attributable to DN, and to overall 
# positive contexts
S2_POS_PosEvnt_DN_Exp_form <- bf(
  POS ~ 1 + c.PosEvnt + prop.PosEvnt + c.DN + (1 + c.PosEvnt|ID)
) + gaussian()

# Running model with priors (see above)
S2_POS_PosEvnt_DN_Exp <- brm_multiple(
  S2_POS_PosEvnt_DN_Exp_form,
  data = dat.study2_list, 
  chains = 3,
  cores = 3,
  iter = 15000,
  warmup = 10000, 
  control = list(adapt_delta = .99, 
                 max_treedepth = 15), 
  seed = 7201969, 
  save_all_pars = TRUE, 
  save_model = "S2_POS_PosEvnt_DN_Exp", 
  open_progress = FALSE,
  refresh = 0, 
  prior = prior_config
)

save(list = c("S2_POS_PosEvnt_DN_Exp", "S2_POS_PosEvnt_DN_Exp_form", "dat.study2_list", "dat.study2_model"), 
     file = '{study2.out}/S2_POS_PosEvnt_DN_Exp.RData' %>% glue()) 

sink(paste0(study2.model, '/S2_POS_PosEvnt_DN_Exp.txt'))
print(summary(S2_POS_PosEvnt_DN_Exp), digits = 5)
sink()

create_diagnostic_plots(S2_POS_PosEvnt_DN_Exp, dat.study2_list, n_samples = 100, dir_path = study2.graphics)

remove(list=c("S2_POS_PosEvnt_DN_Exp", "S2_POS_PosEvnt_DN_Exp_form"))
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# This model actually represents the "final" model for this modeling tree (DN, positive events, 
# positive mood). All variance components obtained in the present set of analyses stem from 
# isolated R2 differences moving from the unconditional model to this model. 
S2_POS_NegEvnt_Rct_form <- bf(
  POS ~ 1 + c.NegEvnt * c.DN + prop.NegEvnt + c.DN + (1 + c.NegEvnt|ID)
) + gaussian()

# Running model with priors (see above)
S2_POS_NegEvnt_Rct <- brm_multiple(
  S2_POS_NegEvnt_Rct_form,
  data = dat.study2_list, 
  chains = 3,
  cores = 3,
  iter = 15000,
  warmup = 10000, 
  control = list(adapt_delta = .99, 
                 max_treedepth = 15), 
  seed = 7201969, 
  save_all_pars = TRUE, 
  save_model = "S2_POS_NegEvnt_Rct", 
  open_progress = FALSE,
  refresh = 0, 
  prior = prior_config
)

save(list = c("S2_POS_NegEvnt_Rct", "S2_POS_NegEvnt_Rct_form", "dat.study2_list", "dat.study2_model"), 
     file = '{study2.out}/S2_POS_NegEvnt_Rct.RData' %>% glue()) 

sink(paste0(study2.model, '/S2_POS_NegEvnt_Rct.txt'))
print(summary(S2_POS_NegEvnt_Rct), digits = 5)
sink()

create_diagnostic_plots(S2_POS_NegEvnt_Rct, dat.study2_list, n_samples = 100, dir_path = study2.graphics)

remove(list=c("S2_POS_NegEvnt_Rct", "S2_POS_NegEvnt_Rct_form"))
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# This model actually represents the "final" model for this modeling tree (DN, positive events, 
# positive mood). All variance components obtained in the present set of analyses stem from 
# isolated R2 differences moving from the unconditional model to this model. 
S2_POS_PosEvnt_Rct_form <- bf(
  POS ~ 1 + c.PosEvnt * c.DN + prop.PosEvnt + c.DN + (1 + c.PosEvnt|ID)
) + gaussian()

# Running model with priors (see above)
S2_POS_PosEvnt_Rct <- brm_multiple(
  S2_POS_PosEvnt_Rct_form,
  data = dat.study2_list, 
  chains = 3,
  cores = 3,
  iter = 15000,
  warmup = 10000, 
  control = list(adapt_delta = .99, 
                 max_treedepth = 15), 
  seed = 7201969, 
  save_all_pars = TRUE, 
  save_model = "S2_POS_PosEvnt_Rct", 
  open_progress = FALSE,
  refresh = 0, 
  prior = prior_config
)

save(list = c("S2_POS_PosEvnt_Rct", "S2_POS_PosEvnt_Rct_form", "dat.study2_list", "dat.study2_model"), 
     file = '{study2.out}/S2_POS_PosEvnt_Rct.RData' %>% glue()) 

sink(paste0(study2.model, '/S2_POS_PosEvnt_Rct.txt'))
print(summary(S2_POS_PosEvnt_Rct), digits = 5)
sink()

create_diagnostic_plots(S2_POS_PosEvnt_Rct, dat.study2_list, n_samples = 100, dir_path = study2.graphics)

remove(list=c("S2_POS_PosEvnt_Rct", "S2_POS_PosEvnt_Rct_form"))
gc()
