###################################################################################################
# Study 1 Modeling Script

# Negative and Positive event MLMs
# Description: 
#   The analyses below examine a bivariate association between individual dispositional negativity 
#   scores and the occurrence of recent negative and positive events.
###################################################################################################

#--------------------------------------------------------------------------------------------------
# Loading relevant packages
library(brms)
library(rstan)
library(rstanarm)
library(bayesplot)
library(tidyverse)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
#--------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------
# Location of repo stored locally
wd<-paste0('~/github/ATNL/shackman-umd-pax-ema-pub')
data.folder<-paste0(wd, '/Data')
study2.model<-paste0(wd, '/Study_2_model_summaries')

# Will save very large posterior files from analyses (not recommended for git repo)
# For anyone attempting to reproduce these analyses be sure to identify a storage location with sufficient memory
posterior_save_dir <- "/media/dr-owner/HDD1"
study2.out<-paste0(posterior_save_dir, '/EMA_S2_Bayesian_Posteriors')

# Also generally not recommended to store image files on GH... 
study2.graphics<-paste0(posterior_save_dir, '/EMA_S2_Graphics')
#--------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------
load(paste0(data.folder, '/study2_data.RData'))

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-
# Simple Event Models:

# Substantive Questions: 
#   1. Is dispositional negativity associated with intensity of recent negative events? 
#   2. Is dispositional negativity associated with intensity of recent positive events? 
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-
dat.study2_model$NegEvnt_r <- as.numeric(dat.study2_model$NegEvnt == "Yes")
dat.study2_model$PosEvnt_r <- as.numeric(dat.study2_model$PosEvnt == "Yes")

#Negative Event Model: 
NegEvnt_DN <- bf(
  NegEvnt_r ~ 1 + c.DN + (1|ID) 
)+bernoulli()

#Creating priors for intercept to help constrain final model in response space

model_priors <- c(
  prior(student_t(3, 0, 2.5), class = "Intercept"),
  prior(student_t(3, 0, 2.5), class = "b")
)

#Running model with priors (see above)
S2_NegEvnt_DN <- brm(NegEvnt_DN,
                     data = dat.study2_model, 
                     chains = 3,
                     prior = model_priors,
                     iter = 15000,
                     warmup = 10000, 
                     control = list(adapt_delta = .99, 
                                    max_treedepth = 15))

save(list = c("S2_NegEvnt_DN", "dat.study2_model"), 
     file = paste0(posterior_save_dir, "/", 
                   '/S2_NegEvnt_DN.RData'))

sink(paste0(study2.model, '/S2_NegEvnt_DN.txt'))
print(summary(S2_NegEvnt_DN), digits = 5)
sink()

#Simple Model Check plotting:
ppc_density <- 
pp_check(S2_NegEvnt_DN, 
         newdata = dat.study2_model[!is.na(dat.study2_model$NegEvnt_r),], 
         nsamples = 100)+
  ggtitle("S2_NegEvnt_DN Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S2_NegEvnt_DN, 
           newdata = dat.study2_model[!is.na(dat.study2_model$NegEvnt_r),], 
           nsamples = 12, 
           type = "error_hist")+
  ggtitle("S2_NegEvnt_DN Model Posterior Residuals")

png(paste0(study2.graphics, '/S2_NegEvnt_DN_ppc.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)

dev.off()
remove("S2_NegEvnt_DN")

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-
#Positive Event Model: 
PosEvnt_DN <- bf(
  PosEvnt ~ 1 + c.DN + (1|ID) 
)+bernoulli()

#Running model with priors (see above)
S2_PosEvnt_DN <- brm(PosEvnt_DN,
                     data = dat.study2_model, 
                     chains = 3,
                     prior = model_priors,
                     iter = 15000,
                     warmup = 10000, 
                     control = list(adapt_delta = .99, 
                                    max_treedepth = 15))

save(list = c("S2_PosEvnt_DN", "dat.study2_model"), 
     file = paste0(posterior_save_dir, "/", 
                   '/S2_PosEvnt_DN.RData'))

sink(paste0(study2.model, '/S2_PosEvnt_DN.txt'))
print(summary(S2_PosEvnt_DN), digits = 5)
sink()

#Simple Model Check plotting:
ppc_density <- 
  pp_check(S2_PosEvnt_DN, 
           newdata = dat.study2_model[!is.na(dat.study2_model$PosEvnt_r),], 
           nsamples = 100)+
  ggtitle("S2_PosEvnt_DN Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S2_PosEvnt_DN, 
           newdata = dat.study2_model[!is.na(dat.study2_model$PosEvnt_r),], 
           nsamples = 12, 
           type = "error_hist")+
  ggtitle("S2_PosEvnt_DN Model Posterior Residuals")

png(paste0(study2.graphics, '/S2_PosEvnt_DN_ppc.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)

dev.off()
remove("S2_PosEvnt_DN")
