###############################################################################
# Study 2 Modeling Script

# Study 2 Negative Mood Models:
###############################################################################

#Imputation will be handled in the modeling stage itself. 

#The same exact imputation model will be used in each case... 

#These models will take longer to run and the output will be denser
#Should be faster than fitting to multiple data sets though

#------------------------------------------------------------------------------
library(brms)
library(rstan)
library(rstanarm)
library(bayesplot)
library(pan)
library(mitml)
library(mice)
library(parallel)
library(RColorBrewer)
library(ggridges)
library(riverplot)
library(ggalluvial)
library(ggpubr)
library(tidyverse)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
user<-ifelse(Sys.getenv("USERPROFILE")=="", "~", Sys.getenv("USERPROFILE"))
wd<-paste0(user, '/Dropbox/UMD/Shackman Lab/EMA_MS')
data.folder<-paste0(wd, '/Data')
study2.out<-paste0(wd, '/Study 2 output')
study2.graphics<-paste0(study2.out, '/Graphics')
study2.model<-paste0(study2.out, '/Model summaries')
stan.code<-paste0(wd, '/Stan_code')
EDA.folder<-paste0(study2.out, '/EDA')
#------------------------------------------------------------------------------
#Loading Study 2 Data with imputed values
load(paste0(data.folder, '/Study2_Clean_w_Impute.RData'))

dat.study2_list_stacked <- data.frame() 
for(m in 1:M) # Assumes M is still defined in the imputed data set loaded from .RData above
  dat.study2_list_stacked <- rbind(dat.study2_list_stacked, dat.study2_list[[i]])

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Creating utility functions to gather priors based on imputed data sets
get_imputed_mean <- function(data, variable, logged=TRUE){
  #browser()
  M <- length(data) # data needs to be a list object
  tmp_vec <- c()
  for(m in 1:M){
    var <- unlist(data[[m]][variable])
    if(logged){
      var <- log(var)
    }
    tmp_vec <- c(tmp_vec, mean(var))
  }
  return(mean(tmp_vec))
}

get_imputed_sd <- function(data, variable, logged=TRUE){
  #browser()
  M <- length(data) # data needs to be a list object
  tmp_vec <- c()
  for(m in 1:M){
    var <- unlist(data[[m]][variable])
    if(logged){
      var <- log(var)
    }
    tmp_vec <- c(tmp_vec, sd(var))
  }
  return(mean(tmp_vec))
}

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Anxiety Null Model - Intercept Only
S2_DEP_ucm_form <- bf(
  DEP ~ 1 + (1|ID)
)+lognormal()

#Creating priors for intercept to help constrain final model in response space
mu <- get_imputed_mean(dat.study2_list, "DEP")
sigma <- get_imputed_sd(dat.study2_list, "DEP")

Int_prior <- set_prior(paste0("normal(", mu, ",", sigma, ")"), 
                       class = "Intercept")

#Running model with priors (see above)
S2_DEP_ucm <- brm_multiple(S2_DEP_ucm_form,
                           data = dat.study2_list,
                           chains = 3,
                           prior = c(Int_prior),
                           iter = 25000,
                           warmup = 20000,
                           control = list(adapt_delta = .99,
                                          max_treedepth = 15))
  
sink(paste0(study2.model, '/S2_DEP_ucm.txt'))
print(summary(S2_DEP_ucm), digits = 5)
sink()
  
#Simple Model Check plotting:
ppc_density <- 
  pp_check(S2_DEP_ucm, 
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_DEP_ucm Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S2_DEP_ucm, 
           newdata = dat.study2_list_stacked,
           nsamples = 12, 
           type = "error_hist")+
  ggtitle("S2_DEP_ucm Model Posterior Residuals")

#Saving Plots:
png(paste0(study2.graphics, '/S2_DEP_ucm_ppc.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)
dev.off()

save(list = c('S2_DEP_ucm', 'S2_DEP_ucm_form'), 
     file = '/media/matthew/My Book/EMA_S2_Bayesian_Posteriors/S2_DEP_ucm.RData')
remove(S2_DEP_ucm)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
beta_prior <- set_prior('normal(0, 2)', class='b')

S2_DEP_DN_form <- bf(
  DEP ~ 1 + c.DN + (1|ID)
)+lognormal()

S2_DEP_DN <- brm_multiple(S2_DEP_DN_form,
                          data = dat.study2_list,
                          chains = 3,
                          prior = c(Int_prior, 
                                    beta_prior),
                          iter = 25000,
                          warmup = 20000,
                          control = list(adapt_delta = .99,
                                         max_treedepth = 15))

sink(paste0(study2.model, '/S2_DEP_DN.txt'))
print(summary(S2_DEP_DN), digits = 5)
sink()

#Simple Model Check plotting:
ppc_density <- 
  pp_check(S2_DEP_DN, 
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_DEP_DN Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S2_DEP_DN, 
           newdata = dat.study2_list_stacked,
           nsamples = 12, 
           type = "error_hist")+
  ggtitle("S2_DEP_DN Model Posterior Residuals")

#Saving Plots:
png(paste0(study2.graphics, '/S2_DEP_DN_ppc.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)
dev.off()

save(list = c('S2_DEP_DN', 'S2_DEP_DN_form'), 
     file = '/media/matthew/My Book/EMA_S2_Bayesian_Posteriors/S2_DEP_DN.RData')
remove(S2_DEP_DN)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
S2_DEP_prop.NegEvnt_form <- bf(
  DEP ~ 1 + prop.NegEvnt + (1|ID)
)+lognormal()

S2_DEP_prop.NegEvnt <- brm_multiple(S2_DEP_prop.NegEvnt_form,
                                    data = dat.study2_list,
                                    chains = 3,
                                    prior = c(Int_prior, 
                                              beta_prior),
                                    iter = 25000,
                                    warmup = 20000,
                                    control = list(adapt_delta = .99,
                                                   max_treedepth = 15))

sink(paste0(study2.model, '/S2_DEP_prop.NegEvnt.txt'))
print(summary(S2_DEP_prop.NegEvnt), digits = 5)
sink()

#Simple Model Check plotting:
ppc_density <- 
  pp_check(S2_DEP_prop.NegEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_DEP_prop.NegEvnt Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S2_DEP_prop.NegEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 12, 
           type = "error_hist")+
  ggtitle("S2_DEP_prop.NegEvnt Model Posterior Residuals")

#Saving Plots:
png(paste0(study2.graphics, '/S2_DEP_prop.NegEvnt_ppc.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)
dev.off()

save(list = c('S2_DEP_prop.NegEvnt', 'S2_DEP_prop.NegEvnt_form'), 
     file = '/media/matthew/My Book/EMA_S2_Bayesian_Posteriors/S2_DEP_prop.NegEvnt.RData')
remove(S2_DEP_prop.NegEvnt)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
S2_DEP_prop.PosEvnt_form <- bf(
  DEP ~ 1 + prop.PosEvnt + (1|ID)
)+lognormal()

S2_DEP_prop.PosEvnt <- brm_multiple(S2_DEP_prop.PosEvnt_form,
                                    data = dat.study2_list,
                                    chains = 3,
                                    prior = c(Int_prior, 
                                              beta_prior),
                                    iter = 25000,
                                    warmup = 20000,
                                    control = list(adapt_delta = .99,
                                                   max_treedepth = 15))

sink(paste0(study2.model, '/S2_DEP_prop.PosEvnt.txt'))
print(summary(S2_DEP_prop.PosEvnt), digits = 5)
sink()

#Simple Model Check plotting:
ppc_density <- 
  pp_check(S2_DEP_prop.PosEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_DEP_prop.PosEvnt Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S2_DEP_prop.PosEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 12, 
           type = "error_hist")+
  ggtitle("S2_DEP_prop.PosEvnt Model Posterior Residuals")

#Saving Plots:
png(paste0(study2.graphics, '/S2_DEP_prop.PosEvnt_ppc.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)
dev.off()

save(list = c('S2_DEP_prop.PosEvnt', 'S2_DEP_prop.PosEvnt_form'), 
     file = '/media/matthew/My Book/EMA_S2_Bayesian_Posteriors/S2_DEP_prop.PosEvnt.RData')
remove(S2_DEP_prop.PosEvnt)
gc()


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
S2_DEP_DN_prop.NegEvnt_form <- bf(
  DEP ~ 1 + c.DN + prop.NegEvnt + (1|ID)
)+lognormal()

S2_DEP_DN_prop.NegEvnt <- brm_multiple(S2_DEP_DN_prop.NegEvnt_form,
                                       data = dat.study2_list,
                                       chains = 3,
                                       prior = c(Int_prior, 
                                                 beta_prior),
                                       iter = 25000,
                                       warmup = 20000,
                                       control = list(adapt_delta = .99,
                                                      max_treedepth = 15))

sink(paste0(study2.model, '/S2_DEP_DN_prop.NegEvnt.txt'))
print(summary(S2_DEP_DN_prop.NegEvnt), digits = 5)
sink()

#Simple Model Check plotting:
ppc_density <- 
  pp_check(S2_DEP_DN_prop.NegEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_DEP_DN_prop.NegEvnt Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S2_DEP_DN_prop.NegEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 12, 
           type = "error_hist")+
  ggtitle("S2_DEP_DN_prop.NegEvnt Model Posterior Residuals")

#Saving Plots:
png(paste0(study2.graphics, '/S2_DEP_DN_prop.NegEvnt_ppc.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)
dev.off()

save(list = c('S2_DEP_DN_prop.NegEvnt', 'S2_DEP_DN_prop.NegEvnt_form'), 
     file = '/media/matthew/My Book/EMA_S2_Bayesian_Posteriors/S2_DEP_DN_prop.NegEvnt.RData')
remove(S2_DEP_DN_prop.NegEvnt)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
S2_DEP_DN_prop.PosEvnt_form <- bf(
  DEP ~ 1 + c.DN + prop.PosEvnt + (1|ID)
)+lognormal()

S2_DEP_DN_prop.PosEvnt <- brm_multiple(S2_DEP_DN_prop.PosEvnt_form,
                                       data = dat.study2_list,
                                       chains = 3,
                                       prior = c(Int_prior, 
                                                 beta_prior),
                                       iter = 25000,
                                       warmup = 20000,
                                       control = list(adapt_delta = .99,
                                                      max_treedepth = 15))

sink(paste0(study2.model, '/S2_DEP_DN_prop.PosEvnt.txt'))
print(summary(S2_DEP_DN_prop.PosEvnt), digits = 5)
sink()

#Simple Model Check plotting:
ppc_density <- 
  pp_check(S2_DEP_DN_prop.PosEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_DEP_DN_prop.PosEvnt Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S2_DEP_DN_prop.PosEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 12, 
           type = "error_hist")+
  ggtitle("S2_DEP_DN_prop.PosEvnt Model Posterior Residuals")

#Saving Plots:
png(paste0(study2.graphics, '/S2_DEP_DN_prop.PosEvnt_ppc.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)
dev.off()

save(list = c('S2_DEP_DN_prop.PosEvnt', 'S2_DEP_DN_prop.PosEvnt_form'), 
     file = '/media/matthew/My Book/EMA_S2_Bayesian_Posteriors/S2_DEP_DN_prop.PosEvnt.RData')
remove(S2_DEP_DN_prop.PosEvnt)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
S2_DEP_NegEvnt_form <- bf(
  DEP ~ 1 + c.NegEvnt + (1 + c.NegEvnt|ID)
)+lognormal()

cor_prior <- set_prior('lkj(2)', class='cor')

#Running model with priors (see above)
S2_DEP_NegEvnt <- brm_multiple(S2_DEP_NegEvnt_form,
                               data = dat.study2_list, 
                               chains = 3,
                               prior = c(Int_prior, 
                                         beta_prior, 
                                         cor_prior),
                               iter = 25000,
                               warmup = 20000,
                               control = list(adapt_delta = .99, 
                                              max_treedepth = 15))

sink(paste0(study2.model, '/S2_DEP_NegEvnt.txt'))
print(summary(S2_DEP_NegEvnt), digits = 5)
sink()

#Simple Model Check plotting:
ppc_density <- 
  pp_check(S2_DEP_NegEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_DEP_NegEvnt Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S2_DEP_NegEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 12, 
           type = "error_hist")+
  ggtitle("S2_DEP_NegEvnt Model Posterior Residuals")

#Saving Plots:
png(paste0(study2.graphics, '/S2_DEP_NegEvnt_ppc.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)
dev.off()

save(list = c('S2_DEP_NegEvnt', 'S2_DEP_NegEvnt_form'), 
     file = '/media/matthew/My Book/EMA_S2_Bayesian_Posteriors/S2_DEP_NegEvnt.RData')
remove(S2_DEP_NegEvnt)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
S2_DEP_PosEvnt_form <- bf(
  DEP ~ 1 + c.PosEvnt + (1 + c.PosEvnt|ID)
)+lognormal()

#Running model with priors (see above)
S2_DEP_PosEvnt <- brm_multiple(S2_DEP_PosEvnt_form,
                               data = dat.study2_list, 
                               chains = 3,
                               prior = c(Int_prior, 
                                         beta_prior, 
                                         cor_prior),
                               iter = 25000,
                               warmup = 20000,
                               control = list(adapt_delta = .99, 
                                              max_treedepth = 15))

sink(paste0(study2.model, '/S2_DEP_PosEvnt.txt'))
print(summary(S2_DEP_PosEvnt), digits = 5)
sink()

#Simple Model Check plotting:
ppc_density <- 
  pp_check(S2_DEP_PosEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_DEP_PosEvnt Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S2_DEP_PosEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 12, 
           type = "error_hist")+
  ggtitle("S2_DEP_PosEvnt Model Posterior Residuals")

#Saving Plots:
png(paste0(study2.graphics, '/S2_DEP_PosEvnt_ppc.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)
dev.off()

save(list = c('S2_DEP_PosEvnt', 'S2_DEP_PosEvnt_form'), 
     file = '/media/matthew/My Book/EMA_S2_Bayesian_Posteriors/S2_DEP_PosEvnt.RData')
remove(S2_DEP_PosEvnt)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
S2_DEP_NegEvnt_DN_form <- bf(
  DEP ~ 1 + c.NegEvnt + c.DN + (1 + c.NegEvnt|ID)
)+lognormal()

#Running model with priors (see above)
S2_DEP_NegEvnt_DN <- brm_multiple(S2_DEP_NegEvnt_DN_form,
                                 data = dat.study2_list, 
                                 chains = 3,
                                 prior = c(Int_prior, 
                                           beta_prior, 
                                           cor_prior),
                                 iter = 25000,
                                 warmup = 20000,
                                 control = list(adapt_delta = .99, 
                                                max_treedepth = 15))

sink(paste0(study2.model, '/S2_DEP_NegEvnt_DN.txt'))
print(summary(S2_DEP_NegEvnt_DN), digits = 5)
sink()

#Simple Model Check plotting:
ppc_density <- 
  pp_check(S2_DEP_NegEvnt_DN, 
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_DEP_NegEvnt_DN Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S2_DEP_NegEvnt_DN, 
           newdata = dat.study2_list_stacked,
           nsamples = 12, 
           type = "error_hist")+
  ggtitle("S2_DEP_NegEvnt_DN Model Posterior Residuals")

#Saving Plots:
png(paste0(study2.graphics, '/S2_DEP_NegEvnt_DN_ppc.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)
dev.off()

save(list = c('S2_DEP_NegEvnt_DN', 'S2_DEP_NegEvnt_DN_form'), 
     file = '/media/matthew/My Book/EMA_S2_Bayesian_Posteriors/S2_DEP_NegEvnt_DN.RData')
remove(S2_DEP_NegEvnt_DN)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
S2_DEP_PosEvnt_DN_form <- bf(
  DEP ~ 1 + c.PosEvnt + c.DN + (1 + c.PosEvnt|ID)
)+lognormal()

#Running model with priors (see above)
S2_DEP_PosEvnt_DN <- brm_multiple(S2_DEP_PosEvnt_DN_form,
                                 data = dat.study2_list, 
                                 chains = 3,
                                 prior = c(Int_prior, 
                                           beta_prior, 
                                           cor_prior),
                                 iter = 25000,
                                 warmup = 20000,
                                 control = list(adapt_delta = .99, 
                                                max_treedepth = 15))

sink(paste0(study2.model, '/S2_DEP_PosEvnt_DN.txt'))
print(summary(S2_DEP_PosEvnt_DN), digits = 5)
sink()

#Simple Model Check plotting:
ppc_density <- 
  pp_check(S2_DEP_PosEvnt_DN, 
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_DEP_PosEvnt_DN Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S2_DEP_PosEvnt_DN, 
           newdata = dat.study2_list_stacked,
           nsamples = 12, 
           type = "error_hist")+
  ggtitle("S2_DEP_PosEvnt_DN Model Posterior Residuals")

#Saving Plots:
png(paste0(study2.graphics, '/S2_DEP_PosEvnt_DN_ppc.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)
dev.off()

save(list = c('S2_DEP_PosEvnt_DN', 'S2_DEP_PosEvnt_DN_form'), 
     file = '/media/matthew/My Book/EMA_S2_Bayesian_Posteriors/S2_DEP_PosEvnt_DN.RData')
remove(S2_DEP_PosEvnt_DN)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
S2_DEP_NegEvnt_prop.NegEvnt_form <- bf(
  DEP ~ 1 + c.NegEvnt + prop.NegEvnt + (1 + c.NegEvnt|ID)
)+lognormal()

#Running model with priors (see above)
S2_DEP_NegEvnt_prop.NegEvnt <- brm_multiple(S2_DEP_NegEvnt_prop.NegEvnt_form,
                                  data = dat.study2_list, 
                                  chains = 3,
                                  prior = c(Int_prior, 
                                            beta_prior, 
                                            cor_prior),
                                  iter = 25000,
                                  warmup = 20000,
                                  control = list(adapt_delta = .99, 
                                                 max_treedepth = 15))

sink(paste0(study2.model, '/S2_DEP_NegEvnt_prop.NegEvnt.txt'))
print(summary(S2_DEP_NegEvnt_prop.NegEvnt), digits = 5)
sink()

#Simple Model Check plotting:
ppc_density <- 
  pp_check(S2_DEP_NegEvnt_prop.NegEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_DEP_NegEvnt_prop.NegEvnt Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S2_DEP_NegEvnt_prop.NegEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 12, 
           type = "error_hist")+
  ggtitle("S2_DEP_NegEvnt_prop.NegEvnt Model Posterior Residuals")

#Saving Plots:
png(paste0(study2.graphics, '/S2_DEP_NegEvnt_prop.NegEvnt_ppc.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)
dev.off()

save(list = c('S2_DEP_NegEvnt_prop.NegEvnt', 'S2_DEP_NegEvnt_prop.NegEvnt_form'), 
     file = '/media/matthew/My Book/EMA_S2_Bayesian_Posteriors/S2_DEP_NegEvnt_prop.NegEvnt.RData')
remove(S2_DEP_NegEvnt_prop.NegEvnt)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
S2_DEP_PosEvnt_prop.PosEvnt_form <- bf(
  DEP ~ 1 + c.PosEvnt + prop.PosEvnt + (1 + c.PosEvnt|ID)
)+lognormal()

#Running model with priors (see above)
S2_DEP_PosEvnt_prop.PosEvnt <- brm_multiple(S2_DEP_PosEvnt_prop.PosEvnt_form,
                                  data = dat.study2_list, 
                                  chains = 3,
                                  prior = c(Int_prior, 
                                            beta_prior, 
                                            cor_prior),
                                  iter = 25000,
                                  warmup = 20000,
                                  control = list(adapt_delta = .99, 
                                                 max_treedepth = 15))

sink(paste0(study2.model, '/S2_DEP_PosEvnt_prop.PosEvnt.txt'))
print(summary(S2_DEP_PosEvnt_prop.PosEvnt), digits = 5)
sink()

#Simple Model Check plotting:
ppc_density <- 
  pp_check(S2_DEP_PosEvnt_prop.PosEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_DEP_PosEvnt_prop.PosEvnt Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S2_DEP_PosEvnt_prop.PosEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 12, 
           type = "error_hist")+
  ggtitle("S2_DEP_PosEvnt_prop.PosEvnt Model Posterior Residuals")

#Saving Plots:
png(paste0(study2.graphics, '/S2_DEP_PosEvnt_prop.PosEvnt_ppc.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)
dev.off()

save(list = c('S2_DEP_PosEvnt_prop.PosEvnt', 'S2_DEP_PosEvnt_prop.PosEvnt_form'), 
     file = '/media/matthew/My Book/EMA_S2_Bayesian_Posteriors/S2_DEP_PosEvnt_prop.PosEvnt.RData')
remove(S2_DEP_PosEvnt_prop.PosEvnt)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
S2_DEP_NegEvnt_DN_prop.NegEvnt_form <- bf(
  DEP ~ 1 + c.NegEvnt + c.DN + prop.NegEvnt + (1 + c.NegEvnt|ID)
)+lognormal()

#Running model with priors (see above)
S2_DEP_NegEvnt_DN_prop.NegEvnt <- brm_multiple(S2_DEP_NegEvnt_DN_prop.NegEvnt_form,
                                            data = dat.study2_list, 
                                            chains = 3,
                                            prior = c(Int_prior, 
                                                      beta_prior, 
                                                      cor_prior),
                                            iter = 25000,
                                            warmup = 20000,
                                            control = list(adapt_delta = .99, 
                                                           max_treedepth = 15))

sink(paste0(study2.model, '/S2_DEP_NegEvnt_DN_prop.NegEvnt.txt'))
print(summary(S2_DEP_NegEvnt_DN_prop.NegEvnt), digits = 5)
sink()

#Simple Model Check plotting:
ppc_density <- 
  pp_check(S2_DEP_NegEvnt_DN_prop.NegEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_DEP_NegEvnt_DN_prop.NegEvnt Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S2_DEP_NegEvnt_DN_prop.NegEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 12, 
           type = "error_hist")+
  ggtitle("S2_DEP_NegEvnt_DN_prop.NegEvnt Model Posterior Residuals")

#Saving Plots:
png(paste0(study2.graphics, '/S2_DEP_NegEvnt_DN_prop.NegEvnt_ppc.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)
dev.off()

save(list = c('S2_DEP_NegEvnt_DN_prop.NegEvnt', 'S2_DEP_NegEvnt_DN_prop.NegEvnt_form'), 
     file = '/media/matthew/My Book/EMA_S2_Bayesian_Posteriors/S2_DEP_NegEvnt_DN_prop.NegEvnt.RData')
remove(S2_DEP_NegEvnt_DN_prop.NegEvnt)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
S2_DEP_PosEvnt_DN_prop.PosEvnt_form <- bf(
  DEP ~ 1 + c.PosEvnt + c.DN + prop.PosEvnt + (1 + c.PosEvnt|ID)
)+lognormal()

#Running model with priors (see above)
S2_DEP_PosEvnt_DN_prop.PosEvnt <- brm_multiple(S2_DEP_PosEvnt_DN_prop.PosEvnt_form,
                                            data = dat.study2_list, 
                                            chains = 3,
                                            prior = c(Int_prior, 
                                                      beta_prior, 
                                                      cor_prior),
                                            iter = 25000,
                                            warmup = 20000,
                                            control = list(adapt_delta = .99, 
                                                           max_treedepth = 15))

sink(paste0(study2.model, '/S2_DEP_PosEvnt_DN_prop.PosEvnt.txt'))
print(summary(S2_DEP_PosEvnt_DN_prop.PosEvnt), digits = 5)
sink()

#Simple Model Check plotting:
ppc_density <- 
  pp_check(S2_DEP_PosEvnt_DN_prop.PosEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_DEP_PosEvnt_DN_prop.PosEvnt Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S2_DEP_PosEvnt_DN_prop.PosEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 12, 
           type = "error_hist")+
  ggtitle("S2_DEP_PosEvnt_DN_prop.PosEvnt Model Posterior Residuals")

#Saving Plots:
png(paste0(study2.graphics, '/S2_DEP_PosEvnt_DN_prop.PosEvnt_ppc.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)
dev.off()

save(list = c('S2_DEP_PosEvnt_DN_prop.PosEvnt', 'S2_DEP_PosEvnt_DN_prop.PosEvnt_form'), 
     file = '/media/matthew/My Book/EMA_S2_Bayesian_Posteriors/S2_DEP_PosEvnt_DN_prop.PosEvnt.RData')
remove(S2_DEP_PosEvnt_DN_prop.PosEvnt)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
S2_DEP_NegEvnt_x_DN_prop.NegEvnt_form <- bf(
  DEP ~ 1 + c.NegEvnt * c.DN + prop.NegEvnt + (1 + c.NegEvnt|ID)
)+lognormal()

#Running model with priors (see above)
S2_DEP_NegEvnt_x_DN_prop.NegEvnt <- brm_multiple(S2_DEP_NegEvnt_x_DN_prop.NegEvnt_form,
                                            data = dat.study2_list, 
                                            chains = 3,
                                            prior = c(Int_prior, 
                                                      beta_prior, 
                                                      cor_prior),
                                            iter = 25000,
                                            warmup = 20000,
                                            control = list(adapt_delta = .99, 
                                                           max_treedepth = 15))

sink(paste0(study2.model, '/S2_DEP_NegEvnt_x_DN_prop.NegEvnt.txt'))
print(summary(S2_DEP_NegEvnt_x_DN_prop.NegEvnt), digits = 5)
sink()

#Simple Model Check plotting:
ppc_density <- 
  pp_check(S2_DEP_NegEvnt_x_DN_prop.NegEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_DEP_NegEvnt_x_DN_prop.NegEvnt Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S2_DEP_NegEvnt_x_DN_prop.NegEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 12, 
           type = "error_hist")+
  ggtitle("S2_DEP_NegEvnt_x_DN_prop.NegEvnt Model Posterior Residuals")

#Saving Plots:
png(paste0(study2.graphics, '/S2_DEP_NegEvnt_x_DN_prop.NegEvnt_ppc.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)
dev.off()

save(list = c('S2_DEP_NegEvnt_x_DN_prop.NegEvnt', 'S2_DEP_NegEvnt_x_DN_prop.NegEvnt_form'), 
     file = '/media/matthew/My Book/EMA_S2_Bayesian_Posteriors/S2_DEP_NegEvnt_x_DN_prop.NegEvnt.RData')
remove(S2_DEP_NegEvnt_DN_prop.NegEvnt)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
S2_DEP_PosEvnt_x_DN_prop.PosEvnt_form <- bf(
  DEP ~ 1 + c.PosEvnt * c.DN + prop.PosEvnt + (1 + c.PosEvnt|ID)
)+lognormal()

#Running model with priors (see above)
S2_DEP_PosEvnt_x_DN_prop.PosEvnt <- brm_multiple(S2_DEP_PosEvnt_x_DN_prop.PosEvnt_form,
                                               data = dat.study2_list, 
                                               chains = 3,
                                               prior = c(Int_prior, 
                                                         beta_prior, 
                                                         cor_prior),
                                               iter = 25000,
                                               warmup = 20000, 
                                               control = list(adapt_delta = .99, 
                                                              max_treedepth = 15))

sink(paste0(study2.model, '/S2_DEP_PosEvnt_x_DN_prop.PosEvnt.txt'))
print(summary(S2_DEP_PosEvnt_x_DN_prop.PosEvnt), digits = 5)
sink()

#Simple Model Check plotting:
ppc_density <- 
  pp_check(S2_DEP_PosEvnt_x_DN_prop.PosEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_DEP_PosEvnt_x_DN_prop.PosEvnt Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S2_DEP_PosEvnt_x_DN_prop.PosEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 12, 
           type = "error_hist")+
  ggtitle("S2_DEP_PosEvnt_x_DN_prop.PosEvnt Model Posterior Residuals")

#Saving Plots:
png(paste0(study2.graphics, '/S2_DEP_PosEvnt_x_DN_prop.PosEvnt_ppc.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)
dev.off()

save(list = c('S2_DEP_PosEvnt_x_DN_prop.PosEvnt', 'S2_DEP_PosEvnt_x_DN_prop.PosEvnt_form'), 
     file = '/media/matthew/My Book/EMA_S2_Bayesian_Posteriors/S2_DEP_PosEvnt_x_DN_prop.PosEvnt.RData')
remove(S2_DEP_PosEvnt_x_DN_prop.PosEvnt)
gc()
