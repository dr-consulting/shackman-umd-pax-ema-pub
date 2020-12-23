####################################################################################################
# Study 1 Modeling Script: Generating Imputation Data Set - Follow-up

# Description: 
#   These analyses were performed for completeness only. To bring into alignment with Study 2 measures,
#   this script was produced after the main set of analyses was performed. This set of analyses can be
#   best thought of as a sensitivity analysis. 
####################################################################################################

#---------------------------------------------------------------------------------------------------
# Package import
library(pan)
library(mitml)
library(mice)
library(tidyverse)
library(corrgram)
library(glue)
library(magrittr)
library(dplyr)
#---------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------
# Folder setup
wd <- '~/github/ATNL/shackman-umd-pax-ema-pub'
data.folder <- '{wd}/Data' %>% glue()
study1.model <- '{wd}/Study1_model_summaries' %>% glue()
#---------------------------------------------------------------------------------------------------
load('{data.folder}/study1_data.RData' %>% glue())

dat.study1_model <- dat.study1_model %>% 
  mutate(
    NegEvnt_dich = ifelse(NegEvnt == 1, "No", "Yes"),
    NegEvnt_dich = forcats::fct_relevel(NegEvnt_dich, "No", "Yes"),
    PosEvnt_dich = ifelse(PosEvnt == 1, "No", "Yes"),
    PosEvnt_dich = forcats::fct_relevel(PosEvnt_dich, "No", "Yes"), 
    Chrfl_f = as.factor(Chrfl),
    Hppy_f = as.factor(Hppy),
    Jyfl_f = as.factor(Jyfl),
    Nrvs_f = as.factor(Nrvs),
    Anxs_f = as.factor(Anxs),
    Unesy_f = as.factor(Unesy)
  )

# The modeling syntax uses all summary data from the EMA mood items to generate posterior 
# distributions for individual missing values with random intercepts allowing for individual 
# differences in posterior distributions. Note that this imputation can take a while... 
fml <- NegEvnt_dich + PosEvnt_dich + Chrfl_f + Hppy_f + Jyfl_f + Nrvs_f + Anxs_f + Unesy_f ~ c.DN + m.NegEvnt + 
  m.PosEvnt + sd.NegEvnt + sd.PosEvnt + m.NEG + m.POS + sd.NEG + sd.POS + (1|ID) 

M <- 10 
imp <- jomoImpute(dat.study1_model, 
                  formula = fml, 
                  n.burn = 10, 
                  n.iter = 250, 
                  m = M, 
                  keep.chains = 'diagonal',
                  seed = 20201108)

# plot(imp)
dat.imp <- mitmlComplete(imp)
dat.study1_list <- list() 

IDs <- unique(dat.study1_model$ID)

for(i in 1:M){
  dat.study1_list[[i]] <- dat.imp[[i]]
  dat.study1_list[[i]]$NegEvnt_dich <- ifelse(dat.study1_list[[i]]$NegEvnt_dich == 'No', 0, 1)
  dat.study1_list[[i]]$PosEvnt_dich <- ifelse(dat.study1_list[[i]]$PosEvnt_dich == 'No', 0, 1)

  for(j in 1:length(IDs)){
    # Force imputed y-vals to remain within range
    dat.study1_list[[i]]$NEG[dat.study1_list[[i]]$ID == IDs[j]] <- ifelse(
      dat.study1_list[[i]]$NEG[dat.study1_list[[i]]$ID == IDs[j]] > 5, 
      sample(3:5, prob = c(sum()))
    )
    dat.study1_list[[i]]$prop.NegEvnt[dat.study1_list[[i]]$ID == IDs[j]] <- mean(dat.study1_list[[i]]$NegEvnt_dich[dat.study1_list[[i]]$ID == IDs[j]])
    dat.study1_list[[i]]$prop.PosEvnt[dat.study1_list[[i]]$ID == IDs[j]] <- mean(dat.study1_list[[i]]$PosEvnt_dich[dat.study1_list[[i]]$ID == IDs[j]])
  }
  dat.study1_list[[i]]$c.NegEvnt <- dat.study1_list[[i]]$NegEvnt_dich - dat.study1_list[[i]]$prop.NegEvnt
  dat.study1_list[[i]]$c.PosEvnt <- dat.study1_list[[i]]$PosEvnt_dich - dat.study1_list[[i]]$prop.PosEvnt
}

save(list = c('M', ls()[grepl('dat.study1_', ls())]), 
     file = '{data.folder}/study1_data.RData' %>% glue())
