####################################################################################################
# Study 1 Modeling Script: Generating Imputation Data Set - Follow-up

# Description: 
#   These analyses were performed for completeness only. To bring into alignment with Study 2 measures,
#   this script was produced after the main set of analyses was performed. This set of analyses can be
#   best thought of as a sensitivity analysis. 
####################################################################################################

#---------------------------------------------------------------------------------------------------
# Package import
library(glue)
library(mitml)
library(tidyverse)
library(parallel)

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
    NegEvnt_f = as.factor(NegEvnt),
    PosEvnt_f = as.factor(PosEvnt),
    Chrfl_f = as.factor(Chrfl),
    Hppy_f = as.factor(Hppy),
    Jyfl_f = as.factor(Jyfl),
    Nrvs_f = as.factor(Nrvs),
    Anxs_f = as.factor(Anxs),
    Unesy_f = as.factor(Unesy)
  )

# basic checks - should see perfect alignment on diagnoal if recodes wored as expected
new_factors <- names(dat.study1_model)[grep('_f', names(dat.study1_model))]
original_vars <- str_replace_all(new_factors, "_f", "")

for(i in seq_along(new_factors)) {
  table(dat.study1_model[[new_factors[i]]], dat.study1_model[[original_vars[i]]]) %>% 
    print()
}

# The modeling syntax uses all summary data from the EMA mood items to generate posterior 
# distributions for individual missing values with random intercepts allowing for individual 
# differences in posterior distributions. Note that this imputation can take a while... 
fml <- NegEvnt_f + PosEvnt_f + Chrfl_f + Hppy_f + Jyfl_f + Nrvs_f + Anxs_f + Unesy_f ~ c.DN + m.NegEvnt + 
  m.PosEvnt + sd.NegEvnt + sd.PosEvnt + m.NEG + m.POS + sd.NEG + sd.POS + (1|ID) 

M <- 10

impute_wrapper <- function(i) {
    imp_df <- jomoImpute(
        dat.study1_model, 
        formula = fml, 
        n.burn = 5000, 
        n.iter = 199, # because it is prime and I like primes :) 
        m = 1, 
        keep.chains = 'diagonal'
    )
    
    save(list = c('imp_df'),
         file = '{data.folder}/study1_data_imp_{i}.RData' %>% glue())
}

mclapply(1:M, mc.cores = M, FUN = impute_wrapper)

print('Imputation complete, checking saved file status')

for(m in 1:M){
	print("Verifying file exists")
	file_path <- '{data.folder}/study1_data_imp_{m}.RData' %>% glue()
	print('{file_path} exists? {file.exists(file_path)}' %>% glue())
}
