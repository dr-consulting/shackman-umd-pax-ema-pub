library(brms)
library(rstan)
library(rstanarm)
library(bayesplot)
library(tidyverse)
library(glue)

rstan_options(auto_write=TRUE)
options(mc.cores = parallel::detectCores())

wd <- '~/github/ATNL/shackman-umd-pax-ema-pub'
data.folder <- '{wd}/Data' %>% glue()
study3.model <- '{wd}/Study_3_model_summaries' %>% glue()

posterior_save_dir <- '/media/dr-owner/HDD1'
study3.out <- '{posterior_save_dir}/EMA_S3_Bayesian_Posteriors/gamma' %>% glue()
study3.graphics <- '{study3.out}/diagnostic_plots' %>% glue()

source('{wd}/utils.R' %>% glue())

# Code was used to bring in updated data for the project - SG performed more recent cleaning of these data
n220_sg_spss <- haven::read_spss('~/Desktop/grogans_pax220_ratings_032521.sav')

n220_model_df <- n220_sg_spss %>%
    select(StudyID, UT_AvgRating, PT_AvgRating, US_AvgRating, PS_AvgRating, ZUSE_THIS_DN) %>%
    mutate(
        DN = ZUSE_THIS_DN,
        ID = StudyID,
    ) %>%
    select(-ZUSE_THIS_DN, -StudyID) %>%
    pivot_longer(-all_of(c('ID', 'DN')), values_to = 'rating', names_to = 'cond') %>%
    mutate(
        valence = ifelse(str_detect(cond, 'T'), 'Threat', 'Safety'),
        valence = forcats::fct_relevel(valence, 'Safety', 'Threat'),
		valence_dich = ifelse(valence == 'Threat', 1, 0),
		c_valence = valence_dich - .5,
        certainty = ifelse(str_detect(cond, 'U'), 'Uncertain', 'Certain'),
        certainty = forcats::fct_relevel(certainty, 'Certain', 'Uncertain'),
        cond = str_split(cond, '_') %>% map_chr(., 1)
    )

dat.study3 <- n220_model_df
remove(list = c('n220_model_df', 'n220_sg_spss'))

data_folder <- "~/github/ATNL/shackman-umd-pax-ema-pub/Data"

###########################################################################################
prior_config <- c(set_prior('lognormal(3, .5)', class='shape'))

S3_ucm_form <- bf(
	rating ~ 1 + (1|ID)
) + Gamma(link='log')

S3_ucm <- brm(
	S3_ucm_form,
	data=dat.study3,
	chains=3, 
	cores=3, 
	iter=15000, 
	warmup=10000, 
	control=list(adapt_delta=.99, 
				 max_treedepth=15),
	seed=8228,
	save_all_pars=TRUE, 
	save_model='S3_ucm',
	open_progress=FALSE, 
	refresh=0, 
	prior=prior_config
)

create_diagnostic_plots(S3_ucm, list(dat.study3), n_samples=100, dir_path=study3.graphics)

sink(paste0(study3.model, '/S3_ucm.txt'))
print(summary(S3_ucm), digits=5)
sink()

save(list=c('S3_ucm', 'S3_ucm_form'), 
     file=paste0(study3.out, '/S3_ucm.RData'))
remove(S3_ucm)
gc()

#######################################################################################
prior_config <- c(prior_config, 
				  set_prior('normal(0, 2)', class='b'))

S3_valence_form <- bf(
	rating ~ 1 + c_valence + (1 + c_valence|ID)
) + Gamma(link='log')

S3_valence <- brm(
	S3_valence_form,
	data=dat.study3,
	chains=3, 
	cores=3, 
	iter=15000, 
	warmup=10000, 
	control=list(adapt_delta=.99, 
				 max_treedepth=15),
	seed=8228,
	save_all_pars=TRUE, 
	save_model='S3_valence',
	open_progress=FALSE, 
	refresh=0, 
	prior=prior_config
)

create_diagnostic_plots(S3_valence, list(dat.study3), n_samples=100, dir_path=study3.graphics)

sink(paste0(study3.model, '/S3_valence.txt'))
print(summary(S3_valence), digits=5)
sink()

save(list=c('S3_valence', 'S3_valence_form'), 
     file=paste0(study3.out, '/S3_valence.RData'))
remove(S3_valence)
gc()

#######################################################################################
S3_DN_form <- bf(
	rating ~ 1 + c_valence + DN + (1 + c_valence|ID)
) + Gamma(link='log')

S3_DN <- brm(
	S3_DN_form,
	data=dat.study3,
	chains=3, 
	cores=3, 
	iter=15000, 
	warmup=10000, 
	control=list(adapt_delta=.99, 
				 max_treedepth=15),
	seed=8228,
	save_all_pars=TRUE, 
	save_model='S3_DN',
	open_progress=FALSE, 
	refresh=0, 
	prior=prior_config
)

create_diagnostic_plots(S3_DN, list(dat.study3), n_samples=100, dir_path=study3.graphics)

sink(paste0(study3.model, '/S3_DN.txt'))
print(summary(S3_DN), digits=5)
sink()

save(list=c('S3_DN', 'S3_DN_form'), 
     file=paste0(study3.out, '/S3_DN.RData'))
remove(S3_DN)
gc()

#######################################################################################
S3_DN_x_valence_form <- bf(
	rating ~ 1 + c_valence * DN + (1 + c_valence|ID)
) + Gamma(link='log')

S3_DN_x_valence <- brm(
	S3_DN_x_valence_form,
	data=dat.study3,
	chains=3, 
	cores=3, 
	iter=15000, 
	warmup=10000, 
	control=list(adapt_delta=.99, 
				 max_treedepth=15),
	seed=8228,
	save_all_pars=TRUE, 
	save_model='S3_DN_x_valence',
	open_progress=FALSE, 
	refresh=0, 
	prior=prior_config
)

create_diagnostic_plots(S3_DN_x_valence, list(dat.study3), n_samples=100, dir_path=study3.graphics)

sink(paste0(study3.model, '/S3_DN_x_valence.txt'))
print(summary(S3_DN_x_valence), digits=5)
sink()

save(list=c('S3_DN_x_valence', 'S3_DN_x_valence_form'), 
     file=paste0(study3.out, '/S3_DN_x_valence.RData'))
remove(S3_DN_x_valence)
gc()

