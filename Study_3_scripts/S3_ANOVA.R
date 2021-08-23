library(tidyverse)
library(brms)
library(gganimate)
library(tidybayes)

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
        certainty = ifelse(str_detect(cond, 'U'), 'Uncertain', 'Certain'),
        certainty = forcats::fct_relevel(certainty, 'Certain', 'Uncertain'),
        cond = str_split(cond, '_') %>% map_chr(., 1)
    )

dat.study3 <- n220_model_df
remove(list = c('n220_model_df', 'n220_sg_spss'))

data_folder <- "~/github/ATNL/shackman-umd-pax-ema-pub/Data"
# load(paste0(data_folder, "/study3_data.RData"))

# Step 1 - Define general MLM formula
anova_form <- bf(
    rating ~ certainty*valence*DN + (1|ID)
)

fit_gaussian <- brm(
    anova_form + gaussian(), 
    data = dat.study3,
    prior = c(set_prior("normal(0, 2)", class = "b")), 
    cores = 3, 
    chains = 3, 
    iter = 20000, 
    warmup = 15000, 
    control = list(adapt_delta = .99, 
                   max_treedepth = 15)
)

loo_fit_gaussian <- loo(fit_gaussian, reloo=TRUE)

# Step 2 fit with a weakly informative Gamma prior on the likelihood
fit_gamma <- brm(
    anova_form + Gamma(link = 'log'), 
    data = dat.study3,
    prior = c(set_prior("normal(0, 2)", class = "b")), 
    cores = 3, 
    chains = 3, 
    iter = 20000, 
    warmup = 15000, 
    control = list(adapt_delta = .99, 
                   max_treedepth = 15)
)

loo_fit_gamma <- loo(fit_gamma, reloo=TRUE)

# Create direct comparison of models 
loocv_gaus_vs_gamma <- loo_compare(loo_fit_gaussian, loo_fit_gamma)
print(loocv_gaus_vs_gamma)

# sink('~/github/ATNL/shackman-umd-pax-ema-pub/Study_3_model_summaries/fit_gamma_3way.txt')
# cat('Guassian vs. Gamma LOOCV\n')
# cat('--------------------------\n')
# cat(print(loocv_gaus_vs_gamma %>% as.data.frame()))
# cat('----------------------------\n')
# Uncomment to save
# save.image(file = "~/github/ATNL/shackman-umd-pax-ema-pub/Data/study3_data.RData")

###########################################################################################
# remove 3-way interaction term 
anova_form_all_2ways <- bf(
    rating ~ certainty*valence + certainty*DN + DN*valence + (1|ID)
)

fit_gamma_all_2ways <- brm(
    anova_form_all_2ways + Gamma(link = 'log'), 
    data = dat.study3,
    prior = c(set_prior("normal(0, 2)", class = "b")), 
    cores = 3, 
    chains = 3, 
    iter = 20000, 
    warmup = 15000, 
    control = list(adapt_delta = .99, 
                   max_treedepth = 15)
)

loo_fit_gamma_all_2ways <- loo(fit_gamma_all_2ways, reloo=TRUE)
loocv_3way_vs_2way <- loo_compare(loo_fit_gamma_all_2ways, loo_fit_gamma)
print(loocv_3way_vs_2way)

# sink('~/github/ATNL/shackman-umd-pax-ema-pub/Study_3_model_summaries/fit_gamma_3way.txt')
# cat('Guassian vs. Gamma LOOCV\n')
# cat('--------------------------\n')
# cat(print(loocv_3way_vs_2way %>% as.data.frame()))
# cat('----------------------------\n')
# cat(print(summary(fit_gamma_all_2ways)))
# sink()

#####################################################
# drop valence_DN
anova_form_drop_valence_DN <- bf(
    rating ~ certainty*valence + certainty*DN + (1|ID)
)

fit_gamma_drop_valence_DN <- brm(
    anova_form_drop_valence_DN + Gamma(link = 'log'), 
    data = dat.study3,
    prior = c(set_prior("normal(0, 2)", class = "b")), 
    cores = 3, 
    chains = 3, 
    iter = 20000, 
    warmup = 15000, 
    control = list(adapt_delta = .99, 
                   max_treedepth = 15)
)

loo_fit_gamma_drop_valence_DN <- loo(fit_gamma_drop_valence_DN, reloo=TRUE)

#####################################################################
# drop certainty x DN
anova_form_drop_certainty_DN <- bf(
    rating ~ certainty*valence + DN*valence + (1|ID)
)

fit_gamma_drop_certainty_DN <- brm(
    anova_form_drop_certainty_DN + Gamma(link = 'log'), 
    data = dat.study3,
    prior = c(set_prior("normal(0, 2)", class = "b")), 
    cores = 3, 
    chains = 3, 
    iter = 20000, 
    warmup = 15000, 
    control = list(adapt_delta = .99, 
                   max_treedepth = 15)
)

loo_fit_gamma_drop_certainty_DN <- loo(fit_gamma_drop_certainty_DN, reloo=TRUE)

#####################################################################
# drop certainty x valence
anova_form_drop_certainty_valence <- bf(
    rating ~ certainty*DN + DN*valence + (1|ID)
)

fit_gamma_drop_certainty_valence <- brm(
    anova_form_drop_certainty_valence + Gamma(link = 'log'), 
    data = dat.study3,
    prior = c(set_prior("normal(0, 2)", class = "b")), 
    cores = 3, 
    chains = 3, 
    iter = 20000, 
    warmup = 15000, 
    control = list(adapt_delta = .99, 
                   max_treedepth = 15)
)
loo_fit_gamma_drop_certainty_valence <- loo(fit_gamma_drop_certainty_valence, reloo=TRUE)

#######################################################################################################################
# Only DN x Valence 
anova_form_DN_valence <- bf(
    rating ~ certainty + DN*valence + (1|ID)
)

fit_gamma_DN_valence <- brm(
    anova_form_DN_valence + Gamma(link = 'log'), 
    data = dat.study3,
    prior = c(set_prior("normal(0, 2)", class = "b")), 
    cores = 3, 
    chains = 3, 
    iter = 20000, 
    warmup = 15000, 
    control = list(adapt_delta = .99, 
                   max_treedepth = 15)
)

loo_fit_gamma_DN_valence <- loo(fit_gamma_DN_valence, reloo=TRUE)



########################################################################################################################
# Only DN x certainty
anova_form_DN_certainty <- bf(
    rating ~ valence + DN*certainty + (1|ID)
)

fit_gamma_DN_certainty <- brm(
    anova_form_DN_certainty + Gamma(link = 'log'), 
    data = dat.study3,
    prior = c(set_prior("normal(0, 2)", class = "b")), 
    cores = 3, 
    chains = 3, 
    iter = 20000, 
    warmup = 15000, 
    control = list(adapt_delta = .99, 
                   max_treedepth = 15)
)

loo_fit_gamma_DN_certainty <- loo(fit_gamma_DN_certainty, reloo=TRUE)
#######################################################################################################################
# Only main effects
anova_form_main <- bf(
    rating ~ valence + DN + certainty + (1|ID)
)

fit_gamma_main <- brm(
    anova_form_main + Gamma(link = 'log'), 
    data = dat.study3,
    prior = c(set_prior("normal(0, 2)", class = "b")), 
    cores = 3, 
    chains = 3, 
    iter = 20000, 
    warmup = 15000, 
    control = list(adapt_delta = .99, 
                   max_treedepth = 15)
)

loo_fit_gamma_main <- loo(fit_gamma_main, reloo=TRUE)

#######################################################################################################################
loocv_all <- loo_compare(loo_fit_gamma_all_2ways, loo_fit_gamma_drop_certainty_valence, loo_fit_gamma, 
                         loo_fit_gamma_drop_certainty_DN, loo_fit_gamma_drop_valence_DN, 
                         loo_fit_gamma_DN_certainty, loo_fit_gamma_DN_valence, 
                         loo_fit_gamma_main)

loocv_all

# sink('~/github/ATNL/shackman-umd-pax-ema-pub/Study_3_model_summaries/fit_gamma_final.txt')
# cat('Guassian vs. Gamma LOOCV\n')
# cat('--------------------------\n')
# cat(print(loocv_3way_vs_2way %>% as.data.frame()))
# cat('----------------------------\n')
# cat('Summary Gamma Model Comparisons LOOCV\n')
# cat('--------------------------\n')
# cat(print(loocv_all %>% as.data.frame()))
# cat('----------------------------\n')
# print(summary(fit_gamma_main))
# sink()

summary(fit_gamma_DN_certainty)
summary(fit_gamma_DN_valence)
