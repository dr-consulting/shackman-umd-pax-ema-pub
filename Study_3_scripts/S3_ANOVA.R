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

# Uncomment to save
save.image(file = "~/github/ATNL/shackman-umd-pax-ema-pub/Data/study3_data.RData")
