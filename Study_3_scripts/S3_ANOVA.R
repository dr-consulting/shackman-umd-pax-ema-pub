library(brms)

data_folder <- "~/dr-consulting_GH/shackman-umd-pax-ema-pub/Data"

load(paste0(data_folder, "/study3_data.RData"))

anova_form <- bf(
	rating ~ certainty*valence*DN + (1|ID)
) + gaussian()

fit <- brm(anova_form, data = dat.study3, 
           prior = c(set_prior("normal(0, 5)", class = "b")), cores=3, 
           chains = 3, iter = 20000, warmup = 15000, 
           control = list(adapt_delta = .99, 
                          max_treedepth = 15)
)

# Initial data checks suggest Gaussian prior may not be appropriate
# Will cross-validate distributional assumptions of the model
loo_fit <- loo(fit, reloo=TRUE)

#######################################################################################################################
# Gaussian distribution is not a great fit of the data. Looking for a model that 
# better predicts the posterior distribution.

anova_form_log <- bf(
	rating ~ certainty*valence*DN+ (1|ID)
) + lognormal()

fit_log <- brm(anova_form_log, data = dat.study3, 
               prior = c(set_prior("normal(0, 2)", class = "b")), cores=3, 
               chains = 3, iter = 20000, warmup = 15000, 
               control = list(adapt_delta = .99, 
                              max_treedepth = 15)
)

# Obtaining loo information criteria for cross-validation of model posterior
loo_fit_log <- loo(fit_log, reloo=TRUE)

# "Preferred" model is in the top row of output
# Mutliply values in elpd_diff column by -2 to place on deviance scale 
loo_compare(loo_fit, loo_fit_log)

