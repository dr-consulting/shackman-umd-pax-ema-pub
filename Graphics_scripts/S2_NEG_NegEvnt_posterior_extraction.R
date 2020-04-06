source("~/dr-consulting_GH/shackman-umd-pax-ema-pub/Graphics_scripts/bayes_variance_riverplots_utils.R")

POSTERIOR_PATH <- "/media/matthew/My Book/EMA_S2_Bayesian_Posteriors"
DATA_FILEPATH <- "~/dr-consulting_GH/shackman-umd-pax-ema-pub/Data/study2_data.RData"

data_loader(POSTERIOR_PATH, "S2_NEG_NegEvnt_x_DN_prop.NegEvnt", DATA_FILEPATH, "S2_NEG_ucm")

# Create interaction variable in dat.study2_list (may make a small helper for this eventually)
for(m in 1:length(dat.study2_list)){
  dat.study2_list[[m]][paste("c.NegEvnt", "c.DN", sep = ":")] <- dat.study2_list[[m]]["c.NegEvnt"]*
    dat.study2_list[[m]]["c.DN"]
}

#----------------------------------------------------------------------------------------------------------------------
within_vars <- c("c.NegEvnt", "c.NegEvnt:c.DN")
between_vars <- c("c.DN", "prop.NegEvnt")
random_vars <- c("c.NegEvnt")

full_model_decomp <- r2MLM_brms_wrapper(dat.study2_list, within_vars, between_vars, random_vars,
                                        focal_model = S2_NEG_NegEvnt_x_DN_prop.NegEvnt, null_model=S2_NEG_ucm, 
                                        has_intercept = TRUE, clustermeancentered = TRUE, link_func = "log")

# Just some memory saving...
remove(list=c("S2_NEG_NegEvnt_x_DN_prop.NegEvnt"))
gc()

#----------------------------------------------------------------------------------------------------------------------
# Next Set of inputs - lv2 no interaction
data_loader(POSTERIOR_PATH, "S2_NEG_NegEvnt_DN_prop.NegEvnt")
within_vars <- c("c.NegEvnt")
between_vars <- c("c.DN", "prop.NegEvnt")
random_vars <- c("c.NegEvnt")

lv2_Exp_DN_decomp <- r2MLM_brms_wrapper(dat.study2_list, within_vars, between_vars, random_vars,
                                        focal_model = S2_NEG_NegEvnt_x_DN_prop.NegEvnt, null_model=S2_NEG_ucm, 
                                        has_intercept = TRUE, clustermeancentered = TRUE, link_func = "log")


# Just some memory saving...
remove(list=c("S2_NEG_NegEvnt_DN_prop.NegEvnt"))
gc()

#----------------------------------------------------------------------------------------------------------------------
# Next Set of inputs - lv2 no interaction
data_loader(POSTERIOR_PATH, "S2_NEG_NegEvnt_DN")
within_vars <- c("c.NegEvnt")
between_vars <- c("c.DN")
random_vars <- c("c.NegEvnt")

lv2_DN_decomp <- r2MLM_brms_wrapper(dat.study2_list, within_vars, between_vars, random_vars,
                                    focal_model = S2_NEG_NegEvnt_x_DN_prop.NegEvnt, null_model=S2_NEG_ucm, 
                                    has_intercept = TRUE, clustermeancentered = TRUE, link_func = "log")


# Just some memory saving...
remove(list=c("S2_NEG_NegEvnt_DN"))
gc()

#----------------------------------------------------------------------------------------------------------------------
# Next Set of inputs - lv2 no interaction
data_loader(POSTERIOR_PATH, "S2_NEG_NegEvnt_prop.NegEvnt")
within_vars <- c("c.NegEvnt")
between_vars <- c("prop.NegEvnt")
random_vars <- c("c.NegEvnt")

lv2_Exp_decomp <- r2MLM_brms_wrapper(dat.study2_list, within_vars, between_vars, random_vars,
                                     focal_model = S2_NEG_NegEvnt_x_DN_prop.NegEvnt, null_model=S2_NEG_ucm, 
                                     has_intercept = TRUE, clustermeancentered = TRUE, link_func = "log")

# Just some memory saving...
remove(list=c("S2_NEG_NegEvnt_prop.NegEvnt"))
gc()

#----------------------------------------------------------------------------------------------------------------------
# Next Set of inputs - lv2 no interaction
data_loader(POSTERIOR_PATH, "S2_NEG_NegEvnt")
within_vars <- c("c.NegEvnt")
between_vars <- NULL
random_vars <- c("c.NegEvnt")

lv1_only_decomp <- r2MLM_brms_wrapper(dat.study2_list, within_vars, between_vars, random_vars,
                                      focal_model = S2_NEG_NegEvnt_x_DN_prop.NegEvnt, null_model=S2_NEG_ucm, 
                                      has_intercept = TRUE, clustermeancentered = TRUE, link_func = "log")
