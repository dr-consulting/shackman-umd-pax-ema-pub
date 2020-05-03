source("~/dr-consulting_GH/shackman-umd-pax-ema-pub/Graphics_scripts/bayes_variance_riverplots_utils.R")

POSTERIOR_PATH <- "/media/matthew/My Book/EMA_S1_Bayesian_Posteriors"
DATA_FILEPATH <- "~/dr-consulting_GH/shackman-umd-pax-ema-pub/Data/study1_data.RData"
SUMMARY_DIRPATH <- "~/dr-consulting_GH/shackman-umd-pax-ema-pub/Study_1_model_summaries/"
data_loader(POSTERIOR_PATH, "S1_POS_NegEvnt_Rct", DATA_FILEPATH, "S1_POS_ucm")

# Create interaction variable in dat.study2_list (may make a small helper for this eventually)
dat.study1_model["c.NegEvnt:c.DN"] <- dat.study1_model["c.NegEvnt"]*dat.study1_model["c.DN"]

#----------------------------------------------------------------------------------------------------------------------
summary_filepath <- "~/dr-consulting_GH/shackman-umd-pax-ema-pub/Study_1_model_summaries/var_decomp_S1_POS_NegEvnt_x_DN.txt"
within_vars <- c("c.NegEvnt", "c.NegEvnt:c.DN")
between_vars <- c("c.DN", "m.NegEvnt")
random_vars <- c("c.NegEvnt")


posterior_df <- posterior_samples_extractor(S1_POS_ucm, S1_POS_NegEvnt_Rct, resp="POS")

full_model_decomp <- posterior_r2mlm_draws(dat.study1_model, posterior_df, between_vars, within_vars, random_vars, 
                                           has_intercept=TRUE, clustermeancentered=TRUE)

sink(summary_filepath)
print(psych::describe(full_model_decomp, skew = FALSE, quant=c(.025, .975)), digits=4)
sink()

# Just some memory saving...
remove(list=c("S1_POS_NegEvnt_Rct"))
gc()

#----------------------------------------------------------------------------------------------------------------------
# Next Set of inputs - lv2 no interaction
data_loader(POSTERIOR_PATH, "S1_POS_NegEvnt_DN_Exp")
summary_filepath <- "~/dr-consulting_GH/shackman-umd-pax-ema-pub/Study_1_model_summaries/var_decomp_S1_POS_NegEvnt_DN_Exp.txt"
within_vars <- c("c.NegEvnt")
between_vars <- c("c.DN", "m.NegEvnt")
random_vars <- c("c.NegEvnt")

posterior_df <- posterior_samples_extractor(S1_POS_ucm, S1_POS_NegEvnt_DN_Exp, resp="POS")

lv2_Exp_DN_decomp <- posterior_r2mlm_draws(dat.study1_model, posterior_df, between_vars, within_vars, random_vars, 
                                           has_intercept=TRUE, clustermeancentered=TRUE)
sink(summary_filepath)
print(psych::describe(lv2_Exp_DN_decomp, skew = FALSE, quant=c(.025, .975)), digits=4)
sink()

# Just some memory saving...
remove(list=c("S1_POS_NegEvnt_DN_Exp"))
gc()

#----------------------------------------------------------------------------------------------------------------------
# Next Set of inputs - only DN at level 2 of the equation
data_loader(POSTERIOR_PATH, "S1_POS_NegEvnt_DN")
summary_filepath <- "~/dr-consulting_GH/shackman-umd-pax-ema-pub/Study_1_model_summaries/var_decomp_S1_POS_NegEvnt_DN.txt"
within_vars <- c("c.NegEvnt")
between_vars <- c("c.DN")
random_vars <- c("c.NegEvnt")

posterior_df <- posterior_samples_extractor(S1_POS_ucm, S1_POS_NegEvnt_DN, resp="POS")

lv2_DN_decomp <- posterior_r2mlm_draws(dat.study1_model, posterior_df, between_vars, within_vars, random_vars, 
                                       has_intercept=TRUE, clustermeancentered=TRUE)
sink(summary_filepath)
print(psych::describe(lv2_DN_decomp, skew = FALSE, quant=c(.025, .975)), digits=4)
sink()

# Just some memory saving...
remove(list=c("S1_POS_NegEvnt_DN"))
gc()

#----------------------------------------------------------------------------------------------------------------------
# Next Set of inputs - lv2 no interaction
data_loader(POSTERIOR_PATH, "S1_POS_NegEvnt_Exp")
summary_filepath <- "~/dr-consulting_GH/shackman-umd-pax-ema-pub/Study_1_model_summaries/var_decomp_S1_POS_NegEvnt_Exp.txt"
within_vars <- c("c.NegEvnt")
between_vars <- c("m.NegEvnt")
random_vars <- c("c.NegEvnt")

posterior_df <- posterior_samples_extractor(S1_POS_ucm, S1_POS_NegEvnt_Exp, resp="POS")

lv2_Exp_decomp <- posterior_r2mlm_draws(dat.study1_model, posterior_df, between_vars, within_vars, random_vars, 
                                        has_intercept=TRUE, clustermeancentered=TRUE)
sink(summary_filepath)
print(psych::describe(lv2_Exp_decomp, skew = FALSE, quant=c(.025, .975)), digits=4)
sink()

# Just some memory saving...
remove(list=c("S1_POS_NegEvnt_Exp"))
gc()
#----------------------------------------------------------------------------------------------------------------------
# Next Set of inputs - lv1 only - no level 2 predictors of any kind
data_loader(POSTERIOR_PATH, "S1_POS_NegEvnt")
summary_filepath <- "~/dr-consulting_GH/shackman-umd-pax-ema-pub/Study_1_model_summaries/var_decomp_S1_POS_NegEvnt.txt"
within_vars <- c("c.NegEvnt")
between_vars <- NULL
random_vars <- c("c.NegEvnt")

posterior_df <- posterior_samples_extractor(S1_POS_ucm, S1_POS_NegEvnt, resp="POS")

lv1_only_decomp <- posterior_r2mlm_draws(dat.study1_model, posterior_df, between_vars, within_vars, random_vars, 
                                           has_intercept=TRUE, clustermeancentered=TRUE)
sink(summary_filepath)
print(psych::describe(lv1_only_decomp, skew = FALSE, quant=c(.025, .975)), digits=4)
sink()

# Just some memory saving...
remove(list=c("S1_POS_NegEvnt"))
gc()
#----------------------------------------------------------------------------------------------------------------------
# Within Formulas: 
event <- "mean(lv1_only[['tot_fix_wthn']])"
reactivity <- "mean(lv2_Exp_DN[['tot_slp_varn']]) + mean(lv2_Exp_DN[['tot_sig_varn']]) - mean(full_model[['tot_slp_varn']]) - mean(full_model[['tot_sig_varn']])"
unmodeled <- "mean(full_model[['tot_slp_varn']]) + mean(full_model[['tot_sig_varn']]) + mean(full_model[['tot_fix_wthn']]) - sum(within_decomp[1:2])"

# Between Formulas: 
tonic_DN <- "mean(lv2_Exp[['tot_int_varn']]) - mean(lv2_Exp_DN[['tot_int_varn']])"
exposure <- "mean(lv2_DN[['tot_int_varn']]) - mean(lv2_Exp_DN[['tot_int_varn']])"
shared_DN_exp <- "mean(lv1_only[['tot_int_varn']]) - mean(lv2_Exp_DN[['tot_int_varn']]) - sum(between_decomp[1:2])"
# Alright what I don't love about this is that it has to be labeled between_decomp under the hood 
# Fairly brittle approach here but going to live with it for now. 
# Could re-visit if I turn this into a more complete plotting package
unmodeled_btwn <- "mean(full_model[['tot_int_varn']]) + mean(full_model[['tot_fix_btwn']]) - sum(between_decomp[1:3])"

within_contrasts <- c("Negative \n Event" = event, 
                      "Reactivity \n" = reactivity, 
                      "Unmod. \n Within" = unmodeled)

between_contrasts <- c("Tonic \n DN" = tonic_DN, 
                       "Negative \n Event \n Exp." = exposure,
                       "DN \n Shared w/ \n Exp." = shared_DN_exp,
                       "Unmod. \n Between" = unmodeled_btwn)

custom_contrasts <- c("Tonic \n DN" = tonic_DN, 
                      "DN \n Shared w/ \n Exp." = shared_DN_exp, 
                      "Reactivity \n" = reactivity)

model_names <- c("full_model", "lv1_only", "lv2_DN", "lv2_Exp", "lv2_Exp_DN")

model_variance_list <- list()
model_variance_list[[1]] <- full_model_decomp
model_variance_list[[2]] <- lv1_only_decomp
model_variance_list[[3]] <- lv2_DN_decomp
model_variance_list[[4]] <- lv2_Exp_decomp
model_variance_list[[5]] <- lv2_Exp_DN_decomp

riverplot_df_helper(model_variance_list, model_names, within_contrasts, between_contrasts, custom_contrasts, 
                    within_color = RColorBrewer::brewer.pal(9, "Reds")[5], 
                    between_color = RColorBrewer::brewer.pal(9, "Blues")[5], 
                    merge_color = RColorBrewer::brewer.pal(9, "Purples")[5], 
                    custom_contrast_name = "Combined \n DN Effect", main_filename = "~/S1_POS_NegEvnt_decomp.eps", 
                    main_title = "Total Variance Decomposition: Positive Mood, DN, and Negative Events", 
                    custom_filename = "~/S1_POS_NegEvnt_DN_combined.eps", 
                    custom_title = "DN Variance Decomposition: Positive Mood, DN, and Negative Events", 
                    combined_plot_filename = "~/S1_POS_NegEvnt_full_decomp.eps")
