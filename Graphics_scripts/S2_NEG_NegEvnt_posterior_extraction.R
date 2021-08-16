source("~/github/ATNL/shackman-umd-pax-ema-pub/Graphics_scripts/bayes_variance_riverplots_utils.R")
library(glue)
library(plyr)
library(tidyverse)

# To ensure stable output in text files up to reasonable "wide-screen" width
options(width = 140)

BASE_DIR <- "/media/dr-owner/HDD1/EMA_S2_Bayesian_Posteriors"
POSTERIOR_DIR <- "{BASE_DIR}/gamma" %>% glue()
DATA_FILEPATH <- "~/github/ATNL/shackman-umd-pax-ema-pub/Data/study2_data.RData"

SUMMARY_DIR <- "{BASE_DIR}/gamma/var_decomp_txt" %>% glue()
if(!dir.exists(SUMMARY_DIR)) dir.create(SUMMARY_DIR, recursive = TRUE)

GRAPHICS_DIR <- "{BASE_DIR}/gamma/riverplots" %>% glue()
if(!dir.exists(GRAPHICS_DIR)) dir.create(GRAPHICS_DIR, recursive = TRUE)

LOAD_PATH <- "{GRAPHICS_DIR}/S2_NEG_NegEvnt_river_base.RData" %>% glue()
EXECUTE_FLAG <- !file.exists(LOAD_PATH)

data_loader(POSTERIOR_DIR, "S2_NEG_NegEvnt_x_DN_prop.NegEvnt", DATA_FILEPATH, "S2_NEG_ucm")

if(EXECUTE_FLAG){
    #----------------------------------------------------------------------------------------------------------------------
    # Pre-processing
    # Create imputed mean data.frame
    # Averages any imputed missing values and returns a single summary data.frame of all required variables
    model_vars <- c('c.DN', 'c.NegEvnt', 'prop.NegEvnt', 'NEG')
    id_vector <- dat.study2_list[[1]][['ID']]
    
    # Extract only the required variables
    for(m in 1:length(dat.study2_list)) {
        dat.study2_list[[m]] <- dat.study2_list[[m]][, model_vars]
    }

    # Calculate element-wise means
    mean_imputed_df <- Reduce(`+`, dat.study2_list) / length(dat.study2_list)
    
    # Add back in the ID variable: 
    mean_imputed_df[['ID']] <- id_vector
    
    # Create a negative event interaction term with the correct label: 
    mean_imputed_df[['c.NegEvnt:c.DN']] <- mean_imputed_df[['c.NegEvnt']] * mean_imputed_df[['c.DN']]
    
    #----------------------------------------------------------------------------------------------------------------------
    summary_filepath <- "{SUMMARY_DIR}/var_decomp_S2_NEG_NegEvnt_x_DN.txt" %>% glue()
    within_vars <- c("c.NegEvnt", "c.NegEvnt:c.DN")
    between_vars <- c("c.DN", "prop.NegEvnt")
    random_vars <- c("c.NegEvnt")
    
    posterior_df <- posterior_samples_extractor(S2_NEG_ucm, S2_NEG_NegEvnt_x_DN_prop.NegEvnt, obs_lvl_var="gamma")
    
    full_model_decomp <- posterior_r2mlm_draws(mean_imputed_df, posterior_df, between_vars, within_vars, random_vars, 
                                               has_intercept=TRUE, clustermeancentered=TRUE, obs_lvl_var="gamma")
    
    sink(summary_filepath)
    print(psych::describe(full_model_decomp, skew = FALSE, quant=c(.025, .975)), digits=4)
    sink()
    
    # Just some memory saving...
    remove(list=c("S2_NEG_NegEvnt_x_DN_prop.NegEvnt"))
    gc()
    
    #----------------------------------------------------------------------------------------------------------------------
    # Next Set of inputs - lv2 no interaction
    data_loader(POSTERIOR_DIR, "S2_NEG_NegEvnt_DN_prop.NegEvnt")
    summary_filepath <- "{SUMMARY_DIR}/var_decomp_S2_NEG_NegEvnt_DN_Exp.txt" %>% glue()
    within_vars <- c("c.NegEvnt")
    between_vars <- c("c.DN", "prop.NegEvnt")
    random_vars <- c("c.NegEvnt")
    
    posterior_df <- posterior_samples_extractor(S2_NEG_ucm, S2_NEG_NegEvnt_DN_prop.NegEvnt, obs_lvl_var="gamma")
    
    lv2_Exp_DN_decomp <- posterior_r2mlm_draws(mean_imputed_df, posterior_df, between_vars, within_vars, random_vars, 
                                               has_intercept=TRUE, clustermeancentered=TRUE, obs_lvl_var="gamma")
    
    sink(summary_filepath)
    print(psych::describe(lv2_Exp_DN_decomp, skew = FALSE, quant=c(.025, .975)), digits=4)
    sink()
    
    # Just some memory saving...
    remove(list=c("S2_NEG_NegEvnt_DN_prop.NegEvnt"))
    gc()
    
    #----------------------------------------------------------------------------------------------------------------------
    # Next Set of inputs - only DN at level 2 of the equation
    data_loader(POSTERIOR_DIR, "S2_NEG_NegEvnt_DN")
    summary_filepath <- "{SUMMARY_DIR}/var_decomp_S2_NEG_NegEvnt_DN.txt" %>% glue()
    within_vars <- c("c.NegEvnt")
    between_vars <- c("c.DN")
    random_vars <- c("c.NegEvnt")
    
    posterior_df <- posterior_samples_extractor(S2_NEG_ucm, S2_NEG_NegEvnt_DN, obs_lvl_var="gamma")
    
    lv2_DN_decomp <- posterior_r2mlm_draws(mean_imputed_df, posterior_df, between_vars, within_vars, random_vars, 
                                           has_intercept=TRUE, clustermeancentered=TRUE, obs_lvl_var="gamma")
    
    sink(summary_filepath)
    print(psych::describe(lv2_DN_decomp, skew = FALSE, quant=c(.025, .975)), digits=4)
    sink()
    
    # Just some memory saving...
    remove(list=c("S2_NEG_NegEvnt_DN"))
    gc()
    
    #----------------------------------------------------------------------------------------------------------------------
    # Next Set of inputs - lv2 no interaction
    data_loader(POSTERIOR_DIR, "S2_NEG_NegEvnt_prop.NegEvnt")
    summary_filepath <- "{SUMMARY_DIR}/var_decomp_S2_NEG_NegEvnt_Exp.txt" %>% glue()
    within_vars <- c("c.NegEvnt")
    between_vars <- c("prop.NegEvnt")
    random_vars <- c("c.NegEvnt")
    
    posterior_df <- posterior_samples_extractor(S2_NEG_ucm, S2_NEG_NegEvnt_prop.NegEvnt, obs_lvl_var="gamma")
    
    lv2_Exp_decomp <- posterior_r2mlm_draws(mean_imputed_df, posterior_df, between_vars, within_vars, random_vars, 
                                            has_intercept=TRUE, clustermeancentered=TRUE, obs_lvl_var="gamma")
    
    sink(summary_filepath)
    print(psych::describe(lv2_Exp_decomp, skew = FALSE, quant=c(.025, .975)), digits=4)
    sink()
    
    # Just some memory saving...
    remove(list=c("S2_NEG_NegEvnt_prop.NegEvnt"))
    gc()
    #----------------------------------------------------------------------------------------------------------------------
    # Next Set of inputs - lv1 only - no level 2 predictors of any kind
    data_loader(POSTERIOR_DIR, "S2_NEG_NegEvnt")
    summary_filepath <- "{SUMMARY_DIR}/var_decomp_S2_NEG_NegEvnt.txt" %>% glue()
    within_vars <- c("c.NegEvnt")
    between_vars <- NULL
    random_vars <- c("c.NegEvnt")
    
    posterior_df <- posterior_samples_extractor(S2_NEG_ucm, S2_NEG_NegEvnt, obs_lvl_var="gamma")
    
    lv1_only_decomp <- posterior_r2mlm_draws(mean_imputed_df, posterior_df, between_vars, within_vars, random_vars, 
                                             has_intercept=TRUE, clustermeancentered=TRUE, obs_lvl_var="gamma")
    sink(summary_filepath)
    print(psych::describe(lv1_only_decomp, skew = FALSE, quant=c(.025, .975)), digits=4)
    sink()
    
    # Just some memory saving...
    remove(list=c("S2_NEG_NegEvnt"))
    gc()
    
    #----------------------------------------------------------------------------------------------------------------------
    model_names <- c("full_model", "lv1_only", "lv2_DN", "lv2_Exp", "lv2_Exp_DN")
    
    model_variance_list <- list()
    model_variance_list[[1]] <- full_model_decomp
    model_variance_list[[2]] <- lv1_only_decomp
    model_variance_list[[3]] <- lv2_DN_decomp
    model_variance_list[[4]] <- lv2_Exp_decomp
    model_variance_list[[5]] <- lv2_Exp_DN_decomp
    
    save.image(LOAD_PATH)
}

# Within Formulas: 
event <- "mean(full_model[['tot_fix_wthn']])"
reactivity <- "(mean(lv2_Exp_DN[['wthn_sig_varn']]) - mean(full_model[['wthn_sig_varn']])) * mean(full_model[['tot_sig_varn']])"
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

color_palette <- c("base"='#ddeeed', 
                   "Tonic \n DN"="#426ebd", 
                   "DN \n Shared w/ \n Exp."="#00878e", 
                   "Reactivity \n"="#ad580b")

riverplot_df_helper(model_variance_list, model_names, within_contrasts, between_contrasts, custom_contrasts, 
                    color_palette = color_palette, 
                    custom_contrast_name = "Combined \n DN Effect", 
                    main_filename = "{GRAPHICS_DIR}/S2_NEG_NegEvnt_decomp.eps" %>% glue(), 
                    main_title = "Total Variance Decomposition: Negative Mood, DN, and Negative Events", 
                    custom_filename = "{GRAPHICS_DIR}/S2_NEG_NegEvnt_DN_combined.eps" %>% glue(), 
                    custom_title = "DN Variance Decomposition: Negative Mood, DN, and Negative Events", 
                    combined_plot_filename = "{GRAPHICS_DIR}/S2_NEG_NegEvnt_full_decomp.eps" %>% glue())
