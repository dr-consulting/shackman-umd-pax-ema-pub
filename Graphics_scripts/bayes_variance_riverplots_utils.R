library(foreach)
library(miceadds)
library(brms)
library(riverplot)
N_CORES <- parallel::detectCores()-2 # for my workstation this gets me 10 cores

#######################################################################################################################
# r2MLM below created by Rights and Sterba
# https://my.vanderbilt.edu/jasonrights/software/r2mlm/

r2MLM <- function(data, within_covs ,between_covs, random_covs, gamma_w, gamma_b, Tau, sigma2, has_intercept=TRUE, 
                  clustermeancentered=TRUE){
  if(has_intercept==T){
    if(length(gamma_b)>1) gamma <- c(1,gamma_w,gamma_b[2:length(gamma_b)])
    if(length(gamma_b)==1) gamma <- c(1,gamma_w)
    if(is.null(within_covs)==T) gamma_w <- 0
  }
  if(has_intercept==F){
    gamma <- c(gamma_w,gamma_b)
    if(is.null(within_covs)==T) gamma_w <- 0
    if(is.null(between_covs)==T) gamma_b <- 0
  }
  if(is.null(gamma)) gamma <- 0
  ##compute phi
  phi <- var(cbind(1,data[,c(within_covs)],data[,c(between_covs)]),na.rm=T)
  if(has_intercept==F) phi <- var(cbind(data[,c(within_covs)],data[,c(between_covs)]),na.rm=T)
  if(is.null(within_covs)==T & is.null(within_covs)==T & has_intercept==F) phi <- 0
  phi_w <- var(data[,within_covs],na.rm=T)
  if(is.null(within_covs)==T) phi_w <- 0
  phi_b <- var(cbind(1,data[,between_covs]),na.rm=T)
  if(is.null(between_covs)==T) phi_b <- 0
  ##compute psi and kappa
  var_randomcovs <- var(cbind(1,data[,c(random_covs)]),na.rm=T)
  if(length(Tau)>1) psi <- matrix(c(diag(Tau)),ncol=1)
  if(length(Tau)==1) psi <- Tau
  if(length(Tau)>1) kappa <- matrix(c(Tau[lower.tri(Tau)==TRUE]),ncol=1)
  if(length(Tau)==1) kappa <- 0
  v <- matrix(c(diag(var_randomcovs)),ncol=1)
  r <- matrix(c(var_randomcovs[lower.tri(var_randomcovs)==TRUE]),ncol=1)
  if(is.null(random_covs)==TRUE){
    v <- 0
    r <- 0
    m <- matrix(1,ncol=1)
  }
  if(length(random_covs)>0) m <- matrix(c(colMeans(cbind(1,data[,c(random_covs)]),na.rm=T)),ncol=1)
  ##total variance
  totalvar_notdecomp <- t(v)%*%psi + 2*(t(r)%*%kappa) + t(gamma)%*%phi%*%gamma + t(m)%*%Tau%*%m + sigma2
  totalwithinvar <- (t(gamma_w)%*%phi_w%*%gamma_w) + (t(v)%*%psi + 2*(t(r)%*%kappa)) + sigma2
  totalbetweenvar <- (t(gamma_b)%*%phi_b%*%gamma_b) + Tau[1]
  totalvar <- totalwithinvar + totalbetweenvar
  ##total decomp
  decomp_fixed_notdecomp <- (t(gamma)%*%phi%*%gamma) / totalvar
  decomp_fixed_within <- (t(gamma_w)%*%phi_w%*%gamma_w) / totalvar
  decomp_fixed_between <- (t(gamma_b)%*%phi_b%*%gamma_b) / totalvar
  decomp_fixed <- decomp_fixed_within + decomp_fixed_between
  decomp_varslopes <- (t(v)%*%psi + 2*(t(r)%*%kappa)) / totalvar
  decomp_varmeans <- (t(m)%*%Tau%*%m) / totalvar
  decomp_sigma <- sigma2/totalvar
  ##within decomp
  decomp_fixed_within_w <- (t(gamma_w)%*%phi_w%*%gamma_w) / totalwithinvar
  decomp_varslopes_w <- (t(v)%*%psi + 2*(t(r)%*%kappa)) / totalwithinvar
  decomp_sigma_w <- sigma2/totalwithinvar
  ##between decomp
  decomp_fixed_between_b <- (t(gamma_b)%*%phi_b%*%gamma_b) / totalbetweenvar
  decomp_varmeans_b <- Tau[1] / totalbetweenvar
  #NEW measures
  if (clustermeancentered==TRUE){
    R2_f <- decomp_fixed
    R2_f1 <- decomp_fixed_within
    R2_f2 <- decomp_fixed_between
    R2_fv <- decomp_fixed + decomp_varslopes
    R2_fvm <- decomp_fixed + decomp_varslopes + decomp_varmeans
    R2_v <- decomp_varslopes
    R2_m <- decomp_varmeans
    R2_f_w <- decomp_fixed_within_w
    R2_f_b <- decomp_fixed_between_b
    R2_fv_w <- decomp_fixed_within_w + decomp_varslopes_w
    R2_v_w <- decomp_varslopes_w
    R2_m_b <- decomp_varmeans_b
  }
  if (clustermeancentered==FALSE){
    R2_f <- decomp_fixed_notdecomp
    R2_fv <- decomp_fixed_notdecomp + decomp_varslopes
    R2_fvm <- decomp_fixed_notdecomp + decomp_varslopes + decomp_varmeans
    R2_v <- decomp_varslopes
    R2_m <- decomp_varmeans
  }
  if(clustermeancentered==TRUE){
    decomp_table <- matrix(c(decomp_fixed_within,decomp_fixed_between,decomp_varslopes,decomp_varmeans,decomp_sigma,
                             decomp_fixed_within_w,"NA",decomp_varslopes_w,"NA",decomp_sigma_w,
                             "NA",decomp_fixed_between_b,"NA",decomp_varmeans_b,"NA"),ncol=3)
    rownames(decomp_table) <- c("fixed, within","fixed, between","slope variation","mean variation","sigma2")
    colnames(decomp_table) <- c("total","within","between")
    R2_table <- matrix(c(R2_f1,R2_f2,R2_v,R2_m,R2_f,R2_fv,R2_fvm,
                         R2_f_w,"NA",R2_v_w,"NA","NA",R2_fv_w,"NA",
                         "NA",R2_f_b,"NA",R2_m_b,"NA","NA","NA")
                       ,ncol=3)
    rownames(R2_table) <- c("f1","f2","v","m","f","fv","fvm")
    colnames(R2_table) <- c("total","within","between")
  }
  ##barchart
  if(clustermeancentered==TRUE){
    contributions_stacked <- matrix(c(decomp_fixed_within,decomp_fixed_between,decomp_varslopes,decomp_varmeans,decomp_sigma,
                                      decomp_fixed_within_w,0,decomp_varslopes_w,0,decomp_sigma_w,
                                      0,decomp_fixed_between_b,0,decomp_varmeans_b,0),5,3)
    colnames(contributions_stacked) <- c("total","within","between")
    rownames(contributions_stacked) <- c("fixed slopes (within)",
                                         "fixed slopes (between)",
                                         "slope variation (within)",
                                         "intercept variation (between)",
                                         "residual (within)")
    #barplot(contributions_stacked, main="Decomposition", horiz=FALSE,
    #        ylim=c(0,1),col=c("darkred","steelblue","darkred","midnightblue","white"),ylab="proportion of variance",
    #        density=c(NA,NA,30,40,NA),angle=c(0,45,0,135,0),xlim=c(0,1),width=c(.3,.3))
    #legend(.30,-.1,legend=rownames(contributions_stacked),fill=c("darkred","steelblue","darkred","midnightblue","white"),
    #       cex=.7, pt.cex = 1,xpd=T,density=c(NA,NA,30,40,NA),angle=c(0,45,0,135,0))
  }
  if(clustermeancentered==FALSE){
    decomp_table <- matrix(c(decomp_fixed_notdecomp,decomp_varslopes,decomp_varmeans,decomp_sigma),ncol=1)
    rownames(decomp_table) <- c("fixed","slope variation","mean variation","sigma2")
    colnames(decomp_table) <- c("total")
    R2_table <- matrix(c(R2_f,R2_v,R2_m,R2_fv,R2_fvm),ncol=1)
    rownames(R2_table) <- c("f","v","m","fv","fvm")
    colnames(R2_table) <- c("total")
    ##barchar
    contributions_stacked <- matrix(c(decomp_fixed_notdecomp,decomp_varslopes,decomp_varmeans,decomp_sigma),4,1)
    colnames(contributions_stacked) <- c("total")
    rownames(contributions_stacked) <- c("fixed slopes",
                                         "slope variation",
                                         "intercept variation",
                                         "residual")
    #barplot(contributions_stacked, main="Decomposition", horiz=FALSE,
    #        ylim=c(0,1),col=c("darkblue","darkblue","darkblue","white"),ylab="proportion of variance",
    #        density=c(NA,30,40,NA),angle=c(0,0,135,0),xlim=c(0,1),width=c(.6))
    #legend(.30,-.1,legend=rownames(contributions_stacked),fill=c("darkblue","darkblue","darkblue","white"),
    #       cex=.7, pt.cex = 1,xpd=TRUE,density=c(NA,30,40,NA),angle=c(0,0,135,0))
  }
  Output <- list(noquote(decomp_table),noquote(R2_table))
  names(Output) <- c("Decompositions","R2s")
  return(Output)
}


#' Utility function that takes in data, a model, level 1 and level 2 variable names, and settings for r2MLM function
#' 
#' @param data is the relevant input data set. Can be either a set list of imputed data sets or a single data set. 
#' 
#' @param within_vars is a character vector of the relevant column names for any within-subjects vars of interest.
#' 
#' @param between_vars is a character vector of the relevant column names for any between-subjects vars of interest.
#' 
#' @param random_vars is a character vector of the within-subject vars that include random variances and covariances.
#' 
#' @param focal_model is a brms model object that contains the primary effects of interest
#' 
#' @param null_model is a brms model object represents a random intercepts model. Required if using lognormal prior for 
#' intercept.
#' 
#' @param has_intercept is a parameter to pass through to the r2mlm function. Defaults to true, and should not be 
#' changed if attempting to replicate this workflow. Currently not robust to changes to this parameter. 
#' 
#' @param clustermeancentered is a paramter to pass through to the r2mlm function. Defaults to true, and should not be 
#' changed if attempting to replicate this workflow. Currently not robust to changes to this parameter. 
#' 
#' @param link_func function for transforming level-1 intercept variance (if link function used). Currently only have 
#' support for lognormal priors

r2MLM_brms_wrapper <- function(df, within_vars, between_vars, random_vars, focal_model, null_model,
                          has_intercept=TRUE, clustermeancentered=TRUE, link_func=NULL){
  
  # First extract the appropriate posterior samples:
  posterior_df <- posterior_samples_extractor(null_model, focal_model, link_func)
  
  if(class(df) == "list"){
    r2mlm_posterior_samples <- data.frame()
    for(l in 1:length(df)){
      start_time <- Sys.time()
      print(paste("Initiating draws of posterior variance decompositions for imputed dataset:", l))
      print(paste("Start date and time:", start_time))
      
      tmp_output <- posterior_r2mlm_draws(df[[l]], posterior_df, between_vars, within_vars, random_vars, has_intercept, 
                                          clustermeancentered)
      tmp_output["imputed_df"] <- l
      r2mlm_posterior_samples <- rbind(r2mlm_posterior_samples, tmp_output)
      
      stop_time <- Sys.time()
      run_time <- round(as.numeric(difftime(stop_time, start_time, units = "min")), 2)
      print(paste("Completed draws for imputed dataset:", l, "out of", length(df))) 
      print(paste("Total runtime:", run_time))
      print("Adjust expectations accordingly")
    }
  }
  
  else{ 
    r2mlm_posterior_samples <- posterior_r2mlm_draws(df, posterior_df, between_vars, within_vars, random_vars, has_intercept, 
                                                     clustermeancentered)
  }
  
  return(r2mlm_posterior_samples)
}


posterior_samples_extractor <- function(null_model, focal_model, link_func=NULL){
  #browser()
  par_vals <- parnames(focal_model)
  # Select out relevant parameters:
  pars_to_select <- grepl("b_.*", par_vals) + 
    grepl("sd_ID__.*", par_vals) +
    grepl("cor_ID__.*", par_vals) +
    grepl("sigma", par_vals)
  
  # Choose only parameters needed for variance calculations
  par_vals <- par_vals[as.logical(pars_to_select)]
  posterior_df <- posterior_samples(focal_model, pars=par_vals, add_chain=TRUE)
  
  # Create variance/covariance columns for random effects: 
  ranef_sd_names <- par_vals[grepl("sd_ID__.*", par_vals)]
  n_ranef <- length(ranef_sd_names)
  if(n_ranef >= 1){
    for(i in 1:n_ranef){
      ranef_var_name <- gsub("sd", "var", ranef_sd_names[i])
      posterior_df[ranef_var_name] <- posterior_df[ranef_sd_names[i]]^2 
    }
  }
  
  # Pull out covariance terms
  ranef_cor_names <- par_vals[grepl("cor_ID__.*", par_vals)]
  n_ranef_cor <- length(ranef_cor_names)
  if(n_ranef_cor >= 1){
    for(i in 1:n_ranef_cor){
      var1_name <- strsplit(ranef_cor_names[i], split = "__")[[i]][2]
      var2_name <- strsplit(ranef_cor_names[i], split = "__")[[i]][3]
      posterior_df[paste("cov_ID", var1_name, var2_name, sep="__")] <- posterior_df[ranef_cor_names[i]]*
        posterior_df[paste0("sd_ID__", var1_name)]*
        posterior_df[paste0("sd_ID__", var2_name)]
    }
  }
  
  # If there is a link function convert as appropriate
  # Current implementation only allows for a lognormal link function
  if(!is.null(link_func)){
    if(link_func == "log"){
      # Grab unconditional or "null" model intercept
      beta_00 <- fixef(null_model)[1,1]
      posterior_df["sigma"] <- lognormal_link_func(beta_00, posterior_df["sigma"])
    }
  }
  
  final_names <- grepl("b_.*", colnames(posterior_df)) + 
    grepl("var_ID__.*", colnames(posterior_df)) +
    grepl("cov_ID__.*", colnames(posterior_df)) +
    grepl("sigma", colnames(posterior_df))
  
  return(posterior_df[as.logical(final_names)])
}


posterior_r2mlm_draws <- function(df, posterior_df, between_vars, within_vars, random_vars, has_intercept, 
                                  clustermeancentered){
  #browser()
  # Note need to create and label interaction function outside of this and pass names in correct locations
  # Added some if/else logic here but honestly this needs a re-factor to catch all edge cases... 
  if(!is.null(within_vars)){
    within_vars_cols <- match(within_vars, colnames(df))
    post_wth_vars <- paste0("b_", within_vars)
  }
  else{
    within_vars_cols <- NULL
  }
  
  if(!is.null(between_vars)){
    between_vars_cols <- match(between_vars, colnames(df))
  }
  else{
    between_vars_cols <- NULL
  }
  
  if(!is.null(random_vars)){
    random_vars_cols <- match(random_vars, colnames(df))
  }
  else{ 
    random_vars_cols <- NULL
  }
  
  # Need to add in the Intercept if present (defuault will be TRUE for top function)
  if(has_intercept){
    random_vars <- c("Intercept", random_vars)
    between_vars <- c("Intercept", between_vars)
    post_btw_vars <- paste0("b_", between_vars)
  }
  
  # Dynamically getting names for Tau matrix variables
  post_tau_vars <- matrix(nrow=length(random_vars), ncol=length(random_vars))
  post_var_names <- colnames(posterior_df)[grepl("var_ID__.*", colnames(posterior_df))]
  post_cov_names <- colnames(posterior_df)[grepl("cov_ID__.*", colnames(posterior_df))]
  
  # Getting all the names in the right places:
  for(i in 1:length(random_vars)){
    for(j in 1:length(random_vars)){
      if(i == j){
        post_tau_vars[i, j] <- post_var_names[i]
      }
      else{
        tmp_cov_name <- paste0("cov_ID__", random_vars[i], "__", random_vars[j])
        if(tmp_cov_name %in% colnames(posterior_df)){
          post_tau_vars[i, j] <- tmp_cov_name
          post_tau_vars[j, i] <- tmp_cov_name
        }
      }
    } 
  }

  if(class(df)  == "data.frame"){
    cl <- parallel::makeCluster(N_CORES)
    doParallel::registerDoParallel(cl)
    post_var_decomp_out <- foreach(r = 1:nrow(posterior_df), .combine = rbind, .export = "r2MLM") %dopar% {
      gamma_w <- unlist(posterior_df[r, post_wth_vars])
      names(gamma_w) <- NULL
      gamma_b <- unlist(posterior_df[r, post_btw_vars])
      names(gamma_b) <-NULL
      sigma <- posterior_df[r, "sigma"]
      tau <- matrix(nrow = nrow(post_tau_vars), ncol = ncol(post_tau_vars))
      
      # Extracting tau matrix values:
      for(i in 1:nrow(post_tau_vars)){
        for(j in 1:nrow(post_tau_vars)){
          tau[i, j] <- posterior_df[r, post_tau_vars[i, j]]
        }
      }
      
      # Written to be re-factored if I want to expand beyond options available for clustermeancentered analyses
      # Currently all analyses include level 1 predictors that have been centered at the individual mean
      r2mlm_out <- r2MLM(data=df, within_covs = within_vars_cols, between_covs = between_vars_cols, 
                         random_covs = random_vars_cols, gamma_w = gamma_w, gamma_b = gamma_b, Tau = tau,
                         sigma2 = sigma, has_intercept = has_intercept, clustermeancentered = clustermeancentered)
      
      decomp <- r2mlm_out$Decompositions
      
      data.frame(tot_fix_wthn = as.numeric(decomp["fixed, within", "total"]),
                 tot_fix_btwn = as.numeric(decomp["fixed, between", "total"]), 
                 tot_slp_varn = as.numeric(decomp["slope variation", "total"]), 
                 tot_int_varn = as.numeric(decomp["mean variation", "total"]), 
                 tot_sig_varn = as.numeric(decomp["sigma2", "total"]), 
                 wthn_fix_wthn = as.numeric(decomp["fixed, within", "within"]), 
                 wthn_slp_varn = as.numeric(decomp["slope variation", "within"]),
                 wthn_sig_varn = as.numeric(decomp["sigma2", "within"]),
                 btwn_fix_btwn = as.numeric(decomp["fixed, between", "between"]), 
                 btwn_int_varn = as.numeric(decomp["mean variation", "between"]))
    }
    parallel::stopCluster(cl)
  }
  return(post_var_decomp_out)
}


lognormal_link_func <- function(beta_00, sigma_samples){
  return(log(1 + sigma_samples/beta_00))
}


data_loader <- function(posterior_path, focal_model, data_filename=NULL, null_model=NULL){
  if(!is.null(data_filename)){
    load(data_filename, envir = .GlobalEnv)
  }
  
  if(!is.null(null_model)){
    null_model_filename <- paste0(posterior_path, "/", null_model, ".RData")
    load(null_model_filename, envir = .GlobalEnv)
  }
  
  focal_model_filename <- paste0(posterior_path, "/", focal_model, ".RData")
  load(focal_model_filename, envir = .GlobalEnv)
}


range01<-function(x){(x-min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T))}


riverplot_df_helper <- function(model_variance_list, model_names, within_constrasts, between_constrasts, 
                                custom_contrasts=NULL, within_color, between_color, merge_color, 
                                custom_contrast_name=NULL, main_filename, main_title, custom_filename=NULL, 
                                custom_title=NULL, combined_plot_filename=NULL){
  
  names(model_variance_list) <- model_names
  
  # Expand posterior data sets: 
  for(i in 1:length(model_names)){
    eval(parse(text = paste(model_names[i], "model_variance_list[[i]]", sep="<-")))
  }
  
  N1 <- c(names(within_constrasts), names(between_constrasts), "Total \n Within", "Total \n Between")
  N2 <- c(rep("Total \n Within", length(within_constrasts)), rep("Total \n Between", length(between_constrasts)), 
          rep("Total \n Variance", 2))
  
  within_decomp <- rep(NA, length(within_constrasts))
  for(i in 1:length(within_constrasts)){
    within_decomp[i] <- eval(parse(text = within_constrasts[i]))
  }
  
  between_decomp <- rep(NA, length(between_contrasts))
  for(i in 1:length(between_contrasts)){
    between_decomp[i] <- eval(parse(text = between_constrasts[i]))
  }
  
  tot_decomp <- c(sum(within_decomp), sum(between_decomp))
  value <- c(within_decomp, between_decomp, tot_decomp)
  
  if(sum(value < 0) > 0){
    print("Detected a negative variance term")
    neg_var_pos <- which(value < 0)
    print(paste(N1[neg_var_pos], "=", value[neg_var_pos]))
    
    within_pos <- which(grepl("*Within", N2))
    if(neg_var_pos %in% within_pos){
      # Reover negative variance from the "Unmodeled" term
      unmod_pos <- which(grepl("Unmod.*Within", N1[within_pos]))
      value[unmod_pos] <- value[unmod_pos] + sum(value[neg_var_pos])
    }
    
    between_pos <- which(grepl("*Between", N2))
    if(neg_var_pos %in% between_pos){
      # Reover negative variance from the "Unmodeled" term
      unmod_pos <- which(grepl("Unmod.*Between", N1))
      value[unmod_pos] <- value[unmod_pos] + sum(value[neg_var_pos])
    }
    value[neg_var_pos] <- 0
  }
  
  # Updating labels (i.e., N1)
  N1 <- paste(N1, "\n", paste0(round(value*100, digits = 2), "%"))
  N2 <- c(N2[grepl("*Within", N2)], 
          N2[grepl("*Between", N2)], 
          rep("Total \n Variance", 2))
  
  N2 <- ifelse(grepl("Total \n Within", N2), N1[grepl("Total \n Within*", N1)], N2)
  N2 <- ifelse(grepl("Total \n Between", N2), N1[grepl("Total \n Between*", N1)], N2)
  
  tot_ind_effects <- length(within_constrasts) + length(between_constrasts)
  effects_y_pos <- 0:(tot_ind_effects-1)*2 + .5
  
  cut_points <- quantile(effects_y_pos, c(.25, .75, .5), names=FALSE)
  
  nodes <- data.frame(ID = c(N1, "Total \n Variance"), 
                      x = c(rep(1, tot_ind_effects), rep(2, 2), 3), 
                      y = c(effects_y_pos, cut_points), 
                      stringsAsFactors = FALSE)
  
  col_pal <- c(rep(within_color, length(within_constrasts)), 
               rep(between_color, length(between_constrasts)),
               within_color, between_color, merge_color)
  
  styles <- lapply(nodes$y, 
                   function(x){
                     list(col = col_pal[x], lty=0, textcol = "black")
                   })
  
  for(i in 1:length(col_pal)){
    styles[[i]]$col <- col_pal[i]
  }
  
  names(styles) <- nodes$ID
  
  river_DF <- data.frame(N1 = N1, 
                         N2 = N2, 
                         Value = value, 
                         stringsAsFactors = FALSE)
  

  main_river_plot <- makeRiver(nodes = nodes, 
                               edges = river_DF, 
                               node_styles = styles)
  
  postscript(main_filename, height = 10, width = 10)
  riverplot(main_river_plot, nodewidth = 3, plot_area = .95)
  title(ylab=main_title)
  dev.off()
  
  if(sum(is.null(c(custom_contrasts, custom_contrast_name, custom_filename, custom_title))) == 0){
    n_custom_constrasts <- length(custom_contrasts)
    cust_river_DF <- data.frame(stringsAsFactors = FALSE)
    
    for(i in 1:n_custom_constrasts){
      cust_river_DF <- rbind(cust_river_DF, river_DF[grepl(paste0(names(custom_contrasts)[i], "*"), river_DF$N1), ])
    }
    
    cust_effects_y_pos <- 0:(n_custom_constrasts-1)*2 + .5
    end_point <- median(cust_effects_y_pos)
    
    nodes <- data.frame(ID = c(cust_river_DF$N1, custom_contrast_name), 
                        x = c(rep(1, n_custom_constrasts), 2), 
                        y = c(cust_effects_y_pos, end_point), 
                        stringsAsFactors = FALSE)
    
    col_pal <- ifelse(grepl("Total \n Within.*", cust_river_DF$N2), within_color, between_color)
    col_pal <- c(col_pal, merge_color)
    
    styles <- lapply(nodes$y, 
                     function(x){
                       list(col = col_pal[x], lty=0, textcol = "black")
                     })
    
    for(i in 1:length(col_pal)){
      styles[[i]]$col <- col_pal[i]
    }
    
    names(styles) <- nodes$ID
    
    # Overwrite this after the colors have been assigned (can leverage previous DF info this way)
    cust_river_DF$N2 <- custom_contrast_name
    
    # Make the custom plot...
    custom_river_plot <- makeRiver(nodes = nodes, 
                                   edges = cust_river_DF, 
                                   node_styles = styles)
    postscript(custom_filename, height = 10, width = 10)
    riverplot(custom_river_plot, nodewidth = 3, plot_area = .95)
    title(ylab=custom_title)
    dev.off()
    
    # Make the combined plot 
    postscript(combined_plot_filename, height = 20, width = 10)
    layout(matrix(c(1,2,2), nrow=3, ncol=1))
    par(mar = c(0,0,0,0), cex = 1)
    riverplot(custom_river_plot, nodewidth = 3, plot_area = .8)
    riverplot(main_river_plot, nodewidth = 3, plot_area = .8)
    dev.off()
  }
}
