#######################################################################################################################
# Function below taken from Rights and Sterba
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

data, within_covs ,between_covs, random_covs, gamma_w, gamma_b, Tau, sigma2, has_intercept=TRUE, 
clustermeancentered=TRUE

#' Utility function that takes in data, a model, level 1 and level 2 variable names, and settings for r2MLM function
#' 
#' @param data is the relevant input data set. Can be either a set list of imputed data sets or a single data set. 
#' 
#' @param within_covs is a character vector of the relevant column names for any within-subjects effects of interest.
#' 
#' @param between_covs is a character vector of the relevant column names for any between-subjects effects of interest.
#' 
#' @param random_covs is a character vector of the within-subject effects that include random variances and covariances.
#' 
#' @param focal_model is a brms model object that contains the primary effects of interest
#' 
#' @param null_model is a brms model object represents a random intercepts model. 
#' 
#' @param has_intercept is a parameter to pass through to the r2mlm function - defaults to true for current models
#' 
#' @param clustermeancentered is a paramter to pass through to the r2mlm function - defaults to true for current models
#' 
#' @param m number of imputations (if applicable)
#' 
#' @param link_func_trans function for transforming level-1 intercept variance (if link function used)

r2MLM_wrapper <- function(data, within_covs, between_covs, random_covs, focal_model, null_model,
                          has_intercept=TRUE, clustermeancentered=TRUE, m=NULL, link_func_trans=NULL){
  
  
}

posterior_samples_extractor <- function(null_model, between_covs=NULL, within_covs=NULL, random_covs=NULL, focal_model, 
                                        link_func_trans=NULL, has_intercept=TRUE){
  # Initialize empty data.frame for storing results 
  posterior_df <- data.frame()
  
  # Generates list of between-subjects effects 
  if(!is.null(between_covs)){
    
    # Add intercept to between-subjects fixed effects if has_intercept = TRUE
    if(has_intercept){
      between_covs <- c("Intercept", between_covs)
    }
    
    # Create names and pull posterior samples for between-subjects effects
    between_covs_names <- paste0("b_", between_covs)
    for(i in 1:length(between_covs_names)){
      posterior_df[[paste0("btw_", between_covs[i])]] <- posterior_samples(focal_model, between_covs_names[i])
    }
  }
  
  # Generates list of within-subjects effects 
  if(!is.null(within_covs)){
    
    # Create names and pull posterior samples for within-subjects effects
    within_covs_names <- paste0("b_", within_covs)
    for(i in 1:length(within_covs_names)){
      posterior_df[[paste0("wthn_", within_covs[i])]] <- posterior_samples(focal_model, within_covs_names[i])
    }
  }
  
  # Generates fields for random effects covariance matrix
  if(!is.null(random_covs)){
    
    # Pull in intercept if present
    if(has_intercept){
      random_covs <- c("Intercept", random_covs)
    }
    
    # Get variances of random effects in model
    sd_random_covs_names <- paste0("sd_ID__", random_covs)
    n_random_covs <- length(random_covs)
    
    for(i in 1:n_random_covs){
      posterior_df[[paste0("var_", random_covs[i])]] <- posterior_samples(focal_model, random_covs[i])
      
      # Using a try statement to more easily find the correlations 
      # Thought was to make things easier in terms of ordering of inputs and actual model names
      tmp <- random_covs[random_covs != random_covs[i]]
      for(j in (i+1):length(tmp)){
        try{
          posterior_df[[paste0("cov_", random_covs[i])]] <- 
        }
      }
    }
  }
  
  posterior_df <- data.frame()
  
  # Storing between-subjects fixed effects: 
  
  
}

lognormal_link_func <- function(beta_00, sigma_samples){
  return(log(1 + sigma_samples/beta_00))
}