source("~/dr-consulting_GH/shackman-umd-pax-ema-pub/Graphics_scripts/bayes_variance_riverplots_utils.R")

POSTERIOR_PATH <- "/media/matthew/My Book/EMA_S1_Bayesian_Posteriors"
DATA_FILEPATH <- "~/dr-consulting_GH/shackman-umd-pax-ema-pub/Data/study1_data.RData"
data_loader(POSTERIOR_PATH, "S1_NEG_NegEvnt_Rct", data_filename = DATA_FILEPATH)

generate_reactivity_data <- function(fixed_effs, intercept_name, m_Evnt, Evnt_cat, Evnt_vals, Evnt_var, Evnt_m_var, Evnt_int_var, DN_cat, 
                                     DN_var, log=FALSE, error_var=NULL){
  # browser()
  intercept_fixef <- fixed_effs$Estimate[rownames(fixed_effs) == intercept_name]
  evnt_fixef <- fixed_effs$Estimate[rownames(fixed_effs) == Evnt_var]
  m_evnt_fixef <- fixed_effs$Estimate[rownames(fixed_effs) == Evnt_m_var]
  DN_fixef <- fixed_effs$Estimate[rownames(fixed_effs) == DN_var]
  interaction_fixef <- fixed_effs$Estimate[rownames(fixed_effs) == Evnt_int_var] 
  
  DN_vals <- c(1, -1)
  value <- c()
  effect <- c()
  DN_cat_out <- c()
  Evnt_cat_out <- c()
  ymin <- c()
  ymax <- c()
  xmin <- c()
  xmax <- c()

  for(i in 1:length(DN_cat)){
    for(j in 1:length(Evnt_cat)){
      tonic <- intercept_fixef + m_Evnt*m_evnt_fixef
      DN_base_effect <- DN_vals[i]*DN_fixef
      event <- evnt_fixef*Evnt_vals[j]
      react <- interaction_fixef*DN_vals[i]*Evnt_vals[j]
      
      tmp_val <- c(tonic, DN_base_effect, event, react)
      tmp_val <- c(tmp_val, sum(tmp_val))
      tmp_effect <- c("Tonic", "DN", "Event", "DN x Event", "Mood Rating")
      
      if(log){ 
        # Finding relative contribution of ech effect above the difference between intercept and final mood rating
        # Semi-arbitrary ordering here, could change the magnitude of contributions slightly
        # Logic for order person --> things about the situation --> things about the person in the situation
        tmp_val_exp <- exp(tmp_val)
        
        exp_base <- exp(sum(tmp_val[1]+error_var))
        DN_plus <- exp_base*exp(tmp_val[2])
        event_plus <- DN_plus*exp(tmp_val[3])
        react_plus <- event_plus*exp(tmp_val[4])
        
        diff_series <- diff(c(exp_base, DN_plus, event_plus, react_plus))
        
        tmp_val <- c(exp_base, diff_series, react_plus)
      }
      
      # Make negative for graphing purposes
      tmp_val[5] <- tmp_val[5]*(-1)
      value <- c(value, tmp_val)
      effect <- c(effect, tmp_effect)
      DN_cat_out <- c(DN_cat_out, rep(DN_cat[i], 5))
      Evnt_cat_out <- c(Evnt_cat_out, rep(Evnt_cat[j], 5))
      ymin <- c(ymin, cumsum(tmp_val[1:4]), 0)
      ymax <- c(ymax, 0, cumsum(tmp_val[1:4]))
      xmin <- c(xmin, 1:5)
      xmax <- c(xmax, 2:6)
    }
  }
  
  df_out <- data.frame(effect = effect, 
                       value = value, 
                       DN_cat = DN_cat_out, 
                       Evnt_cat = Evnt_cat_out, 
                       ymin = ymin, 
                       ymax = ymax, 
                       xmin = xmin, 
                       xmax = xmax)
  
  # Getting a set of standardized plotting variables
  df_out$fill <- ifelse(df_out$value < 0, "Decrease", NA)
  df_out$fill <- ifelse(df_out$value > 0, "Increase", df_out$fill)
  df_out$fill <- ifelse(df_out$effect == "Tonic" | df_out$effect == "Mood Rating", "Composite Mood", df_out$fill)
  
  return(df_out)
}

# ---------------------------------------------------------------------------------------------------------------------
# Working with Negative Mood interaction first (found for DN x NegEvnt)
data_loader(POSTERIOR_PATH, "S1_NEG_PosEvnt_Rct", data_filename = DATA_FILEPATH)

NEG_PosEvnt_fixef <- as.data.frame(fixef(S1_NEG_PosEvnt_Rct))
sd_PosEvnt <- sd(dat.study1_model$c.PosEvnt, na.rm = TRUE)
m_Evnt <- mean(dat.study1_model$m.PosEvnt)

DN_cat <- c("+1 SD DN", "-1 SD DN")
Evnt_cat <- c("Yes", "No")
Evnt_vals <- c(sd_PosEvnt, -sd_PosEvnt)

NEG_PosEvnt_int_df <- generate_reactivity_data(NEG_PosEvnt_fixef, intercept_name="NEG_Intercept", m_Evnt=m_Evnt, 
                                               Evnt_cat=Evnt_cat, Evnt_vals=Evnt_vals, Evnt_var="NEG_mic.PosEvnt", 
                                               Evnt_m_var="NEG_m.PosEvnt", Evnt_int_var="NEG_mic.PosEvnt:c.DN", 
                                               DN_cat=DN_cat, DN_var="NEG_c.DN", log=TRUE, error_var = 0.00849)


w <- .75

g <- ggplot(NEG_PosEvnt_int_df[NEG_PosEvnt_int_df$Evnt_cat == "Yes", ], 
            aes(fill=fill))+
  geom_rect(aes(xmin=xmin - w/2, xmax=xmin + w/2, ymin=ymin, ymax=ymax), color="black")+
  scale_x_discrete(limits = c("Average", "DN", "Event", "Reactivity", "Mood Rating"))+
  scale_fill_manual(values=c("Composite Mood"=RColorBrewer::brewer.pal(9, "Purples")[5], 
                             "Decrease"=RColorBrewer::brewer.pal(9, "Blues")[5], 
                             "Increase"=RColorBrewer::brewer.pal(9, "Reds")[5]))+
  geom_segment(data=NEG_PosEvnt_int_df[NEG_PosEvnt_int_df$Evnt_cat == "Yes" & NEG_PosEvnt_int_df$effect != "Mood Rating", ], 
               aes(x=xmin, xend=xmax, y=ymin, yend=ymin))+
  facet_wrap(.~DN_cat)+
  labs(fill="", x="", y="Negative Mood Rating", title = "Study 1 Decomposing DN x Positive Event Interaction: Negative Mood", 
       caption = "DN = dispositional negativity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60, vjust=.5))

g

png(filename = "~/dr-consulting_GH/shackman-umd-pax-ema-pub/plot_outputs/S1_NEG_PosEvnt_x_DN_waterfall.png", 
    width=9, height=6, units="in", res=600)
g
dev.off()

# ---------------------------------------------------------------------------------------------------------------------
# Working with Positive Mood interaction (found for DN x postive event)
data_loader(POSTERIOR_PATH, "S1_POS_PosEvnt_Rct")

POS_PosEvnt_fixef <- as.data.frame(fixef(S1_POS_PosEvnt_Rct))
sd_PosEvnt <- sd(dat.study1_model$c.PosEvnt, na.rm = TRUE)
m_Evnt <- mean(dat.study1_model$m.PosEvnt)

DN_cat <- c("+1 SD DN", "-1 SD DN")
Evnt_cat <- c("Yes", "No")
Evnt_vals <- c(sd_PosEvnt, -sd_PosEvnt)

POS_PosEvnt_int_df <- generate_reactivity_data(POS_PosEvnt_fixef, intercept_name="POS_Intercept", m_Evnt=m_Evnt, 
                                               Evnt_cat=Evnt_cat, Evnt_vals=Evnt_vals, Evnt_var="POS_mic.PosEvnt", 
                                               Evnt_m_var="POS_m.PosEvnt", Evnt_int_var="POS_mic.PosEvnt:c.DN", 
                                               DN_cat=DN_cat, DN_var="POS_c.DN", log=FALSE)


w <- .75

g <- ggplot(POS_PosEvnt_int_df[POS_PosEvnt_int_df$Evnt_cat == "Yes", ], 
            aes(fill=fill))+
  geom_rect(aes(xmin=xmin - w/2, xmax=xmin + w/2, ymin=ymin, ymax=ymax), color="black")+
  scale_x_discrete(limits = c("Average", "DN", "Event", "Reactivity", "Mood Rating"))+
  scale_fill_manual(values=c("Composite Mood"=RColorBrewer::brewer.pal(9, "Purples")[5], 
                             "Decrease"=RColorBrewer::brewer.pal(9, "Blues")[5], 
                             "Increase"=RColorBrewer::brewer.pal(9, "Reds")[5]))+
  geom_segment(data=POS_PosEvnt_int_df[POS_PosEvnt_int_df$Evnt_cat == "Yes" & POS_PosEvnt_int_df$effect != "Mood Rating", ], 
               aes(x=xmin, xend=xmax, y=ymin, yend=ymin))+
  facet_wrap(.~DN_cat)+
  labs(fill="", x="", y="Negative Mood Rating", title = "Study 1 Decomposing DN x Positive Event Interaction: Positive Mood", 
       caption = "DN = dispositional negativity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60, vjust=.5))

g

png(filename = "~/dr-consulting_GH/shackman-umd-pax-ema-pub/plot_outputs/S1_POS_PosEvnt_x_DN_waterfall.png", 
    width=9, height=6, units="in", res=600)
g
dev.off()

