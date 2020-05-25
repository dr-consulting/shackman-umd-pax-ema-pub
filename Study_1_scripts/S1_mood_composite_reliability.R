# Study 1 - Basic evaluation of construct reliability for momentary mood composites

#######################################################################################################################
library(lavaan)

DATA_FILEPATH <- "~/dr-consulting_GH/shackman-umd-pax-ema-pub/Data/study1_data.RData"
load(DATA_FILEPATH)

CFA_mod<-'
W_PA =~ Jyfl + Chrfl + Hppy
W_NA =~ Nrvs + Anxs + Unesy
'

fit_CFA <- cfa(model = CFA_mod, data=dat.study1_lv1, estimator="MLR")
summary(fit_CFA, standardized=TRUE, fit.measures=TRUE)

mod <- modificationindices(CFA_mod)
mod <- mod[order(mod$mi, decreasing = TRUE),]
mod

# Second model adds in Jyfl ~~ Chrfl covariance
CFA_mod_a<-'
W_PA =~ Jyfl + Chrfl + Hppy
W_NA =~ Nrvs + Anxs + Unesy

# Adding in error covariances - within factor only
Jyfl ~~ Chrfl
'

fit_CFA_a <- cfa(model = CFA_mod_a, data=dat.study1_lv1, estimator="MLR")
summary(fit_CFA_a, standardized=TRUE, fit.measures=TRUE)
anova(fit_CFA, fit_CFA_a)

mod <- modificationindices(fit_CFA_a)
mod <- mod[order(mod$mi, decreasing = TRUE),]
mod

# Third model adds in Nrvs ~~ Anxs covariance
CFA_mod_b<-'
W_PA =~ Jyfl + Chrfl + Hppy
W_NA =~ Nrvs + Anxs + Unesy

# Adding in error covariances - within factor only
Jyfl ~~ Chrfl
Nrvs ~~ Anxs
'

fit_CFA_b <- cfa(model = CFA_mod_b, data=dat.study1_lv1, estimator="MLR")
summary(fit_CFA_b, standardized=TRUE, fit.measures=TRUE)
anova(fit_CFA_a, fit_CFA_b)

mod <- modificationindices(fit_CFA_b)
mod <- mod[order(mod$mi, decreasing = TRUE),]
mod

# Third model adds in Nrvs ~~ Anxs covariance
CFA_mod_c<-'
W_PA =~ Jyfl + Chrfl + Hppy
W_NA =~ Nrvs + Anxs + Unesy

# Adding in error covariances - within factor only
Jyfl ~~ Chrfl
Nrvs ~~ Anxs
Nrvs ~~ Unesy
'

fit_CFA_c <- cfa(model = CFA_mod_c, data=dat.study1_lv1, estimator="MLR")
summary(fit_CFA_c, standardized=TRUE, fit.measures=TRUE)
anova(fit_CFA_b, fit_CFA_c)

mod <- modificationindices(fit_CFA_c)
mod <- mod[order(mod$mi, decreasing = TRUE),]
mod

# Lastly getting alpha and omega from this analysis
semTools::reliability(fit_CFA_c)
