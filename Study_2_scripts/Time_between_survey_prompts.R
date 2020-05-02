install.packages("DataCombine")

# Get time gap... 

dat_time <- dat.study2_lv1[,c("ID", "DAY_Corrected", "SIG_Corrected", "TIME_Corrected")]
dat_All <- dat.study2_model[,c("ID", "DAY_Corrected", "SIG_Corrected")]

dat_All$comp_IDs <- paste(dat_All$ID, dat_All$DAY_Corrected, dat_All$SIG_Corrected, sep="_")
dat_time$comp_IDs <- paste(dat_time$ID, dat_time$DAY_Corrected, dat_time$SIG_Corrected, sep="_")
dat_time <- merge(dat_All, dat_time[,c("comp_IDs", "TIME_Corrected")], by="comp_IDs", all.x = TRUE)

dat_time$ID_and_Day <- paste(dat_time$ID, dat_time$DAY_Corrected, sep="_")

dat_time <- DataCombine::slide(dat_time, Var = "TIME_Corrected", TimeVar = "SIG_Corrected", GroupVar = "ID_and_Day", 
                               NewVar = "TIME_Corrected_lag1", keepInvalid = TRUE)

dat_time$TIME_Corrected <- as.POSIXlt(dat_time$TIME_Corrected, format = "%H:%M")
dat_time$TIME_Corrected_lag1 <- as.POSIXlt(dat_time$TIME_Corrected_lag1, format = "%H:%M")


dat_time$Diff_time <- difftime(dat_time$TIME_Corrected, dat_time$TIME_Corrected_lag1, units = "mins")
dat_time$abs_Diff_time <- abs(difftime(dat_time$TIME_Corrected, dat_time$TIME_Corrected_lag1, units = "mins"))

max(difftime(dat_time$TIME_Corrected, dat_time$TIME_Corrected_lag1, units = "mins"), na.rm = TRUE)

mean(dat_time$Diff_time, na.rm=TRUE)
sd(dat_time$Diff_time, na.rm=TRUE)