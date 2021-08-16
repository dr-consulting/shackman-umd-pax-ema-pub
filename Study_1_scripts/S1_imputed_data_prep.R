DATA_DIR <- "~/github/ATNL/shackman-umd-pax-ema-pub/Data"
load("{DATA_DIR}/study1_data.RData" %>% glue::glue())
rdata_files <- list.files(DATA_DIR, full.names = TRUE)
library(mitml)

impute_files <- rdata_files[grepl("study1_data_imp_", rdata_files)] 

dat.study1_list <- list()

NEG_vars <- c("Nrvs", "Anxs", "Unesy")
POS_vars <- c("Chrfl", "Hppy", "Jyfl")

for(i in seq_along(impute_files)) {
    load(impute_files[i])
    
    # Imputed 1 data set 10x in parallel and thus only have one data.frame per imputation data set 
    tmp_df <- mitmlComplete(imp_df)[[1]]
    tmp_lv1_df <- tmp_df %>% 
        mutate(
            # Negative event / Positive - anything above minimum = 1, else 0
            NegEvnt_dich = ifelse(as.numeric(NegEvnt_f) > 2, 1, 0), 
            PosEvnt_dich = ifelse(as.numeric(PosEvnt_f) > 3, 1, 0),
            # Composite scores for NEG and POS mood
            Nrvs = Nrvs_f %>% as.numeric(), 
            Anxs = Anxs_f %>% as.numeric(),
            Unesy = Unesy_f %>% as.numeric(),
            Chrfl = Chrfl_f %>% as.numeric(), 
            Hppy = Hppy_f %>% as.numeric(),
            Jyfl = Jyfl_f %>% as.numeric()
        ) %>% 
        select(all_of(c('ID', NEG_vars, POS_vars, "NegEvnt_dich", "PosEvnt_dich", 'c.DN', 'NegEvnt_f', 'PosEvnt_f'))) 
    
    tmp_lv1_df[["NEG"]] <- rowMeans(tmp_lv1_df[,NEG_vars])
    tmp_lv1_df[["POS"]] <- rowMeans(tmp_lv1_df[,POS_vars])
    
    tmp_lv2_df <- tmp_lv1_df %>% 
        group_by(ID) %>% 
        summarize(prop_NegEvnt = mean(NegEvnt_dich), 
                  prop_PosEvnt = mean(PosEvnt_dich))
    
    tmp_comb_df <- tmp_lv1_df %>% 
        merge(tmp_lv2_df) %>% 
        mutate(
            c.NegEvnt = NegEvnt_dich - prop_NegEvnt, 
            c.PosEvnt = PosEvnt_dich - prop_PosEvnt
        )
    
    dat.study1_list[[i]] <- tmp_comb_df
}

# Test for non-equivalence
for(m in 1:9){
    testthat::expect_error(
        testthat::expect_equal(
            dat.study1_list[[m]], dat.study1_list[[m+1]]
        )
    )
}

save(list = c("dat.study1_list", "dat.study1_lv1", "dat.study1_lv2", "dat.study1_model"), 
     file = "{DATA_DIR}/study1_data.RData" %>% glue::glue())