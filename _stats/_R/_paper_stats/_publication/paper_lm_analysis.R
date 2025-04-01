## Statistical analysis for paper Märcher-Rørsted et al. "Peripheral and central effects of auditory aging"
# Analysis contains:
# 1. linear models for Y~Age+Sex+PTAlf
# 2. Linear models for Y~Age+Sex
# 3. Linear models for PTA restricted subgroup (Y~Age+Sex)
# 4. post-hoc correlations
# 5. Age tone-number interaction for P2adapt
# 6. Correlations between electrophysiological measures, including residual correlations
# Output in .csv tables

library(lme4)     # subject specific models (GLMM)
library(dplyr)
library(broom)
library(knitr)
library(corrplot)
library(psych)

# load data
setwd('/work3/jonmarc/UHEAL_paper/_stats/') #<-insert data path here
uheal_data <- read.csv("uheal_table.csv", header=TRUE, stringsAsFactors=FALSE)
uheal_data <- filter(uheal_data, CP_new == 0)
setwd('/work3/jonmarc/UHEAL_paper/_stats/_R/_paper_stats/lm_models') #<-set save path here

# Make categorical variables into factors:
uheal_data$subid <- factor(uheal_data$subid)
# gender factor
uheal_data$gender_fact <-factor(uheal_data$gender)
levels(uheal_data$gender_fact) <- c("Female", "Male")
# translate to sex
uheal_data$sex <-uheal_data$gender
uheal_data$sex_fact <-factor(uheal_data$sex)
levels(uheal_data$sex_fact) <- c("Female", "Male")
# NH/HI
uheal_data$CP <- factor(uheal_data$CP_new)

####### FFR_sig ##############
#find significant FFR data points
uheal_data_FFR = filter(uheal_data,FFR_sig == 1)
uheal_data$FFR_SNR[uheal_data$FFR_sig==0] <- NaN
uheal_data$EFR_SNR[uheal_data$EFR_sig==0] <- NaN
uheal_data$FFR_noise[uheal_data$FFR_sig==0] <-NaN
uheal_data$EFR_noise[uheal_data$EFR_sig==0] <- NaN

####### TEOAE_sig ############
#find significant TEOAE data points
uheal_data$teoae_SNR_1[uheal_data$teoae_sig_1==0]<-NaN
uheal_data$teoae_SNR_2[uheal_data$teoae_sig_2==0]<-NaN
uheal_data$teoae_SNR_3[uheal_data$teoae_sig_3==0]<-NaN
uheal_data$teoae_SNR_4[uheal_data$teoae_sig_4==0]<-NaN
uheal_data$teoae_SNR_5[uheal_data$teoae_sig_5==0]<-NaN

summary(uheal_data)
str(uheal_data)

### 1. linear models for Y~Age+Sex+PTAlf ######################

p_val_age = NULL
est_age = NULL
err_age = NULL
t_age = NULL

p_val_sex = NULL
est_sex = NULL
err_sex = NULL
t_sex = NULL

p_val_pta = NULL
est_pta = NULL
err_pta = NULL
t_pta = NULL

df_all = NULL
out=NULL
age_tab=NULL

#list of included measures
for (c in c('PTA_hf',
            'teoae_SNR_1','teoae_SNR_2','teoae_SNR_3','teoae_SNR_4','teoae_SNR_5',
            'SP_amp','AP_amp_pm','WV_amp_pm','abr_IV_ratio',
            'AP_lat','WV_lat',
            'memr_slope',
            'FFR_SNR','FFR_noise','EFR_SNR','EFR_noise',
            'Neg_4Hz','ITPC_ratio',
            'AEP_p2n1_int','AEP_p2n1_slope',
            'acalos_AC_slope_1','acalos_AC_slope_2','acalos_AC_slope_3','acalos_AC_slope_4',
            'nesi','tts','ssq12_mean')){
  temp=print(summary(lm(paste(c,"~ Age+sex+PTA_lf"), data=uheal_data)))
  # age values
  est_age = rbind(est_age,summary(lm(paste(c,"~ Age+sex+PTA_lf"), data=uheal_data))$coefficients[2,1])
  err_age = rbind(err_age,summary(lm(paste(c,"~ Age+sex+PTA_lf"), data=uheal_data))$coefficients[2,2])
  t_age = rbind(t_age,summary(lm(paste(c,"~ Age+sex+PTA_lf"), data=uheal_data))$coefficients[2,3])
  p_val_age = rbind(p_val_age,summary(lm(paste(c,"~ Age+sex+PTA_lf"), data=uheal_data))$coefficients[2,4])
  # sex values
  est_sex = rbind(est_sex,summary(lm(paste(c,"~ Age+sex+PTA_lf"), data=uheal_data))$coefficients[3,1])
  err_sex = rbind(err_sex,summary(lm(paste(c,"~ Age+sex+PTA_lf"), data=uheal_data))$coefficients[3,2])
  t_sex = rbind(t_sex,summary(lm(paste(c,"~ Age+sex+PTA_lf"), data=uheal_data))$coefficients[3,3])
  p_val_sex = rbind(p_val_sex,summary(lm(paste(c,"~ Age+sex+PTA_lf"),data=uheal_data))$coefficients[3,4])
  # PTA values
  est_pta = rbind(est_pta,summary(lm(paste(c,"~ Age+sex+PTA_lf"), data=uheal_data))$coefficients[4,1])
  err_pta = rbind(err_pta,summary(lm(paste(c,"~ Age+sex+PTA_lf"), data=uheal_data))$coefficients[4,2])
  t_pta = rbind(t_pta,summary(lm(paste(c,"~ Age+sex+PTA_lf"), data=uheal_data))$coefficients[4,3])
  p_val_pta = rbind(p_val_pta,summary(lm(paste(c,"~ Age+sex+PTA_lf"),data=uheal_data))$coefficients[4,4])
  
  # df
  df_all = rbind(df_all,temp$df[2])
  
}

# age table
# correction (Benjamini Hochberg)
p_corr_age_hoch = print(p.adjust(p_val_age,method='BH'))
age_tab <-data.frame(df = df_all,
                     est =est_age,
                     err = err_age,
                     t = t_age,
                     pval = p_val_age,
                     p_corr_hoch = p_corr_age_hoch)
rownames(age_tab) <-c('PTA_hf',
                      'teoae_SNR_1','teoae_SNR_2','teoae_SNR_3','teoae_SNR_4','teoae_SNR_5',
                      'SP_amp','AP_amp_pm','WV_amp_pm','abr_IV_ratio',
                      'AP_lat','WV_lat',
                      'memr_slope',
                      'FFR_SNR','FFR_noise','EFR_SNR','EFR_noise',
                      'Neg_4Hz','ITPC_ratio',
                      'AEP_p2n1_int','AEP_p2n1_slope',
                      'acalos_AC_slope_1','acalos_AC_slope_2','acalos_AC_slope_3','acalos_AC_slope_4',
                      'nesi','tts','ssq12_mean')
age_tab

out<-kable(age_tab)
write.csv2(age_tab, "1_age_tab_PTAlf_corrected.csv")

# sex table
p_corr_sex_hoch = print(p.adjust(p_val_sex,method='BH'))
sex_tab <-data.frame(df = df_all,
                     est =est_sex,
                     err = err_sex,
                     t = t_sex,
                     pval = p_val_sex,
                     p_corr_hoch = p_corr_sex_hoch)
rownames(sex_tab) <-c('PTA_hf',
                      'teoae_SNR_1','teoae_SNR_2','teoae_SNR_3','teoae_SNR_4','teoae_SNR_5',
                      'SP_amp','AP_amp_pm','WV_amp_pm','abr_IV_ratio',
                      'AP_lat','WV_lat',
                      'memr_slope',
                      'FFR_SNR','FFR_noise','EFR_SNR','EFR_noise',
                      'Neg_4Hz','ITPC_ratio',
                      'AEP_p2n1_int','AEP_p2n1_slope',
                      'acalos_AC_slope_1','acalos_AC_slope_2','acalos_AC_slope_3','acalos_AC_slope_4',
                      'nesi','tts','ssq12_mean')
sex_tab

out<-kable(sex_tab)
write.csv2(sex_tab, "1_sex_tab_PTAlf_corrected.csv")

# PTA table
p_corr_pta_hoch = print(p.adjust(p_val_pta,method='BH'))
pta_tab <-data.frame(df = df_all,
                     est =est_pta,
                     err = err_pta,
                     t = t_pta,
                     pval = p_val_pta,
                     p_corr_hoch = p_corr_pta_hoch)
rownames(pta_tab) <-c('PTA_hf',
                      'teoae_SNR_1','teoae_SNR_2','teoae_SNR_3','teoae_SNR_4','teoae_SNR_5',
                      'SP_amp','AP_amp_pm','WV_amp_pm','abr_IV_ratio',
                      'AP_lat','WV_lat',
                      'memr_slope',
                      'FFR_SNR','FFR_noise','EFR_SNR','EFR_noise',
                      'Neg_4Hz','ITPC_ratio',
                      'AEP_p2n1_int','AEP_p2n1_slope',
                      'acalos_AC_slope_1','acalos_AC_slope_2','acalos_AC_slope_3','acalos_AC_slope_4',
                      'nesi','tts','ssq12_mean')
pta_tab
out<-kable(pta_tab)
write.csv2(pta_tab, "1_PTAlf_tab_PTAlf_corrected.csv")

### 2. linear models for Y~Age+Sex ######################

p_val_age = NULL
est_age = NULL
err_age = NULL
t_age = NULL

p_val_sex = NULL
est_sex = NULL
err_sex = NULL
t_sex = NULL


df_all = NULL
out=NULL
age_tab=NULL
for (c in c('PTA_lf','PTA_hf',
            'teoae_SNR_1','teoae_SNR_2','teoae_SNR_3','teoae_SNR_4','teoae_SNR_5',
            'SP_amp','AP_amp_pm','WV_amp_pm','abr_IV_ratio',
            'AP_lat','WV_lat',
            'memr_slope',
            'FFR_SNR','FFR_noise','EFR_SNR','EFR_noise',
            'Neg_4Hz','ITPC_ratio',
            'AEP_p2n1_int','AEP_p2n1_slope',
            'acalos_AC_slope_1','acalos_AC_slope_2','acalos_AC_slope_3','acalos_AC_slope_4',
            'nesi','tts','ssq12_mean')){
  temp=print(summary(lm(paste(c,"~ Age+sex"), data=uheal_data)))
  # age values
  est_age = rbind(est_age,summary(lm(paste(c,"~ Age+sex"), data=uheal_data))$coefficients[2,1])
  err_age = rbind(err_age,summary(lm(paste(c,"~ Age+sex"), data=uheal_data))$coefficients[2,2])
  t_age = rbind(t_age,summary(lm(paste(c,"~ Age+sex"), data=uheal_data))$coefficients[2,3])
  p_val_age = rbind(p_val_age,summary(lm(paste(c,"~ Age+sex"), data=uheal_data))$coefficients[2,4])
  # sex values
  est_sex = rbind(est_sex,summary(lm(paste(c,"~ Age+sex"), data=uheal_data))$coefficients[3,1])
  err_sex = rbind(err_sex,summary(lm(paste(c,"~ Age+sex"), data=uheal_data))$coefficients[3,2])
  t_sex = rbind(t_sex,summary(lm(paste(c,"~ Age+sex"), data=uheal_data))$coefficients[3,3])
  p_val_sex = rbind(p_val_sex,summary(lm(paste(c,"~ Age+sex"),data=uheal_data))$coefficients[3,4])
  
  # df
  df_all = rbind(df_all,temp$df[2])
  
}

# age table

p_corr_age_hoch = print(p.adjust(p_val_age,method='BH'))
age_tab <-data.frame(df = df_all,
                     est =est_age,
                     err = err_age,
                     t = t_age,
                     pval = p_val_age,
                     p_corr_hoch = p_corr_age_hoch)
rownames(age_tab) <-c('PTA_lf','PTA_hf',
                      'teoae_SNR_1','teoae_SNR_2','teoae_SNR_3','teoae_SNR_4','teoae_SNR_5',
                      'SP_amp','AP_amp_pm','WV_amp_pm','abr_IV_ratio',
                      'AP_lat','WV_lat',
                      'memr_slope',
                      'FFR_SNR','FFR_noise','EFR_SNR','EFR_noise',
                      'Neg_4Hz','ITPC_ratio',
                      'AEP_p2n1_int','AEP_p2n1_slope',
                      'acalos_AC_slope_1','acalos_AC_slope_2','acalos_AC_slope_3','acalos_AC_slope_4',
                      'nesi','tts','ssq12_mean')
age_tab

out<-kable(age_tab)
write.csv2(age_tab, "2_age_tab.csv")

#sex table
p_corr_sex_hoch = print(p.adjust(p_val_sex,method='BH'))
sex_tab <-data.frame(df = df_all,
                     est =est_sex,
                     err = err_sex,
                     t = t_sex,
                     pval = p_val_sex,
                     p_corr_hoch = p_corr_sex_hoch)
rownames(sex_tab) <-c('PTA_lf','PTA_hf',
                      'teoae_SNR_1','teoae_SNR_2','teoae_SNR_3','teoae_SNR_4','teoae_SNR_5',
                      'SP_amp','AP_amp_pm','WV_amp_pm','abr_IV_ratio',
                      'AP_lat','WV_lat',
                      'memr_slope',
                      'FFR_SNR','FFR_noise','EFR_SNR','EFR_noise',
                      'Neg_4Hz','ITPC_ratio',
                      'AEP_p2n1_int','AEP_p2n1_slope',
                      'acalos_AC_slope_1','acalos_AC_slope_2','acalos_AC_slope_3','acalos_AC_slope_4',
                      'nesi','tts','ssq12_mean')
sex_tab

out<-kable(sex_tab)
write.csv2(sex_tab, "2_sex_tab.csv")

### 3. Linear models for PTA restricted subgroup (Y~Age+Sex) ###############
# find subgroup participants
uheal_data <- filter(uheal_data, CP_SG == 1)

summary(uheal_data)
str(uheal_data)

p_val_age = NULL
est_age = NULL
err_age = NULL
t_age = NULL

p_val_sex = NULL
est_sex = NULL
err_sex = NULL
t_sex = NULL


df_all = NULL
out=NULL
age_tab=NULL
for (c in c('PTA_lf','PTA_hf',
            'teoae_SNR_1','teoae_SNR_2','teoae_SNR_3','teoae_SNR_4','teoae_SNR_5',
            'SP_amp','AP_amp_pm','WV_amp_pm','abr_IV_ratio',
            'AP_lat','WV_lat',
            'memr_slope',
            'FFR_SNR','FFR_noise','EFR_SNR','EFR_noise',
            'Neg_4Hz','ITPC_ratio',
            'AEP_p2n1_int','AEP_p2n1_slope',
            'acalos_AC_slope_1','acalos_AC_slope_2','acalos_AC_slope_3','acalos_AC_slope_4',
            'nesi','tts','ssq12_mean')){
  temp=print(summary(lm(paste(c,"~ Age+sex"), data=uheal_data)))
  # age values
  est_age = rbind(est_age,summary(lm(paste(c,"~ Age+sex"), data=uheal_data))$coefficients[2,1])
  err_age = rbind(err_age,summary(lm(paste(c,"~ Age+sex"), data=uheal_data))$coefficients[2,2])
  t_age = rbind(t_age,summary(lm(paste(c,"~ Age+sex"), data=uheal_data))$coefficients[2,3])
  p_val_age = rbind(p_val_age,summary(lm(paste(c,"~ Age+sex"), data=uheal_data))$coefficients[2,4])
  # sex values
  est_sex = rbind(est_sex,summary(lm(paste(c,"~ Age+sex"), data=uheal_data))$coefficients[3,1])
  err_sex = rbind(err_sex,summary(lm(paste(c,"~ Age+sex"), data=uheal_data))$coefficients[3,2])
  t_sex = rbind(t_sex,summary(lm(paste(c,"~ Age+sex"), data=uheal_data))$coefficients[3,3])
  p_val_sex = rbind(p_val_sex,summary(lm(paste(c,"~ Age+sex"),data=uheal_data))$coefficients[3,4])
  
  # df
  df_all = rbind(df_all,temp$df[2])
  
}

# age table

p_corr_age_hoch = print(p.adjust(p_val_age,method='BH'))
age_tab <-data.frame(df = df_all,
                     est =est_age,
                     err = err_age,
                     t = t_age,
                     pval = p_val_age,
                     p_corr_hoch = p_corr_age_hoch)
rownames(age_tab) <-c('PTA_lf','PTA_hf',
                      'teoae_SNR_1','teoae_SNR_2','teoae_SNR_3','teoae_SNR_4','teoae_SNR_5',
                      'SP_amp','AP_amp_pm','WV_amp_pm','abr_IV_ratio',
                      'AP_lat','WV_lat',
                      'memr_slope',
                      'FFR_SNR','FFR_noise','EFR_SNR','EFR_noise',
                      'Neg_4Hz','ITPC_ratio',
                      'AEP_p2n1_int','AEP_p2n1_slope',
                      'acalos_AC_slope_1','acalos_AC_slope_2','acalos_AC_slope_3','acalos_AC_slope_4',
                      'nesi','tts','ssq12_mean')
age_tab

out<-kable(age_tab)
write.csv2(age_tab, "3_age_tab_subgroup.csv")

#sex table
p_corr_sex_hoch = print(p.adjust(p_val_sex,method='BH'))
sex_tab <-data.frame(df = df_all,
                     est =est_sex,
                     err = err_sex,
                     t = t_sex,
                     pval = p_val_sex,
                     p_corr_hoch = p_corr_sex_hoch)
rownames(sex_tab) <-c('PTA_lf','PTA_hf',
                      'teoae_SNR_1','teoae_SNR_2','teoae_SNR_3','teoae_SNR_4','teoae_SNR_5',
                      'SP_amp','AP_amp_pm','WV_amp_pm','abr_IV_ratio',
                      'AP_lat','WV_lat',
                      'memr_slope',
                      'FFR_SNR','FFR_noise','EFR_SNR','EFR_noise',
                      'Neg_4Hz','ITPC_ratio',
                      'AEP_p2n1_int','AEP_p2n1_slope',
                      'acalos_AC_slope_1','acalos_AC_slope_2','acalos_AC_slope_3','acalos_AC_slope_4',
                      'nesi','tts','ssq12_mean')
sex_tab

out<-kable(sex_tab)
write.csv2(sex_tab, "3_sex_tab_subgroups.csv")


#### post-hoc correlations ####################################################
# reload data
setwd('/work3/jonmarc/UHEAL_paper/_stats/') #<-insert data path here
uheal_data <- read.csv("uheal_table.csv", header=TRUE, stringsAsFactors=FALSE)
uheal_data <- filter(uheal_data, CP_new == 0)
setwd('/work3/jonmarc/UHEAL_paper/_stats/_R/_paper_stats/lm_models') #<-set save path here

# Make categorical variables into factors:
uheal_data$subid <- factor(uheal_data$subid)
# gender factor
uheal_data$gender_fact <-factor(uheal_data$gender)
levels(uheal_data$gender_fact) <- c("Female", "Male")
# translate to sex
uheal_data$sex <-uheal_data$gender
uheal_data$sex_fact <-factor(uheal_data$sex)
levels(uheal_data$sex_fact) <- c("Female", "Male")
# NH/HI
uheal_data$CP <- factor(uheal_data$CP_new)

####### FFR_sig ##############
#find significant FFR data points
uheal_data_FFR = filter(uheal_data,FFR_sig == 1)
uheal_data$FFR_SNR[uheal_data$FFR_sig==0] <- NaN
uheal_data$EFR_SNR[uheal_data$EFR_sig==0] <- NaN
uheal_data$FFR_noise[uheal_data$FFR_sig==0] <-NaN
uheal_data$EFR_noise[uheal_data$EFR_sig==0] <- NaN


# post-hoc correlations for peripheral measures
myvars  = c('AP_amp_pm','FFR_SNR','EFR_SNR','memr_slope')

corr_data <-uheal_data[myvars]

cor_test_p <- corr.test(corr_data, method='spearman', adjust='none')$p 
cor_test_rho <- corr.test(corr_data, method='spearman', adjust='none')$r

# AP correlations
AP_tab <-data.frame(rho = cbind(print(cor_test_rho[2:4])),
                     p =cbind(print(cor_test_p[2:4])))
rownames(AP_tab) <-c('FFR','EFR','MEMR')
# FFR_correlations
FFR_tab <-data.frame(rho = cbind(cor_test_rho[2,c(1,3,4)]),
                     p =cbind(cor_test_p[2,c(1,3,4)]))
rownames(FFR_tab) <-c('AP','EFR','MEMR')

write.csv2(AP_tab, "4_posthoc_AP_corr.csv")
write.csv2(FFR_tab, "4_posthoc_FFR_corr.csv")

# post-hoc correlations for ISI AEPs (N1 @ ISI=2.3s and P2 @ ISI=0.8s)
myvars  = c('Age','AEP_n1_4','AEP_p2_1')

corr_data <-uheal_data[myvars]

cor_test_p <- corr.test(corr_data, method='spearman', adjust='none')$p 
cor_test_rho <- corr.test(corr_data, method='spearman', adjust='none')$r

AEP_tab <-data.frame(rho = cor_test_rho[1,c(2,3)],
                     p =cor_test_p[1,c(2,3)])
rownames(AEP_tab) <-c('AEP N1 (ISI=2.3s)','AEP P2 (ISI=0.8s)')
write.csv2(AEP_tab, "4_posthoc_AEP_age_corr.csv")


### 5. P2 adapt. age x tone number interaction

# generate tone 1 vs rest
uheal_data$tone_p2_first <-uheal_data$tone_p2_1
df <- data.frame(uheal_data$tone_p2_2,uheal_data$tone_p2_3,uheal_data$tone_p2_4,uheal_data$tone_p2_5,uheal_data$tone_p2_6)
uheal_data$tone_p2_rest <- rowMeans(df)
uheal_data$tone_p2_adapt <-cbind(uheal_data$tone_p2_first,uheal_data$tone_p2_rest)
x<-matrix(factor("1"),nrow = length(uheal_data$tone_p2_1),ncol=1)
y<-matrix(factor("2"),nrow = length(uheal_data$tone_p2_1),ncol=1)


tmp<-data.frame(cbind(uheal_data$subid,uheal_data$Age,uheal_data$tone_p2_first,uheal_data$tone_p2_rest))
tmp <-cbind(tmp[1:2],stack(tmp[3:4]))
colnames(tmp) <-c('subid','Age','P2_amp','tone_number')
tmp$tone_number <-factor(tmp$tone_number)
levels(tmp$tone_number) <- c("1", "2")
#tmp$subid <- factor(tmp$subid)
tmp <- tmp[order(tmp$subid),]
# lm model for first vs rest
m.sem <-lmer("P2_amp~Age*tone_number+(1|subid)",data = tmp)
# extract coefficients
coefs <- data.frame(coef(summary(m.sem)))
# use normal distribution to approximate p-value
coefs$p <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

write.csv2(coefs, "5_posthoc_P2adapt_1vsrest.csv")

# 6. Correlations between electrophysiological measures
myvars  = c('AP_amp_pm','FFR_SNR','ITPC_ratio','Neg_4Hz','tone_p2_1vsrest','AEP_p2n1_int')

corr_data <-uheal_data[myvars]

cor_test_p <- corr.test(corr_data, method='spearman', adjust='none')$p 
cor_test_rho <- corr.test(corr_data, method='spearman', adjust='none')$r

#adjust for tests (bonferroni)
cor_test_p = cor_test_p*(length(myvars)^2-length(myvars))/2
rownames(cor_test_p) = c('AP','FFR','ITPC_Q','AEP_{mu}','P2_{adapt}','P2N1_{int}')
colnames(cor_test_p) = c('AP','FFR','ITPC_Q','AEP_{mu}','P2_{adapt}','P2N1_{int}')


write.csv2(cor_test_p, "6_eeg_corr_p.csv")
write.csv2(cor_test_rho,"6_eeg_corr_rho.csv")

#### residual correlations
#generate residual data
uheal_data_res <-uheal_data
uheal_data_res$AP_amp_pm_res <-resid((lm("AP_amp_pm~Age",data=uheal_data_res,na.action=na.exclude)))
uheal_data_res$FFR_SNR_res <-resid(lm("FFR_SNR~Age",data=uheal_data_res,na.action=na.exclude))
uheal_data_res$Neg_4Hz_res <-resid(lm("Neg_4Hz~Age",data=uheal_data_res,na.action=na.exclude))
uheal_data_res$ITPC_ratio_res <-resid(lm("ITPC_ratio~Age",data=uheal_data_res,na.action=na.exclude))
uheal_data_res$AEP_p2n1_int_res <-resid(lm("AEP_p2n1_int~Age",data=uheal_data_res,na.action=na.exclude))
uheal_data_res$tone_p2_1vsrest_res <-resid(lm("tone_p2_1vsrest~Age",data=uheal_data_res,na.action=na.exclude))
         
# Correlations between measures
myvars  = c('AP_amp_pm_res','FFR_SNR_res','ITPC_ratio_res','Neg_4Hz_res','tone_p2_1vsrest_res','AEP_p2n1_int_res')

corr_data <-uheal_data_res[myvars]

cor_test_p <- corr.test(corr_data, method='spearman', adjust='none')$p 
cor_test_rho <- corr.test(corr_data, method='spearman', adjust='none')$r

#adjust for tests
cor_test_p = cor_test_p*(length(myvars)^2-length(myvars))/2
rownames(cor_test_p) = c('AP','FFR','ITPC_Q','AEP_{mu}','P2_{adapt}','P2N1_{int}')
colnames(cor_test_p) = c('AP','FFR','ITPC_Q','AEP_{mu}','P2_{adapt}','P2N1_{int}')


write.csv2(cor_test_p, "6_eeg_residual_corr_p.csv")
write.csv2(cor_test_rho,"6_eeg_residual_corr_rho.csv")