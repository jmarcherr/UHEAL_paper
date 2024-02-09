library(lattice)# graphical display (xyplot)
library(lme4)     # subject specific models (GLMM)
library(psych)    # OPTIONAL: scatterplots with effects
library(ggplot2)  # to make other plots
library(dplyr)
#library(LMMstar)
library(plotrix)               # Load plotrix
#library(xlsx)
#library(broom)
library(knitr)

# load data
setwd('/work1/jonmarc/UHEAL_master/UHEAL_paper/_stats/')
uheal_data <- read.csv("uheal_table.csv", header=TRUE, stringsAsFactors=FALSE)
uheal_data <- filter(uheal_data, CP_new == 0)
setwd('/work1/jonmarc/UHEAL_master/UHEAL_paper/_stats/_R/lm_models')

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
uheal_data_FFR = filter(uheal_data,FFR_sig == 1)
uheal_data$FFR_SNR[uheal_data$FFR_sig==0] <- NaN
uheal_data$EFR_SNR[uheal_data$EFR_sig==0] <- NaN
uheal_data$FFR_noise[uheal_data$FFR_sig==0] <-NaN
uheal_data$EFR_noise[uheal_data$EFR_sig==0] <- NaN

####### TEOAE_sig ############
uheal_data$teoae_SNR_1[uheal_data$teoae_sig_1==0]<-NaN
uheal_data$teoae_SNR_2[uheal_data$teoae_sig_2==0]<-NaN
uheal_data$teoae_SNR_3[uheal_data$teoae_sig_3==0]<-NaN
uheal_data$teoae_SNR_4[uheal_data$teoae_sig_4==0]<-NaN
uheal_data$teoae_SNR_5[uheal_data$teoae_sig_5==0]<-NaN

summary(uheal_data)
str(uheal_data)

### EVERYTHING ######################

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
            'SP_amp','AP_amp_pm','WV_amp_pm',
            'AP_lat','WV_lat',
            'memr_slope',
            'FFR_SNR','FFR_noise','EFR_SNR','EFR_noise',
            'Neg_4Hz','ITPC_ratio',
            'AEP_p2n1_int',
            'acalos_AC_slope_1','acalos_AC_slope_2','acalos_AC_slope_3','acalos_AC_slope_4',
            'rds','nesi','tts','ssq12_mean')){
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

p_corr_age = print(p.adjust(p_val_age,method='bonferroni'))
p_corr_age_hoch = print(p.adjust(p_val_age,method='BH'))
age_tab <-data.frame(df = df_all,
                     est =est_age,
                     err = err_age,
                     t = t_age,
                     pval = p_val_age,
                     p_corr = p_corr_age,
                     p_corr_hoch = p_corr_age_hoch)
rownames(age_tab) <-c('PTA_lf','PTA_hf',
                      'teoae_SNR_1','teoae_SNR_2','teoae_SNR_3','teoae_SNR_4','teoae_SNR_5',
                      'SP_amp','AP_amp_pm','WV_amp_pm',
                      'AP_lat','WV_lat',
                      'memr_slope',
                      'FFR_SNR','FFR_noise','EFR_SNR','EFR_noise',
                      'Neg_4Hz','ITPC_ratio',
                      'AEP_p2n1_int',
                      'acalos_AC_slope_1','acalos_AC_slope_2','acalos_AC_slope_3','acalos_AC_slope_4',
                      'rds','nesi','tts','ssq12_mean')
age_tab

out<-kable(age_tab)
write.csv2(age_tab, "age_tab_EVERYTHING.csv")

#sex table
p_corr_sex = print(p.adjust(p_val_sex,method='bonferroni'))
p_corr_sex_hoch = print(p.adjust(p_val_sex,method='BH'))
sex_tab <-data.frame(df = df_all,
                     est =est_sex,
                     err = err_sex,
                     t = t_sex,
                     pval = p_val_sex,
                     p_corr = p_corr_sex,
                     p_corr_hoch = p_corr_sex_hoch)
rownames(sex_tab) <-c('PTA_lf','PTA_hf',
                      'teoae_SNR_1','teoae_SNR_2','teoae_SNR_3','teoae_SNR_4','teoae_SNR_5',
                      'SP_amp','AP_amp_pm','WV_amp_pm',
                      'AP_lat','WV_lat',
                      'memr_slope',
                      'FFR_SNR','FFR_noise','EFR_SNR','EFR_noise',
                      'Neg_4Hz','ITPC_ratio',
                      'AEP_p2n1_int',
                      'acalos_AC_slope_1','acalos_AC_slope_2','acalos_AC_slope_3','acalos_AC_slope_4',
                      'rds','nesi','tts','ssq12_mean')
sex_tab

out<-kable(sex_tab)
write.csv2(sex_tab, "sex_tab_EVERYTHING.csv")

