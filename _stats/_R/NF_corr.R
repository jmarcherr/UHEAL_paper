## Statistical analysis for paper Märcher-Rørsted et al. "Peripheral and central effects of auditory aging"


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

## SP amp change
uheal_data$SP_amp<-uheal_data$SP_amp_man  

summary(uheal_data)
str(uheal_data)

### 1. linear models for Y~NFest+Sex+PTAlf ######################

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
for (c in c('Age','PTA_lf','PTA_hf',
            'SP_amp','AP_amp_pm','WV_amp_pm','abr_IV_ratio',
            'memr_slope',
            'FFR_SNR','FFR_noise','EFR_SNR','EFR_noise',
            'Neg_4Hz','ITPC_ratio',
            'AEP_p2n1_int','AEP_p2n1_slope',
            'ABR_NF_slope')){
  temp=print(summary(lm(paste(c,"~ ABR_NF_int"), data=uheal_data)))
  # age values
  est_age = rbind(est_age,summary(lm(paste(c,"~ ABR_NF_int"), data=uheal_data))$coefficients[2,1])
  err_age = rbind(err_age,summary(lm(paste(c,"~ ABR_NF_int"), data=uheal_data))$coefficients[2,2])
  t_age = rbind(t_age,summary(lm(paste(c,"~ ABR_NF_int"), data=uheal_data))$coefficients[2,3])
  p_val_age = rbind(p_val_age,summary(lm(paste(c,"~ ABR_NF_int"), data=uheal_data))$coefficients[2,4])
  
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
rownames(age_tab) <-c('Age','PTA_lf','PTA_hf',
                       'SP_amp','AP_amp_pm','WV_amp_pm','abr_IV_ratio',
                       'memr_slope',
                       'FFR_SNR','FFR_noise','EFR_SNR','EFR_noise',
                       'Neg_4Hz','ITPC_ratio',
                       'AEP_p2n1_int','AEP_p2n1_slope',
                       'ABR_NF_slope')
age_tab

out<-kable(age_tab)
#write.csv2(age_tab, "ABR_NF_corr")
