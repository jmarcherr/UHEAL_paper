library(lattice)# graphical display (xyplot)
library(lme4)     # subject specific models (GLMM)
library(psych)    # OPTIONAL: scatterplots with effects
library(ggplot2)  # to make other plots
library(dplyr)
library(LMMstar)
library(plotrix)               # Load plotrix
#library(xlsx)
library(broom)
library(knitr)
library(corrplot)

# load data
# load data
setwd('/work3/jonmarc/UHEAL_paper/_stats/')
uheal_data <- read.csv("uheal_table.csv", header=TRUE, stringsAsFactors=FALSE)
uheal_data <- filter(uheal_data, CP_new == 0)
# path for saving

setwd('/work3/jonmarc/UHEAL_paper/_stats/_R/lm_models/residuals')

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
# FFR sig
uheal_data$FFR_SNR[uheal_data$FFR_sig==0] <- NaN
uheal_data$EFR_SNR[uheal_data$EFR_sig==0] <- NaN
# TTS and AP

summary(lm("tts~ AP_amp_pm", data=uheal_data))
summary(lm("tts~ FFR_SNR", data=uheal_data))

## check for corr with PTA lf and hf
uheal_data_young <-filter(uheal_data,Age<=39)
summary(lm("PTA_lf~ Age", data=uheal_data_young)) # no corr with low freq
summary(lm("PTA_hf~ Age", data=uheal_data_young)) # no corr with high freq

## old
uheal_data_old <-filter(uheal_data,Age>39)
summary(lm("PTA_lf~ Age", data=uheal_data_old)) # no corr with low freq
summary(lm("PTA_hf~ Age", data=uheal_data_old)) # no corr with high freq

## residual data
uheal_data_res <-uheal_data
uheal_data_res$AP_amp_pm_res <-resid((lm("AP_amp_pm~Age",data=uheal_data_res,na.action=na.exclude)))
uheal_data_res$FFR_SNR_res <-resid(lm("FFR_SNR~Age",data=uheal_data_res,na.action=na.exclude))
uheal_data_res$WV_amp_pm_res <-resid(lm("WV_amp_pm~Age",data=uheal_data_res,na.action=na.exclude))
uheal_data_res$EFR_SNR_res <-resid(lm("EFR_SNR~Age",data=uheal_data_res,na.action=na.exclude))
uheal_data_res$neg_4Hz_res <-resid(lm("neg_4Hz~Age",data=uheal_data_res,na.action=na.exclude))
uheal_data_res$itpc_ratio_4Hz_res <-resid(lm("itpc_ratio_4Hz~Age",data=uheal_data_res,na.action=na.exclude))
uheal_data_res$AEP_n1_4_res <-resid(lm("AEP_n1_4~Age",data=uheal_data_res,na.action=na.exclude))
uheal_data_res$AEP_p2_1_res <-resid(lm("AEP_p2_1~Age",data=uheal_data_res,na.action=na.exclude))

#correlations with AP all residuals
p_val_age = NULL
est_age = NULL
err_age = NULL
t_age = NULL

df_all = NULL
out=NULL
age_tab=NULL
for (c in c('WV_amp_pm_res','FFR_SNR_res','EFR_SNR_res',
            'neg_4Hz_res','itpc_ratio_4Hz_res',
            'AEP_n1_4_res','AEP_p2_1_res','Age')){
  temp=print(summary(lm(paste(c,"~ AP_amp_pm_res"), data=uheal_data_res)))
  # age values
  est_age = rbind(est_age,summary(lm(paste(c,"~ AP_amp_pm_res"), data=uheal_data_res))$coefficients[2,1])
  err_age = rbind(err_age,summary(lm(paste(c,"~ AP_amp_pm_res"), data=uheal_data_res))$coefficients[2,2])
  t_age = rbind(t_age,summary(lm(paste(c,"~ AP_amp_pm_res"), data=uheal_data_res))$coefficients[2,3])
  p_val_age = rbind(p_val_age,summary(lm(paste(c,"~ AP_amp_pm_res"), data=uheal_data_res))$coefficients[2,4])
  
  # df
  df_all = rbind(df_all,temp$df[2])
  
}

# age table

p_corr_age = print(p.adjust(p_val_age,method='BH'))
age_tab <-data.frame(df = df_all,
                     est =est_age,
                     err = err_age,
                     t = t_age,
                     pval = p_val_age,
                     p_corr = p_corr_age)
rownames(age_tab) <-c('WV_amp_pm','FFR_SNR','EFR_SNR',
                      'neg_4Hz','itpc_ratio_4Hz',
                      'AEP_n1_4','AEP_p2_4','Age')
age_tab
out<-kable(age_tab)
write.csv2(age_tab, "AP_corr_residuals.csv")

# correlations with FFR for all residuals
p_val_age = NULL
est_age = NULL
err_age = NULL
t_age = NULL

df_all = NULL
out=NULL
age_tab=NULL
for (c in c('AP_amp_pm_res','WV_amp_pm_res','EFR_SNR_res',
            'neg_4Hz_res','itpc_ratio_4Hz_res',
            'AEP_n1_4_res','AEP_p2_4_res','Age')){
  temp=print(summary(lm(paste(c,"~ FFR_SNR"), data=uheal_data_res)))
  # age values
  est_age = rbind(est_age,summary(lm(paste(c,"~ FFR_SNR_res"), data=uheal_data_res))$coefficients[2,1])
  err_age = rbind(err_age,summary(lm(paste(c,"~ FFR_SNR_res"), data=uheal_data_res))$coefficients[2,2])
  t_age = rbind(t_age,summary(lm(paste(c,"~ FFR_SNR_res"), data=uheal_data_res))$coefficients[2,3])
  p_val_age = rbind(p_val_age,summary(lm(paste(c,"~ FFR_SNR_res"), data=uheal_data_res))$coefficients[2,4])
  
  # df
  df_all = rbind(df_all,temp$df[2])
  
}

# age table

p_corr_age = print(p.adjust(p_val_age,method='BH'))
age_tab <-data.frame(df = df_all,
                     est =est_age,
                     err = err_age,
                     t = t_age,
                     pval = p_val_age,
                     p_corr = p_corr_age)
rownames(age_tab) <-c('AP_amp_pm','WV_amp_pm','EFR_SNR',
                      'neg_4Hz','itpc_ratio_4Hz',
                      'AEP_n1_4','AEP_p2_4','Age')
age_tab
setwd('~/Documents/_postdoc/UHEAL_paper/_statsR/Tables')
out<-kable(age_tab)
write.csv2(age_tab, "FFR_corr_residuals.csv")

#correlations with ITPC all residuals
p_val_age = NULL
est_age = NULL
err_age = NULL
t_age = NULL

df_all = NULL
out=NULL
age_tab=NULL
for (c in c('AP_amp_pm_res','WV_amp_pm_res','FFR_SNR_res','EFR_SNR_res',
            'neg_4Hz_res',
            'AEP_n1_4_res','AEP_p2_4_res','Age')){
  temp=print(summary(lm(paste(c,"~ itpc_ratio_4Hz_res"), data=uheal_data_res)))
  # age values
  est_age = rbind(est_age,summary(lm(paste(c,"~ itpc_ratio_4Hz_res"), data=uheal_data_res))$coefficients[2,1])
  err_age = rbind(err_age,summary(lm(paste(c,"~ itpc_ratio_4Hz_res"), data=uheal_data_res))$coefficients[2,2])
  t_age = rbind(t_age,summary(lm(paste(c,"~ itpc_ratio_4Hz_res"), data=uheal_data_res))$coefficients[2,3])
  p_val_age = rbind(p_val_age,summary(lm(paste(c,"~ itpc_ratio_4Hz_res"), data=uheal_data_res))$coefficients[2,4])
  
  # df
  df_all = rbind(df_all,temp$df[2])
  
}

# age table

p_corr_age = print(p.adjust(p_val_age,method='BH'))
age_tab <-data.frame(df = df_all,
                     est =est_age,
                     err = err_age,
                     t = t_age,
                     pval = p_val_age,
                     p_corr = p_corr_age)
rownames(age_tab) <-c('AP_amp_pm','WV_amp_pm','FFR_SNR','EFR_SNR',
                      'neg_4Hz',
                      'AEP_n1_4','AEP_p2_4','Age')
age_tab
out<-kable(age_tab)
write.csv2(age_tab, "ITPC_corr_residuals.csv")

#correlations with AP young
p_val_age = NULL
est_age = NULL
err_age = NULL
t_age = NULL

df_all = NULL
out=NULL
age_tab=NULL
for (c in c('WV_amp_pm','FFR_SNR','EFR_SNR',
            'neg_4Hz','itpc_ratio_4Hz',
            'AEP_n1_4','AEP_p2_4','Age')){
  temp=print(summary(lm(paste(c,"~ AP_amp_pm"), data=uheal_data_young)))
  # age values
  est_age = rbind(est_age,summary(lm(paste(c,"~ AP_amp_pm"), data=uheal_data_young))$coefficients[2,1])
  err_age = rbind(err_age,summary(lm(paste(c,"~ AP_amp_pm"), data=uheal_data_young))$coefficients[2,2])
  t_age = rbind(t_age,summary(lm(paste(c,"~ AP_amp_pm"), data=uheal_data_young))$coefficients[2,3])
  p_val_age = rbind(p_val_age,summary(lm(paste(c,"~ AP_amp_pm"), data=uheal_data_young))$coefficients[2,4])
  
  # df
  df_all = rbind(df_all,temp$df[2])
  
}

# age table

p_corr_age = print(p.adjust(p_val_age,method='BH'))
age_tab <-data.frame(df = df_all,
                     est =est_age,
                     err = err_age,
                     t = t_age,
                     pval = p_val_age,
                     p_corr = p_corr_age)
rownames(age_tab) <-c('WV_amp_pm','FFR_SNR','EFR_SNR',
                      'neg_4Hz','itpc_ratio_4Hz',
                      'AEP_n1_4','AEP_p2_4','Age')
age_tab
out<-kable(age_tab)
write.csv2(age_tab, "AP_corr_young.csv")

#correlations with AP older
p_val_age = NULL
est_age = NULL
err_age = NULL
t_age = NULL

df_all = NULL
out=NULL
age_tab=NULL
for (c in c('WV_amp_pm','FFR_SNR','EFR_SNR',
            'neg_4Hz','itpc_ratio_4Hz',
            'AEP_n1_4','AEP_p2_4','Age')){
  temp=print(summary(lm(paste(c,"~ AP_amp_pm"), data=uheal_data_old)))
  # age values
  est_age = rbind(est_age,summary(lm(paste(c,"~ AP_amp_pm"), data=uheal_data_old))$coefficients[2,1])
  err_age = rbind(err_age,summary(lm(paste(c,"~ AP_amp_pm"), data=uheal_data_old))$coefficients[2,2])
  t_age = rbind(t_age,summary(lm(paste(c,"~ AP_amp_pm"), data=uheal_data_old))$coefficients[2,3])
  p_val_age = rbind(p_val_age,summary(lm(paste(c,"~ AP_amp_pm"), data=uheal_data_old))$coefficients[2,4])
  
  # df
  df_all = rbind(df_all,temp$df[2])
  
}

# age table

p_corr_age = print(p.adjust(p_val_age,method='BH'))
age_tab <-data.frame(df = df_all,
                     est =est_age,
                     err = err_age,
                     t = t_age,
                     pval = p_val_age,
                     p_corr = p_corr_age)
rownames(age_tab) <-c('WV_amp_pm','FFR_SNR','EFR_SNR',
                      'neg_4Hz','itpc_ratio_4Hz',
                      'AEP_n1_4','AEP_p2_4','Age')
age_tab
out<-kable(age_tab)
write.csv2(age_tab, "AP_corr_old.csv")



## corr plot
myvars  = c('FFR_SNR','EFR_SNR','AP_amp_pm','WV_amp_pm',
            'itpc_ratio_4Hz','neg_4Hz','AEP_n1_4','AEP_p2_1')
#tmp_data = filter(uheal_data,Age>=50)
corr_data <-uheal_data[myvars]

cor_test_mat <- corr.test(corr_data)$p    # Apply corr.test function
cor_test_mat 
cordata = cor(corr_data,use="pairwise.complete.obs")
corrplot(cordata,
         p.mat = cor_test_mat,
         insig = "p-value")



## corr plot residuals
myvars  = c('FFR_SNR_res','EFR_SNR_res','AP_amp_pm_res','WV_amp_pm_res',
            'itpc_ratio_4Hz_res','neg_4Hz_res','AEP_n1_4_res','AEP_p2_1_res')
#tmp_data = filter(uheal_data,Age>=50)
corr_data <-uheal_data_res[myvars]

cor_test_mat <- corr.test(corr_data)$p    # Apply corr.test function
cor_test_mat 
cordata = cor(corr_data,use="pairwise.complete.obs")
corrplot(cordata,
         p.mat = cor_test_mat,
         insig = "p-value")

corPlot(corr_data, cex = 1.2)

corrplot.mixed(cor(corr_data,use="pairwise.complete.obs"),
               lower = "number", 
               upper = "circle",
               tl.col = "black")
## just young


## noise tts vs young memr
memr_data <- filter(uheal_data,Age<=30)
for (c in c('tts','nesi')){
  temp=print(summary(lm(paste(c,"~ memr_slope"), data=memr_data)))
}

