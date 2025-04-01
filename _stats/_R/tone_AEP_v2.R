library(lme4)     # subject specific models (GLMM)

# load data
setwd('/work3/jonmarc/UHEAL_paper/_stats/')
uheal_data <- read.csv("uheal_table.csv", header=TRUE, stringsAsFactors=FALSE)
uheal_data <- filter(uheal_data, CP_new == 0)
setwd('/work3/jonmarc/UHEAL_paper/_stats/_R/lm_models')

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

# generate tone 1 vs rest
uheal_data$tone_p2_first <-uheal_data$tone_p2_1
df <- data.frame(uheal_data$tone_p2_2,uheal_data$tone_p2_3,uheal_data$tone_p2_4,uheal_data$tone_p2_5,uheal_data$tone_p2_6)
uheal_data$tone_p2_rest <- rowMeans(df)
uheal_data$tone_p2_adapt <-cbind(uheal_data$tone_p2_first,uheal_data$tone_p2_rest)
x<-matrix(factor("1"),nrow = length(uheal_data$tone_n1_1),ncol=1)
y<-matrix(factor("2"),nrow = length(uheal_data$tone_n1_1),ncol=1)
uheal_data$tone_number <- cbind(x,y)

tmp<-data.frame(cbind(uheal_data$subid,uheal_data$Age,uheal_data$tone_p2_first,uheal_data$tone_p2_rest))
tmp <-cbind(tmp[1:2],stack(tmp[3:4]))
# lm model for first vs rest
summary(lm("values~X2*ind+(1|X1)",data = tmp))


#df2 <-data.frame
#df2$p2<- data.frame(c(uheal_data$tone_p2_first,uheal_data$tone_p2_rest))
#uheal_data$tone_p2_adapt <-cbind(uheal_data$tone_p2_first,uheal_data$tone_p2_rest)
#uheal_data$tone_n2_first <-uheal_data$tone_n2_1
#uheal_data$tone_n2_rest <- mean(c(uheal_data$tone_n2_2,uheal_data$tone_n2_3,uheal_data$tone_n2_4,uheal_data$tone_n2_5,uheal_data$tone_n2_6))
#x<-matrix(factor("1"),nrow = length(uheal_data$tone_n1_1),ncol=1)
#y<-matrix(factor("2"),nrow = length(uheal_data$tone_n1_1),ncol=1)
#uheal_data$tone_number <- cbind(x,y)
