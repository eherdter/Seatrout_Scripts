# ABOUT ####\
#Notes 8/4/17
#2. Must send data to working directory at work
#3. Must clean up script at the end
#4. Must calculate weight at age for ADULTs so that I can use it in the Delta Method script 
#test
# 10/10/2016 This script imports ALK with Bay.xlsx
# Main Objectives of this script: 
# 1. Imports otolith data. Summarizes data, mean age, mean length, sd and se and total sample number. 
# 2. Makes bay-specific observed ALK, calculates some summary statistics.  
# 3. Makes bay specific smoothed (modeled) ALK with multinomial modeling methods in Ogle (87-).
# 4. Likelihood ratio testing to do among group statistical comparisons - Ogle (102-103)
# 5. Plots observed and smoothed ALK for each bay. 
# 6. Calculates proportional age distribution
# 7. Determines mean length -at-age (otolith database) and produces plots
# 8. One way anova to determine if the mean lengths (and mean ages) of males and females are significantly different among estuaries
# 9. T test to determine whether there is a significant difference between male and female age for each estuary
# 10. T test to determine whether there is a significant difference in male and female length for each estuary.
#11.Calculates proportional age distribution for ADULTS to be used in Delta_Method script (different from #6)
#12. Calculates weight at age for ADULTS to be used in Delta_Method script
#tst

# LOAD PACKAGES #####
library(FSA)
library(magrittr)
library(nnet)
library(plotrix)
library(haven)
library(ggplot2)
library(scales)
library(fishmethods) #masking select from dplyr 
library(dplyr)

# SET WORKING DIRECTORY #####
setwd("~/Desktop/PhD project/Projects/Seatrout/Data")
setwd("U:/PhD_projectfiles/Raw_Data/Age_Length_Data")

#1. LOAD DATA ####
#load the csv file
# subset by which bay I want
# subset the FIM program 
# make sure I have just the "aged sample"
# turn mm to cm
# select just a few variables to make it more manageable
# then drop the remaining bay levels still sticking around (droplevels)
# turn tl from mm to cm
# create length categories with FSA package

testNE <- read.csv("ALK_Bay_and_weight.csv", header=T) %>% subset(program == "FIM" & troutreg =="NE")
#JX

testSE <- read.csv("ALK_Bay_and_weight.csv", header=T) %>% subset(program == "FIM" & troutreg =="SE")
#N. Indian River, St. Sebastian River,  Tequesta

testNW <- read.csv("ALK_Bay_and_weight.csv", header=T) %>% subset(program == "FIM" & troutreg =="NW")
#AP, BB, CK, FW, AB
#App, Big Bend, Cedar Key, Santa Rosa Sound, AB

testSW <- read.csv("ALK_Bay_and_weight.csv", header=T) %>% subset(program == "FIM" & troutreg =="SW")
#CH, Sarasota Bay, TB, Florida Bay, Estero Bay 

Agelength_TB<- droplevels(subset(as.data.frame(read.csv("ALK_Bay_and_weight.csv", header=T)), bay=="TB" & tl>14 & final_age >0 & program== 'FIM', select=c(sex,SpecimenNumber, bay, tl,sl, final_age, wt_total, date, program))) %>% mutate(tl=tl/10, sl=sl/10, lcat2 =lencat(tl, w=2)) #, as.fact=TRUE))- can include this when determing ALK below but the smoothed ALK needs to be the nonfactored version of the length categorization variable. 
Agelength_TB$sex[which(Agelength_TB$sex == "m")] = "M"
Agelength_TB$sex <- droplevels(Agelength_TB$sex)

#change date format into a factor so that I can do summary statistics by year later on
Agelength_TB$date=as.character(Agelength_TB$date)
Agelength_TB$DateNew = as.POSIXct(strptime(Agelength_TB$date, format="%d-%B-%y", tz=""))  #B is the selection for when month is spelled out
Agelength_TB = mutate(Agelength_TB, year = strftime(DateNew, format="%Y")) %>% dplyr::select(c(-date, -DateNew))
Agelength_TB$year = as.factor(Agelength_TB$year)

Agelength_AP<- droplevels(subset(as.data.frame(read.csv("ALK_Bay_and_weight.csv", header=T)), bay=="AP" & tl>0 & final_age >0 & program== 'FIM', select=c(sex,SpecimenNumber, bay, tl,sl, final_age, wt_total, date, program))) %>% mutate(tl=tl/10, sl=sl/10, lcat2 =lencat(tl, w=2)) #, as.fact=TRUE))

Agelength_AP$date=as.character(Agelength_AP$date)
Agelength_AP$DateNew = as.POSIXct(strptime(Agelength_AP$date, format="%d-%B-%y", tz="")) 
Agelength_AP = mutate(Agelength_AP, year = strftime(DateNew, format="%Y")) %>% dplyr::select(-date, -DateNew)
Agelength_AP$year = as.factor(Agelength_AP$year) 

Agelength_CK<- droplevels(subset(as.data.frame(read.csv("ALK_Bay_and_weight.csv", header=T)), bay=="CK" & tl>0 & final_age >0 & program== 'FIM', select=c(sex,SpecimenNumber, bay, tl,sl, final_age, wt_total, date, program))) %>% mutate(tl=tl/10, sl=sl/10,lcat2 =lencat(tl, w=2)) # as.fact=TRUE))

Agelength_CK$date=as.character(Agelength_CK$date)
Agelength_CK$DateNew = as.POSIXct(strptime(Agelength_CK$date, format="%d-%B-%y", tz="")) 
Agelength_CK = mutate(Agelength_CK, year = strftime(DateNew, format="%Y")) %>% dplyr::select(-date, -DateNew)
Agelength_CK$year = as.factor(Agelength_CK$year) 

Agelength_CH<- droplevels(subset(as.data.frame(read.csv("ALK_Bay_and_weight.csv", header=T)), bay=="CH" & tl>0 & final_age >0 & program== 'FIM', select=c(sex,SpecimenNumber, bay, tl,sl, final_age, wt_total, date, program))) %>% mutate(tl=tl/10,sl=sl/10, lcat2 =lencat(tl, w=2)) # as.fact=TRUE))
Agelength_CH$sex[which(Agelength_CH$sex == "f")] = "F"
Agelength_CH$sex <- droplevels(Agelength_CH$sex)

Agelength_CH$date=as.character(Agelength_CH$date)
Agelength_CH$DateNew = as.POSIXct(strptime(Agelength_CH$date, format="%d-%B-%y", tz="")) 
Agelength_CH = mutate(Agelength_CH, year = strftime(DateNew, format="%Y")) %>% dplyr::select(-date, -DateNew)
Agelength_CH$year = as.factor(Agelength_CH$year) 

Agelength_IR<- droplevels(subset(as.data.frame(read.csv("ALK_Bay_and_weight.csv", header=T)), bay=="IR" & tl>0 & final_age >0 & program== 'FIM', select=c(sex,SpecimenNumber, bay, tl,sl, final_age, wt_total, date, program))) %>% mutate(tl=tl/10,sl=sl/10, lcat2 =lencat(tl, w=2)) # as.fact=TRUE))

Agelength_IR$date=as.character(Agelength_IR$date)
Agelength_IR$DateNew = as.POSIXct(strptime(Agelength_IR$date, format="%d-%B-%y", tz="")) 
Agelength_IR = mutate(Agelength_IR, year = strftime(DateNew, format="%Y")) %>% dplyr::select(-date, -DateNew)
Agelength_IR$year = as.factor(Agelength_IR$year) 

Agelength_JX<- droplevels(subset(as.data.frame(read.csv("ALK_Bay_and_weight.csv", header=T)), bay=="JX" & tl>0 & final_age >0 & program== 'FIM', select=c(sex,SpecimenNumber, bay, tl,sl, final_age, wt_total, date, program))) %>% mutate(tl=tl/10,sl=sl/10, lcat2 =lencat(tl, w=2)) # as.fact=TRUE))

Agelength_JX$date=as.character(Agelength_JX$date)
Agelength_JX$DateNew = as.POSIXct(strptime(Agelength_JX$date, format="%d-%B-%y", tz="")) 
Agelength_JX = mutate(Agelength_JX, year = strftime(DateNew, format="%Y")) %>% dplyr::select(-date, -DateNew)
Agelength_JX$year = as.factor(Agelength_JX$year) 


#Explore other bays/estuaries ####

#Tequesta (Souther Indian River lagoon)
Agelength_TQ<- droplevels(subset(as.data.frame(read.csv("ALK_Bay_and_weight.csv", header=T)), bay=="TQ" & tl>0 & program== 'FIM', select=c(sex,SpecimenNumber, bay, tl,sl, final_age, wt_total, date, program))) %>% mutate(tl=tl/10,sl=sl/10, lcat2 =lencat(tl, w=2)) # as.fact=TRUE))
Agelength_TQ$date=as.character(Agelength_TQ$date)
Agelength_TQ$DateNew = as.POSIXct(strptime(Agelength_TQ$date, format="%d-%B-%y", tz="")) 
Agelength_TQ = mutate(Agelength_TQ, year = strftime(DateNew, format="%Y")) %>% dplyr::select(-date, -DateNew)
Agelength_TQ$year = as.factor(Agelength_TQ$year) 

#how many observations by year 
table(Agelength_TQ$year)
length_TQ<- droplevels(subset(as.data.frame(read.csv("ALK_Bay_and_weight.csv", header=T)), bay=="TQ" & tl>0 & program== 'FIM', select=c(sex,SpecimenNumber, bay, tl,sl, final_age, wt_total, date, program))) %>% mutate(tl=tl/10,sl=sl/10, lcat2 =lencat(tl, w=2)) # as.fact=TRUE))
age_TQ<- droplevels(subset(as.data.frame(read.csv("ALK_Bay_and_weight.csv", header=T)), bay=="TQ" & final_age >0 & program== 'FIM', select=c(sex,SpecimenNumber, bay, tl,sl, final_age, wt_total, date, program))) %>% mutate(tl=tl/10,sl=sl/10, lcat2 =lencat(tl, w=2)) # as.fact=TRUE))

#Sarasota Bay (part of SW)
Agelength_SB<- droplevels(subset(as.data.frame(read.csv("ALK_Bay_and_weight.csv", header=T)), bay=="SB" & tl>0 & final_age >0 & program== 'FIM', select=c(sex,SpecimenNumber, bay, tl,sl, final_age, wt_total, date, program))) %>% mutate(tl=tl/10,sl=sl/10, lcat2 =lencat(tl, w=2)) # as.fact=TRUE))
Agelength_SB$date=as.character(Agelength_SB$date)
Agelength_SB$DateNew = as.POSIXct(strptime(Agelength_SB$date, format="%d-%B-%y", tz="")) 
Agelength_SB = mutate(Agelength_SB, year = strftime(DateNew, format="%Y")) %>% dplyr::select(-date, -DateNew)
Agelength_SB$year = as.factor(Agelength_SB$year) 

table(Agelength_SB$year)





#Florida Bay (part of SW)
Agelength_KY<- droplevels(subset(as.data.frame(read.csv("ALK_Bay_and_weight.csv", header=T)), bay=="KY" & tl>0 & final_age >0 & program== 'FIM', select=c(sex,SpecimenNumber, bay, tl,sl, final_age, wt_total, date, program))) %>% mutate(tl=tl/10,sl=sl/10, lcat2 =lencat(tl, w=2)) # as.fact=TRUE))
#Estero bay (part of SW)
Agelength_EB<- droplevels(subset(as.data.frame(read.csv("ALK_Bay_and_weight.csv", header=T)), bay=="EB" & tl>0 & final_age >0 & program== 'FIM', select=c(sex,SpecimenNumber, bay, tl,sl, final_age, wt_total, date, program))) %>% mutate(tl=tl/10,sl=sl/10, lcat2 =lencat(tl, w=2)) # as.fact=TRUE))

nrow(Agelength_TQ) #321
nrow(Agelength_SB) #46
nrow(Agelength_KY) #66
nrow(Agelength_EB) #8

# BASIC DATA SUMMARIZATION ####
#total sample number of FIM data
All= rbind(Agelength_AP, Agelength_CK, Agelength_TB, Agelength_CH, Agelength_JX, Agelength_IR) 

str(All <- na.omit(All))


#Age proportion
All_3under =subset(All, final_age <=3)
Proprotion3under <- nrow(All_3under)/nrow(All)

#Min, max length
min(All$tl)
max(All$tl)

#summarize the entire set and group by sex and then bay
Combined_sum <- dplyr::summarize(group_by(All, sex,bay), mean_tl=mean(tl), sd_tl=sd(tl), se_tl= sd_tl/(sqrt(length(final_age))),min_tl=min(tl), max_tl=max(tl), mean_age=mean(final_age), median_age=median(final_age),sd_age= sd(final_age), se_age=sd_age/sqrt(length(final_age)), min_age= min(final_age), max_age=max(final_age))

Combined_sum_year <- summarize(group_by(All,bay), min_year=min(year), max_year=max(year))

#summarize each to get sex ratios etc. 
TB_sum <- summarise(Agelength_TB,mean_age=mean(final_age), median_age=median(final_age),sd_age= sd(final_age), se_age=sd_age/sqrt(length(final_age)), mean_tl = mean(tl), min_tl=min(tl), max_tl=max(tl), sd_tl= sd(tl), se_tl= sd_tl/(sqrt(length(final_age))), prop_F=nrow(subset(Agelength_TB, sex=="F"))/length(final_age),prop_M=nrow(subset(Agelength_TB, sex=="M"))/length(final_age), N_F=nrow(subset(Agelength_TB, sex=="F")), N_M=nrow(subset(Agelength_TB, sex=="M")), totN =N_F+N_M)
AP_sum <- summarise(Agelength_AP,mean_age=mean(final_age), median_age=median(final_age),sd_age= sd(final_age), se_age=sd_age/sqrt(length(final_age)), mean_tl = mean(tl), min_tl=min(tl), max_tl=max(tl), sd_tl= sd(tl), se_tl= sd_tl/(sqrt(length(final_age))), prop_F=nrow(subset(Agelength_AP, sex=="F"))/length(final_age),prop_M=nrow(subset(Agelength_AP, sex=="M"))/length(final_age), N_F=nrow(subset(Agelength_AP, sex=="F")), N_M=nrow(subset(Agelength_AP, sex=="M")), totN =N_F+N_M)
CH_sum <- summarise(Agelength_CH,mean_age=mean(final_age), median_age=median(final_age),sd_age= sd(final_age), se_age=sd_age/sqrt(length(final_age)), mean_tl = mean(tl), min_tl=min(tl), max_tl=max(tl), sd_tl= sd(tl), se_tl= sd_tl/(sqrt(length(final_age))), prop_F=nrow(subset(Agelength_CH, sex=="F"))/length(final_age),prop_M=nrow(subset(Agelength_CH, sex=="M"))/length(final_age), N_F=nrow(subset(Agelength_CH, sex=="F")), N_M=nrow(subset(Agelength_CH, sex=="M")), totN =N_F+N_M)
CK_sum <- summarise(Agelength_CK,mean_age=mean(final_age), median_age=median(final_age),sd_age= sd(final_age), se_age=sd_age/sqrt(length(final_age)), mean_tl = mean(tl), min_tl=min(tl), max_tl=max(tl), sd_tl= sd(tl), se_tl= sd_tl/(sqrt(length(final_age))), prop_F=nrow(subset(Agelength_CK, sex=="F"))/length(final_age),prop_M=nrow(subset(Agelength_CK, sex=="M"))/length(final_age), N_F=nrow(subset(Agelength_CK, sex=="F")), N_M=nrow(subset(Agelength_CK, sex=="M")), totN =N_F+N_M)
IR_sum <- summarise(Agelength_IR,mean_age=mean(final_age), median_age=median(final_age),sd_age= sd(final_age), se_age=sd_age/sqrt(length(final_age)), mean_tl = mean(tl), min_tl=min(tl), max_tl=max(tl), sd_tl= sd(tl), se_tl= sd_tl/(sqrt(length(final_age))), prop_F=nrow(subset(Agelength_IR, sex=="F"))/length(final_age),prop_M=nrow(subset(Agelength_IR, sex=="M"))/length(final_age), N_F=nrow(subset(Agelength_IR, sex=="F")), N_M=nrow(subset(Agelength_IR, sex=="M")), totN =N_F+N_M)
JX_sum <- summarise(Agelength_JX,mean_age=mean(final_age), median_age=median(final_age),sd_age= sd(final_age), se_age=sd_age/sqrt(length(final_age)), mean_tl = mean(tl), min_tl=min(tl), max_tl=max(tl), sd_tl= sd(tl), se_tl= sd_tl/(sqrt(length(final_age))), prop_F=nrow(subset(Agelength_JX, sex=="F"))/length(final_age),prop_M=nrow(subset(Agelength_JX, sex=="M"))/length(final_age), N_F=nrow(subset(Agelength_JX, sex=="F")), N_M=nrow(subset(Agelength_JX, sex=="M")), totN =N_F+N_M)
All_sum <- rbind(TB_sum, AP_sum, CH_sum, CK_sum, IR_sum, JX_sum) %>% mutate(bay=c("TB", "AP", "CH", "CK", "IR", "JX"))
rownames(All_sum) <- c("TB", "AP", "CH", "CK", "IR", "JX")


#limit minimum length to that in the recreational data set and THEN do summarization to compare. 

AL_TB <- subset(Agelength_TB, tl>24)
AL_AP <- subset(Agelength_AP, tl>20)
AL_CH <- subset(Agelength_CH, tl>26)
AL_CK<- subset(Agelength_CK, tl>26)
AL_IR<- subset(Agelength_IR, tl>26)
AL_JX<- subset(Agelength_JX, tl>26)

AL_all <- rbind(AL_TB, AL_AP, AL_CH, AL_CK, AL_IR, AL_JX)
AL_sum <- AL_all %>% group_by(bay) %>% summarise(mean_age=mean(final_age), median_age=median(final_age),sd_age= sd(final_age), se_age=sd_age/sqrt(length(final_age)), mean_tl = mean(tl), min_tl=min(tl), max_tl=max(tl), sd_tl= sd(tl), se_tl= sd_tl/(sqrt(length(final_age))))

# TWO WAY ANOVA ####
# test for differences in age by bay and sex
All_MF <- droplevels(subset(All, sex %in% c("M", "F")))

aov <- aov(final_age ~ bay + sex, data= All_MF)
summary(aov)
plot(aov)

# Tukey test
TukeyHSD(aov)

aov2 <- aov(tl ~ bay + sex, data= All_MF)
TukeyHSD(aov2)

# test for differences in tl by bay and sex

#CHI SQ TEST To test differencs in age distribution between bays #####
#non-parameteric like anova is below so it doesnt assume normality of the age distribution
#can do a test of variance to determine whether can use anova
#anova assumes homogenous variances
#test for variances
#H0= ratio of variances is equal to 1
#Ha = ratio of variances is not equal to 1

var.test(Agelength_TB$final_age, Agelength_AP$final_age) #equal variances
var.test(Agelength_TB$final_age, Agelength_CK$final_age)
var.test(Agelength_TB$final_age, Agelength_CH$final_age)
var.test(Agelength_TB$final_age, Agelength_IR$final_age) 
var.test(Agelength_TB$final_age, Agelength_JX$final_age)
var.test(Agelength_AP$final_age, Agelength_CK$final_age)
var.test(Agelength_AP$final_age, Agelength_CH$final_age) 
var.test(Agelength_AP$final_age, Agelength_IR$final_age)
var.test(Agelength_AP$final_age, Agelength_JX$final_age) #equal variances
var.test(Agelength_CK$final_age, Agelength_CH$final_age)
var.test(Agelength_CK$final_age, Agelength_IR$final_age)
var.test(Agelength_CK$final_age, Agelength_JX$final_age)
var.test(Agelength_CH$final_age, Agelength_IR$final_age)
var.test(Agelength_CH$final_age, Agelength_JX$final_age)
var.test(Agelength_IR$final_age, Agelength_JX$final_age)


#mostly all variances are unequal so use a chi squared test

ALL_to7 <- All %>% subset(final_age<8)
(age_freq <- xtabs(~bay+final_age, data=ALL_to7))

chisq.test(age_freq[1:2,]) #TB to AP, 13.774,6, 0.03
chisq.test(age_freq[1:3,]) #TB to CK, 184.09, 12, 0.001
chisq.test(age_freq[1:4,]) #TB to CH, X-squared = 223.28, df = 18, p-value < 2.2e-16
chisq.test(age_freq[1:5,]) #TB to IR, X-squared = 241.58, df = 24, p-value < 2.2e-16
chisq.test(age_freq[1:6,]) #TB to JX, X-squared = 469.62, df = 30, p-value < 2.2e-16
chisq.test(age_freq[2:3,]) #AP to CK, X-squared = 97.294, df = 6, p-value < 2.2e-16
chisq.test(age_freq[2:4,]) #AP to CH, X-squared = 166.37, df = 12, p-value < 2.2e-16
chisq.test(age_freq[2:5,]) #AP to IR, X-squared = 200.14, df = 18, p-value < 2.2e-16
chisq.test(age_freq[2:6,]) #AP to JX, X-squared = 404.74, df = 24, p-value < 2.2e-16
chisq.test(age_freq[3:4,]) #CK to CH, X-squared = 135.91, df = 6, p-value < 2.2e-16
chisq.test(age_freq[3:5,]) #CK to IR, X-squared = 183.18, df = 12, p-value < 2.2e-16
chisq.test(age_freq[3:6,]) #CK to JX. X-squared = 381.1, df = 18, p-value < 2.2e-16
chisq.test(age_freq[4:5,]) #CH to IR, X-squared = 9.9532, df = 6, p-value = 0.1266
chisq.test(age_freq[4:6,]) #CH to JX, X-squared = 262.89, df = 12, p-value < 2.2e-16
chisq.test(age_freq[5:6,]) #IR to JX, X-squared = 233.88, df = 6, p-value < 2.2e-16

#Adjusting p-values for the multople comparisons
ps_age<- c(chisq.test(age_freq[1:2,])$p.value,
           chisq.test(age_freq[1:3,])$p.value,
           chisq.test(age_freq[1:4,])$p.value,
           chisq.test(age_freq[1:5,])$p.value,
           chisq.test(age_freq[1:6,])$p.value,
           chisq.test(age_freq[2:3,])$p.value,
           chisq.test(age_freq[2:4,])$p.value,
           chisq.test(age_freq[2:5,])$p.value,
           chisq.test(age_freq[2:6,])$p.value,
           chisq.test(age_freq[3:4,])$p.value,
           chisq.test(age_freq[3:5,])$p.value,
           chisq.test(age_freq[3:6,])$p.value,
           chisq.test(age_freq[4:5,])$p.value,
           chisq.test(age_freq[4:6,])$p.value,
           chisq.test(age_freq[5:6,])$p.value)

pdf <- as.data.frame(p.adjust(ps_age))

#all significantly different excpet for TB to AP and CH and IR

#CHI SQ TEST To test differencs in length distribution between bays #####
#test for equal variances
var.test(Agelength_TB$tl, Agelength_AP$tl)
var.test(Agelength_TB$tl, Agelength_CK$tl)
var.test(Agelength_TB$tl, Agelength_CH$tl)
var.test(Agelength_TB$tl, Agelength_IR$tl) #equal variances
var.test(Agelength_TB$tl, Agelength_JX$tl)
var.test(Agelength_AP$tl, Agelength_CK$tl)
var.test(Agelength_AP$tl, Agelength_CH$tl) #equal variances
var.test(Agelength_AP$tl, Agelength_IR$tl)
var.test(Agelength_AP$tl, Agelength_JX$tl)
var.test(Agelength_CK$tl, Agelength_CH$tl)
var.test(Agelength_CK$tl, Agelength_IR$tl)
var.test(Agelength_CK$tl, Agelength_JX$tl)
var.test(Agelength_CH$tl, Agelength_IR$tl)
var.test(Agelength_CH$tl, Agelength_JX$tl)
var.test(Agelength_IR$tl, Agelength_JX$tl)

#unequal variances so will use chi.squared again

All_len_freq_30to58 <- All %>% subset(lcat2>=32 & lcat2<= 58)
(len_freq <- xtabs(~bay+lcat2, data=All_len_freq_30to58))

chisq.test(len_freq[1:2,]) #TB to AP, X-squared = 47.429, df = 13, p-value = 8.177e-06
chisq.test(len_freq[1:3,]) #TB to CK, X-squared = 76.341, df = 26, p-value = 7.589e-07
chisq.test(len_freq[1:4,]) #TB to CH, X-squared = 97.332, df = 39, p-value = 6.842e-07
chisq.test(len_freq[1:5,]) #TB to IR, X-squared = 184.95, df = 52, p-value < 2.2e-16
chisq.test(len_freq[1:6,]) #TB to JX, X-squared = 366.73, df = 65, p-value < 2.2e-16
chisq.test(len_freq[2:3,]) #AP to CK, X-squared = 31.427, df = 13, p-value = 0.002919
chisq.test(len_freq[2:4,]) #AP to CH, X-squared = 55.254, df = 26, p-value = 0.0007036
chisq.test(len_freq[2:5,]) #AP to IR, X-squared = 102.69, df = 39, p-value = 1.231e-07
chisq.test(len_freq[2:6,]) #AP to JX, X-squared = 318.94, df = 52, p-value < 2.2e-16
chisq.test(len_freq[3:4,]) #CK to CH,X-squared = 4.8824, df = 13, p-value = 0.9777
chisq.test(len_freq[3:5,]) #CK to IR, X-squared = 75.013, df = 26, p-value = 1.203e-06
chisq.test(len_freq[3:6,]) #CK to JX. X-squared = 280.31, df = 39, p-value < 2.2e-16
chisq.test(len_freq[4:5,]) #CH to IR, X-squared = 49.094, df = 13, p-value = 4.262e-06
chisq.test(len_freq[4:6,]) #CH to JX, X-squared = 251, df = 26, p-value < 2.2e-16
chisq.test(len_freq[5:6,]) #IR to JX, X-squared = 210.68, df = 13, p-value < 2.2e-16


ps_length <- c(chisq.test(len_freq[1:2,])$p.value, #TB to AP, 80.36, 0.001
chisq.test(len_freq[1:3,])$p.value, #TB to CK, 101.02, 0.001
chisq.test(len_freq[1:4,])$p.value, #TB to CH, 153.5, <0.001
chisq.test(len_freq[1:5,])$p.value, #TB to IR,
chisq.test(len_freq[1:6,])$p.value, #TB to JX, 
chisq.test(len_freq[2:3,])$p.value, #AP to CK,
chisq.test(len_freq[2:4,])$p.value, #AP to CH,
chisq.test(len_freq[2:5,])$p.value, #AP to IR, 
chisq.test(len_freq[2:6,])$p.value, #AP to JX, 
chisq.test(len_freq[3:4,])$p.value, #CK to CH, #NSD
chisq.test(len_freq[3:5,])$p.value, #CK to IR,
chisq.test(len_freq[3:6,])$p.value, #CK to JX. 
chisq.test(len_freq[4:5,])$p.value, #CH to IR, 
chisq.test(len_freq[4:6,])$p.value, #CH to JX,
chisq.test(len_freq[5:6,])$p.value) #IR to JX

plength <- as.data.frame(p.adjust(ps_length))








# MAKE ALKS ####
# Make table with observed total numbers at length by age 
(rawfreq_TB <- xtabs(~lcat2+final_age, data=Agelength_TB)) 
#rawfreq_TB_test_df <- as.data.frame(as.matrix(xtabs(~lcat2+final_age, data=TB_test)))
# there appears to be a fish that was assigned an age of 3 but is in the 0-2 length category. Going to remove this because its probably a typo. Specified in above subsetting step as tl>20mm =(2cm).   
colSums(rawfreq_TB) #number of lengths obs per age
rowSums(rawfreq_TB) #number age obs per length
(rawfreq_AP <- xtabs(~lcat2+final_age, data=Agelength_AP)) 
rowSums(rawfreq_AP)
colSums(rawfreq_AP)
(rawfreq_CK <- xtabs(~lcat2+final_age, data=Agelength_CK)) 
rowSums(rawfreq_CK)
colSums(rawfreq_CK)
(rawfreq_CH <- xtabs(~lcat2+final_age, data=Agelength_CH)) 
rowSums(rawfreq_CH)
colSums(rawfreq_CH)
(rawfreq_IR <- xtabs(~lcat2+final_age, data=Agelength_IR)) 
rowSums(rawfreq_IR)
colSums(rawfreq_IR)
(rawfreq_JX <- xtabs(~lcat2+final_age, data=Agelength_JX))
rowSums(rawfreq_JX)
colSums(rawfreq_JX)



# MAKE OBSERVED ALK FROM ABOVE TABLES #####

#The conditional proportions that form the ALK are calculated by dividing ecah cell of the frequency table by the sum of the corresponding row. 
#These row proportions are constructed by submitting the xtabs() object to prop.table() and including margin=1 to indicate that the proportions are computed by row (page 92). 

#The alkPlot command used for plotting the observed ALK is unable to extend the x axis to the bounds of c(0,80) because xlim is not working. 
# Therefore, in order to produce a plot with an x axis that can span from 0-80 (like what is happening with the length frequency and the smoothed ALK)
# I need to add in "observed" proportions for length categories that were not sampled. 
# I could have added them to the original data frame but I was concerned that in the process of proportion calculations the extra entries would
# affect the proportions or result in proportions that I didn't want. The smoothing process can estimate proportions outside of the range but I wanted 
# to keep the observed plot with just the proportions from the observed data. Therefore, I added in the zero proportion data by writing and editing
# a csv file which I then read below. 

as.data.frame.matrix((prop.table(rawfreq_TB, margin=1))) %>% write.csv("U:/PhD_projectfiles/Exported_R_Datafiles/alk_TB.csv")
alk_TB <- read.csv("U:/PhD_projectfiles/Exported_R_Datafiles/alk_TB_edit.csv", row.names=1)
names(alk_TB) <- c(1,2,3,4,5,6,7,8,9)
round(alk_TB,3)

as.data.frame.matrix((prop.table(rawfreq_CH, margin=1))) %>% write.csv("U:/PhD_projectfiles/Exported_R_Datafiles/alk_CH.csv")
alk_CH <- read.csv("U:/PhD_projectfiles/Exported_R_Datafiles/alk_CH_edit.csv", row.names=1)
names(alk_CH) <- c(1,2,3,4,5,6,7,8,9)
round(alk_CH,3)

as.data.frame.matrix((prop.table(rawfreq_CK, margin=1))) %>% write.csv("U:/PhD_projectfiles/Exported_R_Datafiles/alk_CK.csv")
alk_CK <- read.csv("U:/PhD_projectfiles/Exported_R_Datafiles/alk_CK_edit.csv", row.names=1)
names(alk_CK) <- c(1,2,3,4,5,6,7,8)
round(alk_CK,3)

as.data.frame.matrix((prop.table(rawfreq_AP, margin=1))) %>% write.csv("U:/PhD_projectfiles/Exported_R_Datafiles/alk_AP.csv")
alk_AP <- read.csv("U:/PhD_projectfiles/Exported_R_Datafiles/alk_AP_edit.csv", row.names=1)
names(alk_AP) <- c(1,2,3,4,5,6,7,8,9,10)
round(alk_AP,3)

as.data.frame.matrix((prop.table(rawfreq_JX, margin=1))) %>% write.csv("U:/PhD_projectfiles/Exported_R_Datafiles/alk_JX.csv")
alk_JX <- read.csv("U:/PhD_projectfiles/Exported_R_Datafiles/alk_JX_edit.csv", row.names=1)
names(alk_JX) <- c(1,2,3,4,5,6,7,8)
round(alk_JX,3)

as.data.frame.matrix((prop.table(rawfreq_IR, margin=1))) %>% write.csv("U:/PhD_projectfiles/Exported_R_Datafiles/alk_IR.csv")
alk_IR <- read.csv("U:/PhD_projectfiles/Exported_R_Datafiles/alk_IR_edit.csv", row.names=1)
names(alk_IR) <- c(1,2,3,4,5,6,7,8,9)
round(alk_IR,3)



# Apply Proportional Age Distribution and Mean Length at Age #####
# 1. Proportional age distribution
# 1. Mean length-at-age. (otolith database)

#Proportional age distribution
  # might need to make sure I add some dummy data for bays that don't have equal number of ages
ad_TB <- xtabs(~final_age, data=Agelength_TB)
round(prop.table(ad_TB), 3)

ad_CK <- xtabs(~final_age, data=Agelength_CK)
round(prop.table(ad_CK), 3)

ad_CH <- xtabs(~final_age, data=Agelength_CH)
round(prop.table(ad_CH), 3)

ad_AP <- xtabs(~final_age, data=Agelength_AP)
round(prop.table(ad_AP), 3)

ad_IR <- xtabs(~final_age, data=Agelength_IR)
round(prop.table(ad_IR), 3)

ad_JX <- xtabs(~final_age, data=Agelength_JX)
round(prop.table(ad_JX), 3)

#Mean Length-at-age
TB_sumlen <- Agelength_TB %>% group_by(final_age) %>% summarize(n=validn(tl), mn=mean(tl, na.rm=TRUE),
                                                    sd=sd(tl, na.rm=TRUE), se=se(tl, na.rm=TRUE)) %>%
                                                    as.data.frame()

CK_sumlen <- Agelength_CK %>% group_by(final_age) %>% summarize(n=validn(tl), mn=mean(tl, na.rm=TRUE),
                                                    sd=sd(tl, na.rm=TRUE), se=se(tl, na.rm=TRUE)) %>%
                                                    as.data.frame()                                       

CH_sumlen <- Agelength_CH %>% group_by(final_age) %>% summarize(n=validn(tl), mn=mean(tl, na.rm=TRUE),
                                                    sd=sd(tl, na.rm=TRUE), se=se(tl, na.rm=TRUE)) %>%
                                                    as.data.frame()                                        
                                       
 AP_sumlen <- Agelength_AP %>% group_by(final_age) %>% summarize(n=validn(tl), mn=mean(tl, na.rm=TRUE),
                                                    sd=sd(tl, na.rm=TRUE), se=se(tl, na.rm=TRUE)) %>%
                                                    as.data.frame()
                                        
IR_sumlen <- Agelength_IR %>% group_by(final_age) %>% summarize(n=validn(tl), mn=mean(tl, na.rm=TRUE),
                                                   sd=sd(tl, na.rm=TRUE), se=se(tl, na.rm=TRUE)) %>%
                                                   as.data.frame() 
                                       
JX_sumlen <- Agelength_JX %>% group_by(final_age) %>% summarize(n=validn(tl), mn=mean(tl, na.rm=TRUE),
                                                  sd=sd(tl, na.rm=TRUE), se=se(tl, na.rm=TRUE)) %>%
                                                  as.data.frame()   


# MULTIPLOT AGE HISTOGRAMS ####


age_AP <- ggplot(Agelength_AP, aes(x=final_age))+ 
  geom_histogram(aes(y=(..count..)/sum(..count..)), binwidth=1, boundary=-0.5, fill="white", color="black")+    # plotting in percent frequency
  scale_y_continuous(labels=percent_format(), name="Frequency (%)", limits=c(0,.4))+  ## plotting in percent frequency
  scale_x_continuous(name="Age (years)", limits=c(0,11), breaks=seq(1,10,1), labels=c(1,2,3,4,5,6,7,8,9,"10+"))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_rect(colour="black",fill="white"), 
        axis.text.x=element_text(colour="black", size=24, family="Times New Roman"), #changing  colour of x axisaxis.text.y=element_text(colour="black"), #changing colour of y acis
        axis.text.y=element_text(colour="black", size=24, family="Times New Roman"), #changing colour of y acis
        axis.title.x=element_text(size=24, family="Times New Roman"),
        axis.title.y=element_text(size=24, family="Times New Roman"),
        plot.title=element_text(size=14))+
  annotate("text", x=10, y=.35, label="AP", size=10, family="Times New Roman")

age_CK <- ggplot(Agelength_CK, aes(x=final_age))+ 
  geom_histogram(aes(y=(..count..)/sum(..count..)), binwidth=1, boundary=-.5, fill="white", color="black")+    # plotting in percent frequency
  scale_y_continuous(labels=percent_format(), name="Frequency (%)", limits=c(0,.45))+  ## plotting in percent frequency
  scale_x_continuous(name="Age (years)", limits=c(0,11), breaks=seq(1,10,1), labels=c(1,2,3,4,5,6,7,8,9,"10+"))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_rect(colour="black",fill="white"), 
        axis.text.x=element_text(colour="black", size=24, family="Times New Roman"), #changing  colour of x axisaxis.text.y=element_text(colour="black"), #changing colour of y acis
        axis.text.y=element_text(colour="black", size=24, family="Times New Roman"), #changing colour of y acis
        axis.title.x=element_text(size=24, family="Times New Roman"),
        axis.title.y=element_text(size=24, family="Times New Roman"),
        plot.title=element_text(size=14))+
  annotate("text", x=10, y=.35, label="CK", size=10, family="Times New Roman")

age_TB <- ggplot(Agelength_TB, aes(x=final_age))+ 
  geom_histogram(aes(y=(..count..)/sum(..count..)), binwidth=1, boundary=-.5, fill="white", color="black")+    # plotting in percent frequency
  scale_y_continuous(labels=percent_format(), name="Frequency (%)", limits=c(0,.4))+  ## plotting in percent frequency
  scale_x_continuous(name="Age (years)", limits=c(0,11), breaks=seq(1,10,1), labels=c(1,2,3,4,5,6,7,8,9,"10+"))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_rect(colour="black",fill="white"), 
        axis.text.x=element_text(colour="black", size=24, family="Times New Roman"), #changing  colour of x axisaxis.text.y=element_text(colour="black"), #changing colour of y acis
        axis.text.y=element_text(colour="black", size=24, family="Times New Roman"), #changing colour of y acis
        axis.title.x=element_text(size=24, family="Times New Roman"),
        axis.title.y=element_text(size=24, family="Times New Roman"),
        plot.title=element_text(size=14))+
  annotate("text", x=10, y=.35, label="TB", size=10, family="Times New Roman")

age_CH <- ggplot(Agelength_CH, aes(x=final_age))+ 
  geom_histogram(aes(y=(..count..)/sum(..count..)), binwidth=1, boundary=-.5, fill="white", color="black")+    # plotting in percent frequency
  scale_y_continuous(labels=percent_format(), name="Frequency (%)", limits=c(0,.4))+  ## plotting in percent frequency
  scale_x_continuous(name="Age (years)", limits=c(0,11), breaks=seq(1,10,1), labels=c(1,2,3,4,5,6,7,8,9,"10+"))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_rect(colour="black",fill="white"), 
        axis.text.x=element_text(colour="black", size=24, family="Times New Roman"), #changing  colour of x axisaxis.text.y=element_text(colour="black"), #changing colour of y acis
        axis.text.y=element_text(colour="black", size=24, family="Times New Roman"), #changing colour of y acis
        axis.title.x=element_text(size=24, family="Times New Roman"),
        axis.title.y=element_text(size=24, family="Times New Roman"),
        plot.title=element_text(size=14))+
  annotate("text", x=10, y=.35, label="CH", size=10, family="Times New Roman")

age_JX <- ggplot(Agelength_JX, aes(x=final_age))+ 
  geom_histogram(aes(y=(..count..)/sum(..count..)), binwidth=1, boundary=-.5, fill="white", color="black")+    # plotting in percent frequency
  scale_y_continuous(labels=percent_format(), name="Frequency (%)", limits=c(0,.45))+  ## plotting in percent frequency
  scale_x_continuous(name="Age (years)", limits=c(0,11), breaks=seq(1,10,1), labels=c(1,2,3,4,5,6,7,8,9,"10+"))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_rect(colour="black",fill="white"), 
        axis.text.x=element_text(colour="black", size=24, family="Times New Roman"), #changing  colour of x axisaxis.text.y=element_text(colour="black"), #changing colour of y acis
        axis.text.y=element_text(colour="black", size=24, family="Times New Roman"), #changing colour of y acis
        axis.title.x=element_text(size=24, family="Times New Roman"),
        axis.title.y=element_text(size=24, family="Times New Roman"),
        plot.title=element_text(size=14))+
  annotate("text", x=10, y=.35, label="JX", size=10, family="Times New Roman")

age_IR <- ggplot(Agelength_IR, aes(x=final_age))+ 
  geom_histogram(aes(y=(..count..)/sum(..count..)), binwidth=1, boundary=-.5, fill="white", color="black")+    # plotting in percent frequency
  scale_y_continuous(labels=percent_format(), name="Frequency (%)", limits=c(0,.4))+  ## plotting in percent frequency
  scale_x_continuous(name="Age (years)", limits=c(0,11), breaks=seq(1,10,1), labels=c(1,2,3,4,5,6,7,8,9,"10+"))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_rect(colour="black",fill="white"), 
        axis.text.x=element_text(colour="black", size=24, family="Times New Roman"), #changing  colour of x axisaxis.text.y=element_text(colour="black"), #changing colour of y acis
        axis.text.y=element_text(colour="black", size=24, family="Times New Roman"), #changing colour of y acis
        axis.title.x=element_text(size=24, family="Times New Roman"),
        axis.title.y=element_text(size=24, family="Times New Roman"),
        plot.title=element_text(size=14))+
  annotate("text", x=10, y=.35, label="IR", size=10, family="Times New Roman")

#MULTIPLOT FUNCTION #####
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


multiplot()

## MULTIPLOT LENGTH HISTOGRAMS #####

length_AP <- ggplot(Agelength_AP, aes(x=tl))+ 
  geom_histogram(aes(y=(..count..)/sum(..count..)), binwidth=1, boundary=-0.5, fill="white", color="black")+    # plotting in percent frequency
  scale_y_continuous(labels=percent_format(), name="Frequency (%)", limits=c(0,.10), breaks=seq(0,.10, .02))+  ## plotting in percent frequency
  scale_x_continuous(name="Length (cm)", limits=c(0,80), breaks=seq(0,80,20))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_rect(colour="black",fill="white"), 
        axis.text.x=element_text(colour="black", size=24, family="Times New Roman"), #changing  colour of x axisaxis.text.y=element_text(colour="black"), #changing colour of y acis
        axis.text.y=element_text(colour="black", size=24, family="Times New Roman"), #changing colour of y acis
        axis.title.x=element_text(size=24, family="Times New Roman"),
        axis.title.y=element_text(size=24, family="Times New Roman"),
        plot.title=element_text(size=14))+
  annotate("text", x=65, y= .08, label="AP", size=10, family="Times New Roman")

length_CK <- ggplot(Agelength_CK, aes(x=tl))+ 
  geom_histogram(aes(y=(..count..)/sum(..count..)), binwidth=1, boundary=-0.5, fill="white", color="black")+    # plotting in percent frequency
  scale_y_continuous(labels=percent_format(), name="Frequency (%)", limits=c(0,.10), breaks=seq(0,.10, .02))+  ## plotting in percent frequency
  scale_x_continuous(name="Length (cm)", limits=c(0,80), breaks=seq(0,80,20))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_rect(colour="black",fill="white"), 
        axis.text.x=element_text(colour="black", size=24, family="Times New Roman"), #changing  colour of x axisaxis.text.y=element_text(colour="black"), #changing colour of y acis
        axis.text.y=element_text(colour="black", size=24, family="Times New Roman"), #changing colour of y acis
        axis.title.x=element_text(size=24, family="Times New Roman"),
        axis.title.y=element_text(size=24, family="Times New Roman"),
        plot.title=element_text(size=14))+
  annotate("text", x=65, y= .08, label="CK", size=10, family="Times New Roman")

length_TB <- ggplot(Agelength_TB, aes(x=tl))+ 
  geom_histogram(aes(y=(..count..)/sum(..count..)), binwidth=1, boundary=-0.5, fill="white", color="black")+    # plotting in percent frequency
  scale_y_continuous(labels=percent_format(), name="Frequency (%)", limits=c(0,.10), breaks=seq(0,.10, .02))+  ## plotting in percent frequency
  scale_x_continuous(name="Length (cm)", limits=c(0,80), breaks=seq(0,80,20))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_rect(colour="black",fill="white"), 
        axis.text.x=element_text(colour="black", size=24, family="Times New Roman"), #changing  colour of x axisaxis.text.y=element_text(colour="blaTB"), #changing colour of y acis
        axis.text.y=element_text(colour="black", size=24, family="Times New Roman"), #changing colour of y acis
        axis.title.x=element_text(size=24, family="Times New Roman"),
        axis.title.y=element_text(size=24, family="Times New Roman"),
        plot.title=element_text(size=14))+
  annotate("text", x=65, y= .08, label="TB", size=10, family="Times New Roman")

length_CH <- ggplot(Agelength_CH, aes(x=tl))+ 
  geom_histogram(aes(y=(..count..)/sum(..count..)), binwidth=1, boundary=-0.5, fill="white", color="black")+    # plotting in percent frequency
  scale_y_continuous(labels=percent_format(), name="Frequency (%)", limits=c(0,.10), breaks=seq(0,.10, .02))+  ## plotting in percent frequency
  scale_x_continuous(name="Length (cm)", limits=c(0,80), breaks=seq(0,80,20))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_rect(colour="black",fill="white"), 
        axis.text.x=element_text(colour="black", size=24, family="Times New Roman"), #changing  colour of x axisaxis.text.y=element_text(colour="blaCH"), #changing colour of y acis
        axis.text.y=element_text(colour="black", size=24, family="Times New Roman"), #changing colour of y acis
        axis.title.x=element_text(size=24, family="Times New Roman"),
        axis.title.y=element_text(size=24, family="Times New Roman"),
        plot.title=element_text(size=14))+
  annotate("text", x=65, y= .08, label="CH", size=10, family="Times New Roman")

length_JX <- ggplot(Agelength_JX, aes(x=tl))+ 
  geom_histogram(aes(y=(..count..)/sum(..count..)), binwidth=1, boundary=-0.5, fill="white", color="black")+    # plotting in percent frequency
  scale_y_continuous(labels=percent_format(), name="Frequency (%)", limits=c(0,.10), breaks=seq(0,.10, .02))+  ## plotting in percent frequency
  scale_x_continuous(name="Length (cm)", limits=c(0,80), breaks=seq(0,80,20))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_rect(colour="black",fill="white"), 
        axis.text.x=element_text(colour="black", size=24, family="Times New Roman"), #JXanging  colour of x axisaxis.text.y=element_text(colour="blaJX"), #JXanging colour of y acis
        axis.text.y=element_text(colour="black", size=24, family="Times New Roman"), #JXanging colour of y acis
        axis.title.x=element_text(size=24, family="Times New Roman"),
        axis.title.y=element_text(size=24, family="Times New Roman"),
        plot.title=element_text(size=14))+
  annotate("text", x=65, y= .08, label="JX", size=10, family="Times New Roman")

length_IR <- ggplot(Agelength_IR, aes(x=tl))+ 
  geom_histogram(aes(y=(..count..)/sum(..count..)), binwidth=1, boundary=-0.5, fill="white", color="black")+    # plotting in percent frequency
  scale_y_continuous(labels=percent_format(), name="Frequency (%)", limits=c(0,.10), breaks=seq(0,.10, .02))+  ## plotting in percent frequency
  scale_x_continuous(name="Length (cm)", limits=c(0,80), breaks=seq(0,80,20))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_rect(colour="black",fill="white"), 
        axis.text.x=element_text(colour="black", size=24, family="Times New Roman"), #IRanging  colour of x axisaxis.text.y=element_text(colour="blaIR"), #IRanging colour of y acis
        axis.text.y=element_text(colour="black", size=24, family="Times New Roman"), #IRanging colour of y acis
        axis.title.x=element_text(size=24, family="Times New Roman"),
        axis.title.y=element_text(size=24, family="Times New Roman"),
        plot.title=element_text(size=14))+
  annotate("text", x=65, y= .08, label="IR", size=10, family="Times New Roman")











# Plot of individual lengths at age with mean lengths at age superimpose ##### 

plot(tl~final_age, data=Agelength_TB, pch=19, col=rgb(0,0,0,1/10), xlab="Age", ylab= "Total Length (cm)", ylim=c(0,80), xlim=c(0,9))
lines(mn~final_age, data=TB_sumlen, lwd=2, lty=2)
#change scale of x axis

plot(tl~final_age, data=Agelength_CK, pch=19, col=rgb(0,0,0,1/10), xlab="Age", ylab= "Total Length (cm)", ylim=c(0,80), xlim=c(0,9))
lines(mn~final_age, data=CK_sumlen, lwd=2, lty=2)
#change scale of x axis

plot(tl~final_age, data=Agelength_CH, pch=19, col=rgb(0,0,0,1/10), xlab="Age", ylab= "Total Length (cm)", ylim=c(0,80), xlim=c(0,9))
lines(mn~final_age, data=CH_sumlen, lwd=2, lty=2)
#change scale of x axis

plot(tl~final_age, data=Agelength_AP, pch=19, col=rgb(0,0,0,1/10), xlab="Age", ylab= "Total Length (cm)", ylim=c(0,80), xlim=c(0,9))
lines(mn~final_age, data=AP_sumlen, lwd=2, lty=2)
#change scale of x axis

plot(tl~final_age, data=Agelength_IR, pch=19, col=rgb(0,0,0,1/10), xlab="Age", ylab= "Total Length (cm)", ylim=c(0,80), xlim=c(0,9))
lines(mn~final_age, data=IR_sumlen, lwd=2, lty=2)
#change scale of x axis

plot(tl~final_age, data=Agelength_JX, pch=19, col=rgb(0,0,0,1/10), xlab="Age", ylab= "Total Length (cm)", ylim=c(0,80), xlim=c(0,9))
lines(mn~final_age, data=JX_sumlen, lwd=2, lty=2)
#change scale of x axis

#CROSSBAR PLOT for mean and standard error ####

MinMeanSEMMax <- function(x) {
  v <- c(mean(x) - sd(x)/sqrt(length(x)), mean(x), mean(x) + sd(x)/sqrt(length(x)))
  names(v) <- c("ymin", "y", "ymax")
  v
}


All_min <- subset(All, sex %in% c("F", "M"))
labels <- c(F="Female", M="Male")

File <- ("U:/PhD_projectfiles/Figures/age_crossbar.tiff")
if (file.exists(File)) stop(File, " already exists")
dir.create(dirname(File), showWarnings = FALSE)

tiff(File, units="in", width=5, height=5, res=300)


ggplot(All_min, aes(bay, final_age)) +
  stat_summary(fun.data=MinMeanSEMMax, geom="crossbar", colour="black") + 
  scale_y_continuous(breaks=seq(1.5,3,0.5), labels=seq(1.5,3,0.5))+
  facet_grid(sex ~., labeller=labeller(sex=labels)) +
  xlab("Area")+
  ylab("Age (yrs)")+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), 									
        panel.background=element_rect(fill='white', colour='black'),
        axis.title.y = element_text(colour="black", size=20), # changing font of y axis title
        axis.title.x = element_text(colour="black", size=20),
        axis.text.x=element_text(colour="black", size=16), #changing  colour and font of x axis text
        axis.text.y=element_text(colour="black", size=16),
        strip.text.y = element_text(size=16))  #changing colour and font of y axis
dev.off()

#plot.title=element_text(size=14), # changing size of plot title)



File <- ("U:/PhD_projectfiles/Figures/FIM_length_crossbar.tiff")
#if (file.exists(File)) stop(File, " already exists")
dir.create(dirname(File), showWarnings = FALSE)

tiff(File, units="in", width=5, height=5, res=300)


 ggplot(All_min, aes(bay, tl)) +
  stat_summary(fun.data=MinMeanSEMMax, geom="crossbar", colour="black") + 
  scale_y_continuous(breaks=seq(32.5,42,1.5), labels=seq(32.5,42,1.5))+
  facet_grid(sex ~., labeller=labeller(sex=labels)) +
  xlab("Area")+
  ylab("Total Length (cm) ")+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), 									
        panel.background=element_rect(fill='white', colour='black'),
        axis.title.y = element_text(colour="black", size=20), # changing font of y axis title
        axis.title.x = element_text(colour="black", size=20),
        axis.text.x=element_text(colour="black", size=16), #changing  colour and font of x axis text
        axis.text.y=element_text(colour="black", size=16),
        strip.text.y = element_text(size=16))  #changing colour and font of y axis

#plot.title=element_text(size=14), # changing size of plot title)
dev.off()








# PRODUCE SMOOTHED ALKS #####
#from multinomial modeling exercise which can be used to do Likelihood ratio testing. 
#MODELED AGE LENGTH KEYS (aka SMOOTHED ALK)
#-fixes two common issues with Age Length keys detailed on page 92 Ogle. 
#multinomial logistic regression model -Gerritsen et al. 2006. The response variable has more than two levels

#TB
tb <- multinom(final_age~lcat2, data=Agelength_TB, maxit=500)
lens<- seq(0,80, 1)
alksmo.tb <- predict(tb, data.frame(lcat2=lens), type="probs")
row.names(alksmo.tb) <- lens
round(alksmo.tb, 3)

#CK
ck <- multinom(final_age~lcat2, data=Agelength_CK, maxit=500)
lens<- seq(0,80, 1)
alksmo.ck <- predict(ck, data.frame(lcat2=lens), type="probs")
row.names(alksmo.ck) <- lens
round(alksmo.ck, 3)

#CH
ch <- multinom(final_age~lcat2, data=Agelength_CH, maxit=500)
lens<- seq(0,80, 1)
alksmo.ch <- predict(ch, data.frame(lcat2=lens), type="probs")
row.names(alksmo.ch) <- lens
round(alksmo.ch, 3)

#AP
ap <- multinom(final_age~lcat2, data=Agelength_AP, maxit=500)
lens<- seq(0,80, 1)
alksmo.ap <- predict(ap, data.frame(lcat2=lens), type="probs")
row.names(alksmo.ap) <- lens
round(alksmo.ap, 3)

#JX
jx <- multinom(final_age~lcat2, data=Agelength_JX, maxit=500)
lens<- seq(0,80, 1)
alksmo.jx <- predict(jx, data.frame(lcat2=lens), type="probs")
row.names(alksmo.jx) <- lens
round(alksmo.jx, 3)

#IR
ir <- multinom(final_age~lcat2, data=Agelength_IR, maxit=500)
lens<- seq(0,80, 1)
alksmo.ir <- predict(ir, data.frame(lcat2=lens), type="probs")
row.names(alksmo.ir) <- lens
round(alksmo.ir, 3)




#Plot length frequency by age, observed ALK, and modeled ALK #####
# Multiple options below. 


#Length frequency by age
histStack(lcat2~final_age, data=Agelength_TB, col=gray.colors(9, start=0, end=1), breaks=seq(0,80,1), xlim=c(0,80), ylim=c(0,700), xlab="Total Length (cm)", legend.pos="topright", xaxt="n") #remove col argument if we want it in rainbow format
axis(1, at=seq(0, 80, by=4))


histStack(lcat2~final_age, data=Agelength_CK, col=gray.colors(8, start=0, end=1), breaks=seq(0,80,1), ylim=c(0,80), xlab="Total Length (cm)", legend.pos="topright", xaxt="n") #remove col argument if we want it in rainbow format
axis(1, at=seq(0, 80, by=4))

histStack(lcat2~final_age, data=Agelength_CH, col=gray.colors(9, start=0, end=1), breaks=seq(0,80,1), ylim=c(0,250), xlab="Total Length (cm)", legend.pos="topright", xaxt="n") #remove col argument if we want it in rainbow format
axis(1, at=seq(0, 80, by=4))

histStack(lcat2~final_age, data=Agelength_AP, col=gray.colors(10, start=0, end=1), breaks=seq(0,80,1), ylim=c(0,350), xlab="Total Length (cm)", legend.pos="topright", xaxt="n") #remove col argument if we want it in rainbow format
axis(1, at=seq(0, 80, by=4))

histStack(lcat2~final_age, data=Agelength_IR, col=gray.colors(9, start=0, end=1), breaks=seq(0,80,1), ylim=c(0,500), xlab="Total Length (cm)", legend.pos="topright", xaxt="n") #remove col argument if we want it in rainbow format
axis(1, at=seq(0, 80, by=4))

histStack(lcat2~final_age, data=Agelength_JX, col=gray.colors(8, start=0, end=1), breaks=seq(0,80,1), ylim=c(0,100), xlab="Total Length (cm)", legend.pos="topright", xaxt="n") #remove col argument if we want it in rainbow format
axis(1, at=seq(0, 80, by=4))

####Observed ALK
alkPlot(alk_TB, type="barplot", xlab="Total Length (cm)", showLegend=T, pal="gray") #could remove legend and just reference in figure description

alkPlot(alk_CH, type="barplot", xlab="Total Length (cm)", showLegend=T, pal="gray")     

alkPlot(alk_CK, type="barplot", xlab="Total Length (cm)", showLegend=T, pal="gray")   

alkPlot(alk_AP, type="barplot", xlab="Total Length (cm)", showLegend=T, pal="gray")

alkPlot(alk_IR, type="barplot", xlab="Total Length (cm)", showLegend=T, pal="gray")

alkPlot(alk_JX, type="barplot", xlab="Total Length (cm)", showLegend=T, pal="gray")

 #Other options for plotting
     #obsTB <- alkPlot(alk_TB, type="area", pal="gray", showLegend=TRUE)
     #obsTB <- alkPlot(alk_TB, type="lines", showLegend=TRUE)
     #obsTB <- alkPlot(alk_TB, type="splines", showLegend=TRUE, span=0.1)


####Smoothed ALK####
smoTB <- alkPlot(alksmo.tb, type="barplot", xlab="Total Length (cm)", pal="rainbow", showLegend=TRUE)

smoCK <- alkPlot(alksmo.ck, type="barplot", xlab="Total Length (cm)", pal="rainbow", showLegend=TRUE)

smoCH <- alkPlot(alksmo.ch, type="barplot", xlab="Total Length (cm)", pal="rainbow", showLegend=TRUE)

smoAP <- alkPlot(alksmo.ap, type="barplot", xlab="Total Length (cm)", pal="rainbow", showLegend=TRUE)

smoIR <- alkPlot(alksmo.ir, type="barplot", xlab="Total Length (cm)", pal="rainbow", showLegend=TRUE)

smoJX <- alkPlot(alksmo.jx, type="barplot", xlab="Total Length (cm)", pal="rainbow", showLegend=TRUE)


#AMONG GROUP STATISTICAL COMPARISONS #####
#page 102 in Ogle

Agelength_ALL<- droplevels(subset(as.data.frame(read.csv("ALK_Bay_and_Weight.csv", header=T)),tl>20 & final_age >0 & (bay== "TB" | bay== "CK" | bay== "CH" | bay=="IR" |bay=="AP" | bay=="JX"),select=c(SpecimenNumber, bay, tl, final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1))) # as.fact=TRUE))
mod1 <- multinom(final_age~lcat2, data=Agelength_ALL, maxit=500) #simple model
mod2 <- multinom(final_age~lcat2*bay,data=Agelength_ALL, maxit=500) #more complex model

#likelihood ratio test is computed with anova
test <- anova(mod1, mod2)
# 1       lcat2    219195   66506.76                      
# 2 lcat2 * bay    219105   65007.89 1 vs 2    90 1498.865
# Pr(Chi)

#Drop Bays to test hypothesis of bay influence
      #Null Hypothesis- there is no significant difference in alk between groups

  #removing IR
Agelength_minIR<- droplevels(subset(as.data.frame(read.csv("ALK_Bay_and_Weight.csv", header=T)),tl>20 & final_age >0 & (bay== "TB" | bay== "CK" | bay== "CH" | bay=="AP" | bay=="JX"),select=c(SpecimenNumber, bay, tl, final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1))) # as.fact=TRUE))
mod1 <- multinom(final_age~lcat2, data=Agelength_minIR, maxit=500) #simple model
mod2 <- multinom(final_age~lcat2*bay,data=Agelength_minIR, maxit=500) #more complex model

anova(mod1, mod2)
  #still significantly different

  #now remove JX also
Agelength_minIRJX<- droplevels(subset(as.data.frame(read.csv("ALK_Bay_and_Weight.csv", header=T)),tl>20 & final_age >0 & (bay== "TB" | bay== "CK" | bay== "CH" | bay=="AP"),select=c(SpecimenNumber, bay, tl, final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1))) # as.fact=TRUE))
mod1 <- multinom(final_age~lcat2, data=Agelength_minIRJX, maxit=500) #simple model
mod2 <- multinom(final_age~lcat2*bay,data=Agelength_minIRJX, maxit=500) #more complex model

anova(mod1, mod2)
  #still significantly different

  #now remove AP also
Agelength_minIRJXAP<- droplevels(subset(as.data.frame(read.csv("ALK_Bay_and_Weight.csv", header=T)),tl>20 & final_age >0 & (bay== "TB" | bay== "CK" | bay== "CH"),select=c(SpecimenNumber, bay, tl, final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1))) # as.fact=TRUE))
mod1 <- multinom(final_age~lcat2, data=Agelength_minIRJXAP, maxit=500) #simple model
mod2 <- multinom(final_age~lcat2*bay,data=Agelength_minIRJXAP, maxit=500) #more complex model

anova(mod1, mod2)

#now remove CH also
Agelength_minIRJXAPCH<- droplevels(subset(as.data.frame(read.csv("ALK_Bay_and_Weight.csv", header=T)),tl>20 & final_age >0 & (bay== "TB" | bay== "CK"),select=c(SpecimenNumber, bay, tl, final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1))) # as.fact=TRUE))
mod1 <- multinom(final_age~lcat2, data=Agelength_minIRJXAPCH, maxit=500) #simple model
mod2 <- multinom(final_age~lcat2*bay,data=Agelength_minIRJXAPCH, maxit=500) #more complex model

anova(mod1, mod2)
  #still significantly different

#Bay vs Bay comparison
      #Null Hypothesis- there is no significant difference in alk between groups

#TB vs CK 
Agelength_TBCK<- droplevels(subset(as.data.frame(read.csv("ALK_Bay_and_Weight.csv", header=T)),tl>20 & final_age >0 & (bay== "TB" |  bay== "CK" ),select=c(SpecimenNumber, bay, tl, final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1))) # as.fact=TRUE))
mod1 <- multinom(final_age~lcat2, data=Agelength_TBCK, maxit=500) #simple model
mod2 <- multinom(final_age~lcat2*bay,data=Agelength_TBCK, maxit=500) #more complex model

anova(mod1, mod2)
# The p value for testing the effect of the group in explaining the distribution of lengths within each age is obtained by computing a chi-square test-statistic
#Null Hypothesis- thre is no significant difference in alk between groups
# the likelihood ratio statistic is -2*LL (final value output) of model 1 MINUS
# -2*LL of model 2

#TB vs CH
Agelength_TBCH<- droplevels(subset(as.data.frame(read.csv("ALK_Bay_and_Weight.csv", header=T)),tl>20 & final_age >0 & (bay== "TB" |  bay== "CH" ),select=c(SpecimenNumber, bay, tl, final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1))) # as.fact=TRUE))
mod1 <- multinom(final_age~lcat2, data=Agelength_TBCH, maxit=500) #simple model
mod2 <- multinom(final_age~lcat2*bay,data=Agelength_TBCH, maxit=500) #more complex model

anova(mod1, mod2)

#TB vs AP
Agelength_TBAP<- droplevels(subset(as.data.frame(read.csv("ALK_Bay_and_Weight.csv", header=T)),tl>20 & final_age >0 & (bay== "TB" |  bay== "AP" ),select=c(SpecimenNumber, bay, tl, final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1))) # as.fact=TRUE))
mod1 <- multinom(final_age~lcat2, data=Agelength_TBAP, maxit=500) #simple model
mod2 <- multinom(final_age~lcat2*bay,data=Agelength_TBAP, maxit=500) #more complex model

anova(mod1, mod2)
#An alternative way to do a likelihood ratio test of nested models. 
library(lmtest)
#t <- lrtest(mod1, mod2)

#TB vs JX
Agelength_TBJX<- droplevels(subset(as.data.frame(read.csv("ALK_Bay_and_Weight.csv", header=T)),tl>20 & final_age >0 & (bay== "TB" |  bay== "JX" ),select=c(SpecimenNumber, bay, tl, final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1))) # as.fact=TRUE))
mod1 <- multinom(final_age~lcat2, data=Agelength_TBJX, maxit=500) #simple model
mod2 <- multinom(final_age~lcat2*bay,data=Agelength_TBJX, maxit=500) #more complex model

anova(mod1, mod2)


#TB vs IR
Agelength_TBIR<- droplevels(subset(as.data.frame(read.csv("ALK_Bay_and_Weight.csv", header=T)),tl>20 & final_age >0 & (bay== "TB" |  bay== "IR" ),select=c(SpecimenNumber, bay, tl, final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1))) # as.fact=TRUE))
mod1 <- multinom(final_age~lcat2, data=Agelength_TBIR, maxit=500) #simple model
mod2 <- multinom(final_age~lcat2*bay,data=Agelength_TBIR, maxit=500) #more complex model

anova(mod1, mod2)

#CK vs CH

Agelength_CKCH<- droplevels(subset(as.data.frame(read.csv("ALK_Bay_and_Weight.csv", header=T)),tl>20 & final_age >0 & (bay== "CK" |  bay== "CH" ),select=c(SpecimenNumber, bay, tl, final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1))) # as.fact=TRUE))
mod1 <- multinom(final_age~lcat2, data=Agelength_CKCH, maxit=500) #simple model
mod2 <- multinom(final_age~lcat2*bay,data=Agelength_CKCH, maxit=500) #more complex model

anova(mod1, mod2)

#CK Vs AP
Agelength_CKAP<- droplevels(subset(as.data.frame(read.csv("ALK_Bay_and_Weight.csv", header=T)),tl>20 & final_age >0 & (bay== "CK" |  bay== "AP" ),select=c(SpecimenNumber, bay, tl, final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1))) # as.fact=TRUE))
mod1 <- multinom(final_age~lcat2, data=Agelength_CKAP, maxit=500) #simple model
mod2 <- multinom(final_age~lcat2*bay,data=Agelength_CKAP, maxit=500) #more complex model

anova(mod1, mod2)

#CK Vs JX
Agelength_CKJX<- droplevels(subset(as.data.frame(read.csv("ALK_Bay_and_Weight.csv", header=T)),tl>20 & final_age >0 & (bay== "CK" |  bay== "JX" ),select=c(SpecimenNumber, bay, tl, final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1))) # as.fact=TRUE))
mod1 <- multinom(final_age~lcat2, data=Agelength_CKJX, maxit=500) #simple model
mod2 <- multinom(final_age~lcat2*bay,data=Agelength_CKJX, maxit=500) #more complex model

anova(mod1, mod2)

#CK Vs IR
Agelength_CKIR<- droplevels(subset(as.data.frame(read.csv("ALK_Bay_and_Weight.csv", header=T)),tl>20 & final_age >0 & (bay== "CK" |  bay== "IR" ),select=c(SpecimenNumber, bay, tl, final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1))) # as.fact=TRUE))
mod1 <- multinom(final_age~lcat2, data=Agelength_CKIR, maxit=500) #simple model
mod2 <- multinom(final_age~lcat2*bay,data=Agelength_CKIR, maxit=500) #more complex model

anova(mod1, mod2)


#CH Vs AP
Agelength_CKAP<- droplevels(subset(as.data.frame(read.csv("ALK_Bay_and_Weight.csv", header=T)),tl>20 & final_age >0 & (bay== "CH" |  bay== "AP" ),select=c(SpecimenNumber, bay, tl, final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1))) # as.fact=TRUE))
mod1 <- multinom(final_age~lcat2, data=Agelength_CKAP, maxit=500) #simple model
mod2 <- multinom(final_age~lcat2*bay,data=Agelength_CKAP, maxit=500) #more complex model

anova(mod1, mod2)

#CH Vs JX
Agelength_CKJX<- droplevels(subset(as.data.frame(read.csv("ALK_Bay_and_Weight.csv", header=T)),tl>20 & final_age >0 & (bay== "CH" |  bay== "JX" ),select=c(SpecimenNumber, bay, tl, final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1))) # as.fact=TRUE))
mod1 <- multinom(final_age~lcat2, data=Agelength_CKJX, maxit=500) #simple model
mod2 <- multinom(final_age~lcat2*bay,data=Agelength_CKJX, maxit=500) #more complex model

anova(mod1, mod2)


#CH Vs IR
Agelength_CKIR<- droplevels(subset(as.data.frame(read.csv("ALK_Bay_and_Weight.csv", header=T)),tl>20 & final_age >0 & (bay== "CH" |  bay== "IR" ),select=c(SpecimenNumber, bay, tl, final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1))) # as.fact=TRUE))
mod1 <- multinom(final_age~lcat2, data=Agelength_CKIR, maxit=500) #simple model
mod2 <- multinom(final_age~lcat2*bay,data=Agelength_CKIR, maxit=500) #more complex model

anova(mod1, mod2)

#AP Vs JX
Agelength_APJX<- droplevels(subset(as.data.frame(read.csv("ALK_Bay_and_Weight.csv", header=T)),tl>20 & final_age >0 & (bay== "AP" |  bay== "JX" ),select=c(SpecimenNumber, bay, tl, final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1))) # as.fact=TRUE))
mod1 <- multinom(final_age~lcat2, data=Agelength_APJX, maxit=500) #simple model
mod2 <- multinom(final_age~lcat2*bay,data=Agelength_APJX, maxit=500) #more complex model

anova(mod1, mod2)


#AP Vs IR
Agelength_APIR<- droplevels(subset(as.data.frame(read.csv("ALK_Bay_and_Weight.csv", header=T)),tl>20 & final_age >0 & (bay== "AP" |  bay== "IR" ),select=c(SpecimenNumber, bay, tl, final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1))) # as.fact=TRUE))
mod1 <- multinom(final_age~lcat2, data=Agelength_APIR, maxit=500) #simple model
mod2 <- multinom(final_age~lcat2*bay,data=Agelength_APIR, maxit=500) #more complex model

anova(mod1, mod2)

#JX Vs IR
Agelength_JXIR<- droplevels(subset(as.data.frame(read.csv("ALK_Bay_and_Weight.csv", header=T)),tl>20 & final_age >0 & (bay== "JX" |  bay== "IR" ),select=c(SpecimenNumber, bay, tl, final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1))) # as.fact=TRUE))
mod1 <- multinom(final_age~lcat2, data=Agelength_JXIR, maxit=500) #simple model
mod2 <- multinom(final_age~lcat2*bay,data=Agelength_JXIR, maxit=500) #more complex model

anova(mod1, mod2)

#11. DETERMINE AGE PROPORTION OF ADULTS #####
#of all FIM subset... select sl > = 200 mm which is what was chosen to determine reproductively mature adults
# "Distribution of ages described by either the proportion of fish at age (Pi) or the number of fish at each age (Ni)" Ogle 96
"These calculcations are more easily performed by first computing the Nj values and then dividing each of these by N to compute the Pi values"


TB_adult <- subset(Agelength_TB, sl>= 20)
AP_adult <- subset(Agelength_AP, sl>=20)
CK_adult <- subset(Agelength_CK, sl>=20)
CH_adult <- subset(Agelength_CH, sl>=20)
JX_adult <- subset(Agelength_JX, sl>=20)
IR_adult <- subset(Agelength_IR, sl>=20)

age.n_TB <- xtabs(~final_age, data=TB_adult)
age.n_AP <- xtabs(~final_age, data=AP_adult)
age.n_CK <- xtabs(~final_age, data=CK_adult)
age.n_CH <- xtabs(~final_age, data=CH_adult)
age.n_JX <- xtabs(~final_age, data=JX_adult)
age.n_IR <- xtabs(~final_age, data=IR_adult)

#think I can use prop.table to give the proportion of fish at each age 

prop_TB <- as.data.frame(round(prop.table(age.n_TB),3))
prop_AP <- as.data.frame(round(prop.table(age.n_AP),3))
prop_CK <- as.data.frame(round(prop.table(age.n_CK),3))
prop_CH <- as.data.frame(round(prop.table(age.n_CH),3))
prop_JX <- as.data.frame(round(prop.table(age.n_JX),3))
prop_IR <- as.data.frame(round(prop.table(age.n_IR),3))

#Export proportions at adult age for the FIM catch to use in the Delta_Method script
#For personal computer
write.csv(prop_TB, "~/Desktop/PhD project/Projects/Seatrout/Data/Exported R Dataframes/PropAtAge_TBadult_FIMdata.csv")
write.csv(prop_AP, "~/Desktop/PhD project/Projects/Seatrout/Data/Exported R Dataframes/PropAtAge_APadult_FIMdata.csv")
write.csv(prop_CK, "~/Desktop/PhD project/Projects/Seatrout/Data/Exported R Dataframes/PropAtAge_CKadult_FIMdata.csv")
write.csv(prop_CH, "~/Desktop/PhD project/Projects/Seatrout/Data/Exported R Dataframes/PropAtAge_CHadult_FIMdata.csv")
write.csv(prop_JX, "~/Desktop/PhD project/Projects/Seatrout/Data/Exported R Dataframes/PropAtAge_JXadult_FIMdata.csv")
write.csv(prop_IR, "~/Desktop/PhD project/Projects/Seatrout/Data/Exported R Dataframes/PropAtAge_IRadult_FIMdata.csv")


#For work computer
write.csv(prop_TB, "U:/PhD_projectfiles/Exported_R_Datafiles/PropAtAge_TBadult_FIMdata.csv")
write.csv(prop_AP, "U:/PhD_projectfiles/Exported_R_Datafiles/PropAtAge_APadult_FIMdata.csv")
write.csv(prop_CK, "U:/PhD_projectfiles/Exported_R_Datafiles/PropAtAge_CKadult_FIMdata.csv")
write.csv(prop_CH, "U:/PhD_projectfiles/Exported_R_Datafiles/PropAtAge_CHadult_FIMdata.csv")
write.csv(prop_JX, "U:/PhD_projectfiles/Exported_R_Datafiles/PropAtAge_JXadult_FIMdata.csv")
write.csv(prop_IR, "U:/PhD_projectfiles/Exported_R_Datafiles/PropAtAge_IRadult_FIMdata.csv")

#12. DETERMINE WEIGHT AT AGE SCHEDULE OF ADULTS ######

Weight_at_Age_AP_ad <- Agelength_AP %>% group_by(final_age) %>% summarize(n=length(wt_total), mean_wt=mean(wt_total, na.rm=TRUE), sd=sd(wt_total, na.rm=TRUE), se=se(wt_total, na.rm=TRUE))
Weight_at_Age_CK_ad <- Agelength_CK %>% group_by(final_age) %>% summarize(n=length(wt_total), mean_wt=mean(wt_total, na.rm=TRUE), sd=sd(wt_total, na.rm=TRUE), se=se(wt_total, na.rm=TRUE))
Weight_at_Age_TB_ad <- Agelength_TB %>% group_by(final_age) %>% summarize(n=length(wt_total), mean_wt=mean(wt_total, na.rm=TRUE), sd=sd(wt_total, na.rm=TRUE), se=se(wt_total, na.rm=TRUE))
Weight_at_Age_CH_ad <- Agelength_CH %>% group_by(final_age) %>% summarize(n=length(wt_total), mean_wt=mean(wt_total, na.rm=TRUE), sd=sd(wt_total, na.rm=TRUE), se=se(wt_total, na.rm=TRUE))
Weight_at_Age_JX_ad <- Agelength_JX %>% group_by(final_age) %>% summarize(n=length(wt_total), mean_wt=mean(wt_total, na.rm=TRUE), sd=sd(wt_total, na.rm=TRUE), se=se(wt_total, na.rm=TRUE))
Weight_at_Age_IR_ad <- Agelength_IR %>% group_by(final_age) %>% summarize(n=length(wt_total), mean_wt=mean(wt_total, na.rm=TRUE), sd=sd(wt_total, na.rm=TRUE), se=se(wt_total, na.rm=TRUE))

#Export weight at adult age for the FIM catch to use in the Delta_Method script
#For personal computer
write.csv(Weight_at_Age_AP_ad, "~/Desktop/PhD project/Projects/Seatrout/Data/Exported R Dataframes/Weight_at_Age_APadult_FIMdata.csv")
write.csv(Weight_at_Age_CK_ad, "~/Desktop/PhD project/Projects/Seatrout/Data/Exported R Dataframes/Weight_at_Age_CKadult_FIMdata.csv")
write.csv(Weight_at_Age_TB_ad, "~/Desktop/PhD project/Projects/Seatrout/Data/Exported R Dataframes/Weight_at_Age_TBadult_FIMdata.csv")
write.csv(Weight_at_Age_CH_ad, "~/Desktop/PhD project/Projects/Seatrout/Data/Exported R Dataframes/Weight_at_Age_CHadult_FIMdata.csv")
write.csv(Weight_at_Age_JX_ad, "~/Desktop/PhD project/Projects/Seatrout/Data/Exported R Dataframes/Weight_at_Age_JXadult_FIMdata.csv")
write.csv(Weight_at_Age_IR_ad, "~/Desktop/PhD project/Projects/Seatrout/Data/Exported R Dataframes/Weight_at_Age_IRadult_FIMdata.csv")


#For work computer
write.csv(Weight_at_Age_AP_ad, "U:/PhD_projectfiles/Exported_R_Datafiles/Weight_at_Age_AP_ad_FIMdata.csv")
write.csv(Weight_at_Age_CK_ad, "U:/PhD_projectfiles/Exported_R_Datafiles/Weight_at_Age_CK_ad_FIMdata.csv")
write.csv(Weight_at_Age_TB_ad, "U:/PhD_projectfiles/Exported_R_Datafiles/Weight_at_Age_TB_ad_FIMdata.csv")
write.csv(Weight_at_Age_CH_ad, "U:/PhD_projectfiles/Exported_R_Datafiles/Weight_at_Age_CH_ad_FIMdata.csv")
write.csv(Weight_at_Age_JX_ad, "U:/PhD_projectfiles/Exported_R_Datafiles/Weight_at_Age_JX_ad_FIMdata.csv")
write.csv(Weight_at_Age_IR_ad, "U:/PhD_projectfiles/Exported_R_Datafiles/Weight_at_Age_IR_ad_FIMdata.csv")



# TESTING ####
annuli <- read.csv("U:/PhD_projectfiles/Raw_Data/Age_Length_Data/Access_Tables/annulimeasurements.csv", header=T) %>% dplyr::select(SpecimenNumber, itis, nodc, common_name, Ref)
bio <- read.csv("U:/PhD_projectfiles/Raw_Data/Age_Length_Data/Access_Tables/biology.csv", header=T) %>% rename(SpecimenNumber=specimennumber) %>% dplyr::select(SpecimenNumber, Ref, sl, fl, tl, lengthunits, wt_total, wtunits_total, sex, final_age)
field <- read.csv("U:/PhD_projectfiles/Raw_Data/Age_Length_Data/Access_Tables/field.csv", header=T) %>% dplyr::select(program, date, bay, troutreg, project,Ref)

#annuli and bio link together by specimen number and Ref

#field is linked to bio by Ref

new <- left_join(annuli, bio, by=c("SpecimenNumber", "Ref"))
new2 <- left_join(new, field, by="Ref")

Agelength_TB<- droplevels(subset(new2, bay=="TB" & tl>14 & final_age >0 & program== 'FIM', select=c(sex,SpecimenNumber, bay, tl,sl, final_age, wt_total, date, program))) %>% mutate(tl=tl/10, sl=sl/10, lcat2 =lencat(tl, w=2)) #, as.fact=TRUE))- can include this when determing ALK below but the smoothed ALK needs to be the nonfactored version of the length categorization variable. 
Agelength_TB$sex[which(Agelength_TB$sex == "m")] = "M"
Agelength_TB$sex <- droplevels(Agelength_TB$sex)

#change date format into a factor so that I can do summary statistics by year later on
Agelength_TB$date=as.character(Agelength_TB$date)
Agelength_TB$DateNew = as.POSIXct(strptime(Agelength_TB$date, format="%d-%B-%y", tz=""))  #B is the selection for when month is spelled out
Agelength_TB = mutate(Agelength_TB, year = strftime(DateNew, format="%Y")) %>% dplyr::select(c(-date, -DateNew))
Agelength_TB$year = as.factor(Agelength_TB$year)





MinMeanSEMMax <- function(x) {
  v <- c(mean(x) - sd(x)/sqrt(length(x)), mean(x), mean(x) + sd(x)/sqrt(length(x)))
  names(v) <- c("ymin", "y", "ymax")
  v
}


All_min <- subset(All, sex %in% c("F", "M"))
labels <- c(F="Female", M="Male")


age <- ggplot(All_min, aes(bay, final_age)) +
  stat_summary(fun.data=MinMeanSEMMax, geom="crossbar", colour="black") + 
  scale_y_continuous(breaks=seq(1.5,3,0.5), labels=seq(1.5,3,0.5))+
  facet_grid(sex ~., labeller=labeller(sex=labels)) +
  xlab("Estuary")+
  ylab("Age")+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), 									
        panel.background=element_rect(fill='white', colour='black'),
        axis.title.y = element_text(colour="black", size=20), # changing font of y axis title
        axis.title.x = element_text(colour="black", size=20),
        axis.text.x=element_text(colour="black", size=16), #changing  colour and font of x axis text
        axis.text.y=element_text(colour="black", size=16),
        strip.text.y = element_text(size=16))  #changing colour and font of y axis

















# OLD ANOVAS #####
# ONE WAY ANOVA to determine if there are significant differences among estuary-specific mean lengths ####


plot(tl~bay, data=All)
results=aov(tl~bay, data=All)
summary(results)
TukeyHSD(results, conf.level=0.95)
#deterine if there is significant difference in overall age among estuaries

#determine if there is significnat difference in female length among estuaries
All_F <- subset(All, sex=="F")
plot(tl~bay, data=All_F)
results=aov(tl~bay, data=All_F)
summary(results)

# Df  Sum Sq Mean Sq F value Pr(>F)    
# bay             5   33109    6622   89.58 <2e-16 ***
#   Residuals   19505 1441763      74   

# Multiple comparisons with pairwise.t.test to compute pair wise comparisons between group means with corrections for multiple testing to determine how the means of each differ. 
TukeyHSD(results, conf.level=0.95)

#determine if there is significnat difference in male length among estuaries
All_M <- subset(All, sex=="M")
plot(tl~bay, data=All_M)
results=aov(tl~bay, data=All_M)
summary(results)

# Df Sum Sq Mean Sq F value Pr(>F)    
# bay            5   8881  1776.3   40.77 <2e-16 ***
#   Residuals   7762 338222    43.6                   
# ---
TukeyHSD(results, conf.level=0.95)


#ANOVA determine if there is signficant difference in all age among estuaries ####
plot(final_age~bay, data=All)
results=aov(final_age~bay, data=All)
summary(results)

TukeyHSD(results, conf.level=0.95)


# ANOVA determine if there is significnat difference in female age among estuaries ####
All_F <- subset(All, sex=="F")
plot(final_age~bay, data=All_F)
results=aov(final_age~bay, data=All_F)
summary(results)

# Df Sum Sq Mean Sq F value Pr(>F)    
# bay             5   1223  244.53   150.6 <2e-16 ***
#   Residuals   19505  31664    1.62     

TukeyHSD(results, conf.level=0.95)

# ANOVA determine if there is significnat difference in male age among estuaries####
All_M <- subset(All, sex=="M")
plot(final_age~bay, data=All_M)
results=aov(final_age~bay, data=All_M)
summary(results)

# Df Sum Sq Mean Sq F value Pr(>F)    
# bay            5    400   80.08   45.81 <2e-16 ***
#   Residuals   7762  13568    1.75    

TukeyHSD(results, conf.level=0.95)

# ANOVA determine whether there is a significant difference in male and female tl (not estuary specific) ####
#this should actually be a t-test because thats basically what it is. Comparing two means
All_MF <- droplevels(subset(All, sex=="F" | sex=="M"))
plot(tl~sex, data=All_MF)
result= aov(tl~ sex, data=All_MF)
summary(result)
TukeyHSD(result,conf.level=0.95)

# ANOVA determine whethere there is a significant difference in male and female TOTAL LENGTH for each estuary ####
Agelength_TB_M <- subset(Agelength_TB, sex=="M")
Agelength_TB_F <- subset(Agelength_TB, sex=="F")
chisq.test(Agelength_TB_M$tl, Agelength_TB_F$tl)

Agelength_AP_M <- subset(Agelength_AP, sex=="M")
Agelength_AP_F <- subset(Agelength_AP, sex=="F")
t.test(Agelength_AP_M$tl, Agelength_AP_F$tl)

Agelength_CK_M <- subset(Agelength_CK, sex=="M")
Agelength_CK_F <- subset(Agelength_CK, sex=="F")
t.test(Agelength_CK_M$tl, Agelength_CK_F$tl)

Agelength_CH_M <- subset(Agelength_CH, sex=="M")
Agelength_CH_F <- subset(Agelength_CH, sex=="F")
t.test(Agelength_CH_M$tl, Agelength_CH_F$tl)

Agelength_IR_M <- subset(Agelength_IR, sex=="M")
Agelength_IR_F <- subset(Agelength_IR, sex=="F")
t.test(Agelength_IR_M$tl, Agelength_IR_F$tl)

Agelength_JX_M <- subset(Agelength_JX, sex=="M")
Agelength_JX_F <- subset(Agelength_JX, sex=="F")
t.test(Agelength_JX_M$tl, Agelength_JX_F$tl)


# ANOVA determine whethere there is a significant difference in male and female AGE for each estuary ####
Agelength_TB_M <- subset(Agelength_TB, sex=="M")
Agelength_TB_F <- subset(Agelength_TB, sex=="F")
t.test(Agelength_TB_M$final_age, Agelength_TB_F$final_age)

Agelength_AP_M <- subset(Agelength_AP, sex=="M")
Agelength_AP_F <- subset(Agelength_AP, sex=="F")
t.test(Agelength_AP_M$final_age, Agelength_AP_F$final_age)

Agelength_CK_M <- subset(Agelength_CK, sex=="M")
Agelength_CK_F <- subset(Agelength_CK, sex=="F")
t.test(Agelength_CK_M$final_age, Agelength_CK_F$final_age)

Agelength_CH_M <- subset(Agelength_CH, sex=="M")
Agelength_CH_F <- subset(Agelength_CH, sex=="F")
t.test(Agelength_CH_M$final_age, Agelength_CH_F$final_age)

Agelength_IR_M <- subset(Agelength_IR, sex=="M")
Agelength_IR_F <- subset(Agelength_IR, sex=="F")
t.test(Agelength_IR_M$final_age, Agelength_IR_F$final_age)

Agelength_JX_M <- subset(Agelength_JX, sex=="M")
Agelength_JX_F <- subset(Agelength_JX, sex=="F")
t.test(Agelength_JX_M$final_age, Agelength_JX_F$final_age)



plot(final_age ~ sex +bay, data=All_MF)
result= aov(final_age ~ sex, data=All_MF)
summary(results)
TukeyHSD(results,  conf.level=0.95)


