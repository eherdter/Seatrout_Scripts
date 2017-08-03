# ABOUT ####
# 10/10/2016 This script imports ALK with Bay.xlsx
# Main Objectives of this script: 
# 1. Imports otolith data. Summarizes data, mean age, mean length, sd and se and total sample number. PROJECT TYPE!
# 2. Makes bay-specific observed ALK, calculates some summary statistics.  
# 3. Makes bay specific smoothed (modeled) ALK with multinomial modeling methods in Ogle (87-).
# 4. Likelihood ratio testing to do among group statistical comparisons - Ogle (102-103)
# 5. Plots observed and smoothed ALK for each bay. 
# 6. Calculates proportional age distribution
# 7. Determines mean length -at-age (otolith database) and produces plots
# 8. One way anova to determine if the mean lengths (and mean ages) of males and females are significantly different among estuaries
# 9. T test to determine whether there is a significant difference between male and female age for each estuary
# 10. T test to determine whether there is a significant difference in male and female length for each estuary.
#################################################################
library(FSA)
library(magrittr)
library(dplyr)
library(nnet)
library(plotrix)
library(haven)
library(ggplot2)
library(scales)

# set working directory
setwd("~/Desktop/PhD project/Projects/Seatrout/Data")

# LOAD DATA ####
#load the csv file
# subset by which bay I want
# make sure I have just the "aged sample"
# turn mm to cm
# select just a few variables to make it more manageable
# then drop the remaining bay levels still sticking around (droplevels)
# turn tl from mm to cm
# create length categories with FSA package

test <- read.csv("ALK with Bay.csv", header=T)

Agelength_TB<- droplevels(subset(as.data.frame(read.csv("ALK with Bay.csv", header=T)), bay=="TB" & tl>14 & final_age >0 & program== 'FIM', select=c(sex,specimennumber, bay, tl, final_age, Date))) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1)) #, as.fact=TRUE))- can include this when determing ALK below but the smoothed ALK needs to be the nonfactored version of the length categorization variable. 
Agelength_TB$sex[which(Agelength_TB$sex == "m")] = "M"
Agelength_TB$sex <- droplevels(Agelength_TB$sex)

#change date format into a factor so that I can do summary statistics by year later on
Agelength_TB$Date=as.character(Agelength_TB$Date)
Agelength_TB$DateNew = as.POSIXct(strptime(Agelength_TB$Date, format="%m/%d/%y", tz="")) 
Agelength_TB = mutate(Agelength_TB, year = strftime(DateNew, format="%Y")) %>%select(-Date, -DateNew)
Agelength_TB$year = as.factor(Agelength_TB$year) 

Agelength_AP<- droplevels(subset(as.data.frame(read.csv("ALK with Bay.csv", header=T)), bay=="AP" & tl>0 & final_age >0 & program== 'FIM', select=c(sex,specimennumber, bay, tl, final_age, Date))) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1)) #, as.fact=TRUE))

Agelength_AP$Date=as.character(Agelength_AP$Date)
Agelength_AP$DateNew = as.POSIXct(strptime(Agelength_AP$Date, format="%m/%d/%y", tz="")) 
Agelength_AP = mutate(Agelength_AP, year = strftime(DateNew, format="%Y")) %>%select(-Date, -DateNew)
Agelength_AP$year = as.factor(Agelength_AP$year) 

Agelength_CK<- droplevels(subset(as.data.frame(read.csv("ALK with Bay.csv", header=T)), bay=="CK" & tl>0 & final_age >0 & program== 'FIM', select=c(sex,specimennumber, bay, tl, final_age, Date))) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1)) # as.fact=TRUE))

Agelength_CK$Date=as.character(Agelength_CK$Date)
Agelength_CK$DateNew = as.POSIXct(strptime(Agelength_CK$Date, format="%m/%d/%y", tz="")) 
Agelength_CK = mutate(Agelength_CK, year = strftime(DateNew, format="%Y")) %>%select(-Date, -DateNew)
Agelength_CK$year = as.factor(Agelength_CK$year) 

Agelength_CH<- droplevels(subset(as.data.frame(read.csv("ALK with Bay.csv", header=T)), bay=="CH" & tl>0 & final_age >0 & program== 'FIM', select=c(sex,specimennumber, bay, tl, final_age, Date))) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1)) # as.fact=TRUE))
Agelength_CH$sex[which(Agelength_CH$sex == "f")] = "F"
Agelength_CH$sex <- droplevels(Agelength_CH$sex)

Agelength_CH$Date=as.character(Agelength_CH$Date)
Agelength_CH$DateNew = as.POSIXct(strptime(Agelength_CH$Date, format="%m/%d/%y", tz="")) 
Agelength_CH = mutate(Agelength_CH, year = strftime(DateNew, format="%Y")) %>%select(-Date, -DateNew)
Agelength_CH$year = as.factor(Agelength_CH$year) 

Agelength_IR<- droplevels(subset(as.data.frame(read.csv("ALK with Bay.csv", header=T)), bay=="IR" & tl>0 & final_age >0 & program== 'FIM', select=c(sex,specimennumber, bay, tl, final_age, Date))) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1)) # as.fact=TRUE))

Agelength_IR$Date=as.character(Agelength_IR$Date)
Agelength_IR$DateNew = as.POSIXct(strptime(Agelength_IR$Date, format="%m/%d/%y", tz="")) 
Agelength_IR = mutate(Agelength_IR, year = strftime(DateNew, format="%Y")) %>%select(-Date, -DateNew)
Agelength_IR$year = as.factor(Agelength_IR$year) 

Agelength_JX<- droplevels(subset(as.data.frame(read.csv("ALK with Bay.csv", header=T)), bay=="JX" & tl>0 & final_age >0 & program== 'FIM', select=c(sex,specimennumber, bay, tl, final_age, Date))) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1)) # as.fact=TRUE))

Agelength_JX$Date=as.character(Agelength_JX$Date)
Agelength_JX$DateNew = as.POSIXct(strptime(Agelength_JX$Date, format="%m/%d/%y", tz="")) 
Agelength_JX = mutate(Agelength_JX, year = strftime(DateNew, format="%Y")) %>%select(-Date, -DateNew)
Agelength_JX$year = as.factor(Agelength_JX$year) 


# BASIC SUMMARIZATION ####
#total sample number of FIM data

N= nrow(Agelength_TB) + nrow(Agelength_AP) +nrow(Agelength_CK) +nrow(Agelength_CH) + nrow(Agelength_IR) +nrow(Agelength_JX)

# Added together
All= rbind(Agelength_TB, Agelength_AP, Agelength_CK, Agelength_CH, Agelength_IR, Agelength_JX) 

#Age proportion
All_3under =subset(All, final_age <=3)
Proprotion3under <- nrow(All_3under)/nrow(All)

#Min, max length
min(All$tl)
max(All$tl)

#summarize the entire set and group by sex and then bay
Combined_sum <- summarize(group_by(All, sex,bay), mean_tl=mean(tl), sd_tl=sd(tl), se_tl= sd_tl/(sqrt(length(final_age))),min_tl=min(tl), max_tl=max(tl), mean_age=mean(final_age), median_age=median(final_age),sd_age= sd(final_age), se_age=sd_age/sqrt(length(final_age)), min_age= min(final_age), max_age=max(final_age))

Combined_sum_year <- summarize(group_by(All,bay ), min_year=min(year), max_year=max(year))

#summarize each to get sex ratios etc. 
TB_sum <- summarise(Agelength_TB,mean_age=mean(final_age), median_age=median(final_age),sd_age= sd(final_age), se_age=sd_age/sqrt(length(final_age)), mean_tl = mean(tl),min_tl=min(tl), max_tl=max(tl), sd_tl= sd(tl), se_tl= sd_tl/(sqrt(length(final_age))), prop_F=nrow(subset(Agelength_TB, sex=="F"))/length(final_age),prop_M=nrow(subset(Agelength_TB, sex=="M"))/length(final_age), N_F=nrow(subset(Agelength_TB, sex=="F")), N_M=nrow(subset(Agelength_TB, sex=="M")))
AP_sum <- summarise(Agelength_AP,mean_age=mean(final_age), median_age=median(final_age),sd_age= sd(final_age), se_age=sd_age/sqrt(length(final_age)), mean_tl = mean(tl), sd_tl= sd(tl), se_tl= sd_tl/(sqrt(length(final_age))), prop_F=nrow(subset(Agelength_AP, sex=="F"))/length(final_age),prop_M=nrow(subset(Agelength_AP, sex=="M"))/length(final_age), N_F=nrow(subset(Agelength_AP, sex=="F")), N_M=nrow(subset(Agelength_AP, sex=="M")))
CH_sum <- summarise(Agelength_CH,mean_age=mean(final_age), median_age=median(final_age),sd_age= sd(final_age), se_age=sd_age/sqrt(length(final_age)), mean_tl = mean(tl), sd_tl= sd(tl), se_tl= sd_tl/(sqrt(length(final_age))), prop_F=nrow(subset(Agelength_CH, sex=="F"))/length(final_age),prop_M=nrow(subset(Agelength_CH, sex=="M"))/length(final_age), N_F=nrow(subset(Agelength_CH, sex=="F")), N_M=nrow(subset(Agelength_CH, sex=="M")))
CK_sum <- summarise(Agelength_CK,mean_age=mean(final_age), median_age=median(final_age),sd_age= sd(final_age), se_age=sd_age/sqrt(length(final_age)), mean_tl = mean(tl), sd_tl= sd(tl), se_tl= sd_tl/(sqrt(length(final_age))), prop_F=nrow(subset(Agelength_CK, sex=="F"))/length(final_age),prop_M=nrow(subset(Agelength_CK, sex=="M"))/length(final_age), N_F=nrow(subset(Agelength_CK, sex=="F")), N_M=nrow(subset(Agelength_CK, sex=="M")))
IR_sum <- summarise(Agelength_IR,mean_age=mean(final_age), median_age=median(final_age),sd_age= sd(final_age), se_age=sd_age/sqrt(length(final_age)), mean_tl = mean(tl), sd_tl= sd(tl), se_tl= sd_tl/(sqrt(length(final_age))), prop_F=nrow(subset(Agelength_IR, sex=="F"))/length(final_age),prop_M=nrow(subset(Agelength_IR, sex=="M"))/length(final_age), N_F=nrow(subset(Agelength_IR, sex=="F")), N_M=nrow(subset(Agelength_IR, sex=="M")))
JX_sum <- summarise(Agelength_JX,mean_age=mean(final_age), median_age=median(final_age),sd_age= sd(final_age), se_age=sd_age/sqrt(length(final_age)), mean_tl = mean(tl), sd_tl= sd(tl), se_tl= sd_tl/(sqrt(length(final_age))), prop_F=nrow(subset(Agelength_JX, sex=="F"))/length(final_age),prop_M=nrow(subset(Agelength_JX, sex=="M"))/length(final_age), N_F=nrow(subset(Agelength_JX, sex=="F")), N_M=nrow(subset(Agelength_JX, sex=="M")))
All_sum <- rbind(TB_sum, AP_sum, CH_sum, CK_sum, IR_sum, JX_sum) %>% mutate(bay=c("TB", "AP", "CH", "CK", "IR", "JX"))





# ONE WAY ANOVA to determine if there are significant differences among estuary-specific mean lengths ####
##################################
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
t.test(Agelength_TB_M$tl, Agelength_TB_F$tl)

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

# MAKE ALKS ####
# Make table with observed total numbers at length by age 
###########################################################
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



############################################
# Make observed ALK from the above tables.
############################################
#The conditional proportions that form the ALK are calculated by dividing ecah cell of the frequency table by the sum of the corresponding row. 
#These row proportions are constructed by submitting the xtabs() object to prop.table() and including margin=1 to indicate that the proportions are computed by row (page 92). 

#The alkPlot command used for plotting the observed ALK is unable to extend the x axis to the bounds of c(0,80) because xlim is not working. 
# Therefore, in order to produce a plot with an x axis that can span from 0-80 (like what is happening with the length frequency and the smoothed ALK)
# I need to add in "observed" proportions for length categories that were not sampled. 
# I could have added them to the original data frame but I was concerned that in the process of proportion calculations the extra entries would
# affect the proportions or result in proportions that I didn't want. The smoothing process can estimate proportions outside of the range but I wanted 
# to keep the observed plot with just the proportions from the observed data. Therefore, I added in the zero proportion data by writing and editing
# a csv file which I then read below. 

as.data.frame.matrix((prop.table(rawfreq_TB, margin=1))) %>% write.csv("alk_TB.csv")
alk_TB <- read.csv("alk_TB_edit.csv", row.names=1)
names(alk_TB) <- c(1,2,3,4,5,6,7,8,9)
round(alk_TB,3)

as.data.frame.matrix((prop.table(rawfreq_CH, margin=1))) %>% write.csv("alk_CH.csv")
alk_CH <- read.csv("alk_CH_edit.csv", row.names=1)
names(alk_CH) <- c(1,2,3,4,5,6,7,8,9)
round(alk_CH,3)

as.data.frame.matrix((prop.table(rawfreq_CK, margin=1))) %>% write.csv("alk_CK.csv")
alk_CK <- read.csv("alk_CK_edit.csv", row.names=1)
names(alk_CK) <- c(1,2,3,4,5,6,7,8)
round(alk_CK,3)

as.data.frame.matrix((prop.table(rawfreq_AP, margin=1))) %>% write.csv("alk_AP.csv")
alk_AP <- read.csv("alk_AP_edit.csv", row.names=1)
names(alk_AP) <- c(1,2,3,4,5,6,7,8,9,10)
round(alk_AP,3)

as.data.frame.matrix((prop.table(rawfreq_JX, margin=1))) %>% write.csv("alk_JX.csv")
alk_JX <- read.csv("alk_JX_edit.csv", row.names=1)
names(alk_JX) <- c(1,2,3,4,5,6,7,8)
round(alk_JX,3)

as.data.frame.matrix((prop.table(rawfreq_IR, margin=1))) %>% write.csv("alk_IR.csv")
alk_IR <- read.csv("alk_IR_edit.csv", row.names=1)
names(alk_IR) <- c(1,2,3,4,5,6,7,8,9)
round(alk_IR,3)


######################################################################
# Apply.
# 1. Proportional age distribution
# 1. Mean length-at-age. (otolith database)
######################################################################

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
  scale_y_continuous(labels=percent_format(), name="Frequency (%)", limits=c(0,.4))+  ## plotting in percent frequency
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
  scale_y_continuous(labels=percent_format(), name="Frequency (%)", limits=c(0,.4))+  ## plotting in percent frequency
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












# Plot of individual lengths at age with mean lengths at age superimposed. 

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

###############################################################################################################
# Produce smoothed ALK from multinomial modeling exercise which can be used to do Likelihood ratio testing. 
#MODELED AGE LENGTH KEYS (aka SMOOTHED ALK)
#-fixes two common issues with Age Length keys detailed on page 92 Ogle. 
#multinomial logistic regression model -Gerritsen et al. 2006. The response variable has more than two levels
###############################################################################################################
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



#################################################
#Plot length frequency by age, observed ALK, and modeled ALK
# Multiple options below. 
#################################################

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

####Observed ALK###
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

#########################################
#AMONG GROUP STATISTICAL COMPARISONS
#page 102 in Ogle
##########################################
Agelength_ALL<- droplevels(subset(as.data.frame(read.csv("ALK with Bay.csv", header=T)),tl>20 & final_age >0 & (bay== "TB" | bay== "CK" | bay== "CH" | bay=="IR" |bay=="AP" | bay=="JX"),select=c(specimennumber, bay, tl, final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1))) # as.fact=TRUE))
mod1 <- multinom(final_age~lcat2, data=Agelength_ALL, maxit=500) #simple model
mod2 <- multinom(final_age~lcat2*bay,data=Agelength_ALL, maxit=500) #more complex model

#likelihood ratio test is computed with anova
test <- anova(mod1, mod2)

#Drop Bays to test hypothesis of bay influence
      #Null Hypothesis- there is no significant difference in alk between groups

  #removing IR
Agelength_minIR<- droplevels(subset(as.data.frame(read.csv("ALK with Bay.csv", header=T)),tl>20 & final_age >0 & (bay== "TB" | bay== "CK" | bay== "CH" | bay=="AP" | bay=="JX"),select=c(specimennumber, bay, tl, final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1))) # as.fact=TRUE))
mod1 <- multinom(final_age~lcat2, data=Agelength_minIR, maxit=500) #simple model
mod2 <- multinom(final_age~lcat2*bay,data=Agelength_minIR, maxit=500) #more complex model

anova(mod1, mod2)
  #still significantly different

  #now remove JX also
Agelength_minIRJX<- droplevels(subset(as.data.frame(read.csv("ALK with Bay.csv", header=T)),tl>20 & final_age >0 & (bay== "TB" | bay== "CK" | bay== "CH" | bay=="AP"),select=c(specimennumber, bay, tl, final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1))) # as.fact=TRUE))
mod1 <- multinom(final_age~lcat2, data=Agelength_minIRJX, maxit=500) #simple model
mod2 <- multinom(final_age~lcat2*bay,data=Agelength_minIRJX, maxit=500) #more complex model

anova(mod1, mod2)
  #still significantly different

  #now remove AP also
Agelength_minIRJXAP<- droplevels(subset(as.data.frame(read.csv("ALK with Bay.csv", header=T)),tl>20 & final_age >0 & (bay== "TB" | bay== "CK" | bay== "CH"),select=c(specimennumber, bay, tl, final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1))) # as.fact=TRUE))
mod1 <- multinom(final_age~lcat2, data=Agelength_minIRJXAP, maxit=500) #simple model
mod2 <- multinom(final_age~lcat2*bay,data=Agelength_minIRJXAP, maxit=500) #more complex model

anova(mod1, mod2)

#now remove CH also
Agelength_minIRJXAPCH<- droplevels(subset(as.data.frame(read.csv("ALK with Bay.csv", header=T)),tl>20 & final_age >0 & (bay== "TB" | bay== "CK"),select=c(specimennumber, bay, tl, final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1))) # as.fact=TRUE))
mod1 <- multinom(final_age~lcat2, data=Agelength_minIRJXAPCH, maxit=500) #simple model
mod2 <- multinom(final_age~lcat2*bay,data=Agelength_minIRJXAPCH, maxit=500) #more complex model

anova(mod1, mod2)
  #still significantly different

#Bay vs Bay comparison
      #Null Hypothesis- there is no significant difference in alk between groups

#TB vs CK 
Agelength_TBCK<- droplevels(subset(as.data.frame(read.csv("ALK with Bay.csv", header=T)),tl>20 & final_age >0 & (bay== "TB" |  bay== "CK" ),select=c(specimennumber, bay, tl, final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1))) # as.fact=TRUE))
mod1 <- multinom(final_age~lcat2, data=Agelength_TBCK, maxit=500) #simple model
mod2 <- multinom(final_age~lcat2*bay,data=Agelength_TBCK, maxit=500) #more complex model

anova(mod1, mod2)
# The p value for testing the effect of the group in explaining the distribution of lengths within each age is obtained by computing a chi-square test-statistic
#Null Hypothesis- thre is no significant difference in alk between groups
# the likelihood ratio statistic is -2*LL (final value output) of model 1 MINUS
# -2*LL of model 2

#TB vs CH
Agelength_TBCH<- droplevels(subset(as.data.frame(read.csv("ALK with Bay.csv", header=T)),tl>20 & final_age >0 & (bay== "TB" |  bay== "CH" ),select=c(specimennumber, bay, tl, final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1))) # as.fact=TRUE))
mod1 <- multinom(final_age~lcat2, data=Agelength_TBCH, maxit=500) #simple model
mod2 <- multinom(final_age~lcat2*bay,data=Agelength_TBCH, maxit=500) #more complex model

anova(mod1, mod2)

#TB vs AP
Agelength_TBAP<- droplevels(subset(as.data.frame(read.csv("ALK with Bay.csv", header=T)),tl>20 & final_age >0 & (bay== "TB" |  bay== "AP" ),select=c(specimennumber, bay, tl, final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1))) # as.fact=TRUE))
mod1 <- multinom(final_age~lcat2, data=Agelength_TBAP, maxit=500) #simple model
mod2 <- multinom(final_age~lcat2*bay,data=Agelength_TBAP, maxit=500) #more complex model

anova(mod1, mod2)
#An alternative way to do a likelihood ratio test of nested models. 
library(lmtest)
#t <- lrtest(mod1, mod2)

#TB vs JX
Agelength_TBJX<- droplevels(subset(as.data.frame(read.csv("ALK with Bay.csv", header=T)),tl>20 & final_age >0 & (bay== "TB" |  bay== "JX" ),select=c(specimennumber, bay, tl, final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1))) # as.fact=TRUE))
mod1 <- multinom(final_age~lcat2, data=Agelength_TBJX, maxit=500) #simple model
mod2 <- multinom(final_age~lcat2*bay,data=Agelength_TBJX, maxit=500) #more complex model

lrtest(mod1, mod2)


#TB vs IR
Agelength_TBIR<- droplevels(subset(as.data.frame(read.csv("ALK with Bay.csv", header=T)),tl>20 & final_age >0 & (bay== "TB" |  bay== "IR" ),select=c(specimennumber, bay, tl, final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1))) # as.fact=TRUE))
mod1 <- multinom(final_age~lcat2, data=Agelength_TBIR, maxit=500) #simple model
mod2 <- multinom(final_age~lcat2*bay,data=Agelength_TBIR, maxit=500) #more complex model

anova(mod1, mod2)

#CK vs CH

Agelength_CKCH<- droplevels(subset(as.data.frame(read.csv("ALK with Bay.csv", header=T)),tl>20 & final_age >0 & (bay== "CK" |  bay== "CH" ),select=c(specimennumber, bay, tl, final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1))) # as.fact=TRUE))
mod1 <- multinom(final_age~lcat2, data=Agelength_CKCH, maxit=500) #simple model
mod2 <- multinom(final_age~lcat2*bay,data=Agelength_CKCH, maxit=500) #more complex model

anova(mod1, mod2)

#CK Vs AP
Agelength_CKAP<- droplevels(subset(as.data.frame(read.csv("ALK with Bay.csv", header=T)),tl>20 & final_age >0 & (bay== "CK" |  bay== "AP" ),select=c(specimennumber, bay, tl, final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1))) # as.fact=TRUE))
mod1 <- multinom(final_age~lcat2, data=Agelength_CKAP, maxit=500) #simple model
mod2 <- multinom(final_age~lcat2*bay,data=Agelength_CKAP, maxit=500) #more complex model

anova(mod1, mod2)

#CK Vs JX
Agelength_CKJX<- droplevels(subset(as.data.frame(read.csv("ALK with Bay.csv", header=T)),tl>20 & final_age >0 & (bay== "CK" |  bay== "JX" ),select=c(specimennumber, bay, tl, final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1))) # as.fact=TRUE))
mod1 <- multinom(final_age~lcat2, data=Agelength_CKJX, maxit=500) #simple model
mod2 <- multinom(final_age~lcat2*bay,data=Agelength_CKJX, maxit=500) #more complex model

anova(mod1, mod2)

#CK Vs IR
Agelength_CKIR<- droplevels(subset(as.data.frame(read.csv("ALK with Bay.csv", header=T)),tl>20 & final_age >0 & (bay== "CK" |  bay== "IR" ),select=c(specimennumber, bay, tl, final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1))) # as.fact=TRUE))
mod1 <- multinom(final_age~lcat2, data=Agelength_CKIR, maxit=500) #simple model
mod2 <- multinom(final_age~lcat2*bay,data=Agelength_CKIR, maxit=500) #more complex model

anova(mod1, mod2)


#CH Vs AP
Agelength_CKAP<- droplevels(subset(as.data.frame(read.csv("ALK with Bay.csv", header=T)),tl>20 & final_age >0 & (bay== "CH" |  bay== "AP" ),select=c(specimennumber, bay, tl, final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1))) # as.fact=TRUE))
mod1 <- multinom(final_age~lcat2, data=Agelength_CKAP, maxit=500) #simple model
mod2 <- multinom(final_age~lcat2*bay,data=Agelength_CKAP, maxit=500) #more complex model

anova(mod1, mod2)

#CH Vs JX
Agelength_CKJX<- droplevels(subset(as.data.frame(read.csv("ALK with Bay.csv", header=T)),tl>20 & final_age >0 & (bay== "CH" |  bay== "JX" ),select=c(specimennumber, bay, tl, final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1))) # as.fact=TRUE))
mod1 <- multinom(final_age~lcat2, data=Agelength_CKJX, maxit=500) #simple model
mod2 <- multinom(final_age~lcat2*bay,data=Agelength_CKJX, maxit=500) #more complex model

anova(mod1, mod2)


#CH Vs IR
Agelength_CKIR<- droplevels(subset(as.data.frame(read.csv("ALK with Bay.csv", header=T)),tl>20 & final_age >0 & (bay== "CH" |  bay== "IR" ),select=c(specimennumber, bay, tl, final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1))) # as.fact=TRUE))
mod1 <- multinom(final_age~lcat2, data=Agelength_CKIR, maxit=500) #simple model
mod2 <- multinom(final_age~lcat2*bay,data=Agelength_CKIR, maxit=500) #more complex model

anova(mod1, mod2)

#AP Vs JX
Agelength_APJX<- droplevels(subset(as.data.frame(read.csv("ALK with Bay.csv", header=T)),tl>20 & final_age >0 & (bay== "AP" |  bay== "JX" ),select=c(specimennumber, bay, tl, final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1))) # as.fact=TRUE))
mod1 <- multinom(final_age~lcat2, data=Agelength_APJX, maxit=500) #simple model
mod2 <- multinom(final_age~lcat2*bay,data=Agelength_APJX, maxit=500) #more complex model

anova(mod1, mod2)


#AP Vs IR
Agelength_APIR<- droplevels(subset(as.data.frame(read.csv("ALK with Bay.csv", header=T)),tl>20 & final_age >0 & (bay== "AP" |  bay== "IR" ),select=c(specimennumber, bay, tl, final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1))) # as.fact=TRUE))
mod1 <- multinom(final_age~lcat2, data=Agelength_APIR, maxit=500) #simple model
mod2 <- multinom(final_age~lcat2*bay,data=Agelength_APIR, maxit=500) #more complex model

anova(mod1, mod2)

#JX Vs IR
Agelength_JXIR<- droplevels(subset(as.data.frame(read.csv("ALK with Bay.csv", header=T)),tl>20 & final_age >0 & (bay== "JX" |  bay== "IR" ),select=c(specimennumber, bay, tl, final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1))) # as.fact=TRUE))
mod1 <- multinom(final_age~lcat2, data=Agelength_JXIR, maxit=500) #simple model
mod2 <- multinom(final_age~lcat2*bay,data=Agelength_JXIR, maxit=500) #more complex model

anova(mod1, mod2)

########################################
# ONE WAY ANOVA to determine if there are significant differences among estuary-specific mean lengths
##################################

