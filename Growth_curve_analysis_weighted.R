# ABOUT ####
# 10/13/2016 This script imports ALK with Bay.xlsx
# Main Objectives of this script: 
# 1. Imports otolith data. 
# 2. Determines VBGF parameters for each bay and bootstrapped 95% confidence intervals for each.
# 3. Plots the comparative growth curves
# 4. Compares nested growth models for Bay to Bay comparison to determine which growth parameters are significantly different between bays. 
###############################################################
rm(list=ls())
# load packages####
library(FSA)
library(magrittr)
library(dplyr)
library(nlstools)
library(tidyverse)
library(gridExtra)
library(grid)
library(lattice)

#set working directory ####
setwd("~/Desktop/PhD project/Projects/Seatrout/Data")
setwd("U:/PhD_projectfiles/Raw_Data/Age_Length_Data")

#load the age 0 length data from the McMichael peters study ####
# change some variables and add some dummy variables so the dataframe will match to the Agelength_Bay loaded next 
age0 <- read.csv("Age_length_mcmichael_peters.csv", header=TRUE) %>% mutate(tl=TL..mm./10, lcat2 =lencat(tl, w=1), bay=rep(NA, 149), date=rep(NA, 149)) %>% rename(final_age=Otolith.age..years.) 
age0 <- subset(age0, select=-c(TL..mm.))
sex = sample(c("F", "M"),149,replace=TRUE, prob=c(0.5, 0.5))
age0$sex = sex

#load the csv file #####
# subset by which bay I want
# make sure I have just the "aged sample"
# turn mm to cm
# select just a few variables to make it more manageable
# then drop the remaining bay levels still sticking around (droplevels)
# turn tl from mm to cm
# create length categories with FSA package
#age is >0
# IRpend age0 dataframe with NA bay
# then turn the NA bay into the specific bay that would match for each

test <- read.csv("ALK_Bay_and_weight.csv", header=T)
Agelength_TB<- droplevels(subset(as.data.frame(read.csv("ALK_Bay_and_weight.csv", header=T)), bay=="TB" & tl>14 & final_age >0 & program=='FIM', select=c( bay, tl, final_age, date, sex))) %>% filter(!is.na(final_age)) %>% 
  mutate(tl=tl/10, lcat2 =lencat(tl, w=1)) %>%rbind(age0)
Agelength_TB$bay = "TB"
Agelength_TB <- subset(Agelength_TB, sex %in% c("M", "F"))

Agelength_AP<- droplevels(subset(as.data.frame(read.csv("ALK_Bay_and_weight.csv", header=T)), bay=="AP" & tl>0 & final_age >0 & program == 'FIM', select=c( bay, tl, final_age, date, sex))) %>% filter(!is.na(final_age)) %>% 
  mutate(tl=tl/10, lcat2 =lencat(tl, w=1)) %>%rbind(age0)
Agelength_AP$bay="AP"
Agelength_AP <- subset(Agelength_AP, sex %in% c("M", "F"))

Agelength_CK<- droplevels(subset(as.data.frame(read.csv("ALK_Bay_and_weight.csv", header=T)), bay=="CK" & tl>0 & final_age >0 & program=='FIM', select=c( bay, tl, final_age, date, sex))) %>% filter(!is.na(final_age)) %>% 
  mutate(tl=tl/10, lcat2 =lencat(tl, w=1)) %>%rbind(age0)
Agelength_CK$bay="CK"
Agelength_CK <- subset(Agelength_CK, sex %in% c("M", "F"))

Agelength_CH<- droplevels(subset(as.data.frame(read.csv("ALK_Bay_and_weight.csv", header=T)), bay=="CH" & tl>0 & final_age >0 & program=='FIM', select=c( bay, tl, final_age, date, sex))) %>% filter(!is.na(final_age)) %>% 
  mutate(tl=tl/10, lcat2 =lencat(tl, w=1)) %>%rbind(age0) 
Agelength_CH$bay="CH"
Agelength_CH <- subset(Agelength_CH, sex %in% c("M", "F"))

Agelength_IR<- droplevels(subset(as.data.frame(read.csv("ALK_Bay_and_weight.csv", header=T)), bay=="IR" & tl>0 & final_age >0 & program=='FIM', select=c( bay, tl, final_age, date, sex))) %>% filter(!is.na(final_age)) %>% 
  mutate(tl=tl/10, lcat2 =lencat(tl, w=1)) %>%rbind(age0) 
Agelength_IR$bay="IR"
Agelength_IR <- subset(Agelength_IR, sex %in% c("M", "F"))

Agelength_JX<- droplevels(subset(as.data.frame(read.csv("ALK_Bay_and_weight.csv", header=T)), bay=="JX" & tl>0 & final_age >0 & program=='FIM', select=c( bay, tl, final_age, date, sex))) %>% filter(!is.na(final_age)) %>% 
  mutate(tl=tl/10, lcat2 =lencat(tl, w=1)) %>%rbind(age0) 
Agelength_JX$bay="JX"
Agelength_JX <- subset(Agelength_JX, sex %in% c("M", "F"))



all_raw <- rbind(Agelength_AP, Agelength_CH, Agelength_CK, Agelength_IR, Agelength_JX, Agelength_TB)

all_raw$bay <- as.factor(all_raw$bay)
all_raw$bay <- ordered(all_raw$bay, levels=c("AP", "CK", "TB", "CH", "JX", "IR"))
levels(all_raw$bay) <-  c("W1_AP", "W2_CK", "W3_TB", "W4_CH", "E1_JX", "E2_IR")




# DETERMINE GROWTH PARAMETERS FOR EACH ESTUARY ####


# Define function - more flexible when working with a single group of individuals (one estuary)

vbTyp <- function(age, Linf, K, t0) Linf*(1-exp(-K*(age-t0)))
vbTyp <- function(age, Linf, K, t0) Linf*(1-exp(-K*(age-t0)))
vbTyp <- vbFuns()


#GmTyp <- function(age, Linf, K, t0) Linf*(exp(-K*(age- t0)))



#test = nls(tl ~ GmTyp(final_age, Linf, K, t0), data=AL_TB, start=list(Linf=70, K=1, t0=0.1))

#Fit each bay 
# First define normal starting values

# TB ####
t0=0
AL_TB <-  Agelength_TB %>% left_join(Agelength_TB %>% group_by(final_age) %>% summarize(n= n()), "final_age")
TB_M <- droplevels(Agelength_TB %>% filter(sex == "M")) %>% left_join(Agelength_TB %>% filter(sex == "M") %>% group_by(final_age) %>% summarize(n= n()), "final_age")
TB_F <- droplevels(Agelength_TB %>% filter(sex == "F")) %>% left_join(Agelength_TB %>% filter(sex == "F") %>% group_by(final_age) %>% summarize(n= n()), "final_age")

starting <- list(Linf= max(Agelength_TB$tl, na.rm=TRUE), K=0.1, t0=-3)

#use this as starting if you want to fix t0
#starting <- list(Linf=max(Agelength_TB$tl, na.rm=TRUE), K=0.3)

fitTB <- nls(tl~vbTyp(final_age, Linf, K, t0=t0), data=AL_TB, weights=1/n, start=starting)



# fitTB_a <- nls(tl~vbTyp(final_age, Linf, K, t0=-1.5), data=Agelength_TB, start=starting)
# fitTB_b <- nls(tl~vbTyp(final_age, Linf, K, t0=-1), data=Agelength_TB, start=starting)
# fitTB_c <- nls(tl~vbTyp(final_age, Linf, K, t0=0), data=Agelength_TB, start=starting)
# fitTB_d <- nls(tl~vbTyp(final_age, Linf, K, t0=-2.5), data=Agelength_TB, start=starting)
# fitTB_e <- nls(tl~vbTyp(final_age, Linf, K, t0=-3), data=Agelength_TB, start=starting)
# fitTB_f <- nls(tl~vbTyp(final_age, Linf, K, t0=-3.5), data=Agelength_TB, start=starting)
# fitTB_g <- nls(tl~vbTyp(final_age, Linf, K, t0=-4), data=Agelength_TB, start=starting)


#AIC(fitTB_a, fitTB_b, fitTB_c, fitTB_d, fitTB_e, fitTB_f, fitTB_g)

#fitTB <- nls(tl~vbTyp(final_age, Linf, K), data=Agelength_TB, start=starting)
coef(fitTB)
bootTB <- nlsBoot(fitTB)
confint(bootTB) #, plot=TRUE)

#Linf           K          t0 
#41.97572483  1.32918282  0.05832843 

# 95% LCI    95% UCI
# Linf 41.54687362 42.4521489
# K     1.25258759  1.4040895
# t0    0.03200544  0.0826317

#Visualize the model fit
# - plot the best-fit VBGF with confidence intervals on top of the observed data
# uses a for loop to cycle through all ages

x <- seq(0,9, length.out=30) # ages for prediction
TB_pred <- vbTyp(x, Linf=coef(fitTB)[1],K= coef(fitTB)[2],t0=t0) #predicted lengths
xlmts <- range(c(x, Agelength_TB$final_age)) #set x limits
ylmts <- range(c(TB_pred, Agelength_TB$tl)) #set y limits
plot(tl~final_age, data=Agelength_TB, xlab="Age", ylab="Total Length (cm)", xlim=xlmts, ylim=ylmts, pch=19, col=rgb(0,0,0,1/3))
lines(TB_pred~x, lwd=2)

LCI <- UCI <- numeric(length(x))

for (i in 1:length(x))
{
  tmp <- apply(bootTB$coefboot, MARGIN =1, FUN=vbTyp, t=x[i])
  LCI[i] <- quantile(tmp, 0.025)
  UCI[i] <- quantile(tmp, 0.975)
}

x <- seq(0,9, length.out=30) # ages for prediction
TB_pred <- vbTyp(x, Linf=coef(fitTB)[1],K= coef(fitTB)[2],t0=t0) #predicted lengths
xlmts <- range(c(x, Agelength_TB$final_age)) #set x limits
ylmts <- range(c(TB_pred, Agelength_TB$tl)) #set y limits
plot(tl~final_age, data=Agelength_TB, xlab="Age", ylab="Total Length (cm)", xlim=xlmts, ylim=ylmts, pch=19, col=rgb(0,0,0,1/3))
lines(TB_pred~x, lwd=2)
lines(UCI ~ x, lwd=2, lty="dashed")
lines(LCI ~ x, lwd=2, lty="dashed")


# add confidence bands to each model fit ####

# Now for other estuaries

#Male
#starting <- list(Linf=max(TB_M$tl, na.rm=TRUE), K=0.3, t0=-1)
fitTBM <- nls(tl~vbTyp(final_age, Linf, K, t0=t0), data=TB_M, start=starting, weights=1/n, control=list(maxiter=500))
coef(fitTBM)
bootTBM <- nlsBoot(fitTBM)
confint(bootTBM) #, plot=TRUE)
#Linf           K          t0 
#37.80658295  1.55941940  0.08108728 

# 95% LCI     95% UCI
# Linf 37.32357831 38.24289021
# K     1.45729302  1.71522970
# t0    0.04728832  0.09902628

# Linf           K          t0 
# 37.75719186  1.57668907  0.07459461 

x <- seq(0,9, length.out=30) # ages for prediction
TB_predM <- vbTyp(x, Linf=coef(fitTBM)[1],K= coef(fitTBM)[2],t0=t0) #predicted lengths

# 
#Feamle
#starting <- list(Linf=max(TB_F$tl, na.rm=TRUE), K=0.3, t0=-1)
fitTBF <- nls(tl~vbTyp(final_age, Linf, K, t0=t0), data=TB_F, start=starting, weights=1/n, control=list(maxiter=500))
coef(fitTBF)
bootTBF <- nlsBoot(fitTBF)
confint(bootTBF) #, plot=TRUE)

# 95% LCI     95% UCI
# Linf 45.37935465 46.89330169
# K     0.94476629  1.10744925
# t0   -0.06095907  0.03363978

# Linf           K          t0 
# 46.10298743  1.02473610 -0.01011952 
x <- seq(0,9, length.out=30) # ages for prediction
TB_predF <- vbTyp(x, Linf=coef(fitTBF)[1],K= coef(fitTBF)[2],t0=t0) #predicted lengths


#AP ####
AL_AP <-  Agelength_AP %>% filter(!(final_age == 10)) %>% left_join(Agelength_AP %>% group_by(final_age) %>% summarize(n= n()), "final_age")
AP_M <- droplevels(Agelength_AP %>% filter(!(final_age ==10)) %>% filter(sex == "M")) %>% left_join(Agelength_AP %>% filter(sex == "M") %>% group_by(final_age) %>% summarize(n= n()), "final_age")
AP_F <- droplevels(Agelength_AP %>% filter(!(final_age ==10))  %>% filter(sex == "F")) %>% left_join(Agelength_AP %>% filter(sex == "F") %>% group_by(final_age) %>% summarize(n= n()), "final_age")


starting <- list(Linf=max(AL_AP$tl, na.rm=TRUE), K=0.3, t0=-1)

fitAP <- nls(tl~vbTyp(final_age, Linf, K, t0=t0), data=AL_AP, weights=1/n, start=starting)
coef(fitAP)
bootAP <- nlsBoot(fitAP)
confint(bootAP) #, plot=TRUE)


x <- seq(0,9, length.out=30) # ages for prediction
AP_pred <- vbTyp(x, Linf=coef(fitAP)) #predicted lengths
xlmts <- range(c(x, Agelength_AP$final_age)) #set x limits
ylmts <- range(c(AP_pred, Agelength_AP$tl)) #set y limits
plot(tl~final_age, data=Agelength_AP, xlab="Age", ylab="Total Length (cm)", xlim=xlmts, ylim=ylmts, pch=19, col=rgb(0,0,0,1/3))
lines(AP_pred~x, lwd=2)

#Female 
#starting <- list(Linf=max(AP_F$tl, na.rm=TRUE), K=0.3, t0=-1)
fitAPF <- nls(tl~vbTyp(final_age, Linf, K, t0=t0), data=AP_F, start=starting, weights=1/n, control=list(maxiter=500))
coef(fitAPF)
bootAPF <- nlsBoot(fitAPF)
confint(bootAPF) #, plot=TRUE)

x <- seq(0,9, length.out=30) # ages for prediction
AP_predF <- vbTyp(x, Linf=coef(fitAPF)) #predicted lengths


#Male
#starting <- list(Linf=max(AP_M$tl, na.rm=TRUE), K=0.3, t0=-1)
fitAPM <- nls(tl~vbTyp(final_age, Linf, K, t0=t0), data=AP_M, start=starting, weights=1/n, control=list(maxiter=500))
coef(fitAPM)
bootAPM <- nlsBoot(fitAPM)
confint(bootAPM) #, plot=TRUE)
# 
x <- seq(0,9, length.out=30) # ages for prediction
AP_predM <- vbTyp(x, Linf=coef(fitAPM)) #predicted lengths

#CK ####
starting <- list(Linf=max(AL_CK$tl, na.rm=TRUE), K=0.3, t0=-1)
AL_CK <-  Agelength_CK %>% left_join(Agelength_CK %>% group_by(final_age) %>% summarize(n= n()), "final_age")
CK_M <- droplevels(Agelength_CK %>% filter(sex == "M")) %>% left_join(Agelength_CK %>% filter(sex == "M") %>% group_by(final_age) %>% summarize(n= n()), "final_age")
CK_F <- droplevels(Agelength_CK %>% filter(sex == "F")) %>% left_join(Agelength_CK %>% filter(sex == "F") %>% group_by(final_age) %>% summarize(n= n()), "final_age")



#starting <- list(Linf=max(Agelength_CK$tl, na.rm=TRUE), K=0.3)

fitCK <- nls(tl~vbTyp(final_age, Linf, K, t0=t0), data=AL_CK, weights=1/n, start=starting)
coef(fitCK)
bootCK <- nlsBoot(fitCK)
confint(bootCK) #, plot=TRUE)


x <- seq(0,9, length.out=30) # ages for prediction
CK_pred <- vbTyp(x, Linf=coef(fitCK)) #predicted lengths
xlmts <- range(c(x, Agelength_CK$final_age)) #set x limits
ylmts <- range(c(CK_pred, Agelength_CK$tl)) #set y limits
plot(tl~final_age, data=Agelength_CK, xlab="Age", ylab="Total Length (cm)", xlim=xlmts, ylim=ylmts, pch=19, col=rgb(0,0,0,1/3))
lines(CK_pred~x, lwd=2)

#Female
#starting <- list(Linf=max(CK_F$tl, na.rm=TRUE), K=0.3, t0=-1)
fitCKF <- nls(tl~vbTyp(final_age, Linf, K, t0=t0), data=CK_F, start=starting, weights=1/n, control=list(maxiter=500))
coef(fitCKF)
bootCKF <- nlsBoot(fitCKF)
confint(bootCKF) #, plot=TRUE)

x <- seq(0,9, length.out=30) # ages for prediction
CK_predF <- vbTyp(x, Linf=coef(fitCKF)) #predicted lengths

#Male 
fitCKM <- nls(tl~vbTyp(final_age, Linf, K, t0=t0), data=CK_M, start=starting, weights=1/n,control=list(maxiter=500))
coef(fitCKM)
bootCKM <- nlsBoot(fitCKM)
confint(bootCKM) #, plot=TRUE)

x <- seq(0,9, length.out=30) # ages for prediction
CK_predM <- vbTyp(x, Linf=coef(fitCKM)) #predicted lengths


#CH ####

AL_CH <-  Agelength_CH %>% left_join(Agelength_CH %>% group_by(final_age) %>% summarize(n= n()), "final_age")
CH_M <- droplevels(Agelength_CH %>% filter(sex == "M")) %>% left_join(Agelength_CH %>% filter(sex == "M") %>% group_by(final_age) %>% summarize(n= n()), "final_age")
CH_F <- droplevels(Agelength_CH %>% filter(sex == "F")) %>% left_join(Agelength_CH %>% filter(sex == "F") %>% group_by(final_age) %>% summarize(n= n()), "final_age")


starting <- list(Linf=max(AL_CH$tl, na.rm=TRUE), K=0.3, t0=0)

fitCH <- nls(tl~vbTyp(final_age, Linf, K, t0=t0), data=AL_CH, weights=1/n, start=starting)
coef(fitCH)
bootCH <- nlsBoot(fitCH)
confint(bootCH) #, plot=TRUE)



x <- seq(0,9, length.out=30) # ages for prediction
CH_pred <- vbTyp(x, Linf=coef(fitCH)) #predicted lengths
xlmts <- range(c(x, AL_CH$final_age)) #set x limits
ylmts <- range(c(CH_pred, Agelength_CH$tl)) #set y limits
plot(tl~final_age, data=AL_CH, xlab="Age", ylab="Total Length (cm)", xlim=xlmts, ylim=ylmts, pch=19, col=rgb(0,0,0,1/3))
lines(CH_pred~x, lwd=2)

#Female
#starting <- list(Linf=max(CH_F$tl, na.rm=TRUE), K=0.3, t0=-1)

fitCHF <- nls(tl~vbTyp(final_age, Linf, K, t0=t0), data=CH_F, start=starting,weights=1/n, control=list(maxiter=500))
coef(fitCHF)
bootCHF <- nlsBoot(fitCHF)
confint(bootCHF) #, plot=TRUE)

x <- seq(0,9, length.out=30) # ages for prediction
CH_predF <- vbTyp(x, Linf=coef(fitCHF)) #predicted lengths

#Male 
starting <- list(Linf=max(CH_M$tl, na.rm=TRUE), K=0.3, t0=-1)

fitCHM <- nls(tl~vbTyp(final_age, Linf, K, t0=t0), data=CH_M, start=starting, weights=1/n,control=list(maxiter=500))
coef(fitCHM)
bootCHM <- nlsBoot(fitCHM)
confint(bootCHM) #, plot=TRUE)

x <- seq(0,9, length.out=30) # ages for prediction
CH_predM <- vbTyp(x, Linf=coef(fitCHM)) #predicted lengths


#IR####
AL_IR <-  Agelength_IR  %>% filter(!(final_age == 7 & tl < 40)) %>% left_join(Agelength_IR %>% group_by(final_age) %>% summarize(n= n()), "final_age")
IR_M <- droplevels(Agelength_IR %>% filter(!(final_age == 7 & tl < 40)) %>% filter(sex == "M") %>% filter(!(final_age>6 & tl < 30 ))) %>% left_join(Agelength_IR%>% filter(sex == "M") %>% group_by(final_age) %>% summarize(n= n()), "final_age")
IR_F <- droplevels(Agelength_IR %>% filter(!(final_age == 7 & tl < 40)) %>% filter(sex == "F")) %>% left_join(Agelength_IR %>% filter(sex == "F") %>% group_by(final_age) %>% summarize(n= n()), "final_age")


starting <- list(Linf=max(AL_IR$tl, na.rm=TRUE), K=0.3, t0=-1)

fitIR <- nls(tl~vbTyp(final_age, Linf, K, t0=t0), data=AL_IR, weights=1/n, start=starting)
coef(fitIR)
bootIR <- nlsBoot(fitIR)
confint(bootIR) #, plot=TRUE)


x <- seq(0,9, length.out=30) # ages for prediction
IR_pred <- vbTyp(x, Linf=coef(fitIR)) #predicted lengths
xlmts <- range(c(x, Agelength_IR$final_age)) #set x limits
ylmts <- range(c(IR_pred, Agelength_IR$tl)) #set y limits
plot(tl~final_age, data=AL_IR, xlab="Age", ylab="Total Length (cm)", xlim=xlmts, ylim=ylmts, pch=19, col=rgb(0,0,0,1/3))
lines(IR_pred~x, lwd=2)

#Female
#starting <- list(Linf=max(IR_F$tl, na.rm=TRUE), K=0.3, t0=-1)

fitIRF <- nls(tl~vbTyp(final_age, Linf, K, t0=t0), data=IR_F, start=starting, weights=1/n, control=list(maxiter=500))
coef(fitIRF)
bootIRF <- nlsBoot(fitIRF)
confint(bootIRF) #, plot=TRUE)

x <- seq(0,9, length.out=30) # ages for prediction
IR_predF <- vbTyp(x, Linf=coef(fitIRF)) #predicted lengths

#Male 
#starting <- list(Linf=max(IR_M$tl, na.rm=TRUE), K=0.3, t0=-1)

fitIRM <- nls(tl~vbTyp(final_age, Linf, K, t0=t0), data=IR_M, start=starting,weights=1/n, control=list(maxiter=500))
coef(fitIRM)
bootIRM <- nlsBoot(fitIRM)
confint(bootIRM) #, plot=TRUE)


x <- seq(0,9, length.out=30) # ages for prediction
IR_predM <- vbTyp(x, Linf=coef(fitIRM)) #predicted lengths

#JX####
AL_JX <-  Agelength_JX %>% left_join(Agelength_JX %>% group_by(final_age) %>% summarize(n= n()), "final_age")
JX_M <- droplevels(Agelength_JX %>% filter(sex == "M")) %>% left_join(Agelength_JX%>% filter(sex == "M") %>% group_by(final_age) %>% summarize(n= n()), "final_age")
JX_F <- droplevels(Agelength_JX %>% filter(sex == "F")) %>% left_join(Agelength_JX %>% filter(sex == "F") %>% group_by(final_age) %>% summarize(n= n()), "final_age")


starting <- list(Linf=max(AL_JX$tl, na.rm=TRUE), K=0.3, t0=-1)

fitJX <- nls(tl~vbTyp(final_age, Linf, K, t0), data=AL_JX, weights =1/n, start=starting)
coef(fitJX)
bootJX <- nlsBoot(fitJX)
confint(bootJX) #, plot=TRUE)

x <- seq(0,9, length.out=30) # ages for prediction
JX_pred <- vbTyp(x, Linf=coef(fitJX)) #predicted lengths
xlmts <- range(c(x, Agelength_JX$final_age)) #set x limits
ylmts <- range(c(JX_pred, Agelength_JX$tl)) #set y limits
plot(tl~final_age, data=Agelength_JX, xlab="Age", ylab="Total Length (cm)", xlim=xlmts, ylim=ylmts, pch=19, col=rgb(0,0,0,1/3))
lines(JX_pred~x, lwd=2)

#Female
#starting <- list(Linf=max(JX_F$tl, na.rm=TRUE), K=0.3, t0=-1)

fitJXF <- nls(tl~vbTyp(final_age, Linf, K, t0), data=JX_F, start=starting, weights=1/n, control=list(maxiter=500))
coef(fitJXF)
bootJXF <- nlsBoot(fitJXF)
confint(bootJXF) #, plot=TRUE)

x <- seq(0,9, length.out=30) # ages for prediction
JX_predF <- vbTyp(x, Linf=coef(fitJXF)) #predicted lengths

#Male 
fitJXM <- nls(tl~vbTyp(final_age, Linf, K, t0), data=JX_M, start=starting,weights=1/n, control=list(maxiter=500))
coef(fitJXM)
bootJXM <- nlsBoot(fitJXM)
confint(bootJXM, plot=TRUE)

x <- seq(0,9, length.out=30) # ages for prediction
JX_predM <- vbTyp(x, Linf=coef(fitJXM)) #predicted lengths

#PLOT ALL PREDICTED FITS IN ONE PLOT ####

# combine all predicted matrices into one
x <- seq(0,9, length.out=30)
t= data.frame(cbind(x,TB_pred)) %>% rename(pred=TB_pred) %>% mutate(bay=rep("TB",30))
a= data.frame(cbind(x,AP_pred)) %>% rename(pred=AP_pred) %>% mutate(bay=rep("AP",30))
ck= data.frame(cbind(x,CK_pred)) %>% rename(pred=CK_pred) %>% mutate(bay=rep("CK",30))
c= data.frame(cbind(x,CH_pred)) %>% rename(pred=CH_pred) %>% mutate(bay=rep("CH",30))
j= data.frame(cbind(x,JX_pred))%>% rename(pred=JX_pred) %>% mutate(bay=rep("JX",30))
i= data.frame(cbind(x,IR_pred)) %>% rename(pred=IR_pred) %>% mutate(bay=rep("IR",30))

pred_all <- rbind(t, a, ck, c, j,i)
pred_all$sex <- "All"
pred_all$bay <- as.factor(pred_all$bay)
pred_all$bay <- ordered(pred_all$bay, levels=c("AP", "CK", "TB", "CH", "JX", "IR"))
levels(pred_all$bay) <-  c("W1_AP", "W2_CK", "W3_TB", "W4_CH", "E1_JX", "E2_IR")


tF= data.frame(cbind(x,TB_predF)) %>% rename(pred=TB_predF) %>% mutate(bay=rep("TB",30))
aF= data.frame(cbind(x,AP_predF)) %>% rename(pred=AP_predF) %>% mutate(bay=rep("AP",30))
ckF= data.frame(cbind(x,CK_predF)) %>% rename(pred=CK_predF) %>% mutate(bay=rep("CK",30))
cF= data.frame(cbind(x,CH_predF)) %>% rename(pred=CH_predF) %>% mutate(bay=rep("CH",30))
jF= data.frame(cbind(x,JX_predF))%>% rename(pred=JX_predF) %>% mutate(bay=rep("JX",30))
iF= data.frame(cbind(x,IR_predF)) %>% rename(pred=IR_predF) %>% mutate(bay=rep("IR",30))

pred_allF <- rbind(tF, aF, ckF, cF, jF,iF)
pred_allF$sex <- "F"
pred_allF$bay <- as.factor(pred_allF$bay)
pred_allF$bay <- ordered(pred_allF$bay, levels=c("AP", "CK", "TB", "CH", "JX", "IR"))
levels(pred_allF$bay) <-  c("W1_AP", "W2_CK", "W3_TB", "W4_CH", "E1_JX", "E2_IR")


tM= data.frame(cbind(x,TB_predM)) %>% rename(pred=TB_predM) %>% mutate(bay=rep("TB",30))
aM= data.frame(cbind(x,AP_predM)) %>% rename(pred=AP_predM) %>% mutate(bay=rep("AP",30))
ckM= data.frame(cbind(x,CK_predM)) %>% rename(pred=CK_predM) %>% mutate(bay=rep("CK",30))
cM= data.frame(cbind(x,CH_predM)) %>% rename(pred=CH_predM) %>% mutate(bay=rep("CH",30))
jM= data.frame(cbind(x,JX_predM))%>% rename(pred=JX_predM) %>% mutate(bay=rep("JX",30))
iM= data.frame(cbind(x,IR_predM)) %>% rename(pred=IR_predM) %>% mutate(bay=rep("IR",30))

pred_allM <- rbind(tM, aM, ckM, cM, jM,iM)
pred_allM$sex <- "M"
pred_allM$bay <- as.factor(pred_allM$bay)
pred_allM$bay <- ordered(pred_allM$bay, levels=c("AP", "CK", "TB", "CH", "JX", "IR"))
levels(pred_allM$bay) <-  c("W1_AP", "W2_CK", "W3_TB", "W4_CH", "E1_JX", "E2_IR")


PRED <- rbind(pred_all, pred_allF, pred_allM)
PRED$sex <- as.factor(PRED$sex)


library(grid)
library(ggplot2)

File <- ("U:/PhD_projectfiles/Figures/growth_curves_nosex_WEIGHTED.tiff")
if (file.exists(File)) stop(File, " already exists")
dir.create(dirname(File), showWarnings = FALSE)

tiff(File, units="in", width=8, height=5, res=300)

ggplot(data= pred_all, aes(x=x, y=pred, group=bay))+
  geom_line() + #aes(linetype=bay), size =.5)+ # make line types based on the different labels- this will be our workaround because in a few stps we will specify the first label (obserseved) be a blank line (therefore a scatter plot)
  geom_line(aes(color=bay), size=1.1)+
  #geom_point(aes(shape=bay), size=2) + #, color=bay))+ # groups the points together correctly and then gives them a unique shape them differently based on the line type
  #scale_shape_manual(values=c(0,1,2, 3, 4, 5,6,7,8 ))+
  #scale_linetype_manual(values=c('solid', 'dashed', 'dotted'))+
  scale_y_continuous(limits=c(-10,60), breaks= seq(-10,60,10))+
  scale_x_continuous(limits=c(0,10), breaks=seq(0,10,1))+
  xlab("Age (yrs)")+
  ylab("Total Length (cm)")+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
        panel.background=element_rect(fill='white', colour='black'),
        legend.key=element_blank(), legend.title=element_blank(),
        legend.background=element_rect(fill='white', size=.5),
        legend.position="none",
        legend.key.size= unit(3, "cm"),
        axis.title.y = element_text(colour="black", size=20), # changing font of y axis title
        axis.title.x = element_text(colour="black", size=20),
        axis.text.x=element_text(colour="black", size=16), #changing  colour and font of x axis text
        axis.text.y=element_text(colour="black", size=16), #changing colour and font of y axis
        #plot.title=element_text(size=14), # changing size of plot title)+
        legend.text=element_text( size=18),
        legend.key.height=unit(25,"point")) # changes vertical spacing of the legend text
#ggtitle("Comparative von Bertalanffy Growth Models")
dev.off()


# File <- ("U:/PhD_projectfiles/Figures/growth_curves_coloredversion_WEIGHTED.tiff")
# if (file.exists(File)) stop(File, " already exists")
# dir.create(dirname(File), showWarnings = FALSE)
# 
# tiff(File, units="in", width=8, height=5, res=300)
# 
# ggplot(data= pred_all, aes(x=x, y=pred, group=bay))+  
#   geom_line() + #aes(linetype=bay), size =.5)+ # make line types based on the different labels- this will be our workaround because in a few stps we will specify the first label (obserseved) be a blank line (therefore a scatter plot)
#   geom_line(aes(color=bay), size=1.5)+
#   #geom_point(aes(shape=bay), size=2) + #, color=bay))+ # groups the points together correctly and then gives them a unique shape them differently based on the line type 
#   #scale_shape_manual(values=c(0,1,2, 3, 4, 5,6,7,8 ))+
#   #scale_linetype_manual(values=c('solid', 'dashed', 'dotted'))+
#   #scale_y_continuous(limits=c(-10,55), breaks= seq(-10,55,5))+
#   scale_x_continuous(limits=c(0,10), breaks=seq(0,10,1))+
#   xlab("Age (yrs)")+
#   ylab("Total Length (cm)")+
#   theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), 									
#         panel.background=element_rect(fill='white', colour='black'),
#         legend.key=element_blank(), legend.title=element_blank(),
#         legend.background=element_rect(fill='white', size=.5),
#         legend.position=c(.70,.38),
#         legend.key.size= unit(3, "cm"),
#         axis.title.y = element_text(colour="black", size=20), # changing font of y axis title
#         axis.title.x = element_text(colour="black", size=20),
#         axis.text.x=element_text(colour="black", size=16), #changing  colour and font of x axis text
#         axis.text.y=element_text(colour="black", size=16), #changing colour and font of y axis
#         #plot.title=element_text(size=14), # changing size of plot title)+
#         legend.text=element_text( size=18),
#         legend.key.height=unit(25,"point")) # changes vertical spacing of the legend text
# #ggtitle("Comparative von Bertalanffy Growth Models")
# dev.off()

#CROSS BAR PLOTS FOR ERROR ####

coefs <- data.frame(rbind(coef(fitTB),coef(fitTBF),coef(fitTBM), 
                          coef(fitAP),coef(fitAPF),coef(fitAPM),
                          coef(fitCK),coef(fitCKF), coef(fitCKM),
                          coef(fitCH), coef(fitCHF), coef(fitCHM),
                          coef(fitIR), coef(fitIRF), coef(fitIRM),
                          coef(fitJX), coef(fitJXF), coef(fitJXM)))


coefs$bay <- rep(c("TB", "AP", "CK", "CH", "IR", "JX"), each=3) 
coefs$sex <- rep(c("All", "F", "M"), each = 1)

test <- coefs %>% gather(c(Linf, K, t0), key="params", value="value")


confs <- data.frame(rbind(confint(bootTB),confint(bootTBF),confint(bootTBM),
                          confint(bootAP),confint(bootAPF),confint(bootAPM),
                          confint(bootCK),confint(bootCKF),confint(bootCKM),
                          confint(bootCH),confint(bootCHF),confint(bootCHM),
                          confint(bootIR),confint(bootIRF),confint(bootIRM),
                          confint(bootJX),confint(bootJXF),confint(bootJXM)))

confs$params <- rep(c("Linf", "K", "t0"), each=1)
confs$bay <- rep(c("TB", "AP", "CK", "CH", "IR", "JX"), each=9)
confs$sex <- rep(c("All", "F", "M"), each = 3)

uncertainty_matrix <- left_join(confs, test, by=(c("bay", "params", "sex")))
names(uncertainty_matrix)[1] <- "LCI"
names(uncertainty_matrix)[2] <- "UCI"

uncertainty_matrix$bay <- as.factor(uncertainty_matrix$bay)
uncertainty_matrix$bay <- ordered(uncertainty_matrix$bay, levels=c("AP", "CK", "TB", "CH", "JX", "IR"))
levels(uncertainty_matrix$bay) <-  c("W1_AP", "W2_CK", "W3_TB", "W4_CH", "E1_JX", "E2_IR")




ggplot(uncertainty_matrix %>% filter(sex == "All", params != "t0"), aes(bay, value, color=bay))+ geom_point(position=position_dodge(0.5)) +
  geom_errorbar(aes(ymin= LCI, ymax=UCI, color=bay), position=position_dodge(0.5))+
  facet_wrap(~params, scales = "free") +  theme(panel.grid.minor=element_blank(), 
                                                panel.grid.major=element_blank(), 
                                                panel.background=element_rect(colour="black", fill="white"),
                                                legend.position = "none",
                                               strip.background = element_rect(color= "black"))

#PLOT MEAN SIZE AT AGE & RAW DATA####


# NOTE:::: TO run this the age0 dataframe should not be added in the beginning
all <- rbind(Agelength_AP, Agelength_CH, Agelength_CK, Agelength_IR, Agelength_JX, Agelength_TB)
all <- all%>% ungroup()


sum <- all %>% dplyr::group_by(bay, final_age) %>% dplyr::summarize(N=length(tl), meantl=mean(tl), sdtl=sd(tl), se= sdtl/sqrt(N))
names(sum)[1] <- "Area"

File <- ("U:/PhD_projectfiles/Figures/raw_data_with_fit_PLUS_confints_WEIGHTED.tiff")
if (file.exists(File)) stop(File, " already exists")
dir.create(dirname(File), showWarnings = FALSE)

tiff(File, units="in", width=8, height=6, res=300)


# grid.arrange(ggplot(sum, aes(final_age, meantl, color=Area)) + 
#   geom_errorbar(aes(ymin=meantl-se, ymax=meantl+se), width=0.5, size=1) +
#   geom_point(aes(shape=Area), size=2.5)+
#   geom_line(aes(group=Area), size=0.75) +theme_bw()+
#   #xlab("Age (yrs)")+
#   ylab("Mean total length (cm)")+
#   scale_y_continuous(limits=c(0,80), breaks=seq(0,80, 10))+
#   scale_x_continuous(limits=c(0,10), breaks=seq(0,10, 1))+
#   theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
#         axis.title.y = element_text( colour="black", size=16), # changing font of y axis title
#         axis.title.x = element_blank(),
#         axis.text.x=element_text(colour="black", size=12), #changing  colour and font of x axis text
#         axis.text.y=element_text(colour="black", size=12),
#         legend.text=element_text( size=14),
#         legend.title=element_blank()),

grid.arrange(
  ggplot() + geom_point(data=all_raw, aes(final_age, tl, shape=sex, color=sex),position = position_dodge(0.5), alpha=0.3, size=2 )+
    geom_line(data=PRED, aes(x, pred, linetype=sex), size=1)+
    #geom_line(data=pred_allF, aes(x, pred), size=1, linetype="dashed")+
    #geom_line(data=pred_allM, aes(x, pred), size=1, linetype="dotted")+
    ylab("Total length (cm)") +
    xlab("Age (yrs)")+ 
    scale_color_manual(values=c("#F8766D", "#619CFF"))+
    scale_shape_manual(values=c(17,15))+
    scale_y_continuous(limits=c(-10,80), breaks=seq(-10,80, 10))+
    scale_x_continuous(limits=c(0, 10), breaks=seq(0, 10, 1))+
    theme(panel.grid.minor=element_blank(), 
          panel.grid.major=element_blank(), 
          panel.background=element_rect(colour="black", fill="white"),
          axis.title.x =element_text(colour="black", size=14),
          axis.text.x = element_text(colour="black"),
          strip.background = element_rect(color = "black", size = 1),
          strip.text.x=element_text(size=12),
          axis.title.y =element_text(colour="black", size=16),
          axis.text.y = element_text(colour="black"),
          legend.text=element_text( size=14),
          plot.title=element_text(size=14),
          legend.title=element_blank())+
    facet_wrap(~bay)+
    guides(color = guide_legend(override.aes = list(alpha=1))),
  ggplot(uncertainty_matrix, aes(bay, value, shape=sex, color=sex))+ 
    geom_point(position=position_dodge(0.5), size=2) +
    geom_errorbar(aes(ymin= LCI, ymax=UCI),size=0.6,position=position_dodge(0.5))+
    facet_wrap(~params, scales = "free")+
    xlab("Estuary")+
    ylab("")+
    scale_color_manual(values=c("#00BA38", "#F8766D", "#619CFF"))+
    theme(panel.grid.minor=element_line(color="grey"), 
          panel.grid.major=element_blank(), 
          panel.background=element_rect(colour="black", fill="white"),
          axis.title.x =element_text( colour="black", size=14),
          axis.text.x = element_text(colour="black", angle =25, size=8),
          strip.background = element_rect(color = "black", size = 1),
          strip.text.x=element_text(size=12),
          axis.title.y =element_text(color="black"),
          axis.text.y = element_text(colour="black"),
          legend.text=element_text( size=14),
          plot.title=element_text(size=14),
          legend.title=element_blank()) +
    geom_vline(xintercept=seq(1.5, length(unique(uncertainty_matrix$bay))-0.5, 1), 
               lwd=1,linetype="solid", colour="grey"), nrow =2, heights=c(1.5,0.75))

dev.off()


# COMPARE VGBFs BETWEEN GROUPS ####
##############################
# requires fitting multiple models to determine which parameters between each group is different. 
# For simplicities sake, I will only compare 2 bays at a time. Can refer to ALK_analysis.R for a list of bay to bay comparisons. 
# Refer to page 236 in Ogle for the family of models that must be considered when examining the differences in VBGF among groups. 
# The nested relationships among these models allows use of likelihoo ratio and extra sum of squares tests.

########################################
# Define the generic models to compare
#######################################


vbLKt <- tl~Linf[sex]*(1-exp(-K[sex]*(final_age-t0[sex])))
vbLK <- tl~Linf[sex]*(1-exp(-K[sex]*(final_age-t0)))
vbLt <- tl~Linf[sex]*(1-exp(-K*(final_age-t0[sex])))
vbKt <- tl~Linf*(1-exp(-K[sex]*(final_age-t0[sex])))
vbL <-  tl~Linf[sex]*(1-exp(-K*(final_age-t0)))
vbK <-  tl~Linf*(1-exp(-K[sex]*(final_age-t0)))
vbt <- tl~Linf*(1-exp(-K*(final_age-t0[sex])))
vb0 <- tl~Linf*(1-exp(-K*(final_age-t0)))

###############################################
#Define starting values for each model scenario
###############################################
sv0 <- vbStarts(tl~final_age, data=AL_TB) #using the TB dataframe as a starter
sv0$Linf= 60
sv0$K = 0.3
sv0$t0 = -0.1
svLKt <- Map(rep,sv0,c(2,2,2))
svLK <- Map(rep,sv0,c(2,2,1))
svLt <- Map(rep,sv0,c(2,1,2))
svKt <- Map(rep,sv0,c(1,2,2))
svL <- Map(rep,sv0,c(2,1,1))
svt <- Map(rep,sv0,c(1,1,2))
svK <- Map(rep, sv0,c(1,2,1))

##########################################
#Define dataframes with bay-bay comparison
##########################################
AL_TBCK<- rbind(AL_TB, AL_CK) 
AL_TBCK$bay <- as.factor(AL_TBCK$bay)
AL_TBCH<- rbind(AL_TB, AL_CH)
AL_TBCH$bay <- as.factor(AL_TBCH$bay)
AL_TBAP<- rbind(AL_TB, AL_AP)
AL_TBAP$bay <- as.factor(AL_TBAP$bay)
AL_TBJX<- rbind(AL_TB, AL_JX)
AL_TBJX$bay <- as.factor(AL_TBJX$bay) 
AL_TBIR<- rbind(AL_TB, AL_IR)
AL_TBIR$bay <- as.factor(AL_TBIR$bay)

AL_CKCH<- rbind(AL_CK, AL_CH)
AL_CKCH$bay <- as.factor(AL_CKCH$bay)
AL_CKAP<- rbind(AL_CK, AL_AP) 
AL_CKAP$bay <- as.factor(AL_CKAP$bay)
AL_CKJX<- rbind(AL_CK, AL_JX)
AL_CKJX$bay <- as.factor(AL_CKJX$bay)
AL_CKIR<- rbind(AL_CK, AL_IR)
AL_CKIR$bay <- as.factor(AL_CKIR$bay)

AL_CHAP<- rbind(AL_CH, AL_AP)
AL_CHAP$bay <- as.factor(AL_CHAP$bay)
AL_CHJX<- rbind(AL_CH, AL_JX)
AL_CHJX$bay <- as.factor(AL_CHJX$bay)
AL_CHIR<- rbind(AL_CH, AL_IR)
AL_CHIR$bay <- as.factor(AL_CHIR$bay)

AL_APJX<- rbind(AL_AP, AL_JX)
AL_APJX$bay <- as.factor(AL_APJX$bay)

AL_APIR<- rbind(AL_AP, AL_IR)
AL_APIR$bay <- as.factor(AL_APIR$bay)

AL_JXIR<- rbind(AL_JX, AL_IR)
AL_JXIR$bay <- as.factor(AL_JXIR$bay)

###############################################################
#Define fits to model scenarios- first by sex for each bay 
############################################################## 

#filter missing sex
library(minpack.lm)

tb_mf <- droplevels(subset(AL_TB, sex %in% c("M", "F")))
#Fit to [LKt] model then to the simplest model [vb0]

fitLKt <- nls(vbLKt, data=tb_mf, weights=1/n, start=svLKt)
#residPlot(fitLKt, col=rgb(0,0,0,1/3))
fitLK <- nls(vbLK, data=tb_mf, weights=1/n, start=svLK)
fitLt <- nls(vbLt, data=tb_mf, weights=1/n, start=svLt)
fitKt <- nls(vbKt, data=tb_mf, weights=1/n, start=svKt)
fitL <- nls(vbL, data=tb_mf, weights=1/n, start=svL)
fitt <- nls(vbt, data=tb_mf, weights=1/n, start=svt)
fitK <- nls(vbK, data=tb_mf, weights=1/n, start=svK)
fit0 <- nls(vb0, data=tb_mf, weights=1/n, start=sv0)

cbind(AIC(fitLKt,fitLK,fitLt,fitKt, fitL, fitK, fitt, fit0), BIC(fitLKt, fitLK, fitLt, fitKt, fitL, fitK, fitt, fit0))

# Likelihood ratio and extra for comparing the two models are calculated with lrt() and extraSS().
lrt(fit0, com=fitLKt, com.name="All pars differ", sim.name="No pars differ")
extraSS(fit0, com=fitLKt, com.name="All pars diff", sim.names="No pars diff")

#Likelihood ratio and extra sums of square test both suggest a significant difference
# between the most complex and the simplest. So there is evidence that there is some
# difference in the parameters 

#Compare the nested models

lrt(fitLK, fitLt, fitKt, com=fitLKt, com.name="All pars diff", sim.names=c("Linf, K diff", "Linf,t0 diff", "K, t0, diff"))

# ***** Any nested model that is not statistically different from[L,K,t0] is considered "better" because it fits equally (statistically) well, but is more parsimonious.
#       If two models are better than [L, K, t0] then the one with the greatest LL is chosen to be the best of the nested models. 

ch_mf <- droplevels(subset(AL_CH, sex %in% c("M", "F")))
#Fit to [LKt] model then to the simplest model [vb0]
fitLKt <- nls(vbLKt, data=ch_mf,weights=1/n, start=svLKt)
#residPlot(fitLKt, col=rgb(0,0,0,1/3))
fitLK <- nls(vbLK, data=ch_mf,weights=1/n, start=svLK)
fitLt <- nls(vbLt, data=ch_mf,weights=1/n, start=svLt)
fitKt <- nls(vbKt, data=ch_mf, weights=1/n,start=svKt)
fitL <- nls(vbL, data=ch_mf,weights=1/n, start=svL)
fitt <- nls(vbt, data=ch_mf, weights=1/n,start=svt)
fitK <- nls(vbK, data=ch_mf, weights=1/n,start=svK)
fit0 <- nls(vb0, data=ch_mf, weights=1/n,start=sv0)

cbind(AIC(fitLKt,fitLK,fitLt,fitKt, fitL, fitK, fitt, fit0), BIC(fitLKt, fitLK, fitLt, fitKt, fitL, fitK, fitt, fit0))

ck_mf <- droplevels(subset(AL_CK, sex %in% c("M", "F")))
#Fit to [LKt] model then to the simplest model [vb0]
fitLKt <- nls(vbLKt, data=ck_mf, weights=1/n,start=svLKt)
#residPlot(fitLKt, col=rgb(0,0,0,1/3))
fitLK <- nls(vbLK, data=ck_mf, weights=1/n,start=svLK)
fitLt <- nls(vbLt, data=ck_mf,weights=1/n, start=svLt)
fitKt <- nls(vbKt, data=ck_mf, weights=1/n,start=svKt)
fitL <- nls(vbL, data=ck_mf, weights=1/n,start=svL)
fitt <- nls(vbt, data=ck_mf,weights=1/n,start=svt)
fitK <- nls(vbK, data=ck_mf, weights=1/n,start=svK)
fit0 <- nls(vb0, data=ck_mf, weights=1/n,start=sv0)

cbind(AIC(fitLKt,fitLK,fitLt,fitKt, fitL, fitK, fitt, fit0), BIC(fitLKt, fitLK, fitLt, fitKt, fitL, fitK, fitt, fit0))

ap_mf <- droplevels(subset(AL_AP, sex %in% c("M", "F")))
#Fit to [LKt] model then to the simplest model [vb0]
fitLKt <- nls(vbLKt, data=ap_mf, weights=1/n,start=svLKt)
#residPlot(fitLKt, col=rgb(0,0,0,1/3))
fitLK <- nls(vbLK, data=ap_mf,weights=1/n, start=svLK)
fitLt <- nls(vbLt, data=ap_mf, weights=1/n,start=svLt)
fitKt <- nls(vbKt, data=ap_mf,weights=1/n, start=svKt)
fitL <- nls(vbL, data=ap_mf, weights=1/n,start=svL)
fitt <- nls(vbt, data=ap_mf, weights=1/n,start=svt)
fitK <- nls(vbK, data=ap_mf, weights=1/n,start=svK)
fit0 <- nls(vb0, data=ap_mf, weights=1/n,start=sv0)

cbind(AIC(fitLKt,fitLK,fitLt,fitKt, fitL, fitK, fitt, fit0), BIC(fitLKt, fitLK, fitLt, fitKt, fitL, fitK, fitt, fit0))

jx_mf <- droplevels(subset(AL_JX, sex %in% c("M", "F")))
#Fit to [LKt] model then to the simplest model [vb0]
fitLKt <- nls(vbLKt, data=jx_mf,weights=1/n, start=svLKt)
#residPlot(fitLKt, col=rgb(0,0,0,1/3))
fitLK <- nls(vbLK, data=jx_mf, weights=1/n,start=svLK)
fitLt <- nls(vbLt, data=jx_mf, weights=1/n,start=svLt)
fitKt <- nls(vbKt, data=jx_mf, weights=1/n,start=svKt)
fitL <- nls(vbL, data=jx_mf, weights=1/n,start=svL)
fitt <- nls(vbt, data=jx_mf, weights=1/n,start=svt)
fitK <- nls(vbK, data=jx_mf, weights=1/n,start=svK)
fit0 <- nls(vb0, data=jx_mf,weights=1/n, start=sv0)

cbind(AIC(fitLKt,fitLK,fitLt,fitKt, fitL, fitK, fitt, fit0), BIC(fitLKt, fitLK, fitLt, fitKt, fitL, fitK, fitt, fit0))

ir_mf <- droplevels(subset(AL_IR, sex %in% c("M", "F")))
#Fit to [LKt] model then to the simplest model [vb0]
fitLKt <- nls(vbLKt, data=ir_mf, weights=1/n,start=svLKt)
#residPlot(fitLKt, col=rgb(0,0,0,1/3))
fitLK <- nls(vbLK, data=ir_mf, weights=1/n,start=svLK)
fitLt <- nls(vbLt, data=ir_mf, weights=1/n,start=svLt)
fitKt <- nls(vbKt, data=ir_mf, weights=1/n,start=svKt)
fitL <- nls(vbL, data=ir_mf, weights=1/n,start=svL)
fitt <- nls(vbt, data=ir_mf, weights=1/n,start=svt)
fitK <- nls(vbK, data=ir_mf, weights=1/n,start=svK)
fit0 <- nls(vb0, data=ir_mf, weights=1/n,start=sv0)

cbind(AIC(fitLKt,fitLK,fitLt,fitKt, fitL, fitK, fitt, fit0), BIC(fitLKt, fitLK, fitLt, fitKt, fitL, fitK, fitt, fit0))


###############################################################
#Define fits to model scenarios (one bay at a time) _TBCK here
##############################################################
#vbLKt <- tl~Linf[bay]*(1-exp(-K[bay]*(final_age-t0[bay])))
vbLK <- tl~Linf[bay]*(1-exp(-K[bay]*(final_age-t0)))
#vbLt <- tl~Linf[bay]*(1-exp(-K*(final_age-t0[bay])))
#vbKt <- tl~Linf*(1-exp(-K[bay]*(final_age-t0[bay])))
vbL <-  tl~Linf[bay]*(1-exp(-K*(final_age-t0)))
vbK <-  tl~Linf*(1-exp(-K[bay]*(final_age-t0)))
#vbt <- tl~Linf*(1-exp(-K*(final_age-t0[bay])))
vb0 <- tl~Linf*(1-exp(-K*(final_age-t0)))

#Fit to [LKt] model then to the simplest model [vb0]
fitLKt <- nls(vbLKt, data=AL_TBCK, weights=1/n, start=svLKt)
#residPlot(fitLKt, col=rgb(0,0,0,1/3))
fitLK <- nls(vbLK, data=AL_TBCK, weights=1/n,start=svLK)
fitLt <- nls(vbLt, data=AL_TBCK,weights=1/n, start=svLt)
fitKt <- nls(vbKt, data=AL_TBCK,weights=1/n, start=svKt)
fitL <- nls(vbL, data=AL_TBCK,weights=1/n, start=svL)
fitt <- nls(vbt, data=AL_TBCK,weights=1/n,start=svt)
fitK <- nls(vbK, data=AL_TBCK, weights=1/n,start=svK)
fit0 <- nls(vb0, data=AL_TBCK,weights=1/n, start=sv0)


############################################
#Likelihood ratio testing example with TBCK
###########################################


# Likelihood ratio and extra for comparing the two models are calculated with lrt() and extraSS().
lrt(fit0, com=fitLKt, com.name="All pars differ", sim.name="No pars differ")
extraSS(fit0, com=fitLKt, com.name="All pars diff", sim.names="No pars diff")

#Likelihood ratio and extra sums of square test both suggest a significant difference
# between the most complex and the simplest. So there is evidence that there is some
# difference in the parameters between TB and CK


#Compare the nested models

lrt(fitLK, fitLt, fitKt, com=fitLKt, com.name="All pars diff", sim.names=c("Linf, K diff", "Linf,t0 diff", "K, t0, diff"))

lrt(fit0, com=fitLK, com.name="All pars differ", sim.name="No pars differ")
lrt(fitL, fitK, com=fitLK, com.name="All pars diff", sim.names=c("Linf diff", "K diff"))


lrtest(fitL, fitLK)
lrtest(fitK, fitLK)
# ***** Any nested model that is not statistically different from[L,K,t0] is considered "better" because it fits equally (statistically) well, but is more parsimonious.
#       If two models are better than [L, K, t0] then the one with the greatest LL is chosen to be the best of the nested models. 


##############################################################
# Easier AIC/BIC example that could be used in place of the LRT
###############################################################
cbind(AIC(fitLKt,fitLK,fitLt,fitKt, fitL, fitK, fitt, fit0), BIC(fitLKt, fitLK, fitLt, fitKt, fitL, fitK, fitt, fit0))
#Lowest:
# #       df      AIC df      BIC
# fitLK   6 22568.78  6 22605.44


#################################################################
# Now do it all over again with each Bay to bay comparison
##################################################################

# TBCH
fitLKt <- nls(vbLKt, data=AL_TBCH,weights=1/n,start=svLKt)
fitLK <- nls(vbLK, data=AL_TBCH, weights=1/n,start=svLK)
fitLt <- nls(vbLt, data=AL_TBCH,weights=1/n, start=svLt)
fitKt <- nls(vbKt, data=AL_TBCH,weights=1/n,start=svKt)
fitL <- nls(vbL, data=AL_TBCH, weights=1/n,start=svL)
fitt <- nls(vbt, data=AL_TBCH,weights=1/n, start=svt)
fitK <- nls(vbK, data=AL_TBCH,weights=1/n, start=svK)
fit0 <- nls(vb0, data=AL_TBCH,weights=1/n, start=sv0)

#lrt(fit0, com=fitLKt, com.name="All pars differ", sim.name="No pars differ")
#lrt(fitLK, fitLt, fitKt, com=fitLKt, com.name="All pars diff", sim.names=c("Linf, K diff", "Linf,t0 diff", "K, t0, diff"))

lrt(fit0, com=fitLK, com.name="All pars differ", sim.name="No pars differ")
lrt(fitL, fitK, com=fitLK, com.name="All pars diff", sim.names=c("Linf diff", "K diff"))


cbind(AIC(fitLKt,fitLK,fitLt,fitKt, fitL, fitK, fitt, fit0), BIC(fitLKt, fitLK, fitLt, fitKt, fitL, fitK, fitt, fit0))
#Lowest:
# df      AIC df      BIC
# fitL    5 21809.18  5 21839.78
#fitK    5 21808.97  5 21839.57
#fitt    5 21808.85  5 21839.45
#fit0    4 21808.25  4 21832.73

#TBAP
fitLKt <- nls(vbLKt, data=AL_TBAP,  weights=1/n,start=svLKt)
fitLK <- nls(vbLK, data=AL_TBAP,weights=1/n, start=svLK)
fitLt <- nls(vbLt, data=AL_TBAP,weights=1/n, start=svLt)
fitKt <- nls(vbKt, data=AL_TBAP,weights=1/n, start=svKt)
fitL <- nls(vbL, data=AL_TBAP, weights=1/n,start=svL)
fitt <- nls(vbt, data=AL_TBAP, weights=1/n,start=svt)
fitK <- nls(vbK, data=AL_TBAP,weights=1/n, start=svK)
fit0 <- nls(vb0, data=AL_TBAP,weights=1/n, start=sv0)

#lrt(fit0, com=fitLKt, com.name="All pars differ", sim.name="No pars differ")
#lrt(fitLK, fitLt, fitKt, com=fitLKt, com.name="All pars diff", sim.names=c("Linf, K diff", "Linf,t0 diff", "K, t0, diff"))


lrt(fit0, com=fitLK, com.name="All pars differ", sim.name="No pars differ")
lrt(fitL, fitK, com=fitLK, com.name="All pars diff", sim.names=c("Linf diff", "K diff"))

lrt(fit0, com=fitK, com.name = "K diff")

cbind(AIC(fitLKt,fitLK,fitLt,fitKt, fitL, fitK, fitt, fit0), BIC(fitLKt, fitLK, fitLt, fitKt, fitL, fitK, fitt, fit0))
#Lowest:
# df      AIC df      BIC
# fitLKt  7 22566.44  7 22609.48
# fitLK   6 22566.46  6 22603.35


#TBJX
fitLKt <- nls(vbLKt, data=AL_TBJX,weights=1/n, start=svLKt)
fitLK <- nls(vbLK, data=AL_TBJX, weights=1/n,start=svLK)
fitLt <- nls(vbLt, data=AL_TBJX, weights=1/n,start=svLt)
fitKt <- nls(vbKt, data=AL_TBJX, weights=1/n,start=svKt)
fitL <- nls(vbL, data=AL_TBJX,weights=1/n, start=svL)
fitt <- nls(vbt, data=AL_TBJX, weights=1/n,start=svt)
fitK <- nls(vbK, data=AL_TBJX, weights=1/n,start=svK)
fit0 <- nls(vb0, data=AL_TBJX, weights=1/n,start=sv0)

#lrt(fit0, com=fitLKt, com.name="All pars differ", sim.name="No pars differ")
#lrt(fitLK, fitLt, fitKt, com=fitLKt, com.name="All pars diff", sim.names=c("Linf, K diff", "Linf,t0 diff", "K, t0, diff"))


lrt(fit0, com=fitLK, com.name="All pars differ", sim.name="No pars differ")
lrt(fitL, fitK, com=fitLK, com.name="All pars diff", sim.names=c("Linf diff", "K diff"))

lrt(fit0, com=fitK, com.name="Kdiff")


cbind(AIC(fitLKt,fitLK,fitLt,fitKt, fitL, fitK, fitt, fit0), BIC(fitLKt, fitLK, fitLt, fitKt, fitL, fitK, fitt, fit0))
#Lowest:
# df      AIC df      BIC
# fitLKt  7 19730.15  7 19772.32
# fitLK   6 19729.45  6 19765.60

#TBIR
fitLKt <- nls(vbLKt, data=AL_TBIR, weights=1/n, start=svLKt)
fitLK <- nls(vbLK, data=AL_TBIR,  weights=1/n,start=svLK)
fitLt <- nls(vbLt, data=AL_TBIR,  weights=1/n,start=svLt)
fitKt <- nls(vbKt, data=AL_TBIR,  weights=1/n,start=svKt)
fitL <- nls(vbL, data=AL_TBIR, weights=1/n,start=svL)
fitt <- nls(vbt, data=AL_TBIR, weights=1/n, start=svt)
fitK <- nls(vbK, data=AL_TBIR,  weights=1/n,start=svK)
fit0 <- nls(vb0, data=AL_TBIR,  weights=1/n,start=sv0)

#lrt(fit0, com=fitLKt, com.name="All pars differ", sim.name="No pars differ")
#lrt(fitLK, fitLt, fitKt, com=fitLKt, com.name="All pars diff", sim.names=c("Linf, K diff", "Linf,t0 diff", "K, t0, diff"))

lrt(fit0, com=fitLK, com.name="All pars differ", sim.name="No pars differ")
lrt(fitL, fitK, com=fitLK, com.name="All pars diff", sim.names=c("Linf diff", "K diff"))




cbind(AIC(fitLKt,fitLK,fitLt,fitKt, fitL, fitK, fitt, fit0), BIC(fitLKt, fitLK, fitLt, fitKt, fitL, fitK, fitt, fit0))
#Lowest:
# df      AIC df      BIC
# fitLKt  7 33832.90  7 33878.47

#CKCH
fitLKt <- nls(vbLKt, data=AL_CKCH, weights=1/n,start=svLKt)
fitLK <- nls(vbLK, data=AL_CKCH, weights=1/n,start=svLK)
fitLt <- nls(vbLt, data=AL_CKCH, weights=1/n,start=svLt)
fitKt <- nls(vbKt, data=AL_CKCH,weights=1/n, start=svKt)
fitL <- nls(vbL, data=AL_CKCH, weights=1/n,start=svL)
fitt <- nls(vbt, data=AL_CKCH, weights=1/n,start=svt)
fitK <- nls(vbK, data=AL_CKCH,weights=1/n, start=svK)
fit0 <- nls(vb0, data=AL_CKCH, weights=1/n,start=sv0)

#lrt(fit0, com=fitLKt, com.name="All pars differ", sim.name="No pars differ")
#lrt(fitLK, fitLt, fitKt, com=fitLKt, com.name="All pars diff", sim.names=c("Linf, K diff", "Linf,t0 diff", "K, t0, diff"))

lrt(fit0, com=fitLK, com.name="All pars differ", sim.name="No pars differ")
lrt(fitL, fitK, com=fitLK, com.name="All pars diff", sim.names=c("Linf diff", "K diff"))



cbind(AIC(fitLKt,fitLK,fitLt,fitKt, fitL, fitK, fitt, fit0), BIC(fitLKt, fitLK, fitLt, fitKt, fitL, fitK, fitt, fit0))
#Lowest:
# df      AIC df      BIC
# fitLK   6 15209.74  6 15244.38

#CKAP
fitLKt <- nls(vbLKt, data=AL_CKAP, weights=1/n,start=svLKt)
fitLK <- nls(vbLK, data=AL_CKAP, weights=1/n,start=svLK)
fitLt <- nls(vbLt, data=AL_CKAP,weights=1/n, start=svLt)
fitKt <- nls(vbKt, data=AL_CKAP, weights=1/n,start=svKt)
fitL <- nls(vbL, data=AL_CKAP, weights=1/n,start=svL)
fitt <- nls(vbt, data=AL_CKAP, weights=1/n,start=svt)
fitK <- nls(vbK, data=AL_CKAP, weights=1/n,start=svK)
fit0 <- nls(vb0, data=AL_CKAP, weights=1/n,start=sv0)

#lrt(fit0, com=fitLKt, com.name="All pars differ", sim.name="No pars differ")
#lrt(fitLK, fitLt, fitKt, com=fitLKt, com.name="All pars diff", sim.names=c("Linf, K diff", "Linf,t0 diff", "K, t0, diff"))

lrt(fit0, com=fitLK, com.name="All pars differ", sim.name="No pars differ")
lrt(fitL, fitK, com=fitLK, com.name="All pars diff", sim.names=c("Linf diff", "K diff"))


cbind(AIC(fitLKt,fitLK,fitLt,fitKt, fitL, fitK, fitt, fit0), BIC(fitLKt, fitLK, fitLt, fitKt, fitL, fitK, fitt, fit0))
#Lowest:
# df      AIC df      BIC
# fitLKt  7 15977.29  7 16017.97


#CKJX
fitLKt <- nls(vbLKt, data=AL_CKJX, weights=1/n, start=svLKt)
fitLK <- nls(vbLK, data=AL_CKJX, weights=1/n,start=svLK)
fitLt <- nls(vbLt, data=AL_CKJX, weights=1/n,start=svLt)
fitKt <- nls(vbKt, data=AL_CKJX, weights=1/n,start=svKt)
fitL <- nls(vbL, data=AL_CKJX, weights=1/n,start=svL)
fitt <- nls(vbt, data=AL_CKJX, weights=1/n,start=svt)
fitK <- nls(vbK, data=AL_CKJX, weights=1/n,start=svK)
fit0 <- nls(vb0, data=AL_CKJX,weights=1/n, start=sv0)

#lrt(fit0, com=fitLKt, com.name="All pars differ", sim.name="No pars differ")
#lrt(fitLK, fitLt, fitKt, com=fitLKt, com.name="All pars diff", sim.names=c("Linf, K diff", "Linf,t0 diff", "K, t0, diff"))

lrt(fit0, com=fitLK, com.name="All pars differ", sim.name="No pars differ")
lrt(fitL, fitK, com=fitLK, com.name="All pars diff", sim.names=c("Linf diff", "K diff"))


cbind(AIC(fitLKt,fitLK,fitLt,fitKt, fitL, fitK, fitt, fit0), BIC(fitLKt, fitLK, fitLt, fitKt, fitL, fitK, fitt, fit0))
#lowest:
# df      AIC df      BIC
# fitLK   6 13116.66  6 13150.46


#CKIR
fitLKt <- nls(vbLKt, data=AL_CKIR, weights=1/n,start=svLKt)
fitLK <- nls(vbLK, data=AL_CKIR,weights=1/n, start=svLK)
fitLt <- nls(vbLt, data=AL_CKIR,weights=1/n, start=svLt)
fitKt <- nls(vbKt, data=AL_CKIR,weights=1/n, start=svKt)
fitL <- nls(vbL, data=AL_CKIR,weights=1/n, start=svL)
fitt <- nls(vbt, data=AL_CKIR, weights=1/n,start=svt)
fitK <- nls(vbK, data=AL_CKIR, weights=1/n,start=svK)
fit0 <- nls(vb0, data=AL_CKIR, weights=1/n,start=sv0)

#lrt(fit0, com=fitLKt, com.name="All pars differ", sim.name="No pars differ")
#lrt(fitLK, fitLt, fitKt, com=fitLKt, com.name="All pars diff", sim.names=c("Linf, K diff", "Linf,t0 diff", "K, t0, diff"))

lrt(fit0, com=fitLK, com.name="All pars differ", sim.name="No pars differ")
lrt(fitL, fitK, com=fitLK, com.name="All pars diff", sim.names=c("Linf diff", "K diff"))


cbind(AIC(fitLKt,fitLK,fitLt,fitKt, fitL, fitK, fitt, fit0), BIC(fitLKt, fitLK, fitLt, fitKt, fitL, fitK, fitt, fit0))
#Lowest:
# df      AIC df      BIC
# fitLKt  7 27244.87  7 27288.89

#CHAP
fitLKt <- nls(vbLKt, data=AL_CHAP, weights=1/n,start=svLKt)
fitLK <- nls(vbLK, data=AL_CHAP, weights=1/n,start=svLK)
fitLt <- nls(vbLt, data=AL_CHAP, weights=1/n,start=svLt)
fitKt <- nls(vbKt, data=AL_CHAP, weights=1/n,start=svKt)
fitL <- nls(vbL, data=AL_CHAP, weights=1/n,start=svL)
fitt <- nls(vbt, data=AL_CHAP, weights=1/n,start=svt)
fitK <- nls(vbK, data=AL_CHAP, weights=1/n,start=svK)
fit0 <- nls(vb0, data=AL_CHAP, weights=1/n,start=sv0)

#lrt(fit0, com=fitLKt, com.name="All pars differ", sim.name="No pars differ")
#lrt(fitLK, fitLt, fitKt, com=fitLKt, com.name="All pars diff", sim.names=c("Linf, K diff", "Linf,t0 diff", "K, t0, diff"))

lrt(fit0, com=fitLK, com.name="All pars differ", sim.name="No pars differ")
lrt(fitL, fitK, com=fitLK, com.name="All pars diff", sim.names=c("Linf diff", "K diff"))


cbind(AIC(fitLKt,fitLK,fitLt,fitKt, fitL, fitK, fitt, fit0), BIC(fitLKt, fitLK, fitLt, fitKt, fitL, fitK, fitt, fit0))
#Lowest:
# df      AIC df      BIC
#fitLKt  7 27244.87  7 27288.89


#CHJX
fitLKt <- nls(vbLKt, data=AL_CHJX, weights=1/n,start=svLKt)
fitLK <- nls(vbLK, data=AL_CHJX, weights=1/n,start=svLK)
fitLt <- nls(vbLt, data=AL_CHJX, weights=1/n,start=svLt)
fitKt <- nls(vbKt, data=AL_CHJX, weights=1/n,start=svKt)
fitL <- nls(vbL, data=AL_CHJX, weights=1/n,start=svL)
fitt <- nls(vbt, data=AL_CHJX, weights=1/n,start=svt)
fitK <- nls(vbK, data=AL_CHJX, weights=1/n,start=svK)
fit0 <- nls(vb0, data=AL_CHJX, weights=1/n,start=sv0)

#lrt(fit0, com=fitLKt, com.name="All pars differ", sim.name="No pars differ")
#lrt(fitLK, fitLt, fitKt, com=fitLKt, com.name="All pars diff", sim.names=c("Linf, K diff", "Linf,t0 diff", "K, t0, diff"))

lrt(fit0, com=fitLK, com.name="All pars differ", sim.name="No pars differ")
lrt(fitL, fitK, com=fitLK, com.name="All pars diff", sim.names=c("Linf diff", "K diff"))

cbind(AIC(fitLKt,fitLK,fitLt,fitKt, fitL, fitK, fitt, fit0), BIC(fitLKt, fitLK, fitLt, fitKt, fitL, fitK, fitt, fit0))
#Lowest:
# df      AIC df      BIC
# fitLK   6 14919.18  6 14953.75


#CHIR
fitLKt <- nls(vbLKt, data=AL_CHIR, weights=1/n,start=svLKt)
fitLK <- nls(vbLK, data=AL_CHIR, weights=1/n,start=svLK)
fitLt <- nls(vbLt, data=AL_CHIR, weights=1/n,start=svLt)
fitKt <- nls(vbKt, data=AL_CHIR, weights=1/n,start=svKt)
fitL <- nls(vbL, data=AL_CHIR, weights=1/n,start=svL)
fitt <- nls(vbt, data=AL_CHIR, weights=1/n,start=svt)
fitK <- nls(vbK, data=AL_CHIR, weights=1/n,start=svK)
fit0 <- nls(vb0, data=AL_CHIR, weights=1/n,start=sv0)

#lrt(fit0, com=fitLKt, com.name="All pars differ", sim.name="No pars differ")
#lrt(fitLK, fitLt, fitKt, com=fitLKt, com.name="All pars diff", sim.names=c("Linf, K diff", "Linf,t0 diff", "K, t0, diff"))

lrt(fit0, com=fitLK, com.name="All pars differ", sim.name="No pars differ")
lrt(fitL, fitK, com=fitLK, com.name="All pars diff", sim.names=c("Linf diff", "K diff"))


cbind(AIC(fitLKt,fitLK,fitLt,fitKt, fitL, fitK, fitt, fit0), BIC(fitLKt, fitLK, fitLt, fitKt, fitL, fitK, fitt, fit0))
#Lowest:
# df      AIC df      BIC
# fitLKt  7 29066.00  7 29110.49

#APJX
fitLKt <- nls(vbLKt, data=AL_APJX, weights=1/n,start=svLKt)
fitLK <- nls(vbLK, data=AL_APJX, weights=1/n,start=svLK)
fitLt <- nls(vbLt, data=AL_APJX,weights=1/n, start=svLt)
fitKt <- nls(vbKt, data=AL_APJX, weights=1/n,start=svKt)
fitL <- nls(vbL, data=AL_APJX, weights=1/n,start=svL)
fitt <- nls(vbt, data=AL_APJX, weights=1/n,start=svt)
fitK <- nls(vbK, data=AL_APJX, weights=1/n,start=svK)
fit0 <- nls(vb0, data=AL_APJX, weights=1/n,start=sv0)

#lrt(fit0, com=fitLKt, com.name="All pars differ", sim.name="No pars differ")
#lrt(fitLK, fitLt, fitKt, com=fitLKt, com.name="All pars diff", sim.names=c("Linf, K diff", "Linf,t0 diff", "K, t0, diff"))

lrt(fit0, com=fitLK, com.name="All pars differ", sim.name="No pars differ")
lrt(fitL, fitK, com=fitLK, com.name="All pars diff", sim.names=c("Linf diff", "K diff"))

lrt(fit0, com=fitK)

cbind(AIC(fitLKt,fitLK,fitLt,fitKt, fitL, fitK, fitt, fit0), BIC(fitLKt, fitLK, fitLt, fitKt, fitL, fitK, fitt, fit0))
#Lowest:
# df      AIC df      BIC
# fitLKt  7 15694.29  7 15734.89

#APIR
fitLKt <- nls(vbLKt, data=AL_APIR, weights=1/n,start=svLKt)
fitLK <- nls(vbLK, data=AL_APIR,weights=1/n, start=svLK)
fitLt <- nls(vbLt, data=AL_APIR,weights=1/n, start=svLt)
fitKt <- nls(vbKt, data=AL_APIR,weights=1/n, start=svKt)
fitL <- nls(vbL, data=AL_APIR,weights=1/n, start=svL)
fitt <- nls(vbt, data=AL_APIR,weights=1/n, start=svt)
fitK <- nls(vbK, data=AL_APIR,weights=1/n, start=svK)
fit0 <- nls(vb0, data=AL_APIR,weights=1/n, start=sv0)

#lrt(fit0, com=fitLKt, com.name="All pars differ", sim.name="No pars differ")
#lrt(fitLK, fitLt, fitKt, com=fitLKt, com.name="All pars diff", sim.names=c("Linf, K diff", "Linf,t0 diff", "K, t0, diff"))

lrt(fit0, com=fitLK, com.name="All pars differ", sim.name="No pars differ")
lrt(fitL, fitK, com=fitLK, com.name="All pars diff", sim.names=c("Linf diff", "K diff"))


cbind(AIC(fitLKt,fitLK,fitLt,fitKt, fitL, fitK, fitt, fit0), BIC(fitLKt, fitLK, fitLt, fitKt, fitL, fitK, fitt, fit0))
#Lowest:
# df      AIC df      BIC
# fitLKt  7 29782.81  7 29827.46

#JXIR
fitLKt <- nls(vbLKt, data=AL_JXIR, weights=1/n,start=svLKt)
fitLK <- nls(vbLK, data=AL_JXIR, weights=1/n,start=svLK)
fitLt <- nls(vbLt, data=AL_JXIR, weights=1/n,start=svLt)
fitKt <- nls(vbKt, data=AL_JXIR, weights=1/n,start=svKt)
fitL <- nls(vbL, data=AL_JXIR, weights=1/n,start=svL)
fitt <- nls(vbt, data=AL_JXIR, weights=1/n,start=svt)
fitK <- nls(vbK, data=AL_JXIR, weights=1/n,start=svK)
fit0 <- nls(vb0, data=AL_JXIR, weights=1/n,start=sv0)

#lrt(fit0, com=fitLKt, com.name="All pars differ", sim.name="No pars differ")
#lrt(fitLK, fitLt, fitKt, com=fitLKt, com.name="All pars diff", sim.names=c("Linf, K diff", "Linf,t0 diff", "K, t0, diff"))
lrt(fit0, com=fitLK, com.name="All pars differ", sim.name="No pars differ")
lrt(fitL, fitK, com=fitLK, com.name="All pars diff", sim.names=c("Linf diff", "K diff"))

lrt(fit0, com=fitL)


cbind(AIC(fitLKt,fitLK,fitLt,fitKt, fitL, fitK, fitt, fit0), BIC(fitLKt, fitLK, fitLt, fitKt, fitL, fitK, fitt, fit0))
#Lowest:
# df      AIC df      BIC
# fitLKt  7 26993.29  7 27037.26

###########################################
# PLOT COMPARITIVE CURVES
###################################

library(grid)

VonB_Comparative <- ggplot(data= pred_all, aes(x=x, y=pred, group=bay))+  
  #geom_line(aes(linetype=Label), size =.5)+ # make line types based on the different labels- this will be our workaround because in a few stps we will specify the first label (obserseved) be a blank line (therefore a scatter plot)
  geom_line(aes(color=bay))+
  geom_point(aes(shape=bay, color=bay))+ # groups the points together correctly and then gives them a unique shape them differently based on the line type 
  scale_shape_manual(values=c(0,1,2, 3, 4, 5,6,7,8 ))+
  #scale_linetype_manual(values=c('solid', 'dashed', 'dotted'))+
  #scale_y_continuous(limits=c(0,100), breaks= seq(0,100,10))+
  scale_x_continuous(limits=c(0,10), breaks=seq(0,10,1))+
  xlab("Age (years)")+
  ylab("Fork Length (cm)")+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), 									
        panel.background=element_rect(fill='white', colour='black'),
        legend.key=element_blank(), legend.title=element_blank(),
        legend.background=element_rect(fill='white', size=.5),
        legend.position=c(.70,.38),
        legend.key.size= unit(.85, "cm"),
        axis.title.y = element_text(family= "Times New Roman",colour="black", size=20), # changing font of y axis title
        axis.title.x = element_text(family="Times New Roman", colour="black", size=20),
        axis.text.x=element_text(family= "Times New Roman", colour="black", size=16), #changing  colour and font of x axis text
        axis.text.y=element_text(family= "Times New Roman", colour="black", size=16), #changing colour and font of y axis
        #plot.title=element_text(size=14), # changing size of plot title)+
        legend.text=element_text(family= "Times New Roman", size=12),
        legend.key.height=unit(15,"point")) # changes vertical spacing of the legend text
#ggtitle("Comparative von Bertalanffy Growth Models")



