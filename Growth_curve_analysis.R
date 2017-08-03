# ABOUT ####
# 10/13/2016 This script imports ALK with Bay.xlsx
# Main Objectives of this script: 
# 1. Imports otolith data. 
# 2. Determines VBGF parameters for each bay and bootstrapped 95% confidence intervals for each.
# 3. Plots the comparative growth curves
# 4. Compares nested growth models for Bay to Bay comparison to determine which growth parameters are significantly different between bays. 
###############################################################

# load packages####
library(FSA)
library(magrittr)
library(dplyr)
library(nlstools)

#set working directory ####
setwd("~/Desktop/PhD project/Projects/Seatrout/Data")

#load the age 0 length data from the McMichael peters study ####
# change some variables and add some dummy variables so the dataframe will match to the Agelength_Bay loaded next 
age0 <- read.csv("Age_length_mcmichael_peters.csv", header=TRUE) %>% mutate(tl=TL..mm./10, lcat2 =lencat(tl, w=1), specimennumber = rep(NA, 149), bay=rep(NA, 149), Date=rep(NA, 149)) %>% rename(final_age=Otolith.age..years.) 
age0 <- subset(age0, select=-c(TL..mm.))

#load the csv file
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

Agelength_TB<- droplevels(subset(as.data.frame(read.csv("ALK with Bay.csv", header=T)), bay=="TB" & tl>14, select=c(specimennumber, bay, tl, final_age, Date))) %>% filter(!is.na(final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1)) %>%rbind(age0) %>% mutate(bay=rep("TB", 10151)) #, as.fact=TRUE))- can include this when determing ALK below but the smoothed ALK needs to be the nonfactored version of the length categorization variable. 
Agelength_AP<- droplevels(subset(as.data.frame(read.csv("ALK with Bay.csv", header=T)), bay=="AP" & tl>0, select=c(specimennumber, bay, tl, final_age, Date))) %>% filter(!is.na(final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1)) %>%rbind(age0) %>% mutate(bay=rep("AP", 3306)) #, as.fact=TRUE))
Agelength_CK<- droplevels(subset(as.data.frame(read.csv("ALK with Bay.csv", header=T)), bay=="CK" & tl>0, select=c(specimennumber, bay, tl, final_age, Date))) %>% filter(!is.na(final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1)) %>%rbind(age0) %>% mutate(bay=rep("CK", 1402)) # as.fact=TRUE))
Agelength_CH<- droplevels(subset(as.data.frame(read.csv("ALK with Bay.csv", header=T)), bay=="CH" & tl>0, select=c(specimennumber, bay, tl, final_age, Date))) %>% filter(!is.na(final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1)) %>%rbind(age0) %>% mutate(bay=rep("CH", 3682))# as.fact=TRUE))
Agelength_IR<- droplevels(subset(as.data.frame(read.csv("ALK with Bay.csv", header=T)), bay=="IR" & tl>0, select=c(specimennumber, bay, tl, final_age, Date))) %>% filter(!is.na(final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1)) %>%rbind(age0) %>% mutate(bay=rep("IR", 9345))# as.fact=TRUE))
Agelength_JX<- droplevels(subset(as.data.frame(read.csv("ALK with Bay.csv", header=T)), bay=="JX" & tl>0, select=c(specimennumber, bay, tl, final_age, Date))) %>% filter(!is.na(final_age)) %>% mutate(tl=tl/10, lcat2 =lencat(tl, w=1)) %>%rbind(age0) %>% mutate(bay=rep("JX", 1224))# as.fact=TRUE))




# DETERMINE GROWTH PARAMETERS FOR EACH ESTUARY ####
########################################################

# Define function - more flexible when working with a single group of individuals (one estuary)

vbTyp <- function(age, Linf, K, t0) Linf*(1-exp(-K*(age-t0)))
vbTyp <- vbFuns()
#Fit each bay 
  # First define normal starting values
starting <- list(Linf=max(Agelength_TB$tl, na.rm=TRUE), K=0.3, t0=-1)

fitTB <- nls(tl~vbTyp(final_age, Linf, K, t0), data=Agelength_TB, start=starting)
coef(fitTB)
bootTB <- nlsBoot(fitTB)
confint(bootTB) #, plot=TRUE)

# 95% LCI    95% UCI
# Linf 51.6946058 53.2442601
# K     0.3605172  0.4017859
# t0   -1.2764007 -1.1097925

#Visualize the model fit
# - plot the best-fit VBGF with confidence intervals on top of the observed data
# uses a for loop to cycle through all ages

x <- seq(-1,9, length.out=100) # ages for prediction
TB_pred <- vbTyp(x, Linf=coef(fitTB)) #predicted lengths
xlmts <- range(c(x, Agelength_TB$final_age)) #set x limits
ylmts <- range(c(TB_pred, Agelength_TB$tl)) #set y limits
plot(tl~final_age, data=Agelength_TB, xlab="Age", ylab="Total Length (cm)", xlim=xlmts, ylim=ylmts, pch=19, col=rgb(0,0,0,1/3))
lines(TB_pred~x, lwd=2)

# add confidence bands to each model fit ####

# Now for other estuaries

#AP

fitAP <- nls(tl~vbTyp(final_age, Linf, K, t0), data=Agelength_AP, start=starting)
coef(fitAP)
bootAP <- nlsBoot(fitAP)
confint(bootAP) #, plot=TRUE)

# 95% LCI    95% UCI
# Linf 45.7793531 46.8241912
# K     0.8960905  1.0050156
# t0   -0.3088619 -0.2363552

x <- seq(0,9, length.out=100) # ages for prediction
AP_pred <- vbTyp(x, Linf=coef(fitAP)) #predicted lengths
xlmts <- range(c(x, Agelength_AP$final_age)) #set x limits
ylmts <- range(c(AP_pred, Agelength_AP$tl)) #set y limits
plot(tl~final_age, data=Agelength_AP, xlab="Age", ylab="Total Length (cm)", xlim=xlmts, ylim=ylmts, pch=19, col=rgb(0,0,0,1/3))
lines(AP_pred~x, lwd=2)

#CK

fitCK <- nls(tl~vbTyp(final_age, Linf, K, t0), data=Agelength_CK, start=starting)
coef(fitCK)
bootCK <- nlsBoot(fitCK)
confint(bootCK) #, plot=TRUE)

# 95% LCI    95% UCI
# Linf 50.7272826 54.3839631
# K     0.4637576  0.5833411
# t0   -0.6572329 -0.4930606

x <- seq(0,9, length.out=100) # ages for prediction
CK_pred <- vbTyp(x, Linf=coef(fitCK)) #predicted lengths
xlmts <- range(c(x, Agelength_CK$final_age)) #set x limits
ylmts <- range(c(CK_pred, Agelength_CK$tl)) #set y limits
plot(tl~final_age, data=Agelength_CK, xlab="Age", ylab="Total Length (cm)", xlim=xlmts, ylim=ylmts, pch=19, col=rgb(0,0,0,1/3))
lines(CK_pred~x, lwd=2)

#CH
fitCH <- nls(tl~vbTyp(final_age, Linf, K, t0), data=Agelength_CH, start=starting)
coef(fitCH)
bootCH <- nlsBoot(fitCH)
confint(bootCH) #, plot=TRUE)

# 95% LCI    95% UCI
# Linf 44.9099417 46.0548986
# K     0.7437704  0.8422404
# t0   -0.3842474 -0.2768528

x <- seq(0,9, length.out=100) # ages for prediction
CH_pred <- vbTyp(x, Linf=coef(fitCH)) #predicted lengths
xlmts <- range(c(x, Agelength_CH$final_age)) #set x limits
ylmts <- range(c(CH_pred, Agelength_CH$tl)) #set y limits
plot(tl~final_age, data=Agelength_CH, xlab="Age", ylab="Total Length (cm)", xlim=xlmts, ylim=ylmts, pch=19, col=rgb(0,0,0,1/3))
lines(CH_pred~x, lwd=2)

#IR
fitIR <- nls(tl~vbTyp(final_age, Linf, K, t0), data=Agelength_IR, start=starting)
coef(fitIR)
bootIR <- nlsBoot(fitIR)
confint(bootIR) #, plot=TRUE)

# 95% LCI    95% UCI
# Linf 57.7000010 60.1145785
# K     0.3289741  0.3759782
# t0   -1.0141367 -0.8289982

x <- seq(0,9, length.out=100) # ages for prediction
IR_pred <- vbTyp(x, Linf=coef(fitIR)) #predicted lengths
xlmts <- range(c(x, Agelength_IR$final_age)) #set x limits
ylmts <- range(c(IR_pred, Agelength_IR$tl)) #set y limits
plot(tl~final_age, data=Agelength_IR, xlab="Age", ylab="Total Length (cm)", xlim=xlmts, ylim=ylmts, pch=19, col=rgb(0,0,0,1/3))
lines(IR_pred~x, lwd=2)

#JX
fitJX <- nls(tl~vbTyp(final_age, Linf, K, t0), data=Agelength_JX, start=starting)
coef(fitJX)
bootJX <- nlsBoot(fitJX)
confint(bootJX) #, plot=TRUE)

# 95% LCI    95% UCI
# Linf 43.4891047 47.2441419
# K     0.5845976  0.7564375
# t0   -0.5016745 -0.3470104


x <- seq(0,9, length.out=100) # ages for prediction
JX_pred <- vbTyp(x, Linf=coef(fitJX)) #predicted lengths
xlmts <- range(c(x, Agelength_JX$final_age)) #set x limits
ylmts <- range(c(JX_pred, Agelength_JX$tl)) #set y limits
plot(tl~final_age, data=Agelength_JX, xlab="Age", ylab="Total Length (cm)", xlim=xlmts, ylim=ylmts, pch=19, col=rgb(0,0,0,1/3))
lines(JX_pred~x, lwd=2)

####################################
#PLOT ALL PREDICTED FITS IN ONE PLOT ####
####################################

# combine all predicted matrices into one
x <- seq(0,9, length.out=100)
t= data.frame(cbind(x,TB_pred)) %>% rename(pred=TB_pred) %>% mutate(bay=rep("TB",100))
a= data.frame(cbind(x,AP_pred)) %>% rename(pred=AP_pred) %>% mutate(bay=rep("AP",100))
ck= data.frame(cbind(x,CK_pred)) %>% rename(pred=CK_pred) %>% mutate(bay=rep("CK",100))
c= data.frame(cbind(x,CH_pred)) %>% rename(pred=CH_pred) %>% mutate(bay=rep("CH",100))
j= data.frame(cbind(x,JX_pred))%>% rename(pred=JX_pred) %>% mutate(bay=rep("JX",100))
i= data.frame(cbind(x,IR_pred)) %>% rename(pred=IR_pred) %>% mutate(bay=rep("IR",100))


pred_all = rbind(t,a,ck,c,j,i)

library(grid)
library(ggplot2)

VonB_Comparative <- ggplot(data= pred_all, aes(x=x, y=pred, group=bay))+  
  geom_line() + #aes(linetype=bay), size =.5)+ # make line types based on the different labels- this will be our workaround because in a few stps we will specify the first label (obserseved) be a blank line (therefore a scatter plot)
  geom_line(aes(color=bay))+
  #geom_point(aes(shape=bay), size=2) + #, color=bay))+ # groups the points together correctly and then gives them a unique shape them differently based on the line type 
  #scale_shape_manual(values=c(0,1,2, 3, 4, 5,6,7,8 ))+
  #scale_linetype_manual(values=c('solid', 'dashed', 'dotted'))+
  scale_y_continuous(limits=c(0,60), breaks= seq(0,60,10))+
  scale_x_continuous(limits=c(0,10), breaks=seq(0,9,1))+
  xlab("Age (years)")+
  ylab("Total Length (cm)")+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), 									
        panel.background=element_rect(fill='white', colour='black'),
        legend.key=element_blank(), legend.title=element_blank(),
        legend.background=element_rect(fill='white', size=.5),
        legend.position=c(.70,.38),
        legend.key.size= unit(3, "cm"),
        axis.title.y = element_text(family= "Times New Roman",colour="black", size=20), # changing font of y axis title
        axis.title.x = element_text(family="Times New Roman", colour="black", size=20),
        axis.text.x=element_text(family= "Times New Roman", colour="black", size=16), #changing  colour and font of x axis text
        axis.text.y=element_text(family= "Times New Roman", colour="black", size=16), #changing colour and font of y axis
        #plot.title=element_text(size=14), # changing size of plot title)+
        legend.text=element_text(family= "Times New Roman", size=18),
        legend.key.height=unit(25,"point")) # changes vertical spacing of the legend text
#ggtitle("Comparative von Bertalanffy Growth Models")

##############################
# COMPARE VGBFs BETWEEN GROUPS ####
##############################
# requires fitting multiple models to determine which parameters between each group is different. 
# For simplicities sake, I will only compare 2 bays at a time. Can refer to ALK_analysis.R for a list of bay to bay comparisons. 
# Refer to page 236 in Ogle for the family of models that must be considered when examining the differences in VBGF among groups. 
# The nested relationships among these models allows use of likelihoo ratio and extra sum of squares tests.

########################################
# Define the generic models to compare
#######################################
vbLKt <- tl~Linf[bay]*(1-exp(-K[bay]*(final_age-t0[bay])))
vbLK <- tl~Linf[bay]*(1-exp(-K[bay]*(final_age-t0)))
vbLt <- tl~Linf[bay]*(1-exp(-K*(final_age-t0[bay])))
vbKt <- tl~Linf*(1-exp(-K[bay]*(final_age-t0[bay])))
vbL <-  tl~Linf[bay]*(1-exp(-K*(final_age-t0)))
vbK <-  tl~Linf*(1-exp(-K[bay]*(final_age-t0)))
vbt <- tl~Linf*(1-exp(-K*(final_age-t0[bay])))
vb0 <- tl~Linf*(1-exp(-K*(final_age-t0)))

###############################################
#Define starting values for each model scenario
###############################################
sv0 <- vbStarts(tl~final_age, data=Agelength_TB) #using the TB dataframe as a starter
sv0$Linf= 72
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
Agelength_TBCK<- rbind(Agelength_TB, Agelength_CK) 
  Agelength_TBCK$bay <- as.factor(Agelength_TBCK$bay)
Agelength_TBCH<- rbind(Agelength_TB, Agelength_CH)
  Agelength_TBCH$bay <- as.factor(Agelength_TBCH$bay)
Agelength_TBAP<- rbind(Agelength_TB, Agelength_AP)
  Agelength_TBAP$bay <- as.factor(Agelength_TBAP$bay)
Agelength_TBJX<- rbind(Agelength_TB, Agelength_JX)
  Agelength_TBJX$bay <- as.factor(Agelength_TBJX$bay) 
Agelength_TBIR<- rbind(Agelength_TB, Agelength_IR)
  Agelength_TBIR$bay <- as.factor(Agelength_TBIR$bay)


Agelength_CKCH<- rbind(Agelength_CK, Agelength_CH)
  Agelength_CKCH$bay <- as.factor(Agelength_CKCH$bay)
Agelength_CKAP<- rbind(Agelength_CK, Agelength_AP) 
  Agelength_CKAP$bay <- as.factor(Agelength_CKAP$bay)
Agelength_CKJX<- rbind(Agelength_CK, Agelength_JX)
  Agelength_CKJX$bay <- as.factor(Agelength_CKJX$bay)
Agelength_CKIR<- rbind(Agelength_CK, Agelength_IR)
  Agelength_CKIR$bay <- as.factor(Agelength_CKIR$bay)

Agelength_CHAP<- rbind(Agelength_CH, Agelength_AP)
  Agelength_CHAP$bay <- as.factor(Agelength_CHAP$bay)
Agelength_CHJX<- rbind(Agelength_CH, Agelength_JX)
  Agelength_CHJX$bay <- as.factor(Agelength_CHJX$bay)
Agelength_CHIR<- rbind(Agelength_CH, Agelength_IR)
  Agelength_CHIR$bay <- as.factor(Agelength_CHIR$bay)
  
Agelength_APJX<- rbind(Agelength_AP, Agelength_JX)
  Agelength_APJX$bay <- as.factor(Agelength_APJX$bay)
Agelength_APIR<- rbind(Agelength_AP, Agelength_IR)
  Agelength_APIR$bay <- as.factor(Agelength_APIR$bay)
  
Agelength_JXIR<- rbind(Agelength_JX, Agelength_IR)
  Agelength_JXIR$bay <- as.factor(Agelength_JXIR$bay)

###############################################################
#Define fits to model scenarios (one bay at a time) _TBCK here
##############################################################
#Fit to [LKt] model then to the simplest model [vb0]
fitLKt <- nls(vbLKt, data=Agelength_TBCK, start=svLKt)
  #residPlot(fitLKt, col=rgb(0,0,0,1/3))
fitLK <- nls(vbLK, data=Agelength_TBCK, start=svLK)
fitLt <- nls(vbLt, data=Agelength_TBCK, start=svLt)
fitKt <- nls(vbKt, data=Agelength_TBCK, start=svKt)
fitL <- nls(vbL, data=Agelength_TBCK, start=svL)
fitt <- nls(vbt, data=Agelength_TBCK, start=svt)
fitK <- nls(vbK, data=Agelength_TBCK, start=svK)
fit0 <- nls(vb0, data=Agelength_TBCK, start=sv0)


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

# ***** Any nested model that is not statistically different from[L,K,t0] is considered "better" because it fits equally (statistically) well, but is more parsimonious.
#       If two models are better than [L, K, t0] then the one with the greatest LL is chosen to be the best of the nested models. 


##############################################################
# Easier AIC/BIC example that could be used in place of the LRT
###############################################################
cbind(AIC(fitLKt,fitLK,fitLt,fitKt, fitL, fitK, fitt, fit0), BIC(fitLKt, fitLK, fitLt, fitKt, fitL, fitK, fitt, fit0))
#Lowest:
# #       df      AIC df      BIC
# fitLKt  7 75503.79  7 75555.27


#################################################################
# Now do it all over again with each Bay to bay comparison
##################################################################

# TBCH
fitLKt <- nls(vbLKt, data=Agelength_TBCH, start=svLKt)
fitLK <- nls(vbLK, data=Agelength_TBCH, start=svLK)
fitLt <- nls(vbLt, data=Agelength_TBCH, start=svLt)
fitKt <- nls(vbKt, data=Agelength_TBCH, start=svKt)
fitL <- nls(vbL, data=Agelength_TBCH, start=svL)
fitt <- nls(vbt, data=Agelength_TBCH, start=svt)
fitK <- nls(vbK, data=Agelength_TBCH, start=svK)
fit0 <- nls(vb0, data=Agelength_TBCH, start=sv0)

lrt(fit0, com=fitLKt, com.name="All pars differ", sim.name="No pars differ")

cbind(AIC(fitLKt,fitLK,fitLt,fitKt, fitL, fitK, fitt, fit0), BIC(fitLKt, fitLK, fitLt, fitKt, fitL, fitK, fitt, fit0))
#Lowest:
# df      AIC df      BIC
# fitLKt  7 90412.66  7 90465.41

#TBAP
fitLKt <- nls(vbLKt, data=Agelength_TBAP, start=svLKt)
fitLK <- nls(vbLK, data=Agelength_TBAP, start=svLK)
fitLt <- nls(vbLt, data=Agelength_TBAP, start=svLt)
fitKt <- nls(vbKt, data=Agelength_TBAP, start=svKt)
fitL <- nls(vbL, data=Agelength_TBAP, start=svL)
fitt <- nls(vbt, data=Agelength_TBAP, start=svt)
fitK <- nls(vbK, data=Agelength_TBAP, start=svK)
fit0 <- nls(vb0, data=Agelength_TBAP, start=sv0)

lrt(fit0, com=fitLKt, com.name="All pars differ", sim.name="No pars differ")

cbind(AIC(fitLKt,fitLK,fitLt,fitKt, fitL, fitK, fitt, fit0), BIC(fitLKt, fitLK, fitLt, fitKt, fitL, fitK, fitt, fit0))
#Lowest:
# df      AIC df      BIC
# fitLKt  7 87299.69  7 87352.24


#TBJX
fitLKt <- nls(vbLKt, data=Agelength_TBJX, start=svLKt)
fitLK <- nls(vbLK, data=Agelength_TBJX, start=svLK)
fitLt <- nls(vbLt, data=Agelength_TBJX, start=svLt)
fitKt <- nls(vbKt, data=Agelength_TBJX, start=svKt)
fitL <- nls(vbL, data=Agelength_TBJX, start=svL)
fitt <- nls(vbt, data=Agelength_TBJX, start=svt)
fitK <- nls(vbK, data=Agelength_TBJX, start=svK)
fit0 <- nls(vb0, data=Agelength_TBJX, start=sv0)

lrt(fit0, com=fitLKt, com.name="All pars differ", sim.name="No pars differ")

cbind(AIC(fitLKt,fitLK,fitLt,fitKt, fitL, fitK, fitt, fit0), BIC(fitLKt, fitLK, fitLt, fitKt, fitL, fitK, fitt, fit0))
#Lowest:
# df      AIC df      BIC
# fitLKt  7 74170.11  7 74221.48

#TBIR
fitLKt <- nls(vbLKt, data=Agelength_TBIR, start=svLKt)
fitLK <- nls(vbLK, data=Agelength_TBIR, start=svLK)
fitLt <- nls(vbLt, data=Agelength_TBIR, start=svLt)
fitKt <- nls(vbKt, data=Agelength_TBIR, start=svKt)
fitL <- nls(vbL, data=Agelength_TBIR, start=svL)
fitt <- nls(vbt, data=Agelength_TBIR, start=svt)
fitK <- nls(vbK, data=Agelength_TBIR, start=svK)
fit0 <- nls(vb0, data=Agelength_TBIR, start=sv0)

lrt(fit0, com=fitLKt, com.name="All pars differ", sim.name="No pars differ")

cbind(AIC(fitLKt,fitLK,fitLt,fitKt, fitL, fitK, fitt, fit0), BIC(fitLKt, fitLK, fitLt, fitKt, fitL, fitK, fitt, fit0))
#Lowest:
# df      AIC df      BIC
# fitLKt   6 131645.2  7 131700.4

#CKCH
fitLKt <- nls(vbLKt, data=Agelength_CKCH, start=svLKt)
fitLK <- nls(vbLK, data=Agelength_CKCH, start=svLK)
fitLt <- nls(vbLt, data=Agelength_CKCH, start=svLt)
fitKt <- nls(vbKt, data=Agelength_CKCH, start=svKt)
fitL <- nls(vbL, data=Agelength_CKCH, start=svL)
fitt <- nls(vbt, data=Agelength_CKCH, start=svt)
fitK <- nls(vbK, data=Agelength_CKCH, start=svK)
fit0 <- nls(vb0, data=Agelength_CKCH, start=sv0)

lrt(fit0, com=fitLKt, com.name="All pars differ", sim.name="No pars differ")

cbind(AIC(fitLKt,fitLK,fitLt,fitKt, fitL, fitK, fitt, fit0), BIC(fitLKt, fitLK, fitLt, fitKt, fitL, fitK, fitt, fit0))
#Lowest:
# df      AIC df      BIC
# fitLKt   7 34053.47  7 34099.21

#CKAP
fitLKt <- nls(vbLKt, data=Agelength_CKAP, start=svLKt)
fitLK <- nls(vbLK, data=Agelength_CKAP, start=svLK)
fitLt <- nls(vbLt, data=Agelength_CKAP, start=svLt)
fitKt <- nls(vbKt, data=Agelength_CKAP, start=svKt)
fitL <- nls(vbL, data=Agelength_CKAP, start=svL)
fitt <- nls(vbt, data=Agelength_CKAP, start=svt)
fitK <- nls(vbK, data=Agelength_CKAP, start=svK)
fit0 <- nls(vb0, data=Agelength_CKAP, start=sv0)

lrt(fit0, com=fitLKt, com.name="All pars differ", sim.name="No pars differ")

cbind(AIC(fitLKt,fitLK,fitLt,fitKt, fitL, fitK, fitt, fit0), BIC(fitLKt, fitLK, fitLt, fitKt, fitL, fitK, fitt, fit0))
#Lowest:
# df      AIC df      BIC
# fitLKt  7 31019.14  7 31064.34


#CKJX
fitLKt <- nls(vbLKt, data=Agelength_CKJX, start=svLKt)
fitLK <- nls(vbLK, data=Agelength_CKJX, start=svLK)
fitLt <- nls(vbLt, data=Agelength_CKJX, start=svLt)
fitKt <- nls(vbKt, data=Agelength_CKJX, start=svKt)
fitL <- nls(vbL, data=Agelength_CKJX, start=svL)
fitt <- nls(vbt, data=Agelength_CKJX, start=svt)
fitK <- nls(vbK, data=Agelength_CKJX, start=svK)
fit0 <- nls(vb0, data=Agelength_CKJX, start=sv0)

lrt(fit0, com=fitLKt, com.name="All pars differ", sim.name="No pars differ")

cbind(AIC(fitLKt,fitLK,fitLt,fitKt, fitL, fitK, fitt, fit0), BIC(fitLKt, fitLK, fitLt, fitKt, fitL, fitK, fitt, fit0))
#lowest:
# df      AIC df      BIC
# fitLKt  7 17800.56  7 17841.67
# fitLK   6 17803.66  6 17838.90


#CKIR
fitLKt <- nls(vbLKt, data=Agelength_CKIR, start=svLKt)
fitLK <- nls(vbLK, data=Agelength_CKIR, start=svLK)
fitLt <- nls(vbLt, data=Agelength_CKIR, start=svLt)
fitKt <- nls(vbKt, data=Agelength_CKIR, start=svKt)
fitL <- nls(vbL, data=Agelength_CKIR, start=svL)
fitt <- nls(vbt, data=Agelength_CKIR, start=svt)
fitK <- nls(vbK, data=Agelength_CKIR, start=svK)
fit0 <- nls(vb0, data=Agelength_CKIR, start=sv0)

lrt(fit0, com=fitLKt, com.name="All pars differ", sim.name="No pars differ")

cbind(AIC(fitLKt,fitLK,fitLt,fitKt, fitL, fitK, fitt, fit0), BIC(fitLKt, fitLK, fitLt, fitKt, fitL, fitK, fitt, fit0))
#Lowest:
# df      AIC df      BIC
# fitLKt  7 74754.55  7 74805.53

#CHAP
fitLKt <- nls(vbLKt, data=Agelength_CHAP, start=svLKt)
fitLK <- nls(vbLK, data=Agelength_CHAP, start=svLK)
fitLt <- nls(vbLt, data=Agelength_CHAP, start=svLt)
fitKt <- nls(vbKt, data=Agelength_CHAP, start=svKt)
fitL <- nls(vbL, data=Agelength_CHAP, start=svL)
fitt <- nls(vbt, data=Agelength_CHAP, start=svt)
fitK <- nls(vbK, data=Agelength_CHAP, start=svK)
fit0 <- nls(vb0, data=Agelength_CHAP, start=sv0)

lrt(fit0, com=fitLKt, com.name="All pars differ", sim.name="No pars differ")

cbind(AIC(fitLKt,fitLK,fitLt,fitKt, fitL, fitK, fitt, fit0), BIC(fitLKt, fitLK, fitLt, fitKt, fitL, fitK, fitt, fit0))
#Lowest:
# df      AIC df      BIC
#fitLKt  7 45929.66  7 45977.62
#fitLK   6 45930.26  6 45971.37

#CHJX
fitLKt <- nls(vbLKt, data=Agelength_CHJX, start=svLKt)
fitLK <- nls(vbLK, data=Agelength_CHJX, start=svLK)
fitLt <- nls(vbLt, data=Agelength_CHJX, start=svLt)
fitKt <- nls(vbKt, data=Agelength_CHJX, start=svKt)
fitL <- nls(vbL, data=Agelength_CHJX, start=svL)
fitt <- nls(vbt, data=Agelength_CHJX, start=svt)
fitK <- nls(vbK, data=Agelength_CHJX, start=svK)
fit0 <- nls(vb0, data=Agelength_CHJX, start=sv0)

lrt(fit0, com=fitLKt, com.name="All pars differ", sim.name="No pars differ")

cbind(AIC(fitLKt,fitLK,fitLt,fitKt, fitL, fitK, fitt, fit0), BIC(fitLKt, fitLK, fitLt, fitKt, fitL, fitK, fitt, fit0))
#Lowest:
# df      AIC df      BIC
# fitLKt  7 32742.65  7 32788.13
# fitLK   6 32743.79  6 32782.78

#CHIR
fitLKt <- nls(vbLKt, data=Agelength_CHIR, start=svLKt)
fitLK <- nls(vbLK, data=Agelength_CHIR, start=svLK)
fitLt <- nls(vbLt, data=Agelength_CHIR, start=svLt)
fitKt <- nls(vbKt, data=Agelength_CHIR, start=svKt)
fitL <- nls(vbL, data=Agelength_CHIR, start=svL)
fitt <- nls(vbt, data=Agelength_CHIR, start=svt)
fitK <- nls(vbK, data=Agelength_CHIR, start=svK)
fit0 <- nls(vb0, data=Agelength_CHIR, start=sv0)

lrt(fit0, com=fitLKt, com.name="All pars differ", sim.name="No pars differ")

cbind(AIC(fitLKt,fitLK,fitLt,fitKt, fitL, fitK, fitt, fit0), BIC(fitLKt, fitLK, fitLt, fitKt, fitL, fitK, fitt, fit0))
#Lowest:
# df      AIC df      BIC
# fitLKt  7 89811.60  7 89863.92

#APJX
fitLKt <- nls(vbLKt, data=Agelength_APJX, start=svLKt)
fitLK <- nls(vbLK, data=Agelength_APJX, start=svLK)
fitLt <- nls(vbLt, data=Agelength_APJX, start=svLt)
fitKt <- nls(vbKt, data=Agelength_APJX, start=svKt)
fitL <- nls(vbL, data=Agelength_APJX, start=svL)
fitt <- nls(vbt, data=Agelength_APJX, start=svt)
fitK <- nls(vbK, data=Agelength_APJX, start=svK)
fit0 <- nls(vb0, data=Agelength_APJX, start=sv0)

lrt(fit0, com=fitLKt, com.name="All pars differ", sim.name="No pars differ")

cbind(AIC(fitLKt,fitLK,fitLt,fitKt, fitL, fitK, fitt, fit0), BIC(fitLKt, fitLK, fitLt, fitKt, fitL, fitK, fitt, fit0))
#Lowest:
# df      AIC df      BIC
# fitLKt  7 29692.38  7 29737.31

#APIR
fitLKt <- nls(vbLKt, data=Agelength_APIR, start=svLKt)
fitLK <- nls(vbLK, data=Agelength_APIR, start=svLK)
fitLt <- nls(vbLt, data=Agelength_APIR, start=svLt)
fitKt <- nls(vbKt, data=Agelength_APIR, start=svKt)
fitL <- nls(vbL, data=Agelength_APIR, start=svL)
fitt <- nls(vbt, data=Agelength_APIR, start=svt)
fitK <- nls(vbK, data=Agelength_APIR, start=svK)
fit0 <- nls(vb0, data=Agelength_APIR, start=sv0)

lrt(fit0, com=fitLKt, com.name="All pars differ", sim.name="No pars differ")

cbind(AIC(fitLKt,fitLK,fitLt,fitKt, fitL, fitK, fitt, fit0), BIC(fitLKt, fitLK, fitLt, fitKt, fitL, fitK, fitt, fit0))
#Lowest:
# df      AIC df      BIC
# fitLKt  7 86881.09  7 86933.20

#JXIR
fitLKt <- nls(vbLKt, data=Agelength_JXIR, start=svLKt)
fitLK <- nls(vbLK, data=Agelength_JXIR, start=svLK)
fitLt <- nls(vbLt, data=Agelength_JXIR, start=svLt)
fitKt <- nls(vbKt, data=Agelength_JXIR, start=svKt)
fitL <- nls(vbL, data=Agelength_JXIR, start=svL)
fitt <- nls(vbt, data=Agelength_JXIR, start=svt)
fitK <- nls(vbK, data=Agelength_JXIR, start=svK)
fit0 <- nls(vb0, data=Agelength_JXIR, start=sv0)

lrt(fit0, com=fitLKt, com.name="All pars differ", sim.name="No pars differ")

cbind(AIC(fitLKt,fitLK,fitLt,fitKt, fitL, fitK, fitt, fit0), BIC(fitLKt, fitLK, fitLt, fitKt, fitL, fitK, fitt, fit0))
#Lowest:
# df      AIC df      BIC
# fitLKt  7 73466.40  7 73517.26

###########################################
# PLOT COMPARITIVE CURVES
###################################

library(grid)

VonB_Comparative <- ggplot(data= df_comparative_rearrange, aes(x=x, y=y, group=Label))+  
  #geom_line(aes(linetype=Label), size =.5)+ # make line types based on the different labels- this will be our workaround because in a few stps we will specify the first label (obserseved) be a blank line (therefore a scatter plot)
  geom_line(aes(color=Label))+
  geom_point(aes(shape=Label, color=Label))+ # groups the points together correctly and then gives them a unique shape them differently based on the line type 
  scale_shape_manual(values=c(0,1,2, 3, 4, 5,6,7,8 ))+
  #scale_linetype_manual(values=c('solid', 'dashed', 'dotted'))+
  #scale_y_continuous(limits=c(0,100), breaks= seq(0,100,10))+
  scale_x_continuous(limits=c(0,30), breaks=seq(0,30,5))+
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



