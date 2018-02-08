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

Agelength_AP<- droplevels(subset(as.data.frame(read.csv("ALK_Bay_and_weight.csv", header=T)), bay=="AP" & tl>0 & final_age >0 & program == 'FIM', select=c( bay, tl, final_age, date, sex))) %>% filter(!is.na(final_age)) %>% 
          mutate(tl=tl/10, lcat2 =lencat(tl, w=1)) %>%rbind(age0)
Agelength_AP$bay="AP"

Agelength_CK<- droplevels(subset(as.data.frame(read.csv("ALK_Bay_and_weight.csv", header=T)), bay=="CK" & tl>0 & final_age >0 & program=='FIM', select=c( bay, tl, final_age, date, sex))) %>% filter(!is.na(final_age)) %>% 
          mutate(tl=tl/10, lcat2 =lencat(tl, w=1)) %>%rbind(age0)
Agelength_CK$bay="CK"

Agelength_CH<- droplevels(subset(as.data.frame(read.csv("ALK_Bay_and_weight.csv", header=T)), bay=="CH" & tl>0 & final_age >0 & program=='FIM', select=c( bay, tl, final_age, date, sex))) %>% filter(!is.na(final_age)) %>% 
          mutate(tl=tl/10, lcat2 =lencat(tl, w=1)) %>%rbind(age0) 
Agelength_CH$bay="CH"

Agelength_IR<- droplevels(subset(as.data.frame(read.csv("ALK_Bay_and_weight.csv", header=T)), bay=="IR" & tl>0 & final_age >0 & program=='FIM', select=c( bay, tl, final_age, date, sex))) %>% filter(!is.na(final_age)) %>% 
          mutate(tl=tl/10, lcat2 =lencat(tl, w=1)) %>%rbind(age0) 
Agelength_IR$bay="IR"

Agelength_JX<- droplevels(subset(as.data.frame(read.csv("ALK_Bay_and_weight.csv", header=T)), bay=="JX" & tl>0 & final_age >0 & program=='FIM', select=c( bay, tl, final_age, date, sex))) %>% filter(!is.na(final_age)) %>% 
          mutate(tl=tl/10, lcat2 =lencat(tl, w=1)) %>%rbind(age0) 
Agelength_JX$bay="JX"



# DETERMINE GROWTH PARAMETERS FOR EACH ESTUARY ####


# Define function - more flexible when working with a single group of individuals (one estuary)

vbTyp <- function(age, Linf, K, t0) Linf*(1-exp(-K*(age-t0)))
vbTyp <- vbFuns()
#Fit each bay 
  # First define normal starting values

# TB ####
TB_M <- droplevels(subset(Agelength_TB, sex == "M"))
TB_F <- droplevels(subset(Agelength_TB, sex == "F"))

starting <- list(Linf=max(Agelength_TB$tl, na.rm=TRUE), K=0.3, t0=-1)

fitTB <- nls(tl~vbTyp(final_age, Linf, K, t0), data=Agelength_TB, start=starting)
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
TB_pred <- vbTyp(x, Linf=coef(fitTB)) #predicted lengths
xlmts <- range(c(x, Agelength_TB$final_age)) #set x limits
ylmts <- range(c(TB_pred, Agelength_TB$tl)) #set y limits
plot(tl~final_age, data=Agelength_TB, xlab="Age", ylab="Total Length (cm)", xlim=xlmts, ylim=ylmts, pch=19, col=rgb(0,0,0,1/3))
lines(TB_pred~x, lwd=2)

# add confidence bands to each model fit ####

# Now for other estuaries

#Male
starting <- list(Linf=max(TB_M$tl, na.rm=TRUE), K=0.3, t0=-1)
fitTBM <- nls(tl~vbTyp(final_age, Linf, K, t0), data=TB_M, start=starting, control=list(maxiter=500))
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
# 
#Feamle
starting <- list(Linf=max(TB_F$tl, na.rm=TRUE), K=0.3, t0=-1)
fitTBF <- nls(tl~vbTyp(final_age, Linf, K, t0), data=TB_F, start=starting, control=list(maxiter=500))
coef(fitTBF)
bootTBF <- nlsBoot(fitTBF)
confint(bootTBF) #, plot=TRUE)

# 95% LCI     95% UCI
# Linf 45.37935465 46.89330169
# K     0.94476629  1.10744925
# t0   -0.06095907  0.03363978

# Linf           K          t0 
# 46.10298743  1.02473610 -0.01011952 


#AP ####
AP_M <- droplevels(subset(Agelength_AP, sex == "M"))
AP_F <- droplevels(subset(Agelength_AP, sex == "F"))

starting <- list(Linf=max(Agelength_AP$tl, na.rm=TRUE), K=0.3, t0=-1)

fitAP <- nls(tl~vbTyp(final_age, Linf, K, t0), data=Agelength_AP, start=starting)
coef(fitAP)
bootAP <- nlsBoot(fitAP)
confint(bootAP) #, plot=TRUE)

# 95% LCI     95% UCI
# Linf 46.878365576 48.38449621
# K     0.874009111  0.99885236
# t0   -0.007496384  0.05852243

x <- seq(0,9, length.out=30) # ages for prediction
AP_pred <- vbTyp(x, Linf=coef(fitAP)) #predicted lengths
xlmts <- range(c(x, Agelength_AP$final_age)) #set x limits
ylmts <- range(c(AP_pred, Agelength_AP$tl)) #set y limits
plot(tl~final_age, data=Agelength_AP, xlab="Age", ylab="Total Length (cm)", xlim=xlmts, ylim=ylmts, pch=19, col=rgb(0,0,0,1/3))
lines(AP_pred~x, lwd=2)

#Female 
starting <- list(Linf=max(AP_F$tl, na.rm=TRUE), K=0.3, t0=-1)
fitAPF <- nls(tl~vbTyp(final_age, Linf, K, t0), data=AP_F, start=starting, control=list(maxiter=500))
coef(fitAPF)
bootAPF <- nlsBoot(fitAPF)
confint(bootAPF) #, plot=TRUE)

# Linf           K          t0 
# 50.79605092  0.85143347 -0.01465805 
# 
# 95% LCI     95% UCI
# Linf 49.92990162 51.73683338
# K     0.79000560  0.91375676
# t0   -0.05855143  0.02532161

#Male
starting <- list(Linf=max(AP_M$tl, na.rm=TRUE), K=0.3, t0=-1)
fitAPM <- nls(tl~vbTyp(final_age, Linf, K, t0), data=AP_M, start=starting, control=list(maxiter=500))
coef(fitAPM)
bootAPM <- nlsBoot(fitAPM)
confint(bootAPM) #, plot=TRUE)
# 
# Linf           K          t0 
# 38.90189512  1.20154408  0.05603756 
# 
# 95% LCI     95% UCI
# Linf 38.16540150 39.69453106
# K     1.09273584  1.33316898
# t0    0.02572991  0.08261692

#CK ####
CK_M <- droplevels(subset(Agelength_CK, sex == "M"))
CK_F <- droplevels(subset(Agelength_CK, sex == "F"))

starting <- list(Linf=max(Agelength_CK$tl, na.rm=TRUE), K=0.3, t0=-1)

fitCK <- nls(tl~vbTyp(final_age, Linf, K, t0), data=Agelength_CK, start=starting)
coef(fitCK)
bootCK <- nlsBoot(fitCK)
confint(bootCK) #, plot=TRUE)

# 95% LCI     95% UCI
# Linf 43.99012512 45.64073597
# K     1.12842441  1.29351522
# t0    0.04542966  0.09576885

x <- seq(0,9, length.out=30) # ages for prediction
CK_pred <- vbTyp(x, Linf=coef(fitCK)) #predicted lengths
xlmts <- range(c(x, Agelength_CK$final_age)) #set x limits
ylmts <- range(c(CK_pred, Agelength_CK$tl)) #set y limits
plot(tl~final_age, data=Agelength_CK, xlab="Age", ylab="Total Length (cm)", xlim=xlmts, ylim=ylmts, pch=19, col=rgb(0,0,0,1/3))
lines(CK_pred~x, lwd=2)

#Female
starting <- list(Linf=max(CK_F$tl, na.rm=TRUE), K=0.3, t0=-1)
fitCKF <- nls(tl~vbTyp(final_age, Linf, K, t0), data=CK_F, start=starting, control=list(maxiter=500))
coef(fitCKF)
bootCKF <- nlsBoot(fitCKF)
confint(bootCKF) #, plot=TRUE)

# Linf           K          t0 
# 46.05784259  1.29560014  0.08728294 

# 95% LCI    95% UCI
# Linf 45.282304 46.8522702
# K     1.211701  1.3851679
# t0    0.060295  0.1116069

#Male 
fitCKM <- nls(tl~vbTyp(final_age, Linf, K, t0), data=CK_M, start=starting, control=list(maxiter=500))
coef(fitCKM)
bootCKM <- nlsBoot(fitCKM)
confint(bootCKM) #, plot=TRUE)

# Linf           K          t0 
# 38.57439632  1.46468496  0.07560443 
# 
# 95% LCI     95% UCI
# Linf 37.4733618 39.85064429
# K     1.2839761  1.64244930
# t0    0.0452437  0.09948328


#CH ####

CH_M <- droplevels(subset(Agelength_CH, sex == "M"))
CH_F <- droplevels(subset(Agelength_CH, sex == "F"))

starting <- list(Linf=max(Agelength_CH$tl, na.rm=TRUE), K=0.3, t0=-1)

fitCH <- nls(tl~vbTyp(final_age, Linf, K, t0), data=Agelength_CH, start=starting)
coef(fitCH)
bootCH <- nlsBoot(fitCH)
confint(bootCH) #, plot=TRUE)

# 95% LCI     95% UCI
# Linf 41.23849994 42.32673242
# K     1.25632886  1.43022834
# t0    0.04767128  0.09585322

x <- seq(0,9, length.out=30) # ages for prediction
CH_pred <- vbTyp(x, Linf=coef(fitCH)) #predicted lengths
xlmts <- range(c(x, Agelength_CH$final_age)) #set x limits
ylmts <- range(c(CH_pred, Agelength_CH$tl)) #set y limits
plot(tl~final_age, data=Agelength_CH, xlab="Age", ylab="Total Length (cm)", xlim=xlmts, ylim=ylmts, pch=19, col=rgb(0,0,0,1/3))
lines(CH_pred~x, lwd=2)

#Female
starting <- list(Linf=max(CH_F$tl, na.rm=TRUE), K=0.3, t0=-1)

fitCHF <- nls(tl~vbTyp(final_age, Linf, K, t0), data=CH_F, start=starting, control=list(maxiter=500))
coef(fitCHF)
bootCHF <- nlsBoot(fitCHF)
confint(bootCHF) #, plot=TRUE)

# Linf           K          t0 
# 44.50019994  1.31955946  0.08126805 
# 
# 95% LCI    95% UCI
# Linf 43.83742659 45.2196484
# K     1.22011888  1.4272520
# t0    0.05176346  0.1095023

#Male 
starting <- list(Linf=max(CH_M$tl, na.rm=TRUE), K=0.3, t0=-1)

fitCHM <- nls(tl~vbTyp(final_age, Linf, K, t0), data=CH_M, start=starting, control=list(maxiter=500))
coef(fitCHM)
bootCHM <- nlsBoot(fitCHM)
confint(bootCHM) #, plot=TRUE)

# Linf           K          t0 
# 37.20759294  1.55239495  0.07567399 
# 95% LCI    95% UCI
# Linf 36.64329893 37.8188731
# K     1.42197983  1.6806089
# t0    0.04837503  0.0977519


#IR####
IR_M <- droplevels(subset(Agelength_IR, sex == "M"))
IR_F <- droplevels(subset(Agelength_IR, sex == "F"))

starting <- list(Linf=max(Agelength_IR$tl, na.rm=TRUE), K=0.3, t0=-1)

fitIR <- nls(tl~vbTyp(final_age, Linf, K, t0), data=Agelength_IR, start=starting)
coef(fitIR)
bootIR <- nlsBoot(fitIR)
confint(bootIR) #, plot=TRUE)
#       Linf          K         t0 
#49.8323974  0.6002281 -0.2248435 

# 95% LCI    95% UCI
# Linf 48.8269549 51.1158977
# K     0.5469929  0.6519876
# t0   -0.3015149 -0.1589514

x <- seq(0,9, length.out=30) # ages for prediction
IR_pred <- vbTyp(x, Linf=coef(fitIR)) #predicted lengths
xlmts <- range(c(x, Agelength_IR$final_age)) #set x limits
ylmts <- range(c(IR_pred, Agelength_IR$tl)) #set y limits
plot(tl~final_age, data=Agelength_IR, xlab="Age", ylab="Total Length (cm)", xlim=xlmts, ylim=ylmts, pch=19, col=rgb(0,0,0,1/3))
lines(IR_pred~x, lwd=2)

#Female
starting <- list(Linf=max(IR_F$tl, na.rm=TRUE), K=0.3, t0=-1)

fitIRF <- nls(tl~vbTyp(final_age, Linf, K, t0), data=IR_F, start=starting, control=list(maxiter=500))
coef(fitIRF)
bootIRF <- nlsBoot(fitIRF)
confint(bootIRF) #, plot=TRUE)

# Linf          K         t0 
# 59.4635941  0.3912379 -0.5532394 
# 95% LCI    95% UCI
# Linf 57.1224790 62.2067355
# K     0.3445593  0.4421503
# t0   -0.6870460 -0.4325580

#Male 
starting <- list(Linf=max(IR_M$tl, na.rm=TRUE), K=0.3, t0=-1)

fitIRM <- nls(tl~vbTyp(final_age, Linf, K, t0), data=IR_M, start=starting, control=list(maxiter=500))
coef(fitIRM)
bootIRM <- nlsBoot(fitIRM)
confint(bootIRM) #, plot=TRUE)

# Linf           K          t0 
# 39.91591949  0.84003319 -0.05171365 
# 95% LCI      95% UCI
# Linf 39.0258429 40.985547632
# K     0.7501263  0.941251310
# t0   -0.1136163  0.007920836


#JX####

JX_M <- droplevels(subset(Agelength_JX, sex == "M"))
JX_F <- droplevels(subset(Agelength_JX, sex == "F"))

starting <- list(Linf=max(Agelength_JX$tl, na.rm=TRUE), K=0.3, t0=-1)

fitJX <- nls(tl~vbTyp(final_age, Linf, K, t0), data=Agelength_JX, start=starting)
coef(fitJX)
bootJX <- nlsBoot(fitJX)
confint(bootJX) #, plot=TRUE)
#Linf          K         t0 
#38.8827460  1.5061312  0.0805226 
# 95% LCI   95% UCI
# Linf 38.20193625 39.622566
# K     1.39953830  1.628236
# t0    0.05568952  0.100464

x <- seq(0,9, length.out=30) # ages for prediction
JX_pred <- vbTyp(x, Linf=coef(fitJX)) #predicted lengths
xlmts <- range(c(x, Agelength_JX$final_age)) #set x limits
ylmts <- range(c(JX_pred, Agelength_JX$tl)) #set y limits
plot(tl~final_age, data=Agelength_JX, xlab="Age", ylab="Total Length (cm)", xlim=xlmts, ylim=ylmts, pch=19, col=rgb(0,0,0,1/3))
lines(JX_pred~x, lwd=2)

#Female
starting <- list(Linf=max(JX_F$tl, na.rm=TRUE), K=0.3, t0=-1)

fitJXF <- nls(tl~vbTyp(final_age, Linf, K, t0), data=JX_F, start=starting, control=list(maxiter=500))
coef(fitJXF)
bootJXF <- nlsBoot(fitJXF)
confint(bootJXF) #, plot=TRUE)

# Linf           K          t0 
# 43.16580952  1.20432385  0.05958172 
# 
# 95% LCI     95% UCI
# Linf 41.95049616 44.53443303
# K     1.08234205  1.33971105
# t0    0.01947565  0.09311223

#Male 
fitJXM <- nls(tl~vbTyp(final_age, Linf, K, t0), data=JX_M, start=starting, control=list(maxiter=500))
coef(fitJXM)
bootJXM <- nlsBoot(fitJXM)
confint(bootJXM) #, plot=TRUE)

# Linf           K          t0 
# 35.46350778  1.73433291  0.08071643 
# 
# 95% LCI    95% UCI
# Linf 34.88894772 36.1296114
# K     1.59101979  1.8862930
# t0    0.05707483  0.1007352


####################################
#PLOT ALL PREDICTED FITS IN ONE PLOT ####
####################################

# combine all predicted matrices into one
x <- seq(0,9, length.out=30)
t= data.frame(cbind(x,TB_pred)) %>% rename(pred=TB_pred) %>% mutate(bay=rep("TB",30))
a= data.frame(cbind(x,AP_pred)) %>% rename(pred=AP_pred) %>% mutate(bay=rep("AP",30))
ck= data.frame(cbind(x,CK_pred)) %>% rename(pred=CK_pred) %>% mutate(bay=rep("CK",30))
c= data.frame(cbind(x,CH_pred)) %>% rename(pred=CH_pred) %>% mutate(bay=rep("CH",30))
j= data.frame(cbind(x,JX_pred))%>% rename(pred=JX_pred) %>% mutate(bay=rep("JX",30))
i= data.frame(cbind(x,IR_pred)) %>% rename(pred=IR_pred) %>% mutate(bay=rep("IR",30))


pred_all = rbind(t,a,ck,c,j,i)

library(grid)
library(ggplot2)

File <- ("U:/PhD_projectfiles/Figures/growth_curves.tiff")
if (file.exists(File)) stop(File, " already exists")
dir.create(dirname(File), showWarnings = FALSE)

tiff(File, units="in", width=7, height=5, res=300)

ggplot(data= pred_all, aes(x=x, y=pred, group=bay))+  
  geom_line() + #aes(linetype=bay), size =.5)+ # make line types based on the different labels- this will be our workaround because in a few stps we will specify the first label (obserseved) be a blank line (therefore a scatter plot)
  #geom_line(aes(color=bay), size=1.5)+
  geom_point(aes(shape=bay), size=2) + #, color=bay))+ # groups the points together correctly and then gives them a unique shape them differently based on the line type 
  #scale_shape_manual(values=c(0,1,2, 3, 4, 5,6,7,8 ))+
  #scale_linetype_manual(values=c('solid', 'dashed', 'dotted'))+
  #scale_y_continuous(limits=c(-10,55), breaks= seq(-10,55,5))+
  scale_x_continuous(limits=c(0,10), breaks=seq(0,10,1))+
  xlab("Age (yrs)")+
  ylab("Total Length (cm)")+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), 									
        panel.background=element_rect(fill='white', colour='black'),
        legend.key=element_blank(), legend.title=element_blank(),
        legend.background=element_rect(fill='white', size=.5),
        legend.position=c(.70,.38),
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
# vbLKt <- tl~Linf[bay]*(1-exp(-K[bay]*(final_age-t0[bay])))
# vbLK <- tl~Linf[bay]*(1-exp(-K[bay]*(final_age-t0)))
# vbLt <- tl~Linf[bay]*(1-exp(-K*(final_age-t0[bay])))
# vbKt <- tl~Linf*(1-exp(-K[bay]*(final_age-t0[bay])))
# vbL <-  tl~Linf[bay]*(1-exp(-K*(final_age-t0)))
# vbK <-  tl~Linf*(1-exp(-K[bay]*(final_age-t0)))
# vbt <- tl~Linf*(1-exp(-K*(final_age-t0[bay])))
# vb0 <- tl~Linf*(1-exp(-K*(final_age-t0)))


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
sv0 <- vbStarts(tl~final_age, data=Agelength_CK) #using the TB dataframe as a starter
sv0$Linf= 70
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
#Define fits to model scenarios- first by sex for each bay 
############################################################## 
#filter missing sex
 library(minpack.lm)
  
tb_mf <- droplevels(subset(Agelength_TB, sex %in% c("M", "F")))
#Fit to [LKt] model then to the simplest model [vb0]

fitLKt <- nls(vbLKt, data=tb_mf, start=svLKt)
#residPlot(fitLKt, col=rgb(0,0,0,1/3))
fitLK <- nls(vbLK, data=tb_mf, start=svLK)
fitLt <- nls(vbLt, data=tb_mf, start=svLt)
fitKt <- nls(vbKt, data=tb_mf, start=svKt)
fitL <- nls(vbL, data=tb_mf, start=svL)
fitt <- nls(vbt, data=tb_mf, start=svt)
fitK <- nls(vbK, data=tb_mf, start=svK)
fit0 <- nls(vb0, data=tb_mf, start=sv0)

cbind(AIC(fitLKt,fitLK,fitLt,fitKt, fitL, fitK, fitt, fit0), BIC(fitLKt, fitLK, fitLt, fitKt, fitL, fitK, fitt, fit0))

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

ch_mf <- droplevels(subset(Agelength_CH, sex %in% c("M", "F")))
#Fit to [LKt] model then to the simplest model [vb0]
fitLKt <- nls(vbLKt, data=ch_mf, start=svLKt)
#residPlot(fitLKt, col=rgb(0,0,0,1/3))
fitLK <- nls(vbLK, data=ch_mf, start=svLK)
fitLt <- nls(vbLt, data=ch_mf, start=svLt)
fitKt <- nls(vbKt, data=ch_mf, start=svKt)
fitL <- nls(vbL, data=ch_mf, start=svL)
fitt <- nls(vbt, data=ch_mf, start=svt)
fitK <- nls(vbK, data=ch_mf, start=svK)
fit0 <- nls(vb0, data=ch_mf, start=sv0)

cbind(AIC(fitLKt,fitLK,fitLt,fitKt, fitL, fitK, fitt, fit0), BIC(fitLKt, fitLK, fitLt, fitKt, fitL, fitK, fitt, fit0))

ck_mf <- droplevels(subset(Agelength_CK, sex %in% c("M", "F")))
#Fit to [LKt] model then to the simplest model [vb0]
fitLKt <- nls(vbLKt, data=ck_mf, start=svLKt)
#residPlot(fitLKt, col=rgb(0,0,0,1/3))
fitLK <- nls(vbLK, data=ck_mf, start=svLK)
fitLt <- nls(vbLt, data=ck_mf, start=svLt)
fitKt <- nls(vbKt, data=ck_mf, start=svKt)
fitL <- nls(vbL, data=ck_mf, start=svL)
fitt <- nls(vbt, data=ck_mf, start=svt)
fitK <- nls(vbK, data=ck_mf, start=svK)
fit0 <- nls(vb0, data=ck_mf, start=sv0)

cbind(AIC(fitLKt,fitLK,fitLt,fitKt, fitL, fitK, fitt, fit0), BIC(fitLKt, fitLK, fitLt, fitKt, fitL, fitK, fitt, fit0))

ap_mf <- droplevels(subset(Agelength_AP, sex %in% c("M", "F")))
#Fit to [LKt] model then to the simplest model [vb0]
fitLKt <- nls(vbLKt, data=ap_mf, start=svLKt)
#residPlot(fitLKt, col=rgb(0,0,0,1/3))
fitLK <- nls(vbLK, data=ap_mf, start=svLK)
fitLt <- nls(vbLt, data=ap_mf, start=svLt)
fitKt <- nls(vbKt, data=ap_mf, start=svKt)
fitL <- nls(vbL, data=ap_mf, start=svL)
fitt <- nls(vbt, data=ap_mf, start=svt)
fitK <- nls(vbK, data=ap_mf, start=svK)
fit0 <- nls(vb0, data=ap_mf, start=sv0)

cbind(AIC(fitLKt,fitLK,fitLt,fitKt, fitL, fitK, fitt, fit0), BIC(fitLKt, fitLK, fitLt, fitKt, fitL, fitK, fitt, fit0))

jx_mf <- droplevels(subset(Agelength_JX, sex %in% c("M", "F")))
#Fit to [LKt] model then to the simplest model [vb0]
fitLKt <- nls(vbLKt, data=jx_mf, start=svLKt)
#residPlot(fitLKt, col=rgb(0,0,0,1/3))
fitLK <- nls(vbLK, data=jx_mf, start=svLK)
fitLt <- nls(vbLt, data=jx_mf, start=svLt)
fitKt <- nls(vbKt, data=jx_mf, start=svKt)
fitL <- nls(vbL, data=jx_mf, start=svL)
fitt <- nls(vbt, data=jx_mf, start=svt)
fitK <- nls(vbK, data=jx_mf, start=svK)
fit0 <- nls(vb0, data=jx_mf, start=sv0)

cbind(AIC(fitLKt,fitLK,fitLt,fitKt, fitL, fitK, fitt, fit0), BIC(fitLKt, fitLK, fitLt, fitKt, fitL, fitK, fitt, fit0))

ir_mf <- droplevels(subset(Agelength_IR, sex %in% c("M", "F")))
#Fit to [LKt] model then to the simplest model [vb0]
fitLKt <- nls(vbLKt, data=ir_mf, start=svLKt)
#residPlot(fitLKt, col=rgb(0,0,0,1/3))
fitLK <- nls(vbLK, data=ir_mf, start=svLK)
fitLt <- nls(vbLt, data=ir_mf, start=svLt)
fitKt <- nls(vbKt, data=ir_mf, start=svKt)
fitL <- nls(vbL, data=ir_mf, start=svL)
fitt <- nls(vbt, data=ir_mf, start=svt)
fitK <- nls(vbK, data=ir_mf, start=svK)
fit0 <- nls(vb0, data=ir_mf, start=sv0)

cbind(AIC(fitLKt,fitLK,fitLt,fitKt, fitL, fitK, fitt, fit0), BIC(fitLKt, fitLK, fitLt, fitKt, fitL, fitK, fitt, fit0))


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
# fitLK   6 22568.78  6 22605.44


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
# fitL    5 21809.18  5 21839.78
#fitK    5 21808.97  5 21839.57
#fitt    5 21808.85  5 21839.45
#fit0    4 21808.25  4 21832.73

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
# fitLKt  7 22566.44  7 22609.48
# fitLK   6 22566.46  6 22603.35


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
# fitLKt  7 19730.15  7 19772.32
# fitLK   6 19729.45  6 19765.60

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
# fitLKt  7 33832.90  7 33878.47

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
# fitLK   6 15209.74  6 15244.38

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
# fitLKt  7 15977.29  7 16017.97


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
# fitLK   6 13116.66  6 13150.46


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
# fitLKt  7 27244.87  7 27288.89

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
#fitLKt  7 27244.87  7 27288.89


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
# fitLK   6 14919.18  6 14953.75


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
# fitLKt  7 29066.00  7 29110.49

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
# fitLKt  7 15694.29  7 15734.89

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
# fitLKt  7 29782.81  7 29827.46

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



