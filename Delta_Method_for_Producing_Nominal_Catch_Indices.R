###### ABOUT ##############
# 10/26/2016
# Purpose: To produce adjusted nominal catch indices using the Delta Method
# 1. Use Delta Method to produce YOY indices using predicted positive and proportion positive data sets
# 2. Use Delta Method to produce adult indices using predicted positive and proportion positive data sets
# 3. Fit SR relationships to yoy and adult indices to produce the residuals of Beverton Holt and Ricker
#test
#test4_26
#test 7_30
#test 7_30
#this is a test from work
#newchange
#test

#Delta Method
# This method uses the Delta method for determining the nominal catch rate of Spotted Seatrout in the FIM catch
# Uses methods described in Chyan-huei et al 1992
# Essentially, you can break up the data set into positive hauls and zero hauls where

# positive haul= hauls with Seatrout
# zero haul= hauls without any Seatrout 

# Positive Data = total numbers of Seatrout in a year/ total numbers of positive hauls in year
# Proportion Positive = total number of positive hauls / total numbers of all hauls (positive and zero) 

# Adjusted Index = Positive Data * Proportion Positive

#4/19/17
# Meeting at FWRI suggested that I need to control at least for habitat variables in my estimates of positive and proportion positive
#The FWRI R script DeltaLogNormal_SST_yoy_FIM.R does this for a portion so it is a good guide to use. 

# I will use the Positive Data as input data to a general linear model and habitat grouping variables
# For example     log(positive numbers)= B + B*bottom_type + B*bottom_veg + B*shore_type + B*zone
# This will produce a predicted number of positive hauls.

# Then I will use the Binomial Data set (used to produce the proportion positive) as input data to a general linear model.
# This will produce the predicted number of binary data, from which I can produce the predicted proportion postive. 

# Then I can create the Adjusted Predicted Index = Predicted Positive Data * Predicted Proportion Positive Data

##### SET PACKAGES ######################
###########################
library(haven) #to load sas
library(dplyr) # to do df manipulation
library(lsmeans) #to determine the least squares means
library(AER) #to test for overdispersion
library(propagate) # for propagating error when multiplying predicted positive by predicted proportion positive
##### IMPORT DATA SETS_YOY ######
###########################
# These data sets were produced using the spp_comb_5_13_EG_2bays_yoy_2015_EHedits.sas program which is stored in my scratch folder
# Bay and river stations were denoted by the gear sampling  code. There is a hard copy of this script in my green folder. 
#11/7/16- Realized that the data sets that were produced by the above sas program that I edited was pretty crummy in that
#         a lot of YOY animals were being excluded by my gear selections. Therefore I went back to FWRI and was more general with my gear selections.
#         Chose gears 19,20,22,23 following what gets selected by the spp_comb_5_13_EG_2bays_yoy_2015.sas program 

# select the important recruitment months for each zone and also check on gear codes
# adjust selected gear if it was not chosen in the first place
# _C$month => depends on recruitment window in each estuary
#               => Jax 5<=x<=11
#               => nor. IRL 5<=x<=11
#               => CK  5<=x<=11
#               => TB  4<=x<=10
#               => CH  4<=x<=10
#               => AP  6<=x<=10

#Changed location of saved file on 4/19/17
# Created a new GitHub repository for just my scripts in the Seatrout project
# New repository is called Seatrout_Scripts

#must change working directory for data when working on personal vs work computer
#setwd("~/Desktop/PhD project/Projects/Seatrout/FWRI SCRATCH FOLDER/Elizabeth Herdter/SAS data sets/FIMData/NEWNov7")
#setwd("T:/Elizabeth Herdter/SAS data sets/FIMData/NEWNov7")

<<<<<<< HEAD
=======

>>>>>>> 0b8438883a655df6a9eab7117dca9381163b292a

#load the data, select the peak reproductive months, reorder columns alphabetically so I can combine dataframes (some columns were in different position in other df)
ap = subset(read_sas("ap_yoy_cn_c.sas7bdat"), month %in% c(6,7,8,9,10,11)) %>% mutate(bUnk=bunk) %>% select(-bunk) 
ap <- ap %>% select(noquote(order(colnames(ap))))  #reorders the columns alphabetically 
apl = read_sas("ap_yoy_cn_l.sas7bdat")

#merge with length data and make sure the approporiate lengths are included 
t <- merge(ap, apl, by="bio_reference")

ck = subset(read_sas("ck_yoy_cn_c.sas7bdat"),  month %in% c(5,6,7,8,9,10,11))
ck <- ck %>% select(noquote(order(colnames(ck))))
ckl = read_sas("ck_yoy_cn_l.sas7bdat")

ch = subset(read_sas("ch_yoy_cn_c.sas7bdat"), month %in% c(4,5,6,7,8,9,10)) %>% mutate(bUnk=bunk) %>% select(-bunk) 
ch <- ch %>% select(noquote(order(colnames(ch))))
chl = read_sas("ch_yoy_cn_l.sas7bdat")

tb = subset(read_sas("tb_yoy_cn_c.sas7bdat"), month %in% c(4,5,6,7,8,9,10)) 
tb <- tb %>% select(noquote(order(colnames(tb))))
tbl = read_sas("tb_yoy_cn_l.sas7bdat")

ir = subset(read_sas("ir_yoy_cn_c.sas7bdat"), month %in% c(5,6,7,8,9,10,11)) 
ir <- ir %>% select(noquote(order(colnames(ir))))
irl = read_sas("ir_yoy_cn_l.sas7bdat")

jx = subset(read_sas("jx_yoy_cn_c.sas7bdat") , month %in% c(5,6,7,8,9,10,11)) 
jx <- jx %>% select(noquote(order(colnames(jx))))
jxl = read_sas("jx_yoy_cn_l.sas7bdat")

#also check on zones to make sure the bays and rivers are stratefying correctly => depends on each estuary (for a more in depth description of monthly sampling within bay or riverine see Annual Report)
#               => TB BAY(A-E), RIVERINE (K-N)
#               => CH BAY(A-D), RIVERINE (M-P)
#               => nor.IRL BAY(A-E, H), RIVERINE (F)
#               => AP BAY(A-B), RIVERINE (C)
#               => Jax RIVERINE (A-F)
#               => CK BAY(B-C), RIVERINE (F)

#Join all together so I can perform following steps on a single dataframe instead of individual ones. 

full <- rbind(ap,ck,tb,ch,jx,ir)

##### SELECT CATEGORICAL HABITAT VALUES_YOY ########
# to be used in the model to predict positive numbers
################################################

#Based on FWRI code the three variables were bottom type (bStr, bsan, bmud), bottom vegetation (bveg), and shoreline (Shore)
#There are three different bottom type variables each of them coded in a binary form.
#I want to take bStr, bsan, and bmud and put them into 1 variable so I will make a new variable entirely = 'bottom'
#I also want to turn bveg into a new variable = 'veg' based on the entries. If alg or Sav then turn to SAV because there are only 9 entries for Alg. 
#Same thing for the shore variable = 'shore'. Decided to have only emergent, structure, terrestrial, and mangrove. 
#Removed old variables (bStr, bSan, bMud, bveg, Shore)
#Removed rows when there was no shoreline variable. 

full <-  select(full, c(number,year, month, bio_reference, bStr, bSan, bMud, bveg, Shore, bay, Zone)) %>% 
        mutate(bottom = ifelse(full$bStr ==1, "structure", ifelse(full$bSan>0 | full$bMud>0, "mudsand", "unknown")), 
               veg= ifelse(full$bveg == "SAVAlg", "SAV", ifelse(full$bveg == "Alg", "SAV", ifelse(full$bveg =="SAV", "SAV", "Noveg"))),
               shore = ifelse(substr(full$Shore,1,3)=="Eme", "Emerge", ifelse(substr(full$Shore,1,3) =="Man", "Mangrove", ifelse(substr(full$Shore,1,3)=="Str", "Structure", 
                       ifelse(substr(full$Shore, 1,3)=="Ter", "Terrestrial", "Non"))))) %>% select(-c(bStr, bSan, bMud, bveg, Shore)) %>% subset(!shore=="Non") 
      
#Turn habitat variables into factors so they can be treated as categorical
full[,c(2,5:9)] <- lapply(full[,c(2,5:9)], factor)

###### MAKE POSITIVE & BINARY SET_YOY ##########
# to determine the total number of positive huals and proportion positive
##############################################

full.pos<- full %>% subset(number>0)
full.bin <- full %>% mutate(number=ifelse(number>0,1,0))

ap.pos <- full.pos %>% subset(bay =='AP')
ck.pos <- full.pos %>% subset(bay =='CK')
tb.pos <- full.pos %>% subset(bay =='TB')
ch.pos <- full.pos %>% subset(bay =='CH')
jx.pos <- full.pos %>% subset(bay =='JX')
ir.pos <- full.pos %>% subset(bay =='IR')

ap.bin <- full.bin %>% subset(bay =='AP')
ck.bin <- full.bin %>% subset(bay =='CK')
tb.bin <- full.bin %>% subset(bay =='TB')
ch.bin <- full.bin %>% subset(bay =='CH')
jx.bin <- full.bin %>% subset(bay =='JX')
ir.bin <- full.bin %>% subset(bay =='IR')

##### VISUALIZE THE DATA_YOY ########
################################

#Plot the data

#AP
plot(ap.pos$bottom, ap.pos$number, xlab= "bottom type", ylab="number")
plot(ap.pos$year, ap.pos$number, vlab="year", ylab="number")
plot(ap.pos$shore, ap.pos$number, vlab="shore", ylab="number")
plot(ap.pos$veg, ap.pos$number, vlab="veg", ylab="number")

### BUILD MODELS_YOY #########
# To produce predicted positive numbers and binomial data set
#######################

# 1. Build the full models with all potential variables and a base model with only year as a variable. 
#    Do this for both the positive (Poisson distribution) and binary (Binomial distribution) datasets. 
# 2. Check for overdispersion in the Poisson distribution scenario. 
# 3. If there is overdispersion use quasipoisson 

# 1. Build the full models for the positive and binomial datasets.  
#AP
Full_ap.pos <- glm(number ~ year+month+bottom+veg+shore, data=ap.pos, family=poisson)
Full_ap.bin <- glm(number ~ year+month+bottom+veg+shore, data=ap.bin, family=binomial)

#CK
Full_ck.pos <- glm(number ~ year+month+bottom+veg+shore, data=ck.pos, family=poisson)
Full_ck.bin <- glm(number ~ year+month+bottom+veg+shore, data=ck.bin, family=binomial)

#TB
Full_tb.pos <- glm(number ~ year+month+bottom+veg+shore, data=tb.pos, family=poisson)
Full_tb.bin <- glm(number ~ year+month+bottom+veg+shore, data=tb.bin, family=binomial)

#CH
Full_ch.pos <- glm(number ~ year+month+bottom+veg+shore, data=ch.pos, family=poisson)
Full_ch.bin <- glm(number ~ year+month+bottom+veg+shore, data=ch.bin, family=binomial)

#JX
Full_jx.pos <- glm(number ~ year+month+bottom+veg+shore, data=jx.pos, family=poisson)
Full_jx.bin <- glm(number ~ year+month+bottom+veg+shore, data=jx.bin, family=binomial)

#IR
Full_ir.pos <- glm(number ~ year+month+bottom+veg+shore, data=ir.pos, family=poisson)
Full_ir.bin <- glm(number ~ year+month+bottom+veg+shore, data=ir.bin, family=binomial)

#2. Test the Poisson GLMs for overdispersion
# With the Bernoulli GLM (binomial, response variable is a vector of zeros and ones) overdispersion does not ever occur (Zuur og 253) so I don't need to test for overdispersion in the .bin models. 
dispersiontest(Full_ap.pos, trafo=1)
dispersiontest(Full_ck.pos, trafo=1)
dispersiontest(Full_tb.pos, trafo=1)
dispersiontest(Full_ch.pos, trafo=1)
dispersiontest(Full_jx.pos, trafo=1)
dispersiontest(Full_ir.pos, trafo=1)

# 3. there is evidence of overdispersion for every bay (p values were less than) so use quasipoisson for Positive models (Zuur pg 226)
Full_ap.pos <- glm(number ~ year +month+veg+bottom+shore, data=ap.pos, family=quasipoisson)
Full_ck.pos <- glm(number ~ year +month+veg+bottom+shore, data=ck.pos, family=quasipoisson)
Full_tb.pos <- glm(number ~ year +month+veg+bottom+shore, data=tb.pos, family=quasipoisson)
Full_ch.pos <- glm(number ~ year +month+veg+bottom+shore, data=ch.pos, family=quasipoisson)
Full_jx.pos <- glm(number ~ year +month+veg+bottom+shore, data=jx.pos, family=quasipoisson)
Full_ir.pos <- glm(number ~ year +month+veg+bottom+shore, data=ir.pos, family=quasipoisson)

## MODEL SELECTION POSITIVE w/ DROP1 command_YOY ######
################################
# Pages 220 is to 230 in Zuur are helpful for following methods. 
# The AIC is not defined for quasipoisson models so can't use the step function like what was used in the FWRI code. 
# Instead, use the drop1 function which is applicable for the quassiPoisson GLM and it is more equivalent to hypothesis testing. (Zuur pg 227)
# If just using Poisson or Bernoulli (binomial) can use step command but this gives AIC- not deviance. (Zuur pg 253)
# Explained deviance is nearly the equivalent of R^2 so use this (Zuur pg 218 for equation), "The smaller the residual deviance the better is the model"

###AP_POS (Year, Veg, Shore = significant factors)
summary(Full_ap.pos)
drop1(Full_ap.pos, test="F")  #model selection in quasipoisson is done using F-ratio (Zuur pg 227)
# bottom and month do not appear significant. Drop all sequentially. 

# drop month
M1_ap.pos <- glm(number ~ year+veg+bottom+shore, data=ap.pos, family=quasipoisson)
drop1(M1_ap.pos, test="F")

# drop month, bottom
M2_ap.pos <- glm(number ~ year+veg+shore, data=ap.pos, family=quasipoisson)
drop1(M2_ap.pos, test="F")
# Year, Veg, Shore = significant factors

### CK_POS (Year, Veg = significant factor)
summary(Full_ck.pos)
drop1(Full_ck.pos, test="F")
#month, bottom, and shore do not appear significant. Drop all sequentially. 

#drop month
M1_ck.pos <- glm(number ~ year+veg+bottom+shore, data=ck.pos, family=quasipoisson)
drop1(M1_ck.pos, test="F")

#drop month, bottom
M2_ck.pos <- glm(number ~ year+veg+shore, data=ck.pos, family=quasipoisson)
drop1(M2_ck.pos, test="F")

#drop month, bottom, and shore
M3_ck.pos <- glm(number ~ year+veg, data=ck.pos, family=quasipoisson)
drop1(M3_ck.pos, test="F")
# Year, Veg = significant factors 

### TB_POS (Year, Veg, Shore = significant factors)
summary(Full_tb.pos)
drop1(Full_tb.pos, test="F")
# bottom and month do not appear significant. Drop all sequentially

# drop month
M1_tb.pos <- glm(number ~ year+veg+bottom+shore, data=tb.pos, family=quasipoisson)
drop1(M1_tb.pos, test="F")

# drop month, bottom
M2_tb.pos <- glm(number ~ year+veg+shore, data=tb.pos, family=quasipoisson)
drop1(M2_tb.pos, test="F")
# Year, Veg, Shore = significant factors

### CH_POS (Year, Month, Veg, Bottom, Shore =significnat factors )
summary(Full_ch.pos)
drop1(Full_ch.pos, test="F")
# Year, Month, Veg, Bottom, Shore =significnat factors 

### JX_POS (Year, Veg, Shore = significant factors)
summary(Full_jx.pos)
drop1(Full_jx.pos, test="F")
# bottom and month do not appear significant. Drop all sequentially

#drop month
M1_jx.pos <- glm(number ~ year+veg+bottom+shore, data=jx.pos, family=quasipoisson)
drop1(M1_jx.pos, test="F")

#drop month, bottom
M2_jx.pos <- glm(number ~ year+veg+shore, data=jx.pos, family=quasipoisson)
drop1(M2_jx.pos, test="F")
# Year, Veg, Shore = significant factors

### IR_POS (Year, Veg, Bottom = significant factors)
summary(Full_ir.pos)
drop1(Full_ir.pos, test="F")
#month and shore do not appear significant. Drop them all. 

# drop month
M1_ir.pos <- glm(number ~ year+veg+bottom+shore, data=ir.pos, family=quasipoisson)
drop1(M1_ir.pos, test="F")

# drop shore
M2_ir.pos <- glm(number ~ year+veg+bottom, data=ir.pos, family=quasipoisson)
drop1(M2_ir.pos, test="F")
# Year, Veg, Bottom = significant factors

## MODEL SELECTION BINARY w/ DROP1 command_YOY ######
################################
# pg 253 Zuur
##  AP_BIN (Year, Veg = significant)
summary(Full_ap.bin)
drop1(Full_ap.bin, test ='Chi')
# month, bottom, and shore are not significant so drop them one at a time

#drop month
M1_ap.bin <- glm(number ~ year+bottom+veg+shore, data=ap.bin, family=binomial)
drop1(M1_ap.bin, test ="Chi")

#drop month, bottom
M2_ap.bin <- glm(number ~ year+veg+shore, data=ap.bin, family=binomial)
drop1(M2_ap.bin, test ="Chi")

#drop month, bottom, shore
M3_ap.bin <- glm(number ~ year+veg, data=ap.bin, family=binomial)
drop1(M3_ap.bin, test ="Chi")
# Year, Veg = significant

## CK_BIN (Year, Month, Veg, Shore = significant)
summary(Full_ck.bin)
drop1(Full_ck.bin, test ='Chi')
# bottom is not signficiant

#drop bottom 
M1_ck.bin <- glm(number ~ year+month+veg+shore, data=ck.bin, family=binomial)
drop1(M1_ck.bin, test ="Chi")


## TB_BIN (Year, Month, Veg, Shore = significant)
summary(Full_tb.bin)
drop1(Full_tb.bin, test ='Chi')
# bottom is not signficiant

#drop bottom 
M1_tb.bin <- glm(number ~ year+month+veg+shore, data=tb.bin, family=binomial)
drop1(M1_tb.bin, test ="Chi")

##  CH_BIN (Year, Month, Veg, Shore = significant)
summary(Full_ch.bin)
drop1(Full_ch.bin, test ='Chi')
# bottom is not signficiant

#drop bottom 
M1_ch.bin <- glm(number ~ year+month+veg+shore, data=ch.bin, family=binomial)
drop1(M1_ch.bin, test ="Chi")

## JX_BIN (Year, Bottom, Veg, Shore =signficant)
summary(Full_jx.bin)
drop1(Full_jx.bin, test ='Chi')
# month is not signficiant

#drop month 
M1_jx.bin <- glm(number ~ year+bottom+veg+shore, data=jx.bin, family=binomial)
drop1(M1_jx.bin, test ="Chi")
# Year, Bottom, Veg, Shore = significant

## IR_BIN (Year, Month, Veg, Shore = significant)
summary(Full_ir.bin)
drop1(Full_ir.bin, test ='Chi')
# bottom is not signficiant

#drop bottom 
M1_ir.bin <- glm(number ~ year+month+veg+shore, data=ir.bin, family=binomial)
drop1(M1_ir.bin, test ="Chi")

### ASSIGN FINAL MODELS_YOY ###### 
###############################
final_ap.pos = M2_ap.pos 
final_ck.pos = M3_ck.pos
final_tb.pos = M2_tb.pos
final_ch.pos= Full_ch.pos
final_jx.pos = M2_jx.pos
final_ir.pos = M2_ir.pos

final_ap.bin = M3_ap.bin
final_ck.bin = M1_ck.bin
final_tb.bin = M1_tb.bin
final_ch.bin = M1_ch.bin
final_jx.bin = M1_jx.bin
final_ir.bin = M1_ir.bin

### DETERMINE LEAST SQUARE MEANS_YOY ###########
# DETERMINE LEAST SQUARE MEANS 
#################################
# Same thing as covariate adjusted means. Basically, determine the mean value of total positive numbers 
# of catch per year controlling for covariates (in this case it would be veg and shore variables). 
# Use lsmeans CRAN document. 

# Looking at the reference grid gives a good idea of over what levels the mean is being averaged. 
ap.rf.grid <- ref.grid(final_ap.pos)
ap.bin.rf.grid <- ref.grid(final_ap.bin)

#Can make predictions using the reference grid. It produces a mean value of numbers based on each scenario combination.  
test_ap.pos = summary(ap.rf.grid)
test_ap.bin= summary(ap.bin.rf.grid)

# Use lsmeans to determine the least square mean of positive values. 
# Display the response scale (as opposed to the log scale which is reported for the Poisson, and the logit scale reported for the Binomial)

#POSITIVE
test<- summary(lsmeans(final_ap.pos, 'year', data=ap.pos))
LSM_ap.pos <- summary(lsmeans(final_ap.pos, 'year', data=ap.pos), type="response")
LSM_ck.pos <- summary(lsmeans(final_ck.pos, 'year', data=ck.pos), type="response")
LSM_tb.pos <- summary(lsmeans(final_tb.pos, 'year', data=tb.pos), type="response")
LSM_ch.pos <- summary(lsmeans(final_ch.pos, 'year', data=ch.pos), type="response")
LSM_jx.pos <- summary(lsmeans(final_jx.pos, 'year', data=jx.pos), type="response")
LSM_ir.pos <- summary(lsmeans(final_ir.pos, 'year', data=ir.pos), type="response")

#BINOMIAL (with type="response" this is equal to proportion positive becuase its like a percentage- proportion of 0s to 1s)
LSM_ap.bin <- summary(lsmeans(final_ap.bin, 'year', data=ap.bin), type="response")
LSM_ck.bin <- summary(lsmeans(final_ck.bin, 'year', data=ck.bin), type="response")
LSM_tb.bin <- summary(lsmeans(final_tb.bin, 'year', data=tb.bin), type="response")
LSM_ch.bin <- summary(lsmeans(final_ch.bin, 'year', data=ch.bin), type="response")
LSM_jx.bin <- summary(lsmeans(final_jx.bin, 'year', data=jx.bin), type="response")
LSM_ir.bin <- summary(lsmeans(final_ir.bin, 'year', data=ir.bin), type="response")

### ERROR PROPAGATION TO FIND FINAL VALUE (pos * prop.pos)_YOY #####
# multiply positive lsmean by porportion positive lsmean and use error propagation to determine value and associated error
# using the package Propagate. See example below. 
# https://www.rdocumentation.org/packages/propagate/versions/1.0-4/topics/propagate    

#Must use a for loop to do the error propagation because it goes one row at a time without the loop and its very cumbersome 
#error propagation steps start with the EXPR command where you tell it what the expression is going to be. 
# The the expression setup gets used within the propagation step below with the actual dataframe (DF)


#AP
num.yr = length(LSM_ap.pos$year)  
df <- data.frame(matrix(data=NA, nrow=num.yr, ncol=6)) #make a dataframe for the loop to store results in
for (i in 1:num.yr) {
  x = c(LSM_ap.pos$rate[1], LSM_ap.pos$SE[1])
  y= c(LSM_ap.bin$prob[1], LSM_ap.bin$SE[1])
  EXPR <- expression(x*y)
  DF <- cbind(x,y)
  RES <- propagate(expr=EXPR, data=DF, type='stat', do.sim=TRUE, verbose=TRUE)
  df[i,] <- t(matrix(RES$sim))
}
Mean_AP <- df %>% cbind(LSM_ap.pos$year)
colnames(Mean_AP) <- c("Mean", "SD", "Median", "MAD", "2.5%", "97.5%", "Year")
write.csv(Mean_AP, "~/Desktop/PhD project/Projects/Seatrout/Data/Indices/DeltaMethod Indices/AP_yoy_index.csv")

#CK
num.yr = length(LSM_ck.pos$year)  
df <- data.frame(matrix(data=NA, nrow=num.yr, ncol=6)) #make a dataframe for the loop to store results in
for (i in 1:num.yr) {
  x = c(LSM_ck.pos$rate[i], LSM_ck.pos$SE[i])
  y= c(LSM_ck.bin$prob[i], LSM_ck.bin$SE[i])
  EXPR <- expression(x*y)
  DF <- cbind(x,y)
  RES <- propagate(expr=EXPR, data=DF, type='stat', do.sim=TRUE, verbose=TRUE)
  df[i,] <- t(matrix(RES$sim))
}
Mean_CK <- df %>% cbind(LSM_ck.pos$year)
colnames(Mean_CK) <- c("Mean", "SD", "Median", "MAD", "2.5%", "97.5%", "Year")
write.csv(Mean_CK, "~/Desktop/PhD project/Projects/Seatrout/Data/Indices/DeltaMethod Indices/CK_yoy_index.csv")

#TB
num.yr = length(LSM_tb.pos$year)  
df <- data.frame(matrix(data=NA, nrow=num.yr, ncol=6)) #make a dataframe for the loop to store results in
for (i in 1:num.yr) {
  x = c(LSM_tb.pos$rate[i], LSM_tb.pos$SE[i])
  y= c(LSM_tb.bin$prob[i], LSM_tb.bin$SE[i])
  EXPR <- expression(x*y)
  DF <- cbind(x,y)
  RES <- propagate(expr=EXPR, data=DF, type='stat', do.sim=TRUE, verbose=TRUE)
  df[i,] <- t(matrix(RES$sim))
}
Mean_TB <- df %>% cbind(LSM_tb.pos$year)
colnames(Mean_TB) <- c("Mean", "SD", "Median", "MAD", "2.5%", "97.5%", "Year")
write.csv(Mean_TB, "~/Desktop/PhD project/Projects/Seatrout/Data/Indices/DeltaMethod Indices/TB_yoy_index.csv")

#CH
num.yr = length(LSM_ch.pos$year)  
df <- data.frame(matrix(data=NA, nrow=num.yr, ncol=6)) #make a dataframe for the loop to store results in
for (i in 1:num.yr) {
  x = c(LSM_ch.pos$rate[i], LSM_ch.pos$SE[i])
  y= c(LSM_ch.bin$prob[i], LSM_ch.bin$SE[i])
  EXPR <- expression(x*y)
  DF <- cbind(x,y)
  RES <- propagate(expr=EXPR, data=DF, type='stat', do.sim=TRUE, verbose=TRUE)
  df[i,] <- t(matrix(RES$sim))
}
Mean_CH <- df %>% cbind(LSM_ch.pos$year)
colnames(Mean_CH) <- c("Mean", "SD", "Median", "MAD", "2.5%", "97.5%", "Year")
write.csv(Mean_CH, "~/Desktop/PhD project/Projects/Seatrout/Data/Indices/DeltaMethod Indices/CH_yoy_index.csv")

#JX
num.yr = length(LSM_jx.pos$year)  
df <- data.frame(matrix(data=NA, nrow=num.yr, ncol=6)) #make a dataframe for the loop to store results in
for (i in 1:num.yr) {
  x = c(LSM_jx.pos$rate[i], LSM_jx.pos$SE[i])
  y= c(LSM_jx.bin$prob[i], LSM_jx.bin$SE[i])
  EXPR <- expression(x*y)
  DF <- cbind(x,y)
  RES <- propagate(expr=EXPR, data=DF, type='stat', do.sim=TRUE, verbose=TRUE)
  df[i,] <- t(matrix(RES$sim))
}
Mean_JX <- df %>% cbind(LSM_jx.pos$year)
colnames(Mean_JX) <- c("Mean", "SD", "Median", "MAD", "2.5%", "97.5%", "Year")
write.csv(Mean_JX, "~/Desktop/PhD project/Projects/Seatrout/Data/Indices/DeltaMethod Indices/JX_yoy_index.csv")

#IR
num.yr = length(LSM_ir.pos$year)  
df <- data.frame(matrix(data=NA, nrow=num.yr, ncol=6)) #make a dataframe for the loop to store results in
for (i in 1:num.yr) {
  x = c(LSM_ir.pos$rate[i], LSM_ir.pos$SE[i])
  y= c(LSM_ir.bin$prob[i], LSM_ir.bin$SE[i])
  EXPR <- expression(x*y)
  DF <- cbind(x,y)
  RES <- propagate(expr=EXPR, data=DF, type='stat', do.sim=TRUE, verbose=TRUE)
  df[i,] <- t(matrix(RES$sim))
}
Mean_IR <- df %>% cbind(LSM_ir.pos$year)
colnames(Mean_IR) <- c("Mean", "SD", "Median", "MAD", "2.5%", "97.5%", "Year")
write.csv(Mean_IR, "~/Desktop/PhD project/Projects/Seatrout/Data/Indices/DeltaMethod Indices/IR_yoy_index.csv")


##LOAD ADULT DATA ########
################################################

setwd("~/Desktop/PhD project/Projects/Seatrout/FWRI SCRATCH FOLDER/Elizabeth Herdter/SAS data sets/FIMData")

ap_ad = subset(read_sas("ap_adult_cn_c.sas7bdat"))%>% mutate(bUnk=bunk) %>% select(-bunk) 
ap_ad <- ap_ad %>% select(noquote(order(colnames(ap_ad))))

ap_adl = subset(read_sas("ap_adult_cn_l.sas7bdat"))

ch_ad = subset(read_sas("ch_adult_cn_c.sas7bdat")) %>% mutate(bUnk=bunk) %>% select(-bunk) 
ch_ad <- ch_ad %>% select(noquote(order(colnames(ch_ad))))

ch_adl = subset(read_sas("ch_adult_cn_l.sas7bdat"))

ck_ad = subset(read_sas("ck_adult_cn_c.sas7bdat"))%>% mutate(bUnk=bunk) %>% select(-bunk) 
ck_ad <- ck_ad %>% select(noquote(order(colnames(ck_ad))))

ck_adl = subset(read_sas("ck_adult_cn_l.sas7bdat"))

tb_ad = subset(read_sas("tb_adult_cn_c.sas7bdat"))
tb_ad <- tb_ad %>% select(noquote(order(colnames(tb_ad))))

tb_adl = subset(read_sas("tb_adult_cn_l.sas7bdat"))

jx_ad = subset(read_sas("jx_adult_cn_c.sas7bdat"))
jx_ad <- jx_ad %>% select(noquote(order(colnames(jx_ad))))

jx_adl = subset(read_sas("jx_adult_cn_l.sas7bdat"))

ir_ad = subset(read_sas("ir_adult_cn_c.sas7bdat"))
ir_ad <- ir_ad %>% select(noquote(order(colnames(ir_ad))))

ir_adl = subset(read_sas("ir_adult_cn_l.sas7bdat"))

#join all together so I can perform habitat selection on a complete data set as opposed to on each individual bay
full_ad <- rbind(ap_ad,ck_ad,tb_ad,ch_ad,jx_ad,ir_ad)

##### SELECT CATEGORICAL HABITAT VALUES_ ADULT_ ########
# to be used in the model to predict positive numbers
################################################

#Based on FWRI code the three variables were bottom type (bStr, bsan, bmud), bottom vegetation (bveg), and shoreline (Shore)
#There are three different bottom type variables each of them coded in a binary form.
#I want to take bStr, bsan, and bmud and put them into 1 variable so I will make a new variable entirely = 'bottom'
#I also want to turn bveg into a new variable = 'veg' based on the entries. If alg or Sav then turn to SAV because there are only 9 entries for Alg. 
#Same thing for the shore variable = 'shore'. Decided to have only emergent, structure, terrestrial, and mangrove. 
#Removed old variables (bStr, bSan, bMud, bveg, Shore)
#Removed rows when there was no shoreline variable. 

full_ad <-  select(full_ad, c(number,year, month, bio_reference, bStr, bSan, bMud, bveg, Shore, bay, Zone)) %>% 
  mutate(bottom = ifelse(full_ad$bStr ==1, "structure", ifelse(full_ad$bSan>0 | full_ad$bMud>0, "mudsand", "unknown")), 
         veg= ifelse(full_ad$bveg == "SAVAlg", "SAV", ifelse(full_ad$bveg == "Alg", "SAV", ifelse(full_ad$bveg =="SAV", "SAV", "Noveg"))),
         shore = ifelse(substr(full_ad$Shore,1,3)=="Eme", "Emerge", ifelse(substr(full_ad$Shore,1,3) =="Man", "Mangrove", ifelse(substr(full_ad$Shore,1,3)=="Str", "Structure", 
            ifelse(substr(full_ad$Shore, 1,3)=="Ter", "Terrestrial", "Non"))))) %>% select(-c(bStr, bSan, bMud, bveg, Shore)) %>% subset(!shore=="Non") 

#Turn habitat variables into factors so they can be treated as categorical
full_ad[,c(2,5:9)] <- lapply(full_ad[,c(2,5:9)], factor)

###### MAKE POSITIVE & BINARY SET_ADULT ##########
# to determine the total number of positive huals and proportion positive
##############################################

full_ad.pos<- full_ad %>% subset(number>0)
full_ad.bin <- full_ad %>% mutate(number=ifelse(number>0,1,0))

ap_ad.pos <- full_ad.pos %>% subset(bay =='AP')
ck_ad.pos <- full_ad.pos %>% subset(bay =='CK')
tb_ad.pos <- full_ad.pos %>% subset(bay =='TB')
ch_ad.pos <- full_ad.pos %>% subset(bay =='CH')
jx_ad.pos <- full_ad.pos %>% subset(bay =='JX')
ir_ad.pos <- full_ad.pos %>% subset(bay =='IR')

ap_ad.bin <- full_ad.bin %>% subset(bay =='AP')
ck_ad.bin <- full_ad.bin %>% subset(bay =='CK')
tb_ad.bin <- full_ad.bin %>% subset(bay =='TB')
ch_ad.bin <- full_ad.bin %>% subset(bay =='CH')
jx_ad.bin <- full_ad.bin %>% subset(bay =='JX')
ir_ad.bin <- full_ad.bin %>% subset(bay =='IR')

### BUILD MODELS_ADULT #########
# To produce predicted positive numbers and binomial data set
#######################

# 1. Build the full models with all potential variables and a base model with only year. 
#    Do this for both the positive (Poisson distribution) and binary (Binomial distribution) datasets. 
# 2. Check for overdispersion in the Poisson distribution scenario. 
# 3. If there is overdispersion use quasipoisson 

# 1. Build the full models for the positive and binomial datasets.  
#AP
Full_ap_ad.pos <- glm(number ~ year+month+bottom+veg+shore, data=ap_ad.pos, family=poisson)
Full_ap_ad.bin <- glm(number ~ year+month+bottom+veg+shore, data=ap_ad.bin, family=binomial)

#CK
Full_ck_ad.pos <- glm(number ~ year+month+bottom+veg+shore, data=ck_ad.pos, family=poisson)
Full_ck_ad.bin <- glm(number ~ year+month+bottom+veg+shore, data=ck_ad.bin, family=binomial)

#TB
Full_tb_ad.pos <- glm(number ~ year+month+bottom+veg+shore, data=tb_ad.pos, family=poisson)
Full_tb_ad.bin <- glm(number ~ year+month+bottom+veg+shore, data=tb_ad.bin, family=binomial)

#CH
Full_ch_ad.pos <- glm(number ~ year+month+bottom+veg+shore, data=ch_ad.pos, family=poisson)
Full_ch_ad.bin <- glm(number ~ year+month+bottom+veg+shore, data=ch_ad.bin, family=binomial)

#JX
Full_jx_ad.pos <- glm(number ~ year+month+bottom+veg+shore, data=jx_ad.pos, family=poisson)
Full_jx_ad.bin <- glm(number ~ year+month+bottom+veg+shore, data=jx_ad.bin, family=binomial)

#IR
Full_ir_ad.pos <- glm(number ~ year+month+bottom+veg+shore, data=ir_ad.pos, family=poisson)
Full_ir_ad.bin <- glm(number ~ year+month+bottom+veg+shore, data=ir_ad.bin, family=binomial)

#2. Test the Poisson GLMs for overdispersion
# With the Bernoulli GLM (binomial, response variable is a vector of zeros and ones) overdispersion does not ever occur (Zuur og 253) so I don't need to test for overdispersion in the .bin models. 
dispersiontest(Full_ap_ad.pos, trafo=1)
dispersiontest(Full_ck_ad.pos, trafo=1)
dispersiontest(Full_tb_ad.pos, trafo=1)
dispersiontest(Full_ch_ad.pos, trafo=1)
dispersiontest(Full_jx_ad.pos, trafo=1)
dispersiontest(Full_ir_ad.pos, trafo=1)

# 3. there is evidence of overdispersion for every bay so use quasipoisson for Positive models (Zuur pg 226)
Full_ap_ad.pos <- glm(number ~ year +month+veg+bottom+shore, data=ap_ad.pos, family=quasipoisson)
Full_ck_ad.pos <- glm(number ~ year +month+veg+bottom+shore, data=ck_ad.pos, family=quasipoisson)
Full_tb_ad.pos <- glm(number ~ year +month+veg+bottom+shore, data=tb_ad.pos, family=quasipoisson)
Full_ch_ad.pos <- glm(number ~ year +month+veg+bottom+shore, data=ch_ad.pos, family=quasipoisson)
Full_jx_ad.pos <- glm(number ~ year +month+veg+bottom+shore, data=jx_ad.pos, family=quasipoisson)
Full_ir_ad.pos <- glm(number ~ year +month+veg+bottom+shore, data=ir_ad.pos, family=quasipoisson)

## MODEL SELECTION POSITIVE w/ DROP1 command_ADULT ######
################################
# Pages 220 is to 230 in Zuur are helpful for following methods. 
# The AIC is not defined for quasipoisson models so can't use the step function like what was used in the FWRI code. 
# Instead, use the drop1 function which is applicable for the quassiPoisson GLM and it is more equivalent to hypothesis testing. (Zuur pg 227)
# If just using Poisson or Bernoulli (binomial) can use step command but this gives AIC- not deviance. (Zuur pg 253)
# Explained deviance is nearly the equivalent of R^2 so use this (Zuur pg 218 for equation), "The smaller the residual deviance the better is the model"

###AP_AD_POS (All Variables = significant factors)
summary(Full_ap_ad.pos)
drop1(Full_ap_ad.pos, test="F")  #model selection in quasipoisson is done using F-ratio (Zuur pg 227)
# all seem significant


### CK_POS (Year, Veg = significant factor)
summary(Full_ck_ad.pos)
drop1(Full_ck_ad.pos, test="F")
#month, bottom, and shore do not appear significant. Drop all sequentially. 

#drop month
M1_ck_ad.pos <- glm(number ~ year+veg+bottom+shore, data=ck_ad.pos, family=quasipoisson)
drop1(M1_ck_ad.pos, test="F")

#drop month, bottom
M2_ck_ad.pos <- glm(number ~ year+veg+shore, data=ck_ad.pos, family=quasipoisson)
drop1(M2_ck_ad.pos, test="F")

#drop month, bottom, and shore
M3_ck_ad.pos <- glm(number ~ year+veg, data=ck_ad.pos, family=quasipoisson)
drop1(M3_ck_ad.pos, test="F")
# Year, Veg = significant factors 

### TB_POS (Year, Veg, Shore = significant factors)
summary(Full_tb_ad.pos)
drop1(Full_tb_ad.pos, test="F")
# bottom does not appear significant. 

# drop bottom
M1_tb_ad.pos <- glm(number ~ year+veg+month+shore, data=tb_ad.pos, family=quasipoisson)
drop1(M1_tb_ad.pos, test="F")

### CH_POS (Year, Month, Veg, Bottom, Shore =significnat factors )
summary(Full_ch_ad.pos)
drop1(Full_ch_ad.pos, test="F")
# Month and shore are not significnat factors 

#drop month and shore
M1_ch_ad.pos <- glm(number ~ year+veg+bottom, data=ch_ad.pos, family=quasipoisson)

### JX_POS (Year, Veg, Shore = significant factors)
summary(Full_jx_ad.pos)
drop1(Full_jx_ad.pos, test="F")
# veg, bottom, and shore do not appear significant. Drop all sequentially

#drop veg
M1_jx_ad.pos <- glm(number ~ year+month+bottom+shore, data=jx_ad.pos, family=quasipoisson)
drop1(M1_jx_ad.pos, test="F")

#drop veg and bottom
M2_jx_ad.pos <- glm(number ~ year+month+shore, data=jx_ad.pos, family=quasipoisson)
drop1(M2_jx_ad.pos, test="F")
# Year, Veg, Shore = significant factors

#drop veg, bottom, and shore
M3_jx_ad.pos <- glm(number ~year+month, data=jx_ad.pos, family=quasipoisson)
drop1(M3_jx_ad.pos, test="F")

### IR_POS (Year, Veg, Bottom = significant factors)
summary(Full_ir_ad.pos)
drop1(Full_ir_ad.pos, test="F")
#bottom and shore do not appear significant. Drop them all. 

# drop bottom
M1_ir_ad.pos <- glm(number ~ year+veg+month+shore, data=ir_ad.pos, family=quasipoisson)
drop1(M1_ir_ad.pos, test="F")

# drop shore
M2_ir_ad.pos <- glm(number ~ year+veg+month, data=ir_ad.pos, family=quasipoisson)
drop1(M2_ir_ad.pos, test="F")

## MODEL SELECTION BINARY w/ DROP1 command_ADULT ######
################################
# pg 253 Zuur
##  AP_BIN 
summary(Full_ap_ad.bin)
drop1(Full_ap_ad.bin, test ='Chi')
# bottom not significant

#drop bottom
M1_ap_ad.bin <- glm(number ~ year+month+veg+shore, data=ap_ad.bin, family=binomial)
drop1(M1_ap_ad.bin, test ="Chi")


## CK_BIN 
summary(Full_ck_ad.bin)
drop1(Full_ck_ad.bin, test ='Chi')
#all are significant

## TB_BIN (Year, Month, Veg, Shore = significant)
summary(Full_tb_ad.bin)
drop1(Full_tb_ad.bin, test ='Chi')
# bottom,veg,shore not significant

#drop bottom 
M1_tb_ad.bin <- glm(number ~ year+month+veg+shore, data=tb_ad.bin, family=binomial)
drop1(M1_tb_ad.bin, test ="Chi")

#drop veg, and shore, bottom
M2_tb_ad.bin <- glm(number ~ year+month, data=tb_ad.bin, family=binomial)
drop1(M2_tb_ad.bin, test="Chi")

##  CH_BIN 
summary(Full_ch_ad.bin)
drop1(Full_ch_ad.bin, test ='Chi')
# nothing is significant except for veg but also need to keep year so I will keep year and veg

M1_ch_ad.bin <- glm(number ~ year+veg, data=ch_ad.bin, family=binomial)
drop1(M1_ch_ad.bin, test='Chi')

## JX_BIN 
summary(Full_jx_ad.bin)
drop1(Full_jx_ad.bin, test ='Chi')
# month, year, and bottom (sort of) are not Sign
 
#drop month, bottom
M1_jx_ad.bin <- glm(number ~ year+bottom+veg+shore, data=jx_ad.bin, family=binomial)
drop1(M1_jx_ad.bin, test ="Chi")


## IR_BIN 
summary(Full_ir_ad.bin)
drop1(Full_ir_ad.bin, test ='Chi')
# nothing is significant exceot for year

#drop bottom 
M1_ir_ad.bin <- glm(number ~ year, data=ir_ad.bin, family=binomial)
drop1(M1_ir_ad.bin, test ="Chi")

### ASSIGN FINAL MODELS_ADULT ###### 
###############################
final_ap_ad.pos = Full_ap_ad.pos
final_ck_ad.pos = M3_ck_ad.pos
final_tb_ad.pos = M1_tb_ad.pos
final_ch_ad.pos= M1_ch_ad.pos
final_jx_ad.pos = M3_jx_ad.pos
final_ir_ad.pos = M2_ir_ad.pos

final_ap_ad.bin = M1_ap_ad.bin
final_ck_ad.bin = Full_ck_ad.bin
final_tb_ad.bin = M2_tb_ad.bin
final_ch_ad.bin = M1_ch_ad.bin
final_jx_ad.bin = M1_jx_ad.bin
final_ir_ad.bin = M1_ir_ad.bin

### DETERMINE LEAST SQUARE MEANS_ADULT###########
# DETERMINE LEAST SQUARE MEANS 
#################################
# Same thing as covariate adjusted means. Basically, determine the mean value of total positive numbers 
# of catch per year controlling for covariates (in this case it would be veg and shore variables). 
# Use lsmeans CRAN document. 

# Looking at the reference grid gives a good idea of over what levels the mean is being averaged. 
ap_ad.rf.grid <- ref.grid(final_ap_ad.pos)
ap_ad.bin.rf.grid <- ref.grid(final_ap_ad.bin)

#Can make predictions using the reference grid. It produces a mean value of numbers based on each scenario combination.  
test = summary(ap_ad.rf.grid)
test_ap_ad.bin= summary(ap_ad.bin.rf.grid)

# Use lsmeans to determine the least square mean of positive values. 
# Display the response scale (as opposed to the log scale which is reported for the Poisson, and the logit scale reported for the Binomial)

#POSITIVE
LSM_ap_ad.pos <- summary(lsmeans(final_ap_ad.pos, 'year', data=ap_ad.pos), type="response")
LSM_ck_ad.pos <- summary(lsmeans(final_ck_ad.pos, 'year', data=ck_ad.pos), type="response")
LSM_tb_ad.pos <- summary(lsmeans(final_tb_ad.pos, 'year', data=tb_ad.pos), type="response")
LSM_ch_ad.pos <- summary(lsmeans(final_ch_ad.pos, 'year', data=ch_ad.pos), type="response")
LSM_jx_ad.pos <- summary(lsmeans(final_jx_ad.pos, 'year', data=jx_ad.pos), type="response")
LSM_ir_ad.pos <- summary(lsmeans(final_ir_ad.pos, 'year', data=ir_ad.pos), type="response")

#BINOMIAL (with type="response" this is equal to proportion positive)
LSM_ap_ad.bin <- summary(lsmeans(final_ap_ad.bin, 'year', data=ap_ad.bin), type="response")
LSM_ck_ad.bin <- summary(lsmeans(final_ck_ad.bin, 'year', data=ck_ad.bin), type="response")
LSM_tb_ad.bin <- summary(lsmeans(final_tb_ad.bin, 'year', data=tb_ad.bin), type="response")
LSM_ch_ad.bin <- summary(lsmeans(final_ch_ad.bin, 'year', data=ch_ad.bin), type="response")
LSM_jx_ad.bin <- summary(lsmeans(final_jx_ad.bin, 'year', data=jx_ad.bin), type="response")
LSM_ir_ad.bin <- summary(lsmeans(final_ir_ad.bin, 'year', data=ir_ad.bin), type="response")

### ERROR PROPAGATION TO FIND FINAL VALUE (pos * prop.pos)_ADULT #####
# multiply positive lsmean by porportion positive lsmean and use error propagation to determine value and associated error
# using the package Propagate. See example below. 
# https://www.rdocumentation.org/packages/propagate/versions/1.0-4/topics/propagate    

#Must use a for loop to do the error propagation because it goes one row at a time. 

#AP
num.yr = length(LSM_ap_ad.pos$year)  
df <- data.frame(matrix(data=NA, nrow=num.yr, ncol=6)) #make a dataframe for the loop to store results in
for (i in 1:num.yr) {
  x = c(LSM_ap_ad.pos$rate[i], LSM_ap_ad.pos$SE[i])
  y= c(LSM_ap_ad.bin$prob[i], LSM_ap_ad.bin$SE[i])
  EXPR <- expression(x*y)
  DF <- cbind(x,y)
  RES <- propagate(expr=EXPR, data=DF, type='stat', do.sim=TRUE, verbose=TRUE)
  df[i,] <- t(matrix(RES$sim))
}
Mean_AP_ad <- df %>% cbind(LSM_ap_ad.pos$year)
colnames(Mean_AP_ad) <- c("Mean", "SD", "Median", "MAD", "2.5%", "97.5%", "Year")
write.csv(Mean_AP_ad, "~/Desktop/PhD project/Projects/Seatrout/Data/Indices/DeltaMethod Indices/AP_adult_index.csv")

#CK
num.yr = length(LSM_ck_ad.pos$year)  
df <- data.frame(matrix(data=NA, nrow=num.yr, ncol=6)) #make a dataframe for the loop to store results in
for (i in 1:num.yr) {
  x = c(LSM_ck_ad.pos$rate[i], LSM_ck_ad.pos$SE[i])
  y= c(LSM_ck_ad.bin$prob[i], LSM_ck_ad.bin$SE[i])
  EXPR <- expression(x*y)
  DF <- cbind(x,y)
  RES <- propagate(expr=EXPR, data=DF, type='stat', do.sim=TRUE, verbose=TRUE)
  df[i,] <- t(matrix(RES$sim))
}
Mean_CK_ad <- df %>% cbind(LSM_ck_ad.pos$year)
colnames(Mean_CK_ad) <- c("Mean", "SD", "Median", "MAD", "2.5%", "97.5%", "Year")
write.csv(Mean_CK_ad, "~/Desktop/PhD project/Projects/Seatrout/Data/Indices/DeltaMethod Indices/CK_adult_index.csv")

#TB
num.yr = length(LSM_tb_ad.pos$year)  
df <- data.frame(matrix(data=NA, nrow=num.yr, ncol=6)) #make a dataframe for the loop to store results in
for (i in 1:num.yr) {
  x = c(LSM_tb_ad.pos$rate[i], LSM_tb_ad.pos$SE[i])
  y= c(LSM_tb_ad.bin$prob[i], LSM_tb_ad.bin$SE[i])
  EXPR <- expression(x*y)
  DF <- cbind(x,y)
  RES <- propagate(expr=EXPR, data=DF, type='stat', do.sim=TRUE, verbose=TRUE)
  df[i,] <- t(matrix(RES$sim))
}
Mean_TB_ad <- df %>% cbind(LSM_tb_ad.pos$year)
colnames(Mean_TB_ad) <- c("Mean", "SD", "Median", "MAD", "2.5%", "97.5%", "Year")
write.csv(Mean_TB_ad, "~/Desktop/PhD project/Projects/Seatrout/Data/Indices/DeltaMethod Indices/TB_adult_index.csv")

#CH
num.yr = length(LSM_ch_ad.pos$year)  
df <- data.frame(matrix(data=NA, nrow=num.yr, ncol=6)) #make a dataframe for the loop to store results in
for (i in 1:num.yr) {
  x = c(LSM_ch_ad.pos$rate[i], LSM_ch_ad.pos$SE[i])
  y= c(LSM_ch_ad.bin$prob[i], LSM_ch_ad.bin$SE[i])
  EXPR <- expression(x*y)
  DF <- cbind(x,y)
  RES <- propagate(expr=EXPR, data=DF, type='stat', do.sim=TRUE, verbose=TRUE)
  df[i,] <- t(matrix(RES$sim))
}
Mean_CH_ad <- df %>% cbind(LSM_ch_ad.pos$year)
colnames(Mean_CH_ad) <- c("Mean", "SD", "Median", "MAD", "2.5%", "97.5%", "Year")
write.csv(Mean_CH_ad, "~/Desktop/PhD project/Projects/Seatrout/Data/Indices/DeltaMethod Indices/CH_adult_index.csv")

#JX
num.yr = length(LSM_jx_ad.pos$year)  
df <- data.frame(matrix(data=NA, nrow=num.yr, ncol=6)) #make a dataframe for the loop to store results in
for (i in 1:num.yr) {
  x = c(LSM_jx_ad.pos$rate[i], LSM_jx_ad.pos$SE[i])
  y= c(LSM_jx_ad.bin$prob[i], LSM_jx_ad.bin$SE[i])
  EXPR <- expression(x*y)
  DF <- cbind(x,y)
  RES <- propagate(expr=EXPR, data=DF, type='stat', do.sim=TRUE, verbose=TRUE)
  df[i,] <- t(matrix(RES$sim))
}
Mean_JX_ad <- df %>% cbind(LSM_jx_ad.pos$year)
colnames(Mean_JX_ad) <- c("Mean", "SD", "Median", "MAD", "2.5%", "97.5%", "Year")
write.csv(Mean_JX_ad, "~/Desktop/PhD project/Projects/Seatrout/Data/Indices/DeltaMethod Indices/JX_adult_index.csv")

#IR
num.yr = length(LSM_ir_ad.pos$year)  
df <- data.frame(matrix(data=NA, nrow=num.yr, ncol=6)) #make a dataframe for the loop to store results in
for (i in 1:num.yr) {
  x = c(LSM_ir_ad.pos$rate[i], LSM_ir_ad.pos$SE[i])
  y= c(LSM_ir_ad.bin$prob[i], LSM_ir_ad.bin$SE[i])
  EXPR <- expression(x*y)
  DF <- cbind(x,y)
  RES <- propagate(expr=EXPR, data=DF, type='stat', do.sim=TRUE, verbose=TRUE)
  df[i,] <- t(matrix(RES$sim))
}
Mean_IR_ad <- df %>% cbind(LSM_ir_ad.pos$year)
colnames(Mean_IR_ad) <- c("Mean", "SD", "Median", "MAD", "2.5%", "97.5%", "Year")
write.csv(Mean_IR_ad, "~/Desktop/PhD project/Projects/Seatrout/Data/Indices/DeltaMethod Indices/IR_adult_index.csv")



####OLD #######################
# Make positive dataset

ap_ad.pos <- ap_ad %>% subset(number>0) %>% group_by(year) %>% summarize(totalnumberpositivehauls=length(unique(bio_reference)), TotalNumberOfSeatroutInPosHauls=sum(number))  %>% 
  mutate(positive = TotalNumberOfSeatroutInPosHauls/totalnumberpositivehauls)

ch_ad.pos <- ch_ad %>% subset(number>0) %>% group_by(year) %>% summarize(totalnumberpositivehauls=length(unique(bio_reference)), TotalNumberOfSeatroutInPosHauls=sum(number))  %>% 
  mutate(positive = TotalNumberOfSeatroutInPosHauls/totalnumberpositivehauls)

ck_ad.pos <- ck_ad %>% subset(number>0) %>% group_by(year) %>% summarize(totalnumberpositivehauls=length(unique(bio_reference)), TotalNumberOfSeatroutInPosHauls=sum(number))  %>% 
  mutate(positive = TotalNumberOfSeatroutInPosHauls/totalnumberpositivehauls)

tb_ad.pos <- tb_ad %>% subset(number>0) %>% group_by(year) %>% summarize(totalnumberpositivehauls=length(unique(bio_reference)), TotalNumberOfSeatroutInPosHauls=sum(number))  %>% 
  mutate(positive = TotalNumberOfSeatroutInPosHauls/totalnumberpositivehauls)

ir_ad.pos <- ir_ad %>% subset(number>0) %>% group_by(year) %>% summarize(totalnumberpositivehauls=length(unique(bio_reference)), TotalNumberOfSeatroutInPosHauls=sum(number))  %>% 
  mutate(positive = TotalNumberOfSeatroutInPosHauls/totalnumberpositivehauls)

jx_ad.pos <- jx_ad %>% subset(number>0) %>% group_by(year) %>% summarize(totalnumberpositivehauls=length(unique(bio_reference)), TotalNumberOfSeatroutInPosHauls=sum(number))  %>% 
  mutate(positive = TotalNumberOfSeatroutInPosHauls/totalnumberpositivehauls)

###################
# Make binomial dataset

ap_ad.bin = ap_ad %>% mutate(HaulCategory= ifelse(ap_ad$number>0,1,0)) %>% group_by(year) %>% 
  summarize(TotalPosHauls= sum(HaulCategory), TotalHauls = length(unique(bio_reference))) %>%
  mutate(ProportionPositive = TotalPosHauls/TotalHauls)

ch_ad.bin = ch_ad %>% mutate(HaulCategory= ifelse(ch_ad$number>0,1,0)) %>% group_by(year) %>% 
  summarize(TotalPosHauls= sum(HaulCategory), TotalHauls = length(unique(bio_reference))) %>%
  mutate(ProportionPositive = TotalPosHauls/TotalHauls)

ck_ad.bin = ck_ad %>% mutate(HaulCategory= ifelse(ck_ad$number>0,1,0)) %>% group_by(year) %>% 
  summarize(TotalPosHauls= sum(HaulCategory), TotalHauls = length(unique(bio_reference))) %>%
  mutate(ProportionPositive = TotalPosHauls/TotalHauls)

tb_ad.bin = tb_ad %>% mutate(HaulCategory= ifelse(tb_ad$number>0,1,0)) %>% group_by(year) %>% 
  summarize(TotalPosHauls= sum(HaulCategory), TotalHauls = length(unique(bio_reference))) %>%
  mutate(ProportionPositive = TotalPosHauls/TotalHauls)

ir_ad.bin = ir_ad %>% mutate(HaulCategory= ifelse(ir_ad$number>0,1,0)) %>% group_by(year) %>% 
  summarize(TotalPosHauls= sum(HaulCategory), TotalHauls = length(unique(bio_reference))) %>%
  mutate(ProportionPositive = TotalPosHauls/TotalHauls)

jx_ad.bin = jx_ad %>% mutate(HaulCategory= ifelse(jx_ad$number>0,1,0)) %>% group_by(year) %>% 
  summarize(TotalPosHauls= sum(HaulCategory), TotalHauls = length(unique(bio_reference))) %>%
  mutate(ProportionPositive = TotalPosHauls/TotalHauls)


###################################
# Produce Adjusted Indices for Adult
##################################

AP_ad <- cbind(ap_ad.bin$year, ap_ad.pos$totalnumberpositivehauls, ap_ad.pos$TotalNumberOfSeatroutInPosHauls, ap_ad.bin$TotalHauls, data.frame(ap_ad.pos$positive*ap_ad.bin$ProportionPositive)) 
names(AP_ad) <- c('year','Pos_Hauls', 'Tot_C.Neb_in_Pos', 'All_Hauls', 'index')
write.csv(AP_ad, "~/Desktop/Github Repo/Seatrout/Data/Indices/DeltaMethod Indices/AP_adult_index.csv")

CH_ad <- cbind(ch_ad.bin$year, ch_ad.pos$totalnumberpositivehauls, ch_ad.pos$TotalNumberOfSeatroutInPosHauls, ch_ad.bin$TotalHauls, data.frame(ch_ad.pos$positive*ch_ad.bin$ProportionPositive)) 
names(CH_ad) <- c('year','Pos_Hauls', 'Tot_C.Neb_in_Pos', 'All_Hauls', 'index')
write.csv(CH_ad, "~/Desktop/Github Repo/Seatrout/Data/Indices/DeltaMethod Indices/CH_adult_index.csv")

CK_ad <- cbind(ck_ad.bin$year, ck_ad.pos$totalnumberpositivehauls, ck_ad.pos$TotalNumberOfSeatroutInPosHauls, ck_ad.bin$TotalHauls, data.frame(ck_ad.pos$positive*ck_ad.bin$ProportionPositive)) 
names(CK_ad) <- c('year','Pos_Hauls', 'Tot_C.Neb_in_Pos', 'All_Hauls', 'index')
write.csv(CK_ad, "~/Desktop/Github Repo/Seatrout/Data/Indices/DeltaMethod Indices/CK_adult_index.csv")

TB_ad <- cbind(tb_ad.bin$year, tb_ad.pos$totalnumberpositivehauls, tb_ad.pos$TotalNumberOfSeatroutInPosHauls, tb_ad.bin$TotalHauls, data.frame(tb_ad.pos$positive*tb_ad.bin$ProportionPositive)) 
names(TB_ad) <- c('year','Pos_Hauls', 'Tot_C.Neb_in_Pos', 'All_Hauls', 'index')
write.csv(TB_ad, "~/Desktop/Github Repo/Seatrout/Data/Indices/DeltaMethod Indices/TB_adult_index.csv")

IR_ad <- cbind(ir_ad.bin$year, ir_ad.pos$totalnumberpositivehauls, ir_ad.pos$TotalNumberOfSeatroutInPosHauls, ir_ad.bin$TotalHauls, data.frame(ir_ad.pos$positive*ir_ad.bin$ProportionPositive)) 
names(IR_ad) <- c('year','Pos_Hauls', 'Tot_C.Neb_in_Pos', 'All_Hauls', 'index')
write.csv(IR_ad, "~/Desktop/Github Repo/Seatrout/Data/Indices/DeltaMethod Indices/IR_adult_index.csv")

JX_ad <- cbind(jx_ad.bin$year, jx_ad.pos$totalnumberpositivehauls, jx_ad.pos$TotalNumberOfSeatroutInPosHauls, jx_ad.bin$TotalHauls, data.frame(jx_ad.pos$positive*jx_ad.bin$ProportionPositive)) 
names(JX_ad) <- c('year','Pos_Hauls', 'Tot_C.Neb_in_Pos', 'All_Hauls', 'index')
write.csv(JX_ad, "~/Desktop/Github Repo/Seatrout/Data/Indices/DeltaMethod Indices/JX_adult_index.csv")

###################################################
# Make combined biomass indices for YOY and Adults 
####################################################
#I will trim the yoy data because the adult timeseries are shorter- only for TB, CH, IR, and CK

AP_bio <- data.frame(cbind(AP_ad$index, AP$index)) %>% mutate(logyoy=log(AP$index), logadult=log(AP_ad$index))
names(AP_bio) <- c("adult", "yoy", "logyoy", "logadult")

CH_bio <- data.frame(cbind(CH_ad$index, CH$index[8:27])) %>% mutate(logyoy=log(CH$index[8:27]), logadult=log(CH_ad$index))
names(CH_bio) <- c("adult", "yoy", "logyoy", "logadult")

CK_bio <- data.frame(cbind(CK_ad$index, CK$index[2:20])) %>% mutate(logyoy=log(CK$index[2:20]), logadult=log(CK_ad$index))
names(CK_bio) <- c("adult", "yoy", "logyoy", "logadult")

TB_bio <- data.frame(cbind(TB_ad$index, TB$index[8:27])) %>% mutate(logyoy=log(TB$index[8:27]), logadult=log(TB_ad$index))
names(TB_bio) <- c("adult", "yoy", "logyoy", "logadult")

IR_bio <- data.frame(cbind(IR_ad$index, IR$index[8:26])) %>% mutate(logyoy=log(IR$index[8:26]), logadult=log(IR_ad$index))
names(IR_bio) <- c("adult", "yoy", "logyoy", "logadult")

JX_bio <- data.frame(cbind(JX_ad$index, JX$index)) %>% mutate(logyoy=log(JX$index), logadult=log(JX_ad$index))
names(JX_bio) <- c("adult", "yoy", "logyoy", "logadult")


######################################################
# FIT SR CURVES TO DETERMINE RESIDUALS
#####################################################

library(FSA)
library(nlstools)
#BEVERTON-HOLT, then density idependent, then immediately followed by a Ricker Model

#AP

#bh
bhs <-srStarts(yoy~adult, data=AP_bio, type="BevertonHolt", param=1)
unlist(bhs)
svR_ap <- srStarts(yoy ~ adult, data=AP_bio, type="BevertonHolt")
svR_ap <- list(a=15, b=18)
bh <- srFuns("BevertonHolt")
srBH_ap <- nls(logyoy~log(bh(adult,a,b)), data=AP_bio, start=svR_ap)
overview(srBH_ap)
#density independent
bh0 <- logyoy ~log(a*adult)
bh0s <- bhs[1]
bh0nls <- nls(bh0, data=AP_bio, start=bh0s, algorithm="port") #algorithm helps it not wander into negative territory

anova(bh0nls, srBH_ap) #Beverton holt is better than density independent

#plot bh and density independent
plot(yoy~adult,data=AP_bio)
curve((coef(srBH_ap)[1]*x)/(1+coef(srBH_ap)[2]*x),from=0,to=120,col="red",lwd=2,add=TRUE)
curve(coef(bh0nls)[1]*x,from=0,to=120,col="blue",lwd=2,add=TRUE)
legend("topleft",legend=c("density independent","density dependent"),col=c("blue","red"),lwd=2,cex=0.6)

#plot bh
x=seq(0,5, length.out=999)
pR <- bh(x, a=coef(srBH_ap))
xlmts=range(c(x,AP_bio$adult))
plot(yoy~adult, data=AP_bio, xlim=xlmts)
lines(pR~x, lwd=2)

#ricker
svR_ap <- srStarts(yoy ~ adult, data=AP_bio, type="Ricker")
svR_ap <- list(a=1.64, b=0.45)
rk <- srFuns("Ricker")
srrk_ap <- nls(logyoy~log(rk(adult,a,b)), data=AP_bio, start=svR_ap)
overview(srrk_ap)
x=seq(0,5, length.out=999)
pR <- rk(x, a=coef(srrk_ap))
xlmts=range(c(x,AP_bio$adult))
plot(yoy~adult, data=AP_bio, xlim=xlmts)
lines(pR~x, lwd=2)

AIC(srBH_ap, srrk_ap)

# BH appears to be a better fit 

write.csv(data.frame(residuals(srBH_ap)) %>% mutate(year = c(1998:2015)), "~/Desktop/Github Repo/Seatrout/Data/Exported R Dataframes/ap_resid_BH.csv")
write.csv(data.frame(residuals(srrk_ap)) %>% mutate(year = c(1998:2015)), "~/Desktop/Github Repo/Seatrout/Data/Exported R Dataframes/ap_resid_RK.csv")


#CH
svR_ch <- srStarts(yoy ~ adult, data=CH_bio, type="BevertonHolt")
svR_ch <- list(a=-19.92, b=-21.60)
bh <- srFuns("BevertonHolt")
srBH_ch <- nls(logyoy~log(bh(adult,a,b)), data=CH_bio, start=svR_ch)
overview(srBH_ch)
x=seq(0,2, length.out=999)
pR <- bh(x, a=coef(srBH_ch))
xlmts=range(c(x,CH_bio$adult))
plot(yoy~adult, data=CH_bio, xlim=xlmts)
lines(pR~x, lwd=2)

svR_ch <- srStarts(yoy ~ adult, data=CH_bio, type="Ricker")
svR_ch <- list(a=5, b=2)
rk <- srFuns("Ricker")
srrk_ch <- nls(logyoy~log(rk(adult,a,b)), data=CH_bio, start=svR_ch)
overview(srrk_ch)
x=seq(0,2, length.out=999)
pR <- rk(x, a=coef(srrk_ch))
xlmts=range(c(x,CH_bio$adult))
plot(yoy~adult, data=CH_bio, xlim=xlmts)
lines(pR~x, lwd=2)

AIC(srBH_ch, srrk_ch)

# ricker appears to be a better fit 
write.csv(data.frame(residuals(srrk_ch)) %>% mutate(year = c(1996:2015)), "~/Desktop/Github Repo/Seatrout/Data/Exported R Dataframes/ch_resid_RK.csv")
write.csv(data.frame(residuals(srBH_ch)) %>% mutate(year = c(1996:2015)), "~/Desktop/Github Repo/Seatrout/Data/Exported R Dataframes/ch_resid_BH.csv")

#CK
svR_ck <- srStarts(yoy ~ adult, data=CK_bio, type="BevertonHolt")
svR_ck <- list(a=0.45, b=0.04)
bh <- srFuns("BevertonHolt")
srBH_ck <- nls(logyoy~log(bh(adult,a,b)), data=CK_bio, start=svR_ck)
summary(srBH_ck)
x=seq(0,4, length.out=999)
pR <- bh(x, a=coef(srBH_ck))
xlmts=range(c(x,CK_bio$adult))
plot(yoy~adult, data=CK_bio, xlim=xlmts)
lines(pR~x, lwd=2)

svR_ck <- srStarts(yoy ~ adult, data=CK_bio, type="Ricker")
svR_ck <- list(a=0.52, b=0.04)
rk <- srFuns("Ricker")
srrk_ck <- nls(logyoy~log(rk(adult,a,b)), data=CK_bio, start=svR_ck)
summary(srrk_ck)
x=seq(0,4, length.out=999)
pR <- rk(x, a=coef(srrk_ch))
xlmts=range(c(x,CK_bio$adult))
plot(yoy~adult, data=CK_bio, xlim=xlmts)
lines(pR~x, lwd=2)

#both seem similar at significance level but the Beverton-Holt seems to fit much better
write.csv(data.frame(residuals(srrk_ck)) %>% mutate(year = c(1997:2015)), "~/Desktop/Github Repo/Seatrout/Data/Exported R Dataframes/ck_resid_RK.csv")
write.csv(data.frame(residuals(srBH_ck)) %>% mutate(year = c(1997:2015)), "~/Desktop/Github Repo/Seatrout/Data/Exported R Dataframes/ck_resid_BH.csv")


# TB

svR_tb <- srStarts(yoy ~ adult, data=TB_bio, type="BevertonHolt")
svR_tb <- list(a=10.415, b=6.192)
bh <- srFuns("BevertonHolt")
srBH_tb <- nls(logyoy~log(bh(adult,a,b)), data=TB_bio, start=svR_tb)
summary(srBH_tb)
x=seq(0,4, length.out=999)
pR <- bh(x, a=coef(srBH_tb))
xlmts=range(c(x,TB_bio$adult))
plot(yoy~adult, data=TB_bio, xlim=xlmts)
lines(pR~x, lwd=2)

svR_tb <- srStarts(yoy ~ adult, data=TB_bio, type="Ricker")
svR_tb <- list(a=5, b=1)
rk <- srFuns("Ricker")
srrk_tb <- nls(logyoy~log(rk(adult,a,b)), data=TB_bio, start=svR_tb)
summary(srrk_tb)
x=seq(0,4, length.out=999)
pR <- rk(x, a=coef(srrk_tb))
xlmts=range(c(x,TB_bio$adult))
plot(yoy~adult, data=TB_bio, xlim=xlmts)
lines(pR~x, lwd=2)

#ricker is better fit
write.csv(data.frame(residuals(srBH_tb)) %>% mutate(year= c(1996:2015)), "~/Desktop/Github Repo/Seatrout/Data/Exported R Dataframes/tb_resid_BH.csv")
write.csv(data.frame(residuals(srrk_tb)) %>% mutate(year= c(1996:2015)), "~/Desktop/Github Repo/Seatrout/Data/Exported R Dataframes/tb_resid_RK.csv")


# IR

svR_ir <- srStarts(yoy ~ adult, data=IR_bio, type="BevertonHolt")
svR_ir <- list(a=-56.92, b=-57.60)
bh <- srFuns("BevertonHolt")
srBH_ir <- nls(logyoy~log(bh(adult,a,b)), data=IR_bio, start=svR_ir)
summary(srBH_ir)
x=seq(0,4, length.out=999)
pR <- bh(x, a=coef(srBH_ir))
xlmts=range(c(x,IR_bio$adult))
plot(yoy~adult, data=IR_bio, xlim=xlmts)
lines(pR~x, lwd=2)

svR_ir <- srStarts(yoy ~ adult, data=IR_bio, type="Ricker")
svR_ir <- list(a=4, b=1)
rk <- srFuns("Ricker")
srrk_ir <- nls(logyoy~log(rk(adult,a,b)), data=IR_bio, start=svR_ir)
summary(srrk_ir)
x=seq(0,4, length.out=999)
pR <- rk(x, a=coef(srrk_ir))
xlmts=range(c(x, IR_bio$adult))
plot(yoy~adult, data=IR_bio, xlim=xlmts)
lines(pR~x, lwd=2)

#ricker is better fit
write.csv(data.frame(residuals(srBH_ir))  %>% mutate(year = c(1997:2015)), "~/Desktop/Github Repo/Seatrout/Data/Exported R Dataframes/ir_resid_BH.csv")
write.csv(data.frame(residuals(srrk_ir))  %>% mutate(year = c(1997:2015)), "~/Desktop/Github Repo/Seatrout/Data/Exported R Dataframes/ir_resid_RK.csv")

#JX

svR_jx <- srStarts(yoy ~ adult, data=JX_bio, type="BevertonHolt")
svR_jx <- list(a=0.40, b=0.34)
bh <- srFuns("BevertonHolt")
srBH_jx <- nls(logyoy~log(bh(adult,a,b)), data=JX_bio, start=svR_jx)
summary(srBH_jx)
x=seq(0,4, length.out=999)
pR <- bh(x, a=coef(srBH_jx))
xlmts=range(c(x,JX_bio$adult))
plot(yoy~adult, data=JX_bio, xlim=xlmts)
lines(pR~x, lwd=2)

svR_jx <- srStarts(yoy ~ adult, data=JX_bio, type="Ricker")
svR_jx <- list(a=0.75, b=0.79)
rk <- srFuns("Ricker")
srrk_jx <- nls(logyoy~log(rk(adult,a,b)), data=JX_bio, start=svR_jx)
summary(srrk_jx)
x=seq(0,4, length.out=999)
pR <- rk(x, a=coef(srrk_jx))
xlmts=range(c(x, JX_bio$adult))
plot(yoy~adult, data=JX_bio, xlim=xlmts)
lines(pR~x, lwd=2)

#ricker marinally better
write.csv(data.frame(residuals(srBH_jx)) %>% mutate(year = c(2001:2015)), "~/Desktop/Github Repo/Seatrout/Data/Exported R Dataframes/jx_resid_BH.csv")
write.csv(data.frame(residuals(srrk_jx)) %>% mutate(year = c(2001:2015)), "~/Desktop/Github Repo/Seatrout/Data/Exported R Dataframes/jx_resid_RK.csv")
