###### ABOUT ##############
# 10/26/2016
# Purpose: To produce adjusted nominal catch indices using the Delta Method
# 1. Use Delta Method to produce YOY indices using predicted positive and proportion positive data sets
# 2. Use Delta Method to produce adult indices using predicted positive and proportion positive data sets
# 3. Fit SR relationships to yoy and adult indices to produce the residuals of Beverton Holt and Ricker

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

##### IMPORT DATA SETS ######
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

setwd("~/Desktop/PhD project/Projects/Seatrout/FWRI SCRATCH FOLDER/Elizabeth Herdter/SAS data sets/FIMData/NEWNov7")

#load the data, select the peak reproductive months, reorder columns alphabetically so I can combine dataframes (some columns were in different position in other df)
ap = subset(read_sas("ap_yoy_cn_c.sas7bdat"), month %in% c(6,7,8,9,10,11)) %>% mutate(bUnk=bunk) %>% select(-bunk) 
ap <- ap %>% select(noquote(order(colnames(ap))))
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

##### SELECT CATEGORICAL HABITAT VALUES ########
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
full[,c(2,5:8)] <- lapply(full[,c(2,5:8)], factor)

###### MAKE POSITIVE & BINARY SET ##########
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

##### VISUALIZE THE DATA ########
################################

#Plot the data

#AP
plot(ap.pos$bottom, ap.pos$number, xlab= "bottom type", ylab="number")
plot(ap.pos$year, ap.pos$number, vlab="year", ylab="number")
plot(ap.pos$shore, ap.pos$number, vlab="shore", ylab="number")
plot(ap.pos$veg, ap.pos$number, vlab="veg", ylab="number")

### BUILD MODELS #########
# To produce predicted positive numbers and binomial data set
#######################

# 1. Build the full models with all potential variables and a base model with only year. 
#    Do this for both the positive (Poisson distribution) and binary (Binomial distribution) datasets. 
# 2. Check for overdispersion in the Poisson distribution scenario. 
# 3. If there is overdispersion use quasipoisson 



# 1. Build the full and base models for the positive and binomial datasets.  
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
# With the Bernoulli GLM (binomial) overdispersion does not ever occur so I don't need to test for overdispersion in the .bin models. 

dispersiontest(Base_ap.pos,trafo=1)
dispersiontest(Full_ap.pos, trafo=1)
dispersiontest(Base_ck.pos,trafo=1)
dispersiontest(Full_ck.pos, trafo=1)
dispersiontest(Base_tb.pos,trafo=1)
dispersiontest(Full_tb.pos, trafo=1)
dispersiontest(Base_ch.pos,trafo=1)
dispersiontest(Full_ch.pos, trafo=1)
dispersiontest(Base_jx.pos,trafo=1)
dispersiontest(Full_jx.pos, trafo=1)
dispersiontest(Base_ir.pos,trafo=1)
dispersiontest(Full_ir.pos, trafo=1)

# 3. there is evidence of overdispersion for every bay so use quasipoisson for all Positive models (Zuur pg 226)

M_full_ap.pos <- glm(number ~ year +month+veg+bottom+shore, data=ap.pos, family=quasipoisson)
M_base_ap.pos <- glm(number ~ year, data=ap.pos, family=quasipoisson)

M_full_ck.pos <- glm(number ~ year +month+veg+bottom+shore, data=ck.pos, family=quasipoisson)
M_base_ck.pos <- glm(number ~ year, data=ck.pos, family=quasipoisson)

M_full_tb.pos <- glm(number ~ year +month+veg+bottom+shore, data=tb.pos, family=quasipoisson)
M_base_tb.pos <- glm(number ~ year, data=tb.pos, family=quasipoisson)

M_full_ch.pos <- glm(number ~ year +month+veg+bottom+shore, data=ch.pos, family=quasipoisson)
M_base_ch.pos <- glm(number ~ year, data=ch.pos, family=quasipoisson)

M_full_jx.pos <- glm(number ~ year +month+veg+bottom+shore, data=jx.pos, family=quasipoisson)
M_base_jx.pos <- glm(number ~ year, data=jx.pos, family=quasipoisson)

M_full_ir.pos <- glm(number ~ year +month+veg+bottom+shore, data=ir.pos, family=quasipoisson)
M_base_ir.pos <- glm(number ~ year, data=ir.pos, family=quasipoisson)

## MODEL SELECTION w/ DROP1 command ######
################################

# The AIC is not defined for quasipoisson models so can't use the step function like what was used in the FWRI code. 
# Instead, use the drop1 function which is applicable for the quassiPoisson GLM and it is more equivalent to hypothesis testing. (Zuur pg 227)
# If just using Poisson or Bernoulli (binomial) can use step command but this gives AIC- not deviance. (Zuur pg 253)
# Explained deviance is nearly the equivalent of R^2 so use this (Zuur pg 218 for equation)

#AP_POS
drop1(M_full_ap.pos, test="F")
#The only questionable variables might be 'bottom' and 'month'. With the full model the deviance is 3044.9. 
#If we drop month or bottom the deviance only moves to 3067.9 and it appears that the Pr(>F) value is 0.069823 and 0.06999 for bottom and month, respectively 
#Let's try to drop 'month' and see if it results in a smaller deviance

M1_ap.pos <- glm(number ~ year+veg+bottom+shore, data=ap.pos, family=quasipoisson)
drop1(M1_ap.pos, test="F")
#Now the p value for bottom increases to 0.07. Lets drop the 'bottom' variable 

M2_ap.pos <- glm(number ~ year+veg+shore, data=ap.pos, family=quasipoisson)
drop1(M2_ap.pos, test='F')
#Year, veg, and shore are all significant so I can stop now with model selection. 

# AP_BIN
drop1(M)

### ASSIGN FINAL MODELS ###### 
###############################
final_ap.pos = M2_ap.pos 
final_ck.pos = 
final_tb.pos =
final_ch.pos=
final_jx.pos =
final_ir.pos = 
  



### DETERMINE LEAST SQUARE MEANS ###########
# DETERMINE LEAST SQUARE MEANS 
#################################

# Same thing as covariate adjusted means. Basically, determine the mean value of total positive numbers 
# of catch per year controlling for covariates (in this case it would be veg and shore variables). 
# Use lsmeans CRAN document. 

#Look at the reference grid
ap.rf.grid <- ref.grid(glm.final.pos)

#Can make predictions using the reference grid. It produces a mean value of numbers based on each scenario combination.  
summary(ap.rf.grid)

# Use lsmeans to determine the least square mean of positive values. 
lsmeans(glm.final.pos, "year", data=ap.pos)

# Display the response scale (as opposed to the log scale above)
summary(lsmeans(glm.final.pos, 'year', data=ap.pos), type="response")




#ap.pos<- ap %>% subset(number>0) %>% group_by(year) %>% summarize(totalnumberpositivehauls=length(unique(bio_reference)), TotalNumberOfSeatroutInPosHauls=sum(number))  %>% 
#mutate(positive = TotalNumberOfSeatroutInPosHauls/totalnumberpositivehauls)

#apr.pos <- apr %>% subset(number>0) %>% group_by(year) %>% summarize(totalnumberpositivehauls=length(unique(bio_reference)), TotalNumberOfSeatroutInPosHauls=sum(number))  %>% 
#mutate(positive = TotalNumberOfSeatroutInPosHauls/totalnumberpositivehauls)

ap_all.pos <- ap %>% subset(number>0) %>% group_by(year) %>% summarize(totalnumberpositivehauls=length(unique(bio_reference)), TotalNumberOfSeatroutInPosHauls=sum(number))  %>% 
  mutate(positive = TotalNumberOfSeatroutInPosHauls/totalnumberpositivehauls)


# apr.pos doesnt have any positives for 2003, 2005 or 2013 which causes a problem when trying to make the indices because the years are of unequal length

#ck.pos<- ck %>% subset(number>0) %>% group_by(year) %>% summarize(totalnumberpositivehauls=length(unique(bio_reference)), TotalNumberOfSeatroutInPosHauls=sum(number))  %>% 
#mutate(positive = TotalNumberOfSeatroutInPosHauls/totalnumberpositivehauls)

#ckr.pos<- ckr %>% subset(number>0) %>% group_by(year) %>% summarize(totalnumberpositivehauls=length(unique(bio_reference)), TotalNumberOfSeatroutInPosHauls=sum(number))  %>% 
#mutate(positive = TotalNumberOfSeatroutInPosHauls/totalnumberpositivehauls)

ck_all.pos<- ck %>% subset(number>0) %>% group_by(year) %>% summarize(totalnumberpositivehauls=length(unique(bio_reference)), TotalNumberOfSeatroutInPosHauls=sum(number))  %>% 
  mutate(positive = TotalNumberOfSeatroutInPosHauls/totalnumberpositivehauls)

#ch.pos<- ch %>% subset(number>0) %>% group_by(year) %>% summarize(totalnumberpositivehauls=length(unique(bio_reference)), TotalNumberOfSeatroutInPosHauls=sum(number))  %>% 
#mutate(positive = TotalNumberOfSeatroutInPosHauls/totalnumberpositivehauls)

#chr.pos<- chr %>% subset(number>0) %>% group_by(year) %>% summarize(totalnumberpositivehauls=length(unique(bio_reference)), TotalNumberOfSeatroutInPosHauls=sum(number))  %>% 
#mutate(positive = TotalNumberOfSeatroutInPosHauls/totalnumberpositivehauls)

ch_all.pos<- ch %>% subset(number>0) %>% group_by(year) %>% summarize(totalnumberpositivehauls=length(unique(bio_reference)), TotalNumberOfSeatroutInPosHauls=sum(number))  %>% 
  mutate(positive = TotalNumberOfSeatroutInPosHauls/totalnumberpositivehauls)

#tb.pos<- tb %>% subset(number>0) %>% group_by(year) %>% summarize(totalnumberpositivehauls=length(unique(bio_reference)), TotalNumberOfSeatroutInPosHauls=sum(number))  %>% 
#mutate(positive = TotalNumberOfSeatroutInPosHauls/totalnumberpositivehauls)

#tbr.pos<- tbr %>% subset(number>0) %>% group_by(year) %>% summarize(totalnumberpositivehauls=length(unique(bio_reference)), TotalNumberOfSeatroutInPosHauls=sum(number))  %>% 
#mutate(positive = TotalNumberOfSeatroutInPosHauls/totalnumberpositivehauls)

tb_all.pos <-  tb %>% subset(number>0) %>% group_by(year) %>% summarize(totalnumberpositivehauls=length(unique(bio_reference)), TotalNumberOfSeatroutInPosHauls=sum(number))  %>% 
  mutate(positive = TotalNumberOfSeatroutInPosHauls/totalnumberpositivehauls)


#ir.pos<- ir %>% subset(number>0) %>% group_by(year) %>% summarize(totalnumberpositivehauls=length(unique(bio_reference)), TotalNumberOfSeatroutInPosHauls=sum(number))  %>% 
#mutate(positive = TotalNumberOfSeatroutInPosHauls/totalnumberpositivehauls)
# missing positive years of data for 98,99,00,01,03,14,15

#irb.pos<- irb %>% subset(number>0) %>% group_by(year) %>% summarize(totalnumberpositivehauls=length(unique(bio_reference)), TotalNumberOfSeatroutInPosHauls=sum(number))  %>% 
#mutate(positive = TotalNumberOfSeatroutInPosHauls/totalnumberpositivehauls)

ir_all.pos <-  ir %>% subset(number>0) %>% group_by(year) %>% summarize(totalnumberpositivehauls=length(unique(bio_reference)), TotalNumberOfSeatroutInPosHauls=sum(number))  %>% 
  mutate(positive = TotalNumberOfSeatroutInPosHauls/totalnumberpositivehauls)


jx_all.pos<- jx %>% subset(number>0) %>% group_by(year) %>% summarize(totalnumberpositivehauls=length(unique(bio_reference)), TotalNumberOfSeatroutInPosHauls=sum(number))  %>% 
  mutate(positive = TotalNumberOfSeatroutInPosHauls/totalnumberpositivehauls)

#########################################################
# Make Bionomial set for Identifying Positive/Zeros Hauls
#########################################################

# ap.bin = ap %>% mutate(HaulCategory= ifelse(ap$number>0,1,0)) %>% group_by(year) %>%
#           summarize(TotalPosHauls= sum(HaulCategory), TotalHauls = length(unique(bio_reference))) %>%
#           mutate(ProportionPositive = TotalPosHauls/TotalHauls)

# apr.bin = apr %>% mutate(HaulCategory= ifelse(apr$number>0,1,0)) %>% group_by(year) %>%
#           summarize(TotalPosHauls= sum(HaulCategory), TotalHauls = length(unique(bio_reference))) %>%
#           mutate(ProportionPositive = TotalPosHauls/TotalHauls)


ap_all.bin = ap %>% mutate(HaulCategory= ifelse(ap$number>0,1,0)) %>% group_by(year) %>% 
  summarize(TotalPosHauls= sum(HaulCategory), TotalHauls = length(unique(bio_reference))) %>%
  mutate(ProportionPositive = TotalPosHauls/TotalHauls)

# ck.bin = ck %>% mutate(HaulCategory= ifelse(ck$number>0,1,0)) %>% group_by(year) %>% 
#         summarize(TotalPosHauls= sum(HaulCategory), TotalHauls = length(unique(bio_reference))) %>%
#         mutate(ProportionPositive = TotalPosHauls/TotalHauls)
# 
# ckr.bin = ckr %>% mutate(HaulCategory= ifelse(ckr$number>0,1,0)) %>% group_by(year) %>% 
#         summarize(TotalPosHauls= sum(HaulCategory), TotalHauls = length(unique(bio_reference))) %>%
#         mutate(ProportionPositive = TotalPosHauls/TotalHauls)

ck_all.bin <- ck %>% mutate(HaulCategory= ifelse(ck$number>0,1,0)) %>% group_by(year) %>% 
  summarize(TotalPosHauls= sum(HaulCategory), TotalHauls = length(unique(bio_reference))) %>%
  mutate(ProportionPositive = TotalPosHauls/TotalHauls)


# tb.bin = tb %>% mutate(HaulCategory= ifelse(tb$number>0,1,0)) %>% group_by(year) %>% 
#         summarize(TotalPosHauls= sum(HaulCategory), TotalHauls = length(unique(bio_reference))) %>%
#         mutate(ProportionPositive = TotalPosHauls/TotalHauls)
# 
# tbr.bin = tbr %>% mutate(HaulCategory= ifelse(tbr$number>0,1,0)) %>% group_by(year) %>% 
#         summarize(TotalPosHauls= sum(HaulCategory), TotalHauls = length(unique(bio_reference))) %>%
#         mutate(ProportionPositive = TotalPosHauls/TotalHauls)

tb_all.bin <- tb %>% mutate(HaulCategory= ifelse(tb$number>0,1,0)) %>% group_by(year) %>% 
  summarize(TotalPosHauls= sum(HaulCategory), TotalHauls = length(unique(bio_reference))) %>%
  mutate(ProportionPositive = TotalPosHauls/TotalHauls)

# ch.bin = ch %>% mutate(HaulCategory= ifelse(ch$number>0,1,0)) %>% group_by(year) %>% 
#         summarize(TotalPosHauls= sum(HaulCategory), TotalHauls = length(unique(bio_reference))) %>%
#         mutate(ProportionPositive = TotalPosHauls/TotalHauls)
# 
# chr.bin = chr %>% mutate(HaulCategory= ifelse(chr$number>0,1,0)) %>% group_by(year) %>% 
#         summarize(TotalPosHauls= sum(HaulCategory), TotalHauls = length(unique(bio_reference))) %>%
#         mutate(ProportionPositive = TotalPosHauls/TotalHauls)

ch_all.bin <- ch %>% mutate(HaulCategory= ifelse(ch$number>0,1,0)) %>% group_by(year) %>% 
  summarize(TotalPosHauls= sum(HaulCategory), TotalHauls = length(unique(bio_reference))) %>%
  mutate(ProportionPositive = TotalPosHauls/TotalHauls)

# ir.bin = ir %>% mutate(HaulCategory= ifelse(ir$number>0,1,0)) %>% group_by(year) %>% 
#           summarize(TotalPosHauls= sum(HaulCategory), TotalHauls = length(unique(bio_reference))) %>%
#           mutate(ProportionPositive = TotalPosHauls/TotalHauls)
# 
# irb.bin = irb %>% mutate(HaulCategory= ifelse(irb$number>0,1,0)) %>% group_by(year) %>% 
#   summarize(TotalPosHauls= sum(HaulCategory), TotalHauls = length(unique(bio_reference))) %>%
#   mutate(ProportionPositive = TotalPosHauls/TotalHauls)

ir_all.bin <- ir %>% mutate(HaulCategory= ifelse(ir$number>0,1,0)) %>% group_by(year) %>% 
  summarize(TotalPosHauls= sum(HaulCategory), TotalHauls = length(unique(bio_reference))) %>%
  mutate(ProportionPositive = TotalPosHauls/TotalHauls)

# jx.bin = jx %>% mutate(HaulCategory= ifelse(jx$number>0,1,0)) %>% group_by(year) %>% 
#           summarize(TotalPosHauls= sum(HaulCategory), TotalHauls = length(unique(bio_reference))) %>%
#            mutate(ProportionPositive = TotalPosHauls/TotalHauls)

jx_all.bin <- jx %>% mutate(HaulCategory= ifelse(jx$number>0,1,0)) %>% group_by(year) %>% 
  summarize(TotalPosHauls= sum(HaulCategory), TotalHauls = length(unique(bio_reference))) %>%
  mutate(ProportionPositive = TotalPosHauls/TotalHauls)

############################
# Produce Adjusted Indices 
############################

AP <- cbind(ap_all.bin$year, ap_all.pos$totalnumberpositivehauls, ap_all.pos$TotalNumberOfSeatroutInPosHauls, ap_all.bin$TotalHauls, data.frame(ap_all.pos$positive*ap_all.bin$ProportionPositive)) 
names(AP) <- c('year','Pos_Hauls', 'Tot_C.Neb_in_Pos', 'All_Hauls', 'index')
write.csv(AP, "~/Desktop/Github Repo/Seatrout/Data/Indices/DeltaMethod Indices/AP_yoy_index.csv")

#APR <- cbind(apr.bin$year, data.frame(apr_sum$positive*apr.bin$ProportionPositive)) 
#names(APR) <- c('year', 'index')

CH <- cbind(ch_all.bin$year, ch_all.pos$totalnumberpositivehauls,ch_all.pos$TotalNumberOfSeatroutInPosHauls, ch_all.bin$TotalHauls, data.frame(ch_all.pos$positive*ch_all.bin$ProportionPositive))
#CHR <- cbind(chr.bin$year, data.frame(chr.pos$positive*chr.bin$ProportionPositive))
names(CH) <-  c('year','Pos_Hauls', 'Tot_C.Neb_in_Pos', 'All_Hauls', 'index')
write.csv(CH,"~/Desktop/Github Repo/Seatrout/Data/Indices/DeltaMethod Indices/CH_yoy_index.csv" )

CK <- cbind(ck_all.bin$year, ck_all.pos$totalnumberpositivehauls,ck_all.pos$TotalNumberOfSeatroutInPosHauls, ck_all.bin$TotalHauls, data.frame(ck_all.pos$positive*ck_all.bin$ProportionPositive))
#CKR <- cbind(ckr.bin$year, data.frame(ckr.pos$positive*ckr.bin$ProportionPositive))
names(CK) <- c('year','Pos_Hauls', 'Tot_C.Neb_in_Pos', 'All_Hauls', 'index')
write.csv(CK,"~/Desktop/Github Repo/Seatrout/Data/Indices/DeltaMethod Indices/CK_yoy_index.csv" )


TB <- cbind(tb_all.bin$year, tb_all.pos$totalnumberpositivehauls,tb_all.pos$TotalNumberOfSeatroutInPosHauls, tb_all.bin$TotalHauls, data.frame(tb_all.pos$positive*tb_all.bin$ProportionPositive))
#TBR <- cbind(tbr.bin$year, data.frame(tbr.pos$positive*tbr.bin$ProportionPositive))
names(TB) <- c('year','Pos_Hauls', 'Tot_C.Neb_in_Pos', 'All_Hauls', 'index')
write.csv(TB,"~/Desktop/Github Repo/Seatrout/Data/Indices/DeltaMethod Indices/TB_yoy_index.csv" )


IR <- cbind(ir_all.bin$year, ir_all.pos$totalnumberpositivehauls,ir_all.pos$TotalNumberOfSeatroutInPosHauls, ir_all.bin$TotalHauls, data.frame(ir_all.pos$positive*ir_all.bin$ProportionPositive))
#IR <- cbind(ir.bin$year, data.frame(ir.pos$positive*ir.bin$ProportionPositive))
names(IR) <- c('year','Pos_Hauls', 'Tot_C.Neb_in_Pos', 'All_Hauls', 'index')
write.csv(IR,"~/Desktop/Github Repo/Seatrout/Data/Indices/DeltaMethod Indices/IR_yoy_index.csv" )


JX <- cbind(jx_all.bin$year, jx_all.pos$totalnumberpositivehauls,jx_all.pos$TotalNumberOfSeatroutInPosHauls, jx_all.bin$TotalHauls, data.frame(jx_all.pos$positive*jx_all.bin$ProportionPositive))
#JXB <- cbind(jxb.bin$year, data.frame(jxb.pos$positive*jxb.bin$ProportionPositive))
names(JX) <- c('year','Pos_Hauls', 'Tot_C.Neb_in_Pos', 'All_Hauls', 'index')
write.csv(JX,"~/Desktop/Github Repo/Seatrout/Data/Indices/DeltaMethod Indices/JX_yoy_index.csv" )


################################################

###########################
# LOAD ADULT DATA
###########################

################################################
# adult data are only 

setwd("~/Desktop/Github Repo/Seatrout/FWRI SCRATCH FOLDER/Elizabeth Herdter/SAS data sets/FIMData")
ap_ad = subset(read_sas("ap_adult_cn_c.sas7bdat"))
ap_adl = subset(read_sas("ap_adult_cn_l.sas7bdat"))


ch_ad = subset(read_sas("ch_adult_cn_c.sas7bdat")) # *******
ch_adl = subset(read_sas("ch_adult_cn_l.sas7bdat"))


ck_ad = subset(read_sas("ck_adult_cn_c.sas7bdat"))
ck_adl = subset(read_sas("ck_adult_cn_l.sas7bdat"))


tb_ad = subset(read_sas("tb_adult_cn_c.sas7bdat"))
tb_adl = subset(read_sas("tb_adult_cn_l.sas7bdat"))


jx_ad = subset(read_sas("jx_adult_cn_c.sas7bdat"))
jx_adl = subset(read_sas("jx_adult_cn_l.sas7bdat"))


ir_ad = subset(read_sas("ir_adult_cn_c.sas7bdat"))
ir_adl = subset(read_sas("ir_adult_cn_l.sas7bdat"))


###########################
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
