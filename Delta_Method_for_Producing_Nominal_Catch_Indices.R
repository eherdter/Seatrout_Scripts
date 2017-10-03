###### ABOUT ##############
# 10/26/2016
# Purpose: To produce adjusted nominal catch indices using the Delta Method
# 1. Use Delta Method to produce YOY indices using predicted positive and proportion positive data sets
# 2. Use Delta Method to produce adult indices using predicted positive and proportion positive data sets
# 3. Fit SR relationships to yoy and adult indices to produce the residuals of Beverton Holt and Ricker
# 4. Apply proportion of adult numbers and a weight schedule to the predicted index to determine the total predicted adult SSB that is also used
#  when applying the stock recruitment curves 


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

# 9/18/17
# Added in different model formula for fitting abundnace indices. 
# See edits below. 

##### SET WORKING DIRECTORY_YOY ####
#must change working directory for data when working on personal vs work computer
setwd("~/Desktop/PhD project/Projects/Seatrout/FWRI SCRATCH FOLDER/Elizabeth Herdter/SAS data sets/FIMData/NEWNov7")
setwd("U:/PhD_projectfiles/Raw_Data/Seatrout_FIM_Data/FIMData/NEWNov7")

##### SET PACKAGES ######################
library(propagate) #error propagation.. MUST load this first or else it will mask other good functions in dplyr
library(haven) #to load sas
library(dplyr) # to do df manipulation
library(lsmeans) #to determine the least squares means
library(AER) #to test for overdispersion
library(pscl) # to run zeroinflated
library(lmtest) #to test ZINB and ZIP
library(boot) #to diagnose model fits for generalized linear models 
library(MASS) #to use negative binomial

##### IMPORT DATA SETS_YOY ######
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

##### MAKE POSITIVE & BINARY SET_YOY ##########
# to determine the total number of positive huals and proportion positive
# make complete dataset also for use in the mixture model

# full.pos<- full %>% subset(number>0)
# full.bin <- full %>% mutate(number=ifelse(number>0,1,0))
# full.bin$month <- as.factor(as.character(full.bin$month))

# ap.pos <- droplevels(full.pos %>% subset(bay =='AP'))
# ck.pos <- droplevels(full.pos %>% subset(bay =='CK')) 
# tb.pos <- droplevels(full.pos %>% subset(bay =='TB')) 
# ch.pos <- droplevels(full.pos %>% subset(bay =='CH')) 
# jx.pos <- droplevels(full.pos %>% subset(bay =='JX')) 
# ir.pos <- droplevels(full.pos %>% subset(bay =='IR')) 
# 
# ap.bin <- droplevels(full.bin %>% subset(bay =='AP'))
# ck.bin <- droplevels(full.bin %>% subset(bay =='CK'))
# tb.bin <- droplevels(full.bin %>% subset(bay =='TB'))
# ch.bin <- droplevels(full.bin %>% subset(bay =='CH'))
# jx.bin <- droplevels(full.bin %>% subset(bay =='JX'))
# ir.bin <- droplevels(full.bin %>% subset(bay =='IR'))

ap.fl <- droplevels(full %>% subset(bay =='AP'))
ck.fl <- droplevels(full %>% subset(bay =='CK'))
tb.fl <- droplevels(full %>% subset(bay =='TB'))
ch.fl <- droplevels(full %>% subset(bay =='CH'))
jx.fl <- droplevels(full %>% subset(bay =='JX'))
ir.fl <- droplevels(full %>% subset(bay =='IR'))


#check histograms to determine mean 
hist(ap.pos$number)
hist(ck.pos$number)
hist(tb.pos$number)
hist(ch.pos$number)
hist(jx.pos$number)
hist(ir.pos$number)

# CHECK FOR AGGREGATIONS #####
#AP
#month
with(ap.pos,tapply(number, list(year,month),sum))
ap.pos$month[ap.pos$month<7]=7
ap.pos$month[ap.pos$month>10]=10
ap.pos$month <- as.factor(as.character(ap.pos$month))
#veg
with(ap.pos,tapply(number, list(year,veg),sum))
#bottom
with(ap.pos,tapply(number, list(year,bottom),sum))
#shore
with(ap.pos,tapply(number, list(year,shore),sum))
ap.pos <- subset(ap.pos, shore !=  "Mangrove") %>% droplevels(ap.pos$shore)

with(ap.fl,tapply(number, list(year,month),sum))
with(ap.fl,tapply(number, list(year,veg),sum))
with(ap.fl,tapply(number, list(year,bottom),sum))
with(ap.fl,tapply(number, list(year,shore),sum))
ap.fl <- subset(ap.fl, shore != "Mangrove") %>% droplevels(ap.fl$shore)
ap.fl$month <- as.factor(as.character(ap.fl$month))


#CK
#month
with(ck.pos,tapply(number, list(year,month),sum))
ck.pos$month[ck.pos$month<6]=6
ck.pos$month <- as.factor(as.character(ck.pos$month))
#veg
with(ck.pos,tapply(number, list(year,veg),sum))
#bottom
with(ck.pos,tapply(number, list(year,bottom),sum))
ck.pos <- droplevels(subset(ck.pos, bottom != "unknown"))
#shore 
with(ck.pos,tapply(number, list(year,shore),sum))
#drop terrestrial and join structure 
ck.pos <- droplevels(subset(ck.pos, shore != "Terrestrial"))
ck.pos$shore[ck.pos$shore== "Mangrove"] = "Structure"
ck.pos <- droplevels(subset(ck.pos, shore != "Mangrove"))

with(ck.fl,tapply(number, list(year,month),sum))
with(ck.fl,tapply(number, list(year,veg),sum))
with(ck.fl,tapply(number, list(year,bottom),sum))
ck.fl = droplevels(subset(ck.fl, bottom != "unknown"))
with(ck.fl,tapply(number, list(year,shore),sum))
#drop terrestrial and join structure 
ck.fl <- droplevels(subset(ck.fl, shore != "Terrestrial"))
ck.fl$shore[ck.fl$shore== "Mangrove"] = "Structure"
ck.fl <- droplevels(subset(ck.fl, shore != "Mangrove"))
ck.fl$month <- as.factor(as.character(ck.fl$month))

#TB
with(tb.pos,tapply(number, list(year,month),sum))
tb.pos$month[tb.pos$month<5]=5
tb.pos$month <- as.factor(as.character(tb.pos$month))
with(tb.pos,tapply(number, list(year,veg),sum))
with(tb.pos,tapply(number, list(year,bottom),sum))
tb.pos <- droplevels(subset(tb.pos, bottom != "unknown"))
with(tb.pos,tapply(number, list(year,shore),sum))

with(tb.fl,tapply(number, list(year,month),sum))
with(tb.fl,tapply(number, list(year,veg),sum))
with(tb.fl,tapply(number, list(year,bottom),sum))
tb.fl <- droplevels(subset(tb.fl, bottom != "unknown"))
with(tb.fl,tapply(number, list(year,shore),sum))
tb.fl$month <- as.factor(as.character(tb.fl$month))

#CH - no aggregation needed
with(ch.pos,tapply(number, list(year,month),sum))
with(ch.pos,tapply(number, list(year,veg),sum))
with(ch.pos,tapply(number, list(year,bottom),sum))
with(ch.pos,tapply(number, list(year,shore),sum))
ch.pos$month <- as.factor(as.character(ch.pos$month))

with(ch.fl,tapply(number, list(year,month),sum))
with(ch.fl,tapply(number, list(year,veg),sum))
with(ch.fl,tapply(number, list(year,bottom),sum))
with(ch.fl,tapply(number, list(year,shore),sum))
ch.fl$month <- as.factor(as.character(ch.fl$month))

#JX 
with(jx.pos,tapply(number, list(year,month),sum))
with(jx.pos,tapply(number, list(year,veg),sum))
with(jx.pos,tapply(number, list(year,bottom),sum))
with(jx.pos,tapply(number, list(year,shore),sum))
jx.pos$shore[jx.pos$shore== "Terrestrial"] = "Emerge"
jx.pos <- na.omit(jx.pos) #randomly an NA present
jx.pos$month <- as.factor(as.character(jx.pos$month))

with(jx.fl,tapply(number, list(year,month),sum))
with(jx.fl,tapply(number, list(year,veg),sum))
with(jx.fl,tapply(number, list(year,bottom),sum))
jx.fl <- droplevels(subset(jx.fl, bottom != "unknown"))
with(jx.fl,tapply(number, list(year,shore),sum))
jx.fl <- na.omit(jx.fl)
jx.fl$month <- as.factor(as.character(jx.fl$month))

#IR - no aggregation needed 
with(ir.pos,tapply(number, list(year,month),sum))
with(ir.pos,tapply(number, list(year,veg),sum))
with(ir.pos,tapply(number, list(year,bottom),sum))
with(ir.pos,tapply(number, list(year,shore),sum))
ir.pos <- na.omit(ir.pos)
ir.fl <- na.omit(ir.fl)
ir.pos$month <- as.factor(as.character(ir.pos$month))
ir.fl$month <- as.factor(as.character(ir.fl$month))

with(ir.fl,tapply(number, list(year,month),sum))
with(ir.fl,tapply(number, list(year,veg),sum))
with(ir.fl,tapply(number, list(year,bottom),sum))
with(ir.fl,tapply(number, list(year,shore),sum))
ir.fl <- na.omit(ir.fl)
ir.fl$month <- as.factor(as.character(ir.fl$month))

##### VISUALIZE THE DATA_YOY ########
#Plot the data

#AP
plot(ap.pos$bottom, ap.pos$number, xlab= "bottom type", ylab="number")
plot(ap.pos$year, ap.pos$number, vlab="year", ylab="number")
plot(ap.pos$shore, ap.pos$number, vlab="shore", ylab="number")
plot(ap.pos$veg, ap.pos$number, vlab="veg", ylab="number")

##### BUILD MODELS_YOY #########
# There are a few options  for modeling these data. 
# First three options: Two-part models: build a manual two part model where presence/absence is modeled with binomial and positive catch is 
# modeled with a poisson, quasi poisson, or a lognormal. 
# Fourth option: Mixture models: model all data with either a poisson or negative binomial. if that doesnt work try zero inflated negative binomial
# Fifth option: zero truncated model for positive. See note below. Do this last case scenario because it's quite involved. 
    # Because the positive data is forced into zero-truncated (the zeros were removed) they may need to be treated differently.
    # Zuur chapter 11 outlines use of zero-truncated models. You can either use a zero-truncated Poissan model which is good 
    # for count data or you can use a zero trucated negative binomial which will deal with overdispersion if the data are 
    # overdispersed. There is no such thing as a zero-truncated, quasi poisson model (where quasi poisson deals with overdispersion). 
    # Must check for overdispersion. If there isn't then use zero-truncated Poisson. If there is then use zero-truncated NB. 
    # (Chapter 11 Zuur). If you can't decide which models to use in terms of a zero truncated or not (i.e. sometimes even if the data are zero truncated the mean will be large so results will be unaffected by model choice, page 269)
    # you can run different model options and compare model validation plots to decide which fit best. 


# For first three options: 
# 1. Build the binary datasets.
# 2. Build the poisson, and lognormal. 
# 3. Check for overdispersion. Make adjustments for each. i.e. quasi Poisson for overdispersed poisson
# For fourth option: 
# 1. Build Poisson and the negative binomial model mixture models. Test. and If necessary build the mixture zero inflated poisson (ZIP) and zero inflated negative binomial (ZINB)
# For fifth option, last case
# 1. Build zero truncated model- last case scenario. 

# 1. Build the binary data sets ####
# With the Bernoulli GLM (binomial, response variable is a vector of zeros and ones) overdispersion does not ever occur (Zuur og 253) so I don't need to test for overdispersion in the .bin models. 
# Full_ap.bin <- glm(number ~ year+month+bottom+veg+shore, data=ap.bin, family=binomial, na.action = na.exclude)
# Full_ck.bin <- glm(number ~ year+month+bottom+veg+shore, data=ck.bin, family=binomial, na.action = na.exclude)
# Full_tb.bin <- glm(number ~ year+month+bottom+veg+shore, data=tb.bin, family=binomial, na.action = na.exclude)
# Full_ch.bin <- glm(number ~ year+month+bottom+veg+shore, data=ch.bin, family=binomial, na.action = na.exclude)
# Full_jx.bin <- glm(number ~ year+month+bottom+veg+shore, data=jx.bin, family=binomial, na.action = na.exclude)
# Full_ir.bin <- glm(number ~ year+month+bottom+veg+shore, data=ir.bin, family=binomial, na.action = na.exclude)
# 
# # 2. Build the poisson and lognormal ####
# ap.pos.P <- glm(number ~ year+month+bottom+veg+shore, data=ap.pos, family=poisson, na.action = na.exclude)
# ck.pos.P <- glm(number ~ year+month+bottom+veg+shore, data=ck.pos, family=poisson, na.action = na.exclude)
# tb.pos.P <- glm(number ~ year+month+bottom+veg+shore, data=tb.pos, family=poisson, na.action = na.exclude)
# ch.pos.P <- glm(number ~ year+month+bottom+veg+shore, data=ch.pos, family=poisson, na.action = na.exclude)
# jx.pos.P <- glm(number ~ year+month+bottom+veg+shore, data=jx.pos, family=poisson, na.action = na.exclude)
# ir.pos.P <- glm(number ~ year+month+bottom+veg+shore, data=ir.pos, family=poisson, na.action = na.exclude)

# For lognormal. first need to log transform the data so that they can be modeled (lognormal isnt a family definition for genearlized linear) 
# ap.pos$lnum <- log(ap.pos$number)
# ck.pos$lnum <- log(ck.pos$number)
# tb.pos$lnum <- log(tb.pos$number)
# ch.pos$lnum <- log(ch.pos$number)
# jx.pos$lnum <- log(jx.pos$number)
# ir.pos$lnum <- log(ir.pos$number)
# 
# #now model the transformed numbers with the guassian
# ap.pos.L <- glm(lnum ~ year+month+bottom+veg+shore, data=ap.pos, family=gaussian, na.action = na.exclude)
# ck.pos.L <- glm(lnum ~ year+month+bottom+veg+shore, data=ck.pos, family=gaussian, na.action = na.exclude)
# tb.pos.L <- glm(lnum ~ year+month+bottom+veg+shore, data=tb.pos, family=gaussian, na.action = na.exclude)
# ch.pos.L <- glm(lnum ~ year+month+bottom+veg+shore, data=ch.pos, family=gaussian, na.action = na.exclude)
# jx.pos.L <- glm(lnum ~ year+month+bottom+veg+shore, data=jx.pos, family=gaussian, na.action = na.exclude)
# ir.pos.L <- glm(lnum ~ year+month+bottom+veg+shore, data=ir.pos, family=gaussian, na.action = na.exclude)
# 
# #3. Test the Poisson GLMs for overdispersion ####
# #manual way for count models 
# # Pearson chi1/residual deviance
# sum((resid(ap.pos.P,type="pearson"))^2)/ap.pos.P$df.resid #9.37
# sum((resid(ck.pos.P,type="pearson"))^2)/ck.pos.P$df.resid #6.30
# sum((resid(tb.pos.P,type="pearson"))^2)/tb.pos.P$df.resid #14.26
# sum((resid(ch.pos.P,type="pearson"))^2)/ch.pos.P$df.resid #8.12
# sum((resid(jx.pos.P,type="pearson"))^2)/jx.pos.P$df.resid #5.79
# sum((resid(ir.pos.P,type="pearson"))^2)/ir.pos.P$df.resid #18.46

#3A. Address overdispersion with different model types ####
#there is evidence of overdispersion for every bay (p values were less than) so use quasipoisson for Positive models (Zuur pg 226)
# #decided to avoid QP and opted just for Negative binomial below
# ap.pos.QP <- glm(number ~ year +month+veg+bottom+shore, data=ap.pos, family=quasipoisson, na.action = na.exclude)
# ck.pos.QP <- glm(number ~ year +month+veg+bottom+shore, data=ck.pos, family=quasipoisson, na.action = na.exclude)
# tb.pos.QP <- glm(number ~ year +month+veg+bottom+shore, data=tb.pos, family=quasipoisson, na.action = na.exclude)
# ch.pos.QP <- glm(number ~ year +month+veg+bottom+shore, data=ch.pos, family=quasipoisson, na.action = na.exclude)
# jx.pos.QP <- glm(number ~ year +month+veg+bottom+shore, data=jx.pos, family=quasipoisson, na.action = na.exclude)
# ir.pos.QP <- glm(number ~ year +month+veg+bottom+shore, data=ir.pos, family=quasipoisson, na.action = na.exclude)

# OR, use negative binomial for the positives 
# ap.pos.NB <- glm.nb(number ~ year +month+veg+bottom+shore, data=ap.pos, na.action = na.exclude)
# ck.pos.NB <- glm.nb(number ~ year +month+veg+bottom+shore, data=ck.pos, na.action = na.exclude)
# tb.pos.NB <- glm.nb(number ~ year +month+veg+bottom+shore, data=tb.pos, na.action = na.exclude)
# ch.pos.NB <- glm.nb(number ~ year +month+veg+bottom+shore, data=ch.pos, na.action = na.exclude)
# jx.pos.NB <- glm.nb(number ~ year +month+veg+bottom+shore, data=jx.pos, na.action = na.exclude)
# ir.pos.NB <- glm.nb(number ~ year +month+veg+bottom+shore, data=ir.pos, na.action = na.exclude)


#4. Mixture models. Build the poisson, negative binomial with the full datasets ####
apP <- glm(number ~ year +month+veg+bottom+shore, data=ap.fl, family=poisson, na.action = na.exclude)
ckP <- glm(number ~ year +month+veg+bottom+shore, data=ck.fl, family=poisson,  na.action = na.exclude)
tbP <- glm(number ~ year +month+veg+bottom+shore, data=tb.fl, family=poisson,  na.action = na.exclude)
chP <- glm(number ~ year +month+veg+bottom+shore, data=ch.fl, family=poisson, na.action = na.exclude)
jxP <- glm(number ~ year +month+veg+bottom+shore, data=jx.fl, family=poisson, na.action = na.exclude)
irP <- glm(number ~ year +month+veg+bottom+shore, data=ir.fl, family=poisson, na.action = na.exclude)

apNB <- glm.nb(number ~ year +month+veg+bottom+shore, data=ap.fl, na.action = na.exclude)
ckNB <- glm.nb(number ~ year +month+veg+bottom+shore, data=ck.fl,  na.action = na.exclude)
tbNB <- glm.nb(number ~ year +month+veg+bottom+shore, data=tb.fl,  na.action = na.exclude)
chNB <- glm.nb(number ~ year +month+veg+bottom+shore, data=ch.fl,  na.action = na.exclude)
jxNB <- glm.nb(number ~ year +month+veg+bottom+shore, data=jx.fl,  na.action = na.exclude)
irNB <- glm.nb(number ~ year +month+veg+bottom+shore, data=ir.fl,  na.action = na.exclude)

summary(apNB)
summary(ckNB)
summary(tbNB)
summary(chNB)
summary(jxNB)
summary(irNB)

#odds test shows negative binomial better than poisson
odTest(apNB)
odTest(ckNB)
odTest(tbNB)
odTest(chNB)
odTest(jxNB)
odTest(irNB)


#check for basic model fit if it is more that around 1 then there may be excess zeros 
#see  page 293 of Analysis of Categorical data with R for equation
# if returns FALSE then we are all good
# apNB$deviance/apNB$df.residual> 1+2*(sqrt(2/(apNB$df.residual)))
# ckNB$deviance/ckNB$df.residual > 1+2*(sqrt(2/(ckNB$df.residual)))
# tbNB$deviance/tbNB$df.residual >1+2*(sqrt(2/(tbNB$df.residual)))
# chNB$deviance/chNB$df.residual >1+2*(sqrt(2/(chNB$df.residual)))
# jxNB$deviance/jxNB$df.residual >1+2*(sqrt(2/(jxNB$df.residual)))
# irNB$deviance/irNB$df.residual >1+2*(sqrt(2/(irNB$df.residual)))
# * Hilbe says we cant use the deviance statistic on count data and its more appropriate to use the pearson stat. 

#check for over-dispersion for full poisson and full NB

sum((resid(apP,type="pearson"))^2)/apP$df.resid #9.80
sum((resid(apNB,type="pearson"))^2)/apNB$df.resid #1.39
sum((resid(apNB,type="pearson"))^2)/apNB$df.resid> 1+2*(sqrt(2/(apNB$df.residual)))
sum((resid(apNB,type="pearson"))^2)/apNB$df.resid> 1+3*(sqrt(2/(apNB$df.residual)))

sum((resid(ckP,type="pearson"))^2)/ckP$df.resid #7.36
sum((resid(ckNB,type="pearson"))^2)/ckNB$df.resid #1.28
sum((resid(ckNB,type="pearson"))^2)/ckNB$df.resid> 1+2*(sqrt(2/(ckNB$df.residual)))
sum((resid(ckNB,type="pearson"))^2)/ckNB$df.resid> 1+3*(sqrt(2/(ckNB$df.residual)))

sum((resid(tbP,type="pearson"))^2)/tbP$df.resid #35.14
sum((resid(tbNB,type="pearson"))^2)/tbNB$df.resid #1.68
sum((resid(tbNB,type="pearson"))^2)/tbNB$df.resid> 1+2*(sqrt(2/(tbNB$df.residual)))
sum((resid(tbNB,type="pearson"))^2)/tbNB$df.resid> 1+3*(sqrt(2/(tbNB$df.residual)))

sum((resid(chP,type="pearson"))^2)/chP$df.resid #16.98
sum((resid(chNB,type="pearson"))^2)/chNB$df.resid #1.52
sum((resid(chNB,type="pearson"))^2)/chNB$df.resid> 1+2*(sqrt(2/(chNB$df.residual)))
sum((resid(chNB,type="pearson"))^2)/chNB$df.resid> 1+3*(sqrt(2/(chNB$df.residual)))

sum((resid(jxP,type="pearson"))^2)/jxP$df.resid #3.95
sum((resid(jxNB,type="pearson"))^2)/jxNB$df.resid #1.34
sum((resid(jxNB,type="pearson"))^2)/jxNB$df.resid> 1+2*(sqrt(2/(jxNB$df.residual)))
sum((resid(jxNB,type="pearson"))^2)/jxNB$df.resid> 1+3*(sqrt(2/(jxNB$df.residual)))

sum((resid(irP,type="pearson"))^2)/irP$df.resid #46
sum((resid(irNB,type="pearson"))^2)/irNB$df.resid #1.93

#Negative binomial reduces dispersion quite a bit but some guidelines indicate poor fit. Lets try model selection and see if that improves things
# before we foray into the ZIP and ZINB models 

# 5. Build the mixture (full) zero inflated poisson (ZIP) and zero inflated negative binomial (ZINB) ####
library(pscl)
#page 279 Zuur

#build additional models used for model selection below 
f1 = formula(number ~ year+month+veg+bottom+shore | year+month+veg+bottom+shore)
# f1C = formula(number ~ year+month+veg+bottom | year+month+veg+bottom)       #drop shore

ap.ZP <- zeroinfl(f1, dist='poisson', link="logit", data=ap.fl, na.action = na.exclude)
ck.ZP <- zeroinfl(f1, dist='poisson', link="logit", data=ck.fl, na.action = na.exclude)
tb.ZP <- zeroinfl(f1, dist='poisson', link="logit", data=tb.fl, na.action = na.exclude)
ch.ZP <- zeroinfl(f1, dist='poisson', link="logit", data=ch.fl, na.action = na.exclude)
jx.ZP <- zeroinfl(f1, dist='poisson', link="logit", data=jx.fl, na.action = na.exclude)
ir.ZP <- zeroinfl(f1, dist='poisson', link="logit", data=ir.fl, na.action = na.exclude)

ap.ZNB <- zeroinfl(f1, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
ck.ZNB <- zeroinfl(f1, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
tb.ZNB <- zeroinfl(f1, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude)

#with f1 I was getting an error abotu the system being computationally singular.
# stack exchange says this is due to the fact that my design matrix is non invertible which is due to
# linearly dependent columns (ie.strongly correlated variables)
#https://stats.stackexchange.com/questions/76488/error-system-is-computationally-singular-when-running-a-glm
#problem fitting ch and jx. address this in model validation section below

ch.ZNB <- zeroinfl(f1C, dist='negbin', link="logit", data=ch.fl, na.action = na.exclude) 
jx.ZNB <- zeroinfl(f1C, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude)
ir.ZNB <- zeroinfl(f1, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude)

lrtest(ap.ZNB, ap.ZP) #ZINB is better for AP than is ZIP
lrtest(ck.ZNB, ck.ZP) #NB is better than ZIP
lrtest(tb.ZNB, tb.ZP) #NB is better than ZIP
lrtest(ch.ZNB, ch.ZP) #NB is better than ZIP
lrtest(jx.ZNB, jx.ZP) #NB is better than ZIP
lrtest(ir.ZNB, ir.ZP) #NB is better than ZIP

# 6. Build the mixture (full) zero altered truncated poisson (ZAP) and zero altered truncated (hurdle) negative binomial(ZANB) #####
#pg 288 Zuur

apZAP <- hurdle(number ~ year+month+veg+bottom+shore, dist="poisson", lin="logit", data=ap.fl)
ckZAP <- hurdle(number ~ year+month+veg+bottom+shore, dist="poisson", lin="logit", data=ck.fl)
tbZAP <- hurdle(number ~ year+month+veg+bottom+shore, dist="poisson", lin="logit", data=tb.fl)
chZAP <- hurdle(number ~ year+month+veg+bottom+shore, dist="poisson", lin="logit", data=ch.fl)
jxZAP <- hurdle(number ~ year+month+veg+bottom+shore, dist="poisson", lin="logit", data=jx.fl)
irZAP <- hurdle(number ~ year+month+veg+bottom+shore, dist="poisson", lin="logit", data=ir.fl)

apZANB <- hurdle(number ~ year+month+veg+bottom+shore, dist="negbin", lin="logit", data=ap.fl)
ckZANB <- hurdle(number ~ year+month+veg+bottom+shore, dist="negbin", lin="logit", data=ck.fl)
tbZANB <- hurdle(number ~ year+month+veg+bottom+shore, dist="negbin", lin="logit", data=tb.fl)
chZANB <- hurdle(number ~ year+month+veg+bottom+shore, dist="negbin", lin="logit", data=ch.fl)
jxZANB <- hurdle(number ~ year+month+veg+bottom+shore, dist="negbin", lin="logit", data=jx.fl)
irZANB <- hurdle(number ~ year+month+veg+bottom+shore, dist="negbin", lin="logit", data=ir.fl)

lrtest(apZAP, apZANB)
lrtest(ckZAP, ckZANB)
lrtest(tbZAP, tbZANB)
lrtest(chZAP, chZANB)
lrtest(jxZAP, jxZANB)
lrtest(irZAP, irZANB)
#chi square says that zero altered NB is better than zero altered Poisson

AIC(ap.ZNB, apZANB) #likely zanb is best
AIC(ck.ZNB, ckZANB) #likely zinb is best
AIC(tb.ZNB, tbZANB) #likely zinb is best
#AIC(ch.ZNB, chZANB) #likely zinb is best
#AIC(jx.ZNB, jxZANB) #likely zinb is best
AIC(ir.ZNB, irZANB) #likely zinb is best

##### MODEL SELECTION_ YOY ######
# Model selection for positives: quasipoisson and lognormal
# Model selection for mixture: poisson and NegativeBi (NB), ZINB

  # For quassipoisson model can use the drop command. 
  # Pages 220 is to 230 in Zuur are helpful for following methods. 
  # The AIC is not defined for quasipoisson models so can't use the step function like what was used in the FWRI code. 
  # Instead, use the drop1 function which is applicable for the quassiPoisson GLM and it is more equivalent to hypothesis testing. (Zuur pg 227)
  # If just using Poisson or Bernoulli (binomial) can use step command but this gives AIC- not deviance. (Zuur pg 253)
  # Explained deviance is nearly the equivalent of R^2 so use this (Zuur pg 218 for equation), "The smaller the residual deviance the better is the model"

# 1. Model selection and validation plots for Delta quasipoisson (QP) ####
# see page 227 in Zuur for full explanation of drop1 and F test for quaispoissona

#AP_POS (Year, Veg, Shore = significant factors)
    summary(ap.pos.QP)
    drop1(ap.pos.QP, test="F")  #model selection in quasipoisson is done using F-ratio (Zuur pg 227)
    # bottom  NS
    
    # drop bottom
    M1_ap.pos <- glm(number ~ year+veg+month+shore, data=ap.pos, family=quasipoisson)
    drop1(M1_ap.pos, test="F")
    summary(M1_ap.pos)
    
# Model validation, page 231 Zuur
    plot(residuals(M1_ap.pos, type="pearson") ~ fitted(M1_ap.pos))
    
    EP <- resid(M1_ap.pos, type="pearson")
    ED <- resid(M1_ap.pos, type ="deviance")
    mu <- predict(M1_ap.pos, type="response")
    E <- ap.pos$number-mu
    EP2 <- E/sqrt(9.422133*mu)
    op <- par(mfrow=c(2,2))
    plot(x=mu, y=E, main="Response residuals")
    plot(x=mu, y=EP, main="Pearson residuals")
    plot(x=mu, Y=EP2, main = "Pearson residual scaled")
    plot(x=mu, y=ED, main = "Deviance residuals")
    par(op) 

#Model validation with boot package
    glm.diag.plots(M1_ap.pos)

#CK_POS (Year, Veg = significant factor)
    summary(ck.pos.QP)
    drop1(ck.pos.QP, test="F")
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
    summary(M3_ck.pos)
    # Year, Veg = significant factors 
    
# Model validation, page 231 Zuur
    EP <- resid(M3_ck.pos, type="pearson")
    ED <- resid(M3_ck.pos, type ="deviance")
    mu <- predict(M3_ck.pos, type="response")
    E <- ck.pos$number-mu
    EP2 <- E/sqrt(6.868047*mu)
    op <- par(mfrow=c(2,2))
    plot(x=mu, y=E, main="Response residuals")
    plot(x=mu, y=EP, main="Pearson residuals")
    plot(x=mu, Y=EP2, main = "Pearson residual scaled")
    plot(x=mu, y=ED, main = "Deviance residuals")
    par(op) 
      
# TB_POS (Year, Veg, Shore = significant factors)
    summary(tb.pos.QP)
    drop1(tb.pos.QP, test="F")
    # bottom and month do not appear significant. Drop all sequentially
    
    # drop month
    M1_tb.pos <- glm(number ~ year+veg+bottom+shore, data=tb.pos, family=quasipoisson)
    drop1(M1_tb.pos, test="F")
    
    # drop month, bottom
    M2_tb.pos <- glm(number ~ year+veg+shore, data=tb.pos, family=quasipoisson)
    drop1(M2_tb.pos, test="F")
    summary(M2_tb.pos)
    # Year, Veg, Shore = significant factors
    
# Model validation, page 231 Zuur
    EP <- resid(M2_tb.pos, type="pearson")
    ED <- resid(M2_tb.pos, type ="deviance")
    mu <- predict(M2_tb.pos, type="response")
    E <- tb.pos$number-mu
    EP2 <- E/sqrt(14.82482*mu)
    op <- par(mfrow=c(2,2))
    plot(x=mu, y=E, main="Response residuals")
    plot(x=mu, y=EP, main="Pearson residuals")
    plot(x=mu, Y=EP2, main = "Pearson residual scaled")
    plot(x=mu, y=ED, main = "Deviance residuals")
    par(op) 

### CH_POS (Year, Month, Veg, Bottom, Shore =significnat factors )
    summary(ch.pos.QP)
    drop1(ch.pos.QP, test="F")
    summary(ch.pos.QP)
    # Year, Month, Veg, Bottom, Shore =significnat factors 
    
# Model validation, page 231 Zuur
    EP <- resid(ch.pos.QP, type="pearson")
    ED <- resid(ch.pos.QP, type ="deviance")
    mu <- predict(ch.pos.QP, type="response")
    E <- ch.pos$number-mu
    EP2 <- E/sqrt(8.48237*mu)
    op <- par(mfrow=c(2,2))
    plot(x=mu, y=E, main="Response residuals")
    plot(x=mu, y=EP, main="Pearson residuals")
    plot(x=mu, Y=EP2, main = "Pearson residual scaled")
    plot(x=mu, y=ED, main = "Deviance residuals")
    par(op)    
    
  
### JX_POS (Year, Veg, Shore = significant factors)
    summary(jx.pos.QP)
    drop1(jx.pos.QP, test="F")
    # bottom and month do not appear significant. Drop all sequentially
    
    #drop month
    M1_jx.pos <- glm(number ~ year+veg+bottom+shore, data=jx.pos, family=quasipoisson)
    drop1(M1_jx.pos, test="F")
    
    #drop month, bottom
    M2_jx.pos <- glm(number ~ year+veg+shore, data=jx.pos, family=quasipoisson)
    drop1(M2_jx.pos, test="F")
    summary(M2_jx.pos)
    # Year, Veg, Shore = significant factors
    
# Model validation, page 231 Zuur
    EP <- resid(M2_jx.pos, type="pearson")
    ED <- resid(M2_jx.pos, type ="deviance")
    mu <- predict(M2_jx.pos, type="response")
    E <- jx.pos$number-mu
    EP2 <- E/sqrt(6.352203*mu)
    op <- par(mfrow=c(2,2))
    plot(x=mu, y=E, main="Response residuals")
    plot(x=mu, y=EP, main="Pearson residuals")
    plot(x=mu, Y=EP2, main = "Pearson residual scaled")
    plot(x=mu, y=ED, main = "Deviance residuals")
    par(op)    
      

### IR_POS (Year, Veg, Bottom = significant factors)
    summary(ir.pos.QP)
    drop1(ir.pos.QP, test="F")
    #month and shore do not appear significant. Drop them all. 
    
    # drop month
    M1_ir.pos <- glm(number ~ year+veg+bottom+shore, data=ir.pos, family=quasipoisson)
    drop1(M1_ir.pos, test="F")
    
    # drop shore
    M2_ir.pos <- glm(number ~ year+veg+bottom, data=ir.pos, family=quasipoisson)
    drop1(M2_ir.pos, test="F")
    summary(M2_ir.pos)
    # Year, Veg, Bottom = significant factors

# Model validation, page 231 Zuur
    EP <- resid(M2_ir.pos, type="pearson")
    ED <- resid(M2_ir.pos, type ="deviance")
    mu <- predict(M2_ir.pos, type="response")
    E <- ir.pos$number-mu
    EP2 <- E/sqrt(21.39846*mu)
    op <- par(mfrow=c(2,2))
    plot(x=mu, y=E, main="Response residuals")
    plot(x=mu, y=EP, main="Pearson residuals")
    plot(x=mu, Y=EP2, main = "Pearson residual scaled")
    plot(x=mu, y=ED, main = "Deviance residuals")
    par(op)  
    
    glm.diag.plots(M2_ir.pos)
    
# 2. Model selection and validation plots for Delta lognormal ####
    #decided this is not a valid method because lognormal is not for count data and it is hard to compare log(number)
    #as a dependent variable as oppsed to number whihc is how all of the other models are modeling the dependent variable
# F test with drop1 command
# Model validation to be done with boot package- glm.diag.plots because I couldnt
# find a good example of model validation plots done in the Zuur book. 
# # AP
#     summary(ap.pos.L)
#     drop1(ap.pos.L, test = "F")
#     #bottom NS
#     
#     M1_ap.pos.L <- glm(lnum ~ year+veg+month+shore, data=ap.pos, family=gaussian)
#     drop1(M1_ap.pos.L, test="F")
#     summary(M1_ap.pos.L)
#     glm.diag.plots(M1_ap.pos.L)
#     
#     plot(residuals(M1_ap.pos.L) ~ fitted(M1_ap.pos.L))
#     
# #CK
#     summary(ck.pos.L)
#     drop1(ck.pos.L, test="F")
#     #month and shore NS
#     
#     M1_ck.pos.L <- glm(lnum ~ year +bottom+veg, data=ck.pos, family=gaussian)
#     drop1(M1_ck.pos.L, test="F")
#     glm.diag.plots(M1_ck.pos.L)
#     
# #TB 
#     summary(tb.pos.L)
#     drop1(tb.pos.L, test="F")
#     #bottom NS
#     
#     M1_tb.pos.L <- glm(lnum ~ year+month+veg+shore, data=tb.pos, family=gaussian)
#     drop1(M1_tb.pos.L, test="F")
#     glm.diag.plots(M1_tb.pos.L)
#     
# #CH
#     summary(ch.pos.L)
#     drop1(ch.pos.L, test="F")
#     glm.diag.plots(ch.pos.L)
#     
# #JX
#     summary(jx.pos.L)
#     drop1(jx.pos.L, test="F")
#     # year, month, bottom, shore NS
#     
#     M1_jx.pos.L <- glm(lnum ~ year+veg, data=jx.pos, family=gaussian)
#     glm.diag.plots(M1_jx.pos.L)
#     
# #IR
#     summary(ir.pos.L)
#     drop1(ir.pos.L, test="F")
#     #year,month,bottom NS
#    
#     M1_ir.pos.L <- glm(lnum ~ year+veg+shore, data=ir.pos, family=gaussian)
#     glm.diag.plots(M1_ir.pos.L) 
# 3. Model selection and validation plots for Delta NB ####
#AP
drop1(ap.pos.NB, test="Chi")
    #bottom NS
ap.pos.NB1 <- glm.nb(number ~ year +month+veg+shore, data=ap.pos, na.action = na.exclude)
lrtest(ap.pos.NB, ap.pos.NB1)

sum((resid(ap.pos.NB1,type="pearson"))^2)/ap.pos.NB1$df.resid #1.59
sum((resid(ap.pos.NB1,type="pearson"))^2)/ap.pos.NB1$df.resid> 1+2*(sqrt(2/(ap.pos.NB1$df.residual)))

#CK
drop1(ck.pos.NB, test="Chi")
#month and shore NS
ck.pos.NB1 <- glm.nb(number ~ year+veg+bottom, data=ck.pos, na.action = na.exclude)
lrtest(ck.pos.NB, ck.pos.NB1)
sum((resid(ck.pos.NB1,type="pearson"))^2)/ck.pos.NB1$df.resid #1.97

#TB
drop1(tb.pos.NB, test="Chi")
#drop bottom
tb.pos.NB1 <- glm.nb(number ~ year +month+veg+shore, data=tb.pos, na.action = na.exclude)
lrtest(tb.pos.NB, tb.pos.NB1)
sum((resid(tb.pos.NB1,type="pearson"))^2)/tb.pos.NB1$df.resid #1.95

#CH
drop1(ch.pos.NB, test="Chi")
sum((resid(ch.pos.NB,type="pearson"))^2)/ch.pos.NB$df.resid #1.83

#JX
drop1(jx.pos.NB, test="Chi")
#drop bottom
jx.pos.NB1 <- glm.nb(number ~ year +month+veg+shore, data=jx.pos, na.action = na.exclude)
lrtest(jx.pos.NB, jx.pos.NB1)
sum((resid(jx.pos.NB1,type="pearson"))^2)/jx.pos.NB1$df.resid #1.71

#IR
drop1(ir.pos.NB, test="Chi")
#drop bottom
ir.pos.NB1 <- glm.nb(number ~ year +month+veg+shore, data=ir.pos, na.action = na.exclude)
drop1(ir.pos.NB1, test="Chi")
lrtest(jx.pos.NB, jx.pos.NB1)
sum((resid(ir.pos.NB1,type="pearson"))^2)/ir.pos.NB1$df.resid #1.71

# 4. Model selection and validation plots for  NB ####
# page 235 Zuur    
#AP
summary(apNB) #Mean Deviance = 0.412
drop1(apNB, test= "Chi")
    #bottom NS 

apNB1 <- glm.nb(number ~ year+month+ veg+shore, data=ap.fl, na.action = na.exclude)
summary(apNB1) 
drop1(apNB1, test="Chi")
#all remaining covariates are significant

lrtest(apNB, apNB1) #equal so drop the bottom
AIC(apNB, apNB1) #AIC doesnt support
#go with AIC because there is really no huge difference in mean deviance
#final apNB
sum((resid(apNB,type="pearson"))^2)/apNB$df.resid #1.39
sum((resid(apNB,type="pearson"))^2)/apNB$df.resid> 1+2*(sqrt(2/(apNB$df.residual)))
sum((resid(apNB,type="pearson"))^2)/apNB$df.resid> 1+3*(sqrt(2/(apNB$df.residual)))
#model selection did not improve fit of NB model to the ap.fl data 

plot(residuals(apNB) ~ fitted(apNB))

# CK
summary(ckNB) #Mean Deviance = 0.43
drop1(ckNB, test= "Chi")
# bottom NS 
ckNB1 <- glm.nb(number ~ year+veg+month+shore, data=ck.fl, na.action = na.exclude)
summary(ckNB1) 

lrtest(ckNB, ckNB1) #equal so go with simpler model ckNB1
AIC(ckNB, ckNB1) #NB1 is better

sum((resid(ckNB1,type="pearson"))^2)/ckNB1$df.resid #1.27
sum((resid(ckNB1,type="pearson"))^2)/ckNB1$df.resid> 1+2*(sqrt(2/(ckNB1$df.residual)))
sum((resid(ckNB1,type="pearson"))^2)/ckNB1$df.resid> 1+3*(sqrt(2/(ckNB1$df.residual)))
#model selection did not improve fit of NB model to the ck.fl data 

#TB
summary(tbNB)
drop1(tbNB, test="Chi")
#bottom NS

tbNB1 <- glm.nb(number ~ year+veg+shore+month, data=tb.fl, na.action = na.exclude)
summary(tbNB1)

lrtest(tbNB, tbNB1) #equal so go with simpler model tbNB1
AIC(tbNB, tbNB1) #AIC says NB1 is better

sum((resid(tbNB1,type="pearson"))^2)/tbNB1$df.resid #1.67
sum((resid(tbNB1,type="pearson"))^2)/tbNB1$df.resid> 1+2*(sqrt(2/(tbNB1$df.residual)))
sum((resid(tbNB1,type="pearson"))^2)/tbNB1$df.resid> 1+3*(sqrt(2/(tbNB1$df.residual)))
#model selection did not improve fit of NB model to the tb.fl data 

#CH
summary(chNB)
drop1(ckNB, test="Chi")
#bottom NS

chNB1 <- glm.nb(number ~ year+veg+shore+month, data=ch.fl, na.action = na.exclude)

lrtest(chNB, chNB1) #different, say chNB is better
AIC(chNB, chNB1) #AIC says chNB is better

sum((resid(chNB,type="pearson"))^2)/chNB$df.resid #1.52
sum((resid(chNB,type="pearson"))^2)/chNB$df.resid> 1+2*(sqrt(2/(chNB$df.residual)))
sum((resid(chNB,type="pearson"))^2)/chNB$df.resid> 1+3*(sqrt(2/(chNB$df.residual)))
#model selection did not improve fit of NB model to the ch.fl data 

#JX
summary(jxNB)
#month, bottom NS
drop1(jxNB, test="Chi")
#drop bottom

jxNB1 <- glm.nb(number ~ year+month+veg+shore, data=jx.fl, na.action = na.exclude)
summary(jxNB1)
lrtest(jxNB, jxNB1) #equal so go with simpler one, jxNB1
AIC(jxNB, jxNB1) # doesnt agree. both look similar

sum((resid(jxNB,type="pearson"))^2)/jxNB$df.resid #1.34
sum((resid(jxNB1,type="pearson"))^2)/jxNB1$df.resid #1.33
sum((resid(jxNB1,type="pearson"))^2)/jxNB1$df.resid> 1+2*(sqrt(2/(jxNB1$df.residual)))
sum((resid(jxNB1,type="pearson"))^2)/jxNB1$df.resid> 1+3*(sqrt(2/(jxNB1$df.residual)))
#model selection did not improve fit of NB model to the jx.fl data 

#IR
summary(irNB)
drop1(irNB, test="Chi")
#bttom NS

irNB1 <- glm.nb(number ~ year+veg+shore+month, data=ir.fl, na.action = na.exclude)
lrtest(irNB, irNB1) #equal so go with simpler
AIC(irNB, irNB1) #AIC agrees that irNB1 is better

sum((resid(irNB1,type="pearson"))^2)/irNB1$df.resid #1.93
sum((resid(irNB1,type="pearson"))^2)/irNB1$df.resid> 1+2*(sqrt(2/(irNB1$df.residual)))
sum((resid(irNB1,type="pearson"))^2)/irNB1$df.resid> 1+3*(sqrt(2/(irNB1$df.residual)))
#model selection did not improve fit of NB model to the ir.fl data 

#because none seemed to really improve fit lets try the zero inflated NB

# 5. Model selection and validation plots for ZINB ####
# Page 280 model selection must be done by hand with ZINB
# Most comprehensive, but takes the most work, is to just drop each term 
    #(count| binomial) 
# Using LRtest defined in code from Shanae
#process to drop each variable is from:
    #StandardizeIndexGLM_adult_6ottertrawl.R which was put together with help from Zuur.

#Function for likelihood ratio test to condense coding
  LLRatioTest <- function(L1,L2,df){
  #L1 is full model, L2 is nested model
  #df is difference in paramters
  L <- 2 * (L1 - L2)
  L <- abs(as.numeric(L))
  p <- 1 - pchisq(L,df)
  round(c(L,p), digits = 4)}

#5A. AP model selection ZINB ####
f1 = formula(number ~ year+month+veg+bottom+shore | year+month+veg+bottom+shore)

#drop year from count
A1= zeroinfl(number ~ month+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop month from count
A2= zeroinfl(number ~ year+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop veg from count
A3= zeroinfl(number ~ year+month+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop bottom from count
A4= zeroinfl(number ~ year+month+veg+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop shore from count
A5= zeroinfl(number ~ year+month+veg+bottom | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)

#drop year from binary
A6= zeroinfl(number ~ year+month+veg+bottom+shore |month+veg+bottom+shore , dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop month from binary
A7= zeroinfl(number ~ year+month+veg+bottom+shore | year+veg+bottom+shore, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop veg from binary
A8= zeroinfl(number ~ year+month+veg+bottom+shore | year+month+bottom+shore, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop bottom from binary
A9= zeroinfl(number ~ year+month+veg+bottom+shore | year+month+veg+shore, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop shore from binary
A10= zeroinfl(number ~ year+month+veg+bottom+shore | year+month+veg+bottom, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)

#Run LL ratio test for each sub-model
#Compare each model to base model; last argument is degrees of freedom, which is the number of levels minus 1 for
#categorical covariates and is 1 for continuous covariates for those that are dropped.
# for example, year is categorical and has 18 levels (For AP) so degrees of freedom for models that drop year from
#either count or binary part of model is 17
Z1 <- LLRatioTest(logLik(ap.ZNB),logLik(A1),17)		#drop year from count part
Z2 <- LLRatioTest(logLik(ap.ZNB),logLik(A2),5)		#drop month from count part
Z3 <- LLRatioTest(logLik(ap.ZNB),logLik(A3),1) #drop veg from count part
Z4 <- LLRatioTest(logLik(ap.ZNB),logLik(A4),1) #drop bottom from count part
Z5 <- LLRatioTest(logLik(ap.ZNB),logLik(A5),2) #drop shore from count part

Z6 <- LLRatioTest(logLik(ap.ZNB),logLik(A6),17)		#drop year from binary
Z7 <- LLRatioTest(logLik(ap.ZNB),logLik(A7),5)		#drop month from binary
Z8 <- LLRatioTest(logLik(ap.ZNB),logLik(A8),1) #drop veg from binary
Z9 <- LLRatioTest(logLik(ap.ZNB),logLik(A9),1) #drop bottom from binary
Z10 <- LLRatioTest(logLik(ap.ZNB),logLik(A10),2) #drop shore from binary

#bottom in log link
#bottom in logistic link
#shore from logistic link

#check with AIC

AIC(ap.ZNB, A1, A2, A3, A4, A5, A6, A8, A9, A10)

#starting model with bottom removed from log link
A <- zeroinfl(number ~ year +month+veg+shore | year+month+bottom+veg+shore, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)

#drop year from count
A1 <- zeroinfl(number ~ month+veg+shore | year+month+bottom+veg+shore, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop month from count
A2 <- zeroinfl(number ~ year+veg+shore | year+month+bottom+veg+shore, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop veg from count
A3 <- zeroinfl(number ~ year +month+shore | year+month+bottom+veg+shore, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop shore from count
A4 <- zeroinfl(number ~ year +month+veg | year+month+bottom+veg+shore, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop year from binary
A5 <- zeroinfl(number ~ year +month+veg+shore |month+bottom+veg+shore, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop month from binary
#A6 <- zeroinfl(number ~ year +month+veg+shore | year+bottom+veg+shore, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop veg from binary
#A7 <- zeroinfl(number ~ year +month+veg+shore | year+month+bottom+shore, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop bottom from binary
A8 <-  zeroinfl(number ~ year +month+veg+shore | year+month+veg+shore, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop shore from binary 
A9 <- zeroinfl(number ~ year +month+veg+shore | year+month+bottom+veg, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)


Z1 <- LLRatioTest(logLik(A),logLik(A1),17)		#drop year from count part
Z2 <- LLRatioTest(logLik(A),logLik(A2),5)		#drop month from count part
Z3 <- LLRatioTest(logLik(A),logLik(A3),1) #drop veg from count part
Z4 <- LLRatioTest(logLik(A),logLik(A4),2) #drop shore from count part

Z5 <- LLRatioTest(logLik(A),logLik(A5),17)		#drop year from binary
#Z6 <- LLRatioTest(logLik(A),logLik(A6),5)		#drop month from binary
#Z7 <- LLRatioTest(logLik(A),logLik(A7),1) #drop veg from binary
Z8 <- LLRatioTest(logLik(A), logLik(A8),1) #drop bottom
Z9 <- LLRatioTest(logLik(A), logLik(A9), 2) #drop shore

#month in log link
#shore in logistic link

#check with AIC
AIC(A, A1, A2, A3, A4, A5, A8, A9)


#starting model with shore in logistic link removed because both upper steps indicate removal
A <- zeroinfl(number ~ year +month+veg+shore | year+month+bottom+veg, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop year
A1 <- zeroinfl(number ~ month+veg+shore | year+month+bottom+veg, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop month
A2 <- zeroinfl(number ~ year+veg+shore | year+month+bottom+veg, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop veg
A3 <- zeroinfl(number ~ year +month+shore | year+month+bottom+veg, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop shore
A4 <- zeroinfl(number ~ year +month+veg | year+month+bottom+veg, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop year binary
A5 <- zeroinfl(number ~ year +month+veg+shore | month+bottom+veg, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop month binary
A6 <- zeroinfl(number ~ year +month+veg+shore | year+bottom+veg, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop bottom binary
A7 <- zeroinfl(number ~ year +month+veg+shore | year+month+veg, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop veg binary
A8 <- zeroinfl(number ~ year +month+veg+shore | year+month+bottom, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)

Z1 <- LLRatioTest(logLik(A),logLik(A1),summary(A1)$df.residual - summary(A)$df.residual)		#drop year from count part
Z2 <- LLRatioTest(logLik(A),logLik(A2),5)		#drop month from count part
Z3 <- LLRatioTest(logLik(A),logLik(A3),1) #drop veg from count part
Z4 <- LLRatioTest(logLik(A),logLik(A4),2) #drop shore from count part
Z5 <- LLRatioTest(logLik(A),logLik(A5), 17) #drop year binary 
Z6 <- LLRatioTest(logLik(A),logLik(A6),5)		#drop month binary
Z7 <- LLRatioTest(logLik(A),logLik(A7),1)		#drop bottom binary 
Z8 <- LLRatioTest(logLik(A), logLik(A8),1) #drop veg from binary

AIC(A, A1, A2, A3, A4, A5, A6, A7, A8)

#all remaining covariates significant; A best model; can be confirmed using AIC

summary(A)	#==>examine standard errors of covariates; unusually large values indicate problem with model

N<-nrow(ap.fl)
EZIP <- resid(A,type="pearson")
Dispersion <- sum(EZIP^2)/(N-52)	#==>be sure to change value for degrees of freedom based on summary output
Dispersion #1.12

#==>Dispersion is 1.12 for model which is an improvement from  1.39 in the NB 
sum((resid(A,type="pearson"))^2)/(N-52) #1.12 
sum((resid(A,type="pearson"))^2)/(N-52) > 1+3*(sqrt(2/(N-52)))
BestModelAP_zinb <- A

#5B. CK model selection ZINB ####

#drop year from count
C1= zeroinfl(number ~ month+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop month from count
#C2= zeroinfl(number ~ year+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop veg from count
#C3= zeroinfl(number ~ year+month+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop bottom from count
C4= zeroinfl(number ~ year+month+veg+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop shore from count
C5= zeroinfl(number ~ year+month+veg+bottom | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)

#drop year from binary
C6= zeroinfl(number ~ year+month+veg+bottom+shore |month+veg+bottom+shore , dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop month from binary
#C7= zeroinfl(number ~ year+month+veg+bottom+shore | year+veg+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop veg from binary
#C8= zeroinfl(number ~ year+month+veg+bottom+shore | year+month+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop bottom from binary
C9= zeroinfl(number ~ year+month+veg+bottom+shore | year+month+veg+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop shore from binary
#C10= zeroinfl(number ~ year+month+veg+bottom+shore | year+month+veg+bottom, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)


#Run LL ratio test for each sub-model
#Compare each model to base model; last argument is degrees of freedom, which is the number of levels minus 1 for
#categorical covariates and is 1 for continuous covariates for those that are dropped.
# for example, year is categorical and has 18 levels (For AP) so degrees of freedom for models that drop year from
#either count or binary part of model is 17
Z1 <- LLRatioTest(logLik(ck.ZNB),logLik(C1),19)		#drop year from count part
Z4 <- LLRatioTest(logLik(ck.ZNB),logLik(C4),1) #drop bottom from count part
Z5 <- LLRatioTest(logLik(ck.ZNB),logLik(C5),1) #drop shore from count part

Z6 <- LLRatioTest(logLik(ck.ZNB),logLik(C6),19)		#drop year from binary
Z9 <- LLRatioTest(logLik(ck.ZNB),logLik(C9),1) #drop bottom from binary


#bottom from count NS
#shore from count NS
#bottom from binary

#first drop bottom from count
C= zeroinfl(number ~ year+ month+veg+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop year
#C1= zeroinfl(number ~  month+veg+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop month
C2= zeroinfl(number ~ year+ veg+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop veg
#C3= zeroinfl(number ~ year+ month+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop shore
#C4= zeroinfl(number ~ year+ month+veg | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop year
C5= zeroinfl(number ~ year+ month+veg+shore | month+veg+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop month
C6= zeroinfl(number ~ year+ month+veg+shore | year+veg+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop veg
#C7= zeroinfl(number ~ year+ month+veg+shore | year+month+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop bottom
#C8= zeroinfl(number ~ year+ month+veg+shore | year+month+veg+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop shore
#C9= zeroinfl(number ~ year+ month+veg+shore | year+month+veg+bottom, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)

Z2 <- LLRatioTest(logLik(C),logLik(C2),6)		#drop month from count part
Z5 <- LLRatioTest(logLik(C),logLik(C5),19)		#drop year from binary
Z6 <- LLRatioTest(logLik(C),logLik(C6),6) #drop month from binary
#no change

#try dropping shore from count
C= zeroinfl(number ~ year+ month+veg +bottom | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop year
C1= zeroinfl(number ~  month+veg +bottom | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop month
C2= zeroinfl(number ~ year+ veg+bottom | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop veg
#C3= zeroinfl(number ~ year+ month+bottom | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop bottom
#C4= zeroinfl(number ~ year+ month+veg | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop year
C5= zeroinfl(number ~ year+ month+veg+bottom | month+veg+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop month
C6= zeroinfl(number ~ year+ month+veg+bottom | year+veg+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop veg
#C7= zeroinfl(number ~ year+ month+veg+bottom | year+month+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop bottom
C8= zeroinfl(number ~ year+ month+veg+bottom | year+month+veg+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop shore 
#C9= zeroinfl(number ~ year+ month+veg+bottom | year+month+veg+bottom, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)

lrtest(C, C1)
lrtest(C, C2)
lrtest(C, C5)
lrtest(C, C6)
lrtest(C, C8) #equal so drop bottom from binary 

#shore from count dropped, bottom from binary dropped
C= zeroinfl(number ~ year+ month+veg +bottom | year+month+veg+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop year
C1= zeroinfl(number ~  month+veg +bottom | year+month+veg+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop month
C2= zeroinfl(number ~ year+ veg +bottom | year+month+veg+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop veg
#C3= zeroinfl(number ~ year+ month +bottom | year+month+veg+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop bottom
#C4= zeroinfl(number ~ year+ month+veg  | year+month+veg+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop year
C5= zeroinfl(number ~ year+ month+veg +bottom | month+veg+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop month
#C6= zeroinfl(number ~ year+ month+veg +bottom | year+veg+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop veg
#C7= zeroinfl(number ~ year+ month+veg +bottom | year+month+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop shore 
#C8= zeroinfl(number ~ year+ month+veg +bottom | year+month+veg, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)

lrtest(C, C1)
lrtest(C, C2)
lrtest(C, C5)

AIC(C, C1, C2, C5, C5A)
# model C appears best

N<-nrow(ck.fl)
EZIP <- resid(C,type="pearson")
Dispersion <- sum(EZIP^2)/(N-57)	#==>be sure to change value for degrees of freedom based on summary output
Dispersion #1.17.. improvement from 1.27 with negative binomial 

BestModelCK_zinb <- C

#5C. TB model selection ZINB ####
#full
t <- zeroinfl(number ~ year+month+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop yea
t1 <- zeroinfl(number ~month+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop month
t2 <- zeroinfl(number ~ year+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop veg
t3 <- zeroinfl(number ~ year+month+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop bottom
t4 <- zeroinfl(number ~ year+month+veg+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop shore
t5 <- zeroinfl(number ~ year+month+veg+bottom | year+month+veg+bottom+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop year binary
t6 <- zeroinfl(number ~ year+month+veg+bottom+shore | month+veg+bottom+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop month binary
#t7 <- zeroinfl(number ~ year+month+veg+bottom+shore | year+veg+bottom+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop veg binary
t8 <- zeroinfl(number ~ year+month+veg+bottom+shore | year+month+bottom+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop bottom binary
t9 <- zeroinfl(number ~ year+month+veg+bottom+shore | year+month+veg+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop shore binary
t10 <- zeroinfl(number ~ year+month+veg+bottom+shore | year+month+veg+bottom, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 

lrtest(t, t1)
lrtest(t, t2)
lrtest(t, t3)
lrtest(t, t4) #equal so drop bottom in count
lrtest(t, t5)
lrtest(t, t6)
lrtest(t, t8) #equal so drop veg from binary
lrtest(t, t9) #equal so drop shore from binary

#drop bottom from count
t <- zeroinfl(number ~ year+month+veg+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop year
t1 <- zeroinfl(number ~ month+veg+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop month
t2 <- zeroinfl(number ~ year+veg+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#dropveg
t3 <- zeroinfl(number ~ year+month+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop shore
t4 <- zeroinfl(number ~ year+month+veg | year+month+veg+bottom+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop year from binary
t5 <- zeroinfl(number ~ year+month+veg+shore | month+veg+bottom+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop month from binary
#t6 <- zeroinfl(number ~ year+month+veg+shore | year+veg+bottom+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop veg from binary
t7 <- zeroinfl(number ~ year+month+veg+shore | year+month+bottom+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop bottom from binary
t8 <- zeroinfl(number ~ year+month+veg+shore | year+month+veg+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop shore from binary
t9 <- zeroinfl(number ~ year+month+veg+shore | year+month+veg+bottom, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 

lrtest(t, t1)
lrtest(t, t2)
lrtest(t, t3)
lrtest(t, t4) 
lrtest(t, t5)
#lrtest(t, t6)
lrtest(t, t7) #equal so drop veg from binary
lrtest(t, t8) #equal so drop shore form binary

#drop veg from binary
t <- zeroinfl(number ~ year+month+veg+shore | year+month+bottom+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop year
t1 <- zeroinfl(number ~ month+veg+shore | year+month+bottom+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop month
t2 <- zeroinfl(number ~ year+veg+shore | year+month+bottom+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop veg
t3 <- zeroinfl(number ~ year+month+shore | year+month+bottom+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop shore
t4 <- zeroinfl(number ~ year+month+veg | year+month+bottom+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop year
t5 <- zeroinfl(number ~ year+month+veg+shore | month+bottom+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop month
t6 <- zeroinfl(number ~ year+month+veg+shore | year+bottom+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop bottom
t7 <- zeroinfl(number ~ year+month+veg+shore | year+month+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop shore 
t8 <- zeroinfl(number ~ year+month+veg+shore | year+month+bottom, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 

lrtest(t, t1)
lrtest(t, t2)
lrtest(t, t3)
lrtest(t, t4) 
lrtest(t, t5)
lrtest(t, t6)
lrtest(t, t7) #equal so drop bottom from binary
lrtest(t, t8) #equal so drop shore from binary

#drop bottom from binary
t <- zeroinfl(number ~ year+month+veg+shore | year+month+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop year
t1 <- zeroinfl(number ~ month+veg+shore | year+month+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop month
t2 <- zeroinfl(number ~ year+veg+shore | year+month+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop veg
t3 <- zeroinfl(number ~ year+month+shore | year+month+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop shore
t4 <- zeroinfl(number ~ year+month+veg | year+month+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop year
t5 <- zeroinfl(number ~ year+month+veg+shore | month+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop month
#t6 <- zeroinfl(number ~ year+month+veg+shore | year+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop shore
t7 <- zeroinfl(number ~ year+month+veg+shore | year+month, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 

lrtest(t, t1)
lrtest(t, t2)
lrtest(t, t3)
lrtest(t, t4) 
lrtest(t, t5)
lrtest(t, t7) #equal so drop shore from binary

#drop shore from binary
t <- zeroinfl(number ~ year+month+veg+shore | year+month, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop year
t1 <- zeroinfl(number ~ month+veg+shore | year+month, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop month
t2 <- zeroinfl(number ~ year+veg+shore | year+month, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop veg
t3 <- zeroinfl(number ~ year+month+shore | year+month, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#srop shore
t4 <- zeroinfl(number ~ year+month+veg | year+month, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop year
t5 <- zeroinfl(number ~ year+month+veg+shore | month, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop month
#t6 <- zeroinfl(number ~ year+month+veg+shore | year, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 

lrtest(t, t1)
lrtest(t, t2)
lrtest(t, t3)
lrtest(t, t4) 
lrtest(t, t5)

AIC(t, t1, t2, t3, t4)
# model t appears best
#All covariates remain significant

N<-nrow(tb.fl)
EZIP <- resid(t,type="pearson")
Dispersion <- sum(EZIP^2)/(N-71)	#==>be sure to change value for degrees of freedom based on summary output
Dispersion #1.19.. improvement from 1.67 with negative binomial 

BestModelTB_zinb <- t

#5D. CH model selection ZINB ####
#ch <- zeroinfl(number ~ year+month+veg+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ch.fl, na.action = na.exclude) 
#doesnt solve

# start with month dropped
ch <- zeroinfl(number ~ year+veg+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ch.fl, na.action = na.exclude) 

ch1 <- zeroinfl(number ~ veg+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ch.fl, na.action = na.exclude) 
ch2 <- zeroinfl(number ~ year+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ch.fl, na.action = na.exclude) 
ch3 <- zeroinfl(number ~ year+veg+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ch.fl, na.action = na.exclude) 
ch4 <- zeroinfl(number ~ year+veg+shore | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ch.fl, na.action = na.exclude) 
ch5 <- zeroinfl(number ~ year+veg+shore+bottom | month+veg+shore+bottom, dist='negbin', link="logit", data=ch.fl, na.action = na.exclude) 
ch6 <- zeroinfl(number ~ year+veg+shore+bottom | year+veg+shore+bottom, dist='negbin', link="logit", data=ch.fl, na.action = na.exclude) 
ch7 <- zeroinfl(number ~ year+veg+shore+bottom | year+month+shore+bottom, dist='negbin', link="logit", data=ch.fl, na.action = na.exclude) 
ch8 <- zeroinfl(number ~ year+veg+shore+bottom | year+month+veg+bottom, dist='negbin', link="logit", data=ch.fl, na.action = na.exclude) 
ch9 <- zeroinfl(number ~ year+veg+shore+bottom | year+month+veg+shore, dist='negbin', link="logit", data=ch.fl, na.action = na.exclude) 

lrtest(ch, ch1)
lrtest(ch, ch2)
lrtest(ch, ch3)
lrtest(ch, ch4)
lrtest(ch, ch5)
lrtest(ch, ch6)
lrtest(ch, ch7) #equal so drop veg from binary
lrtest(ch, ch8)
lrtest(ch, ch9)

# drop veg from binary
ch <- zeroinfl(number ~ year+veg+shore+bottom | year+month+shore+bottom, dist='negbin', link="logit", data=ch.fl, na.action = na.exclude) 

ch1 <- zeroinfl(number ~ veg+shore+bottom | year+month+shore+bottom, dist='negbin', link="logit", data=ch.fl, na.action = na.exclude) 
ch2 <- zeroinfl(number ~ year+shore+bottom | year+month+shore+bottom, dist='negbin', link="logit", data=ch.fl, na.action = na.exclude) 
ch3 <- zeroinfl(number ~ year+veg+bottom | year+month+shore+bottom, dist='negbin', link="logit", data=ch.fl, na.action = na.exclude) 
ch4 <- zeroinfl(number ~ year+veg+shore| year+month+shore+bottom, dist='negbin', link="logit", data=ch.fl, na.action = na.exclude) 
ch5 <- zeroinfl(number ~ year+veg+shore+bottom | month+shore+bottom, dist='negbin', link="logit", data=ch.fl, na.action = na.exclude) 
ch6 <- zeroinfl(number ~ year+veg+shore+bottom | year+shore+bottom, dist='negbin', link="logit", data=ch.fl, na.action = na.exclude) 
ch7 <- zeroinfl(number ~ year+veg+shore+bottom | year+month+bottom, dist='negbin', link="logit", data=ch.fl, na.action = na.exclude) 
ch8 <- zeroinfl(number ~ year+veg+shore+bottom | year+month+shore, dist='negbin', link="logit", data=ch.fl, na.action = na.exclude) 

lrtest(ch, ch1)
lrtest(ch, ch2)
lrtest(ch, ch3)
lrtest(ch, ch4)
lrtest(ch, ch5)
lrtest(ch, ch6)
lrtest(ch, ch7) 
lrtest(ch, ch8)

#model Ch seems to be the best
#All covariates are significant
AIC(ch, ch1, ch2, ch3, ch4, ch5, ch6, ch7, ch8)

N<-nrow(ch.fl)
EZIP <- resid(ch,type="pearson")
Dispersion <- sum(EZIP^2)/(N-65)	#==>be sure to change value for degrees of freedom based on summary output
Dispersion #1.24.. improvement from 1.52 with negative binomial 

BestModelCH_zinb <- ch

#5E. JX model selection ZINB####
j <- zeroinfl(number ~ year+month+veg+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 

#drop year
j1 <- zeroinfl(number ~ month+veg+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop month
#j2 <- zeroinfl(number ~ year+veg+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop veg
#j3 <- zeroinfl(number ~ year+month+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop shore
#j4 <- zeroinfl(number ~ year+month+veg+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop bottom
j5 <- zeroinfl(number ~ year+month+veg+shore | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop year binary
j6 <- zeroinfl(number ~ year+month+veg+shore+bottom | month+veg+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop month binary
j7 <- zeroinfl(number ~ year+month+veg+shore+bottom | year+veg+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop veg binary
j8 <- zeroinfl(number ~ year+month+veg+shore+bottom | year+month+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop shore binary
j9 <- zeroinfl(number ~ year+month+veg+shore+bottom | year+month+veg+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop bottom binary
j10 <- zeroinfl(number ~ year+month+veg+shore+bottom | year+month+veg+shore, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 

lrtest(j, j1)
lrtest(j, j5) #equal so drop bottom from count
lrtest(j, j6)
lrtest(j, j7)
lrtest(j, j8)
lrtest(j, j9) #equal so drop shore from binary
lrtest(j, j10) #equal so drop bottom from binary

#drop bottom from count
j <- zeroinfl(number ~ year+month+veg+shore | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop year
j1 <- zeroinfl(number ~ month+veg+shore | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop month
#j2<- zeroinfl(number ~ year+veg+shore | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop veg
j3 <- zeroinfl(number ~ year+month+shore | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop shore
j4 <- zeroinfl(number ~ year+month+veg | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop year
j5 <- zeroinfl(number ~ year+month+veg+shore | month+veg+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop month
j6 <- zeroinfl(number ~ year+month+veg+shore | year+veg+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop veg
#j7 <- zeroinfl(number ~ year+month+veg+shore | year+month+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop shore
j8 <- zeroinfl(number ~ year+month+veg+shore | year+month+veg+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop bottom
#j9 <- zeroinfl(number ~ year+month+veg+shore | year+month+veg+shore, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 

lrtest(j, j1)
lrtest(j, j3)
lrtest(j, j4)
lrtest(j, j5)
lrtest(j, j6) 
lrtest(j, j8) #equal so drop shore from binary

#drop shore from binary 
j <- zeroinfl(number ~ year+month+veg+shore | year+month+veg+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop year
j1 <- zeroinfl(number ~ month+veg+shore | year+month+veg+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop month
j2 <- zeroinfl(number ~ year+veg+shore | year+month+veg+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop veg
j3 <- zeroinfl(number ~ year+month+shore | year+month+veg+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop shore
j4 <- zeroinfl(number ~ year+month+veg | year+month+veg+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop year
j5 <- zeroinfl(number ~ year+month+veg+shore | month+veg+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop month
#j6 <- zeroinfl(number ~ year+month+veg+shore | year+veg+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop veg
j7 <- zeroinfl(number ~ year+month+veg+shore | year+month+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop bottom
j8 <- zeroinfl(number ~ year+month+veg+shore | year+month+veg, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 

lrtest(j, j1)
lrtest(j, j3)
lrtest(j, j4)
lrtest(j, j5)
lrtest(j, j7) 
lrtest(j, j8) #equal so drop bottom from binary 

#drop bottom from binary
j <- zeroinfl(number ~ year+month+veg+shore | year+month+veg, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop year
j1 <- zeroinfl(number ~ month+veg+shore | year+month+veg, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop month
#j2 <- zeroinfl(number ~ year+veg+shore | year+month+veg, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop veg
j3 <- zeroinfl(number ~ year+month+shore | year+month+veg, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop shore
j4 <- zeroinfl(number ~ year+month+veg | year+month+veg, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop year
#j5 <- zeroinfl(number ~ year+month+veg+shore | month+veg, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop month
j6 <- zeroinfl(number ~ year+month+veg+shore | year+veg, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop veg
j7 <- zeroinfl(number ~ year+month+veg+shore | year+month, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 

lrtest(j, j1)
lrtest(j, j3)
lrtest(j, j4)
lrtest(j, j6)
lrtest(j, j7) 

#all remaining covariates are significant

AIC(j, j1, j2, j3,j4, j5, j6, j7)
#best model is j

N<-nrow(jx.fl)
EZIP <- resid(j,type="pearson")
Dispersion <- sum(EZIP^2)/(N-47)	#==>be sure to change value for degrees of freedom based on summary output
Dispersion #1.07.. improvement from 1.34 with negative binomial 

BestModelJX_zinb <- j

#5F. IR model selection ZINB ####

i <- zeroinfl(number ~ year+month+veg+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
#drop year
i1 <- zeroinfl(number ~ month+veg+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
#drop month
i2 <- zeroinfl(number ~ year+veg+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
#drop veg
i3 <- zeroinfl(number ~ year+month+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
#drop shore
i4 <- zeroinfl(number ~ year+month+veg+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
#drop bottom
i5 <- zeroinfl(number ~ year+month+veg+shore | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
#drop year
i6 <- zeroinfl(number ~ year+month+veg+shore+bottom | month+veg+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
#drop month
i7 <- zeroinfl(number ~ year+month+veg+shore+bottom | year+veg+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
#drop veg
i8 <- zeroinfl(number ~ year+month+veg+shore+bottom | year+month+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
#drop shore
i9 <- zeroinfl(number ~ year+month+veg+shore+bottom | year+month+veg+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
#drop bottom
i10 <- zeroinfl(number ~ year+month+veg+shore+bottom | year+month+veg+shore, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 

lrtest(i, i1)
lrtest(i, i2)
lrtest(i, i3)
lrtest(i, i4) #equal so drop shore in count
lrtest(i, i5) #equal so drop bottom i count
lrtest(i, i6)
lrtest(i, i7)
lrtest(i, i8)
lrtest(i, i9)
lrtest(i, i10) #equal so drop bottom in binary

#drop shore in count
i <- zeroinfl(number ~ year+month+veg+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
#drop year
i1 <- zeroinfl(number ~ month+veg+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
#drop month
i2 <- zeroinfl(number ~ year+veg+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
#drop veg
i3 <- zeroinfl(number ~ year+month+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
#drop bottom
i4 <- zeroinfl(number ~ year+month+veg | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
#drop year
i5 <- zeroinfl(number ~ year+month+veg+bottom | month+veg+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
#drop month
i6 <- zeroinfl(number ~ year+month+veg+bottom | year+veg+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
#drop veg
i7 <- zeroinfl(number ~ year+month+veg+bottom | year+month+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
#drop shore
i8 <- zeroinfl(number ~ year+month+veg+bottom | year+month+veg+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
#drop bottom
i9 <- zeroinfl(number ~ year+month+veg+bottom | year+month+veg+shore, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 

lrtest(i, i1)
lrtest(i, i2)
lrtest(i, i3)
lrtest(i, i4) #equal so drop bottom from count
lrtest(i, i5) 
lrtest(i, i6)
lrtest(i, i7)
lrtest(i, i8)
lrtest(i, i9) #equal so drop bottom from binary

#drop bottom from count
i <- zeroinfl(number ~ year+month+veg | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 

i1 <- zeroinfl(number ~ month+veg | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
i2 <- zeroinfl(number ~ year+veg | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
i3 <- zeroinfl(number ~ year+month | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
i4 <- zeroinfl(number ~ year+month+veg | month+veg+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
i5 <- zeroinfl(number ~ year+month+veg | year+veg+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
i6 <- zeroinfl(number ~ year+month+veg | year+month+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
i7 <- zeroinfl(number ~ year+month+veg | year+month+veg+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
i8 <- zeroinfl(number ~ year+month+veg | year+month+veg+shore, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 

lrtest(i, i1)
lrtest(i, i2)
lrtest(i, i3)
lrtest(i, i4) 
lrtest(i, i5) 
lrtest(i, i6)
lrtest(i, i7)
lrtest(i, i8) #equal so drop bottom from binary

#drop bttom from binary
i <- zeroinfl(number ~ year+month+veg | year+month+veg+shore, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 

i1 <- zeroinfl(number ~ month+veg | year+month+veg+shore, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
i2 <- zeroinfl(number ~ year+veg | year+month+veg+shore, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
i3 <- zeroinfl(number ~ year+month | year+month+veg+shore, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
i4 <- zeroinfl(number ~ year+month+veg | month+veg+shore, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
i5 <- zeroinfl(number ~ year+month+veg | year+veg+shore, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
i6 <- zeroinfl(number ~ year+month+veg | year+month+shore, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
i7 <- zeroinfl(number ~ year+month+veg | year+month+veg, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 

lrtest(i, i1)
lrtest(i, i2)
lrtest(i, i3)
lrtest(i, i4) 
lrtest(i, i5) 
lrtest(i, i6)
lrtest(i, i7)

#all covariates remain significant

AIC(i, i1, i2, i3, i4, i5, i6, i7)
#model i is best

N<-nrow(ir.fl)
EZIP <- resid(i,type="pearson")
Dispersion <- sum(EZIP^2)/(N-70)	#==>be sure to change value for degrees of freedom based on summary output
Dispersion #1.46.. improvement from 1.93 with negative binomial 

BestModelIR_zinb <- i

#https://stats.stackexchange.com/questions/121661/how-to-interpret-anova-output-when-comparing-two-nested-mixed-effect-models
#The log-likelihoods of the two models are almost exactly equal 
#indicating the two models do a similar job of fitting the data. 
#The LRT is telling you that you'd be very likely to observe a test 
#statistic (Chisq) as large as the one reported if the two models 
#provided the same fit. Hence you fail to reject the null 
#hypothesis that the likelihoods of the two models are equivalent.
#therefore, chose the simpler model over the more complex model if the p value is NS and proceed 

#5G. Model validation plots for ZINB ####

#BestModelAP
#BestModelCK
#BestModelTB
#BestModelCH
#BestModelJX
#BestModelIR 

#create index plot of pearson residuals for ZINB models 
plot(residuals(BestModelAP, type="pearson"), ylab="Pearson Residuals")
plot(residuals(BestModelAP, type="pearson") ~ fitted(BestModelAP))

plot(residuals(BestModelCK, type="pearson"), ylab="Pearson Residuals")
plot(residuals(BestModelCK, type="pearson") ~ fitted(BestModelCK))

plot(residuals(BestModelTB, type="pearson"), ylab="Pearson Residuals")
plot(residuals(BestModelTB, type="pearson") ~ fitted(BestModelTB))

plot(residuals(BestModelCH, type="pearson"), ylab="Pearson Residuals")
plot(residuals(BestModelCH, type="pearson") ~ fitted(BestModelCH))

plot(residuals(BestModelJX, type="pearson"), ylab="Pearson Residuals")
plot(residuals(BestModelJX, type="pearson") ~ fitted(BestModelJX))

plot(residuals(BestModelIR, type="pearson"), ylab="Pearson Residuals")
plot(residuals(BestModelIR, type="pearson") ~ fitted(BestModelIR))

AIC(apNB, BestModelAP)

EP <- residuals(ap.ZNB6, type="pearson")
mu <- predict(ap.ZNB6, type="response")
E <- ap.pos$number-mu
#EP2 <- E/sqrt(9.422133*mu)
op <- par(mfrow=c(1,2))
plot(x=mu, y=E, main="Response residuals")
plot(x=mu, y=EP, main="Pearson residuals")
#plot(x=mu, Y=EP2, main = "Pearson residual scaled")
#plot(x=mu, y=ED, main = "Deviance residuals")
par(op) 

plot(predict(ap.ZNB6, type="response"), ap.fl$number)


EP <- residuals(ck.ZNB11, type="pearson")
mu <- predict(ck.ZNB11, type="response")
E <- ck.fl$number-mu
#EP2 <- E/sqrt(9.422133*mu)
op <- par(mfrow=c(1,2))
plot(x=mu, y=E, main="Response residuals")
plot(x=mu, y=EP, main="Pearson residuals")
#plot(x=mu, Y=EP2, main = "Pearson residual scaled")
#plot(x=mu, y=ED, main = "Deviance residuals")
par(op) 



EP <- residuals(tb.ZNB11, type="pearson")
mu <- residuals(tb.ZNB11, type="response")
E <- tb.fl$number-mu
#EP2 <- E/sqrt(9.422133*mu)
op <- par(mfrow=c(1,2))
plot(x=mu, y=E, main="Response residuals")
plot(x=mu, y=EP, main="Pearson residuals")
#plot(x=mu, Y=EP2, main = "Pearson residual scaled")
#plot(x=mu, y=ED, main = "Deviance residuals")
par(op) 

plot(predict(tb.ZNB11), tb.fl$number)

    
#validation plots page 285 
EstPar <- coef(ZNB, type="zero")
Z <- model.matrix(ZNB, model="zero")
g <- Z %*% EstPar
p <- exp(g)/(1+exp(g))

EstPar2 <- coef(ZNB, model="count")
X <- model.matrix(ZNB, model="count")
g <- X %*% EstPar2
mu1 <- exp(g)

mu <- (1-p)*mu1

#variance and pearson residual

K <- ZNB$theta 
VarY <- ((mu^2)/ k+mu)*(1-p) + (mu^2) * (p^2+p)
EP <- (ZNB$number - mu) /sqrt(VarY) #Pearson residuals
plot(x=mu, y=EP, main="Pearson residuals")   


# 6. Model selection and validation plots for ZANB ####
# Like ZINB model selection must be done by hand with ZANB
# Most comprehensive, but takes the most work, is to just drop each term 
#(count| binomial) 
# Using LRtest defined in code from Shanae
#process to drop each variable is from:
#StandardizeIndexGLM_adult_6ottertrawl.R which was put together with help from Zuur.

#6A. AP model selection ZANB ####
A = hurdle(number ~ year+month+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ap.fl, na.action=na.exclude)

#drop year from count
A1= hurdle(number ~ month+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop month from count
A2= hurdle(number ~ year+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop veg from count
A3= hurdle(number ~ year+month+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop bottom from count
A4= hurdle(number ~ year+month+veg+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop shore from count
A5= hurdle(number ~ year+month+veg+bottom | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)

#drop year from binary
A6= hurdle(number ~ year+month+veg+bottom+shore |month+veg+bottom+shore , dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop month from binary
A7= hurdle(number ~ year+month+veg+bottom+shore | year+veg+bottom+shore, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop veg from binary
A8= hurdle(number ~ year+month+veg+bottom+shore | year+month+bottom+shore, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop bottom from binary
A9= hurdle(number ~ year+month+veg+bottom+shore | year+month+veg+shore, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop shore from binary
A10= hurdle(number ~ year+month+veg+bottom+shore | year+month+veg+bottom, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)

#Run LL ratio test for each sub-model
#Compare each model to base model; last argument is degrees of freedom, which is the number of levels minus 1 for
#categorical covariates and is 1 for continuous covariates for those that are dropped.
# for example, year is categorical and has 18 levels (For AP) so degrees of freedom for models that drop year from
#either count or binary part of model is 17
lrtest(A, A1)
lrtest(A, A2)
lrtest(A, A3)
lrtest(A, A4) #equal so drop bottom from count
lrtest(A, A5)
lrtest(A, A6)
lrtest(A, A7)
lrtest(A, A8)
lrtest(A, A9) #equal so drop bottom from binary
lrtest(A, A10) #equal so drop shore from binary

#check with AIC
AIC(ap.ZNB, A1, A2, A3, A4, A5, A6, A8, A9, A10)

#starting model with bottom removed from count
A <- hurdle(number ~ year +month+veg+shore | year+month+bottom+veg+shore, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)

#drop year from count
A1 <- hurdle(number ~ month+veg+shore | year+month+bottom+veg+shore, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop month from count
A2 <- hurdle(number ~ year+veg+shore | year+month+bottom+veg+shore, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop veg from count
A3 <- hurdle(number ~ year +month+shore | year+month+bottom+veg+shore, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop shore from count
A4 <- hurdle(number ~ year +month+veg | year+month+bottom+veg+shore, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop year from binary
A5 <- hurdle(number ~ year +month+veg+shore |month+bottom+veg+shore, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop month from binary
A6 <- hurdle(number ~ year +month+veg+shore | year+bottom+veg+shore, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop veg from binary
A7 <- hurdle(number ~ year +month+veg+shore | year+month+bottom+shore, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop bottom from binary
A8 <-  hurdle(number ~ year +month+veg+shore | year+month+veg+shore, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop shore from binary 
A9 <- hurdle(number ~ year +month+veg+shore | year+month+bottom+veg, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)

lrtest(A, A1)
lrtest(A, A2)
lrtest(A, A3)
lrtest(A, A4) 
lrtest(A, A5)
lrtest(A, A6)
lrtest(A, A7)
lrtest(A, A8) #equal so drop bottom from binary
lrtest(A, A9) #equal so drop shore from binary

#starting model with shore in logistic link removed because both upper steps indicate removal
A <- hurdle(number ~ year +month+veg+shore | year+month+bottom+veg, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop year
A1 <- hurdle(number ~ month+veg+shore | year+month+bottom+veg, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop month
A2 <- hurdle(number ~ year+veg+shore | year+month+bottom+veg, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop veg
A3 <- hurdle(number ~ year +month+shore | year+month+bottom+veg, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop shore
A4 <- hurdle(number ~ year +month+veg | year+month+bottom+veg, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop year binary
A5 <- hurdle(number ~ year +month+veg+shore | month+bottom+veg, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop month binary
A6 <- hurdle(number ~ year +month+veg+shore | year+bottom+veg, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop bottom binary
A7 <- hurdle(number ~ year +month+veg+shore | year+month+veg, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
#drop veg binary
A8 <- hurdle(number ~ year +month+veg+shore | year+month+bottom, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)

lrtest(A, A1)
lrtest(A, A2)
lrtest(A, A3)
lrtest(A, A4) 
lrtest(A, A5)
lrtest(A, A6)
lrtest(A, A7) #equal so drop bottom from binary
lrtest(A, A8)

#starting model with bottom in logistic link removed 
A <- hurdle(number ~ year +month+veg+shore | year+month+veg, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)

A1 <- hurdle(number ~ month+veg+shore | year+month+veg, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
A2 <- hurdle(number ~ year +veg+shore | year+month+veg, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
A3 <- hurdle(number ~ year +month+shore | year+month+veg, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
A4 <- hurdle(number ~ year +month+veg | year+month+veg, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
A5 <- hurdle(number ~ year +month+veg+shore | month+veg, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
A6 <- hurdle(number ~ year +month+veg+shore | year+veg, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
A7 <- hurdle(number ~ year +month+veg+shore | year+month, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)

lrtest(A, A1)
lrtest(A, A2)
lrtest(A, A3)
lrtest(A, A4) 
lrtest(A, A5)
lrtest(A, A6)
lrtest(A, A7) 

#all remaining covariates are significant
AIC(A, A1, A2, A3, A4, A5, A6, A7) #AIC agrees
#model A is best
BestModelAP_zanb<- A
summary(BestModelAP_zanb)	#==>examine standard errors of covariates; unusually large values indicate problem with model

AIC(BestModelAPzinb, BestModelAP_zanb)

N<-nrow(ap.fl)
EZANB <- resid(BestModelAP_zanb,type="pearson")
Dispersion <- sum(EZANB^2)/(N-51)	#==>be sure to change value for degrees of freedom based on summary output
Dispersion #1.02

sum((resid(BestModelAP_zanb,type="pearson"))^2)/(N-51) > 1+2*(sqrt(2/(N-51)))
sum((resid(BestModelAP_zanb,type="pearson"))^2)/(N-51) > 1+3*(sqrt(2/(N-51)))
#good fit on both accounts!!! 

#6B. CK model selection ZANB ####
C= hurdle(number ~ year+ month+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)

#drop year from count
C1= hurdle(number ~ month+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop month from count
C2= hurdle(number ~ year+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop veg from count
C3= hurdle(number ~ year+month+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop bottom from count
C4= hurdle(number ~ year+month+veg+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop shore from count
C5= hurdle(number ~ year+month+veg+bottom | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop year from binary
C6= hurdle(number ~ year+month+veg+bottom+shore |month+veg+bottom+shore , dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop month from binary
C7= hurdle(number ~ year+month+veg+bottom+shore | year+veg+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop veg from binary
C8= hurdle(number ~ year+month+veg+bottom+shore | year+month+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop bottom from binary
C9= hurdle(number ~ year+month+veg+bottom+shore | year+month+veg+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop shore from binary
C10= hurdle(number ~ year+month+veg+bottom+shore | year+month+veg+bottom, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)

lrtest(C, C1)
lrtest(C, C2)
lrtest(C, C3) #equal so drop veg from count
lrtest(C, C4) #equal so drop bottom from count
lrtest(C, C5) #equal so drop shore from count
lrtest(C, C6)
lrtest(C, C7)
lrtest(C, C8)
lrtest(C, C9) #equal so drop bottom from binary 
lrtest(C, C10)

#first drop bottom from count
C= hurdle(number ~ year+ month+veg+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop year
C1= hurdle(number ~  month+veg+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop month
C2= hurdle(number ~ year+ veg+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop veg
C3= hurdle(number ~ year+ month+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop shore
C4= hurdle(number ~ year+ month+veg | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop year
C5= hurdle(number ~ year+ month+veg+shore | month+veg+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop month
C6= hurdle(number ~ year+ month+veg+shore | year+veg+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop veg
C7= hurdle(number ~ year+ month+veg+shore | year+month+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop bottom
C8= hurdle(number ~ year+ month+veg+shore | year+month+veg+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop shore
C9= hurdle(number ~ year+ month+veg+shore | year+month+veg+bottom, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)

lrtest(C, C1)
lrtest(C, C2)
lrtest(C, C3) 
lrtest(C, C4) #equal so drop shore from count
lrtest(C, C5) 
lrtest(C, C6)
lrtest(C, C7)
lrtest(C, C8) #equal so drop bottom from binary
lrtest(C, C9) 

#dropping shore from count
C= hurdle(number ~ year+ month+veg | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop year
C1= hurdle(number ~  month+veg  | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop month
C2= hurdle(number ~ year+ veg| year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop veg
C3= hurdle(number ~ year+ month| year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop year
C4= hurdle(number ~ year+ month+veg | month+veg+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop month
C5= hurdle(number ~ year+ month+veg | year+veg+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop veg
C6= hurdle(number ~ year+ month+veg | year+month+bottom+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop bottom
C7= hurdle(number ~ year+ month+veg| year+month+veg+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop shore 
C8= hurdle(number ~ year+ month+veg | year+month+veg+bottom, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)

lrtest(C, C1)
lrtest(C, C2)
lrtest(C, C3) 
lrtest(C, C4) 
lrtest(C, C5) 
lrtest(C, C6)
lrtest(C, C7) #equal so drop bottom from binary
lrtest(C, C8) 

#dropping bottom from binary
C= hurdle(number ~ year+ month+veg | year+month+veg+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop year
C1= hurdle(number ~  month+veg  | year+month+veg+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop month
C2= hurdle(number ~ year+ veg| year+month+veg+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop veg
C3= hurdle(number ~ year+ month| year+month+veg+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop year
C4= hurdle(number ~ year+ month+veg | month+veg+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop month
C5= hurdle(number ~ year+ month+veg | year+veg+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop veg
C6= hurdle(number ~ year+ month+veg | year+month+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
#drop shore 
C7= hurdle(number ~ year+ month+veg | year+month+veg, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)

lrtest(C, C1)
lrtest(C, C2)
lrtest(C, C3) 
lrtest(C, C4) 
lrtest(C, C5) 
lrtest(C, C6)
lrtest(C, C7) 

#all remaining covariates are significant
AIC(C, C1, C2, C3, C4, C5, C6, C7) #best model is C
BestModelCK_zanb <- C

N<-nrow(ck.fl)
EZANB <- resid(BestModelCK_zanb,type="pearson")
Dispersion <- sum(EZANB^2)/(N-56)	#==>be sure to change value for degrees of freedom based on summary output
Dispersion #1.04

sum((resid(BestModelCK_zanb,type="pearson"))^2)/(N-56) > 1+2*(sqrt(2/(N-56)))
sum((resid(BestModelCK_zanb,type="pearson"))^2)/(N-56) > 1+3*(sqrt(2/(N-56)))
#good fit on both accounts!!! 

#6C. TB model selection ZANB ####

#full
t <- hurdle(number ~ year+month+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop yea
t1 <- hurdle(number ~month+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop month
t2 <- hurdle(number ~ year+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop veg
t3 <- hurdle(number ~ year+month+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop bottom
t4 <- hurdle(number ~ year+month+veg+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop shore
t5 <- hurdle(number ~ year+month+veg+bottom | year+month+veg+bottom+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop year binary
t6 <- hurdle(number ~ year+month+veg+bottom+shore | month+veg+bottom+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop month binary
t7 <- hurdle(number ~ year+month+veg+bottom+shore | year+veg+bottom+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop veg binary
t8 <- hurdle(number ~ year+month+veg+bottom+shore | year+month+bottom+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop bottom binary
t9 <- hurdle(number ~ year+month+veg+bottom+shore | year+month+veg+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop shore binary
t10 <- hurdle(number ~ year+month+veg+bottom+shore | year+month+veg+bottom, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 

lrtest(t, t1)
lrtest(t, t2)
lrtest(t, t3)
lrtest(t, t4) #equal so drop bottom from count
lrtest(t, t5)
lrtest(t, t6)
lrtest(t, t7)
lrtest(t, t8)
lrtest(t, t9) #equal so drop bottom from binary
lrtest(t, t10)

#drop bottom from count
t <- hurdle(number ~ year+month+veg+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop year
t1 <- hurdle(number ~ month+veg+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop month
t2 <- hurdle(number ~ year+veg+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#dropveg
t3 <- hurdle(number ~ year+month+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop shore
t4 <- hurdle(number ~ year+month+veg | year+month+veg+bottom+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop year from binary
t5 <- hurdle(number ~ year+month+veg+shore | month+veg+bottom+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop month from binary
t6 <- hurdle(number ~ year+month+veg+shore | year+veg+bottom+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop veg from binary
t7 <- hurdle(number ~ year+month+veg+shore | year+month+bottom+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop bottom from binary
t8 <- hurdle(number ~ year+month+veg+shore | year+month+veg+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop shore from binary
t9 <- hurdle(number ~ year+month+veg+shore | year+month+veg+bottom, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 

lrtest(t, t1)
lrtest(t, t2)
lrtest(t, t3)
lrtest(t, t4)
lrtest(t, t5)
lrtest(t, t6)
lrtest(t, t7)
lrtest(t, t8) #equal so drop bottom from binary
lrtest(t, t9) 

#drop bottom from binary
t <- hurdle(number ~ year+month+veg+shore | year+month+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop year
t1 <- hurdle(number ~ month+veg+shore | year+month+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop month
t2 <- hurdle(number ~ year+veg+shore | year+month+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop veg
t3 <- hurdle(number ~ year+month+shore | year+month+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop shore
t4 <- hurdle(number ~ year+month+veg | year+month+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop year
t5 <- hurdle(number ~ year+month+veg+shore | month+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop month
#t6 <- hurdle(number ~ year+month+veg+shore | year+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
#drop shore
t7 <- hurdle(number ~ year+month+veg+shore | year+month, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 

lrtest(t, t1)
lrtest(t, t2)
lrtest(t, t3)
lrtest(t, t4)
lrtest(t, t5)
lrtest(t, t6)
lrtest(t, t7)

#all remaining covariates are significant
AIC(t, t1, t2, t3, t4, t5, t6, t7) #AIC says best model is t
BestModelTB_zanb <- t

N<-nrow(tb.fl)
EZANB <- resid(BestModelTB_zanb,type="pearson")
Dispersion <- sum(EZANB^2)/(N-74)	#==>be sure to change value for degrees of freedom based on summary output
Dispersion #1.12

sum((resid(BestModelTB_zanb,type="pearson"))^2)/(N-74) > 1+2*(sqrt(2/(N-74)))
sum((resid(BestModelTB_zanb,type="pearson"))^2)/(N-74) > 1+3*(sqrt(2/(N-74)))
#not great

#6D. CH model selection ZANB ####

ch <- hurdle(number ~ year+month+veg+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ch.fl, na.action = na.exclude) 
#drop year
ch1 <- hurdle(number ~ month+veg+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ch.fl, na.action = na.exclude) 
#drop month
ch2 <- hurdle(number ~ year+veg+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ch.fl, na.action = na.exclude) 
#drop veg
ch3 <- hurdle(number ~ year+month+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ch.fl, na.action = na.exclude) 
#drop shore
ch4 <- hurdle(number ~ year+month+veg+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ch.fl, na.action = na.exclude) 
#drop bottom
ch5 <- hurdle(number ~ year+month+veg+shore | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ch.fl, na.action = na.exclude) 
#drop year binary
ch6 <- hurdle(number ~ year+month+veg+shore+bottom | month+veg+shore+bottom, dist='negbin', link="logit", data=ch.fl, na.action = na.exclude) 
#drop month binary
ch7 <- hurdle(number ~ year+month+veg+shore+bottom | year+veg+shore+bottom, dist='negbin', link="logit", data=ch.fl, na.action = na.exclude) 
#drop veg binary
ch8 <- hurdle(number ~ year+month+veg+shore+bottom | year+month+shore+bottom, dist='negbin', link="logit", data=ch.fl, na.action = na.exclude) 
#drop shore binary
ch9 <- hurdle(number ~ year+month+veg+shore+bottom | year+month+veg+bottom, dist='negbin', link="logit", data=ch.fl, na.action = na.exclude) 
#drop bottom binary
ch10 <- hurdle(number ~ year+month+veg+shore+bottom | year+month+veg+shore, dist='negbin', link="logit", data=ch.fl, na.action = na.exclude) 

lrtest(ch, ch1)
lrtest(ch, ch2)
lrtest(ch, ch3)
lrtest(ch, ch4)
lrtest(ch, ch5)
lrtest(ch, ch6)
lrtest(ch, ch7)
lrtest(ch, ch8)
lrtest(ch, ch9)
lrtest(ch, ch10) #drop bottom from binary

#drop bottom from binary
ch <- hurdle(number ~ year+month+veg+shore+bottom | year+month+veg+shore, dist='negbin', link="logit", data=ch.fl, na.action = na.exclude) 
#drop year
ch1 <- hurdle(number ~ month+veg+shore+bottom | year+month+veg+shore, dist='negbin', link="logit", data=ch.fl, na.action = na.exclude) 
#drop month
ch2 <- hurdle(number ~ year+veg+shore+bottom | year+month+veg+shore, dist='negbin', link="logit", data=ch.fl, na.action = na.exclude) 
#drop veg
ch3 <- hurdle(number ~ year+month+shore+bottom | year+month+veg+shore, dist='negbin', link="logit", data=ch.fl, na.action = na.exclude) 
#drop shore
ch4 <- hurdle(number ~ year+month+veg+bottom | year+month+veg+shore, dist='negbin', link="logit", data=ch.fl, na.action = na.exclude) 
#drop bottom
ch5 <- hurdle(number ~ year+month+veg+shore | year+month+veg+shore, dist='negbin', link="logit", data=ch.fl, na.action = na.exclude) 
#drop year binary
ch6 <- hurdle(number ~ year+month+veg+shore+bottom | month+veg+shore, dist='negbin', link="logit", data=ch.fl, na.action = na.exclude) 
#drop month binary
ch7 <- hurdle(number ~ year+month+veg+shore+bottom | year+veg+shore, dist='negbin', link="logit", data=ch.fl, na.action = na.exclude) 
#drop veg binary
ch8 <- hurdle(number ~ year+month+veg+shore+bottom | year+month+shore, dist='negbin', link="logit", data=ch.fl, na.action = na.exclude) 
#drop shore binary
ch9 <- hurdle(number ~ year+month+veg+shore+bottom | year+month+veg, dist='negbin', link="logit", data=ch.fl, na.action = na.exclude) 

lrtest(ch, ch1)
lrtest(ch, ch2)
lrtest(ch, ch3)
lrtest(ch, ch4)
lrtest(ch, ch5)
lrtest(ch, ch6)
lrtest(ch, ch7)
lrtest(ch, ch8)
lrtest(ch, ch9)

#all remaining covariates are significant
AIC(ch, ch1, ch2, ch3, ch4, ch5, ch6, ch7, ch8, ch9) #AIC says ch is best

BestModelCH_zanb <- ch

N<-nrow(ch.fl)
EZANB <- resid(BestModelCH_zanb,type="pearson")
Dispersion <- sum(EZANB^2)/(N-76)	#==>be sure to change value for degrees of freedom based on summary output
Dispersion #1.18

sum((resid(BestModelCH_zanb,type="pearson"))^2)/(N-76) > 1+2*(sqrt(2/(N-76)))
sum((resid(BestModelCH_zanb,type="pearson"))^2)/(N-76) > 1+3*(sqrt(2/(N-76)))
#not great

#6E. JX model selection ZANB ####
j <- hurdle(number ~ year+month+veg+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop year
j1 <- hurdle(number ~ month+veg+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop month
j2 <- hurdle(number ~ year+veg+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop veg
j3 <- hurdle(number ~ year+month+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop shore
j4 <- hurdle(number ~ year+month+veg+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop bottom
j5 <- hurdle(number ~ year+month+veg+shore | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop year binary
j6 <- hurdle(number ~ year+month+veg+shore+bottom | month+veg+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop month binary
j7 <- hurdle(number ~ year+month+veg+shore+bottom | year+veg+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop veg binary
j8 <- hurdle(number ~ year+month+veg+shore+bottom | year+month+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop shore binary
j9 <- hurdle(number ~ year+month+veg+shore+bottom | year+month+veg+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop bottom binary
j10 <- hurdle(number ~ year+month+veg+shore+bottom | year+month+veg+shore, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 

lrtest(j, j1) #equal so drop year in count
lrtest(j, j2) 
lrtest(j, j3)
lrtest(j, j4) #equal so drop shore in count
lrtest(j, j5) #equal so frop bottom in count
lrtest(j, j6) 
lrtest(j, j7)
lrtest(j, j8)
lrtest(j, j9)
lrtest(j, j10)

#dropped bottom in count
j <- hurdle(number ~ year+month+veg+shore | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop year
j1 <- hurdle(number ~ month+veg+shore | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop month
j2<- hurdle(number ~ year+veg+shore | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop veg
j3 <- hurdle(number ~ year+month+shore | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop shore
j4 <- hurdle(number ~ year+month+veg | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop year
j5 <- hurdle(number ~ year+month+veg+shore | month+veg+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop month
j6 <- hurdle(number ~ year+month+veg+shore | year+veg+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop veg
j7 <- hurdle(number ~ year+month+veg+shore | year+month+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop shore
j8 <- hurdle(number ~ year+month+veg+shore | year+month+veg+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
#drop bottom
j9 <- hurdle(number ~ year+month+veg+shore | year+month+veg+shore, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 

lrtest(j, j1)#equal so drop year count
lrtest(j, j2)
lrtest(j, j3)
lrtest(j, j4) #equal so drop shore from count
lrtest(j, j5)
lrtest(j, j6)
lrtest(j, j7)
lrtest(j, j8) 
lrtest(j, j9) 

#drop year from count
j <- hurdle(number ~ month+veg+shore | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
j1 <- hurdle(number ~ veg+shore | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
j2 <- hurdle(number ~ month+shore | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
j3 <- hurdle(number ~ month+veg | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
j4 <- hurdle(number ~ month+veg+shore | month+veg+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
j5 <- hurdle(number ~ month+veg+shore | year+veg+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
j6 <- hurdle(number ~ month+veg+shore | year+month+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
j7 <- hurdle(number ~ month+veg+shore | year+month+veg+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
j8 <- hurdle(number ~ month+veg+shore | year+month+veg+shore, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 

lrtest(j, j1)
lrtest(j, j2)
lrtest(j, j3) #equal so drop shore from count
lrtest(j, j4) 
lrtest(j, j5)
lrtest(j, j6)
lrtest(j, j7)
lrtest(j, j8) 

#drop shore from count
j <- hurdle(number ~ veg +month | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
j1 <- hurdle(number ~ month | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
j2 <- hurdle(number ~ veg | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
j3 <- hurdle(number ~ month+veg | month+veg+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
j4 <- hurdle(number ~ month+veg | year+veg+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
j5 <- hurdle(number ~ month+veg | year+month+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
j6 <- hurdle(number ~ month+veg | year+month+veg+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
j7 <- hurdle(number ~ month+veg | year+month+veg+shore, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 

lrtest(j, j1)
lrtest(j, j2)
lrtest(j, j3) 
lrtest(j, j4) 
lrtest(j, j5)
lrtest(j, j6)
lrtest(j, j7)

#all remaining covariates are significant
AIC(j, j1, j2, j3, j4, j5,j6, j7, j8, j9) #AIC says j is best

BestModelJX_zanb <- j

N<-nrow(jx.fl)
EZANB <- resid(BestModelJX_zanb,type="pearson")
Dispersion <- sum(EZANB^2)/(N-34)	#==>be sure to change value for degrees of freedom based on summary output
Dispersion #0.97

sum((resid(BestModelJX_zanb,type="pearson"))^2)/(N-34) > 1+2*(sqrt(2/(N-34)))
sum((resid(BestModelJX_zanb,type="pearson"))^2)/(N-34) > 1+3*(sqrt(2/(N-34)))
#good fit! 

#6F. IR model selection ZANB ####
i <- hurdle(number ~ year+month+veg+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
#drop year
i1 <- hurdle(number ~ month+veg+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
#drop month
i2 <- hurdle(number ~ year+veg+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
#drop veg
i3 <- hurdle(number ~ year+month+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
#drop shore
i4 <- hurdle(number ~ year+month+veg+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
#drop bottom
i5 <- hurdle(number ~ year+month+veg+shore | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
#drop year
i6 <- hurdle(number ~ year+month+veg+shore+bottom | month+veg+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
#drop month
i7 <- hurdle(number ~ year+month+veg+shore+bottom | year+veg+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
#drop veg
i8 <- hurdle(number ~ year+month+veg+shore+bottom | year+month+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
#drop shore
i9 <- hurdle(number ~ year+month+veg+shore+bottom | year+month+veg+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
#drop bottom
i10 <- hurdle(number ~ year+month+veg+shore+bottom | year+month+veg+shore, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 

lrtest(i, i1)
lrtest(i, i2)
lrtest(i, i3)
lrtest(i, i4) #equal so drop shore in count
lrtest(i, i5) #equal so drop bottom i count
lrtest(i, i6)
lrtest(i, i7)
lrtest(i, i8)
lrtest(i, i9)
lrtest(i, i10) 

#drop shore in count
i <- hurdle(number ~ year+month+veg+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
#drop year
i1 <- hurdle(number ~ month+veg+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
#drop month
i2 <- hurdle(number ~ year+veg+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
#drop veg
i3 <- hurdle(number ~ year+month+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
#drop bottom
i4 <- hurdle(number ~ year+month+veg | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
#drop year
i5 <- hurdle(number ~ year+month+veg+bottom | month+veg+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
#drop month
i6 <- hurdle(number ~ year+month+veg+bottom | year+veg+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
#drop veg
i7 <- hurdle(number ~ year+month+veg+bottom | year+month+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
#drop shore
i8 <- hurdle(number ~ year+month+veg+bottom | year+month+veg+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
#drop bottom
i9 <- hurdle(number ~ year+month+veg+bottom | year+month+veg+shore, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 

lrtest(i, i1)
lrtest(i, i2)
lrtest(i, i3)
lrtest(i, i4) #equal so drop bottom from count
lrtest(i, i5) 
lrtest(i, i6)
lrtest(i, i7)
lrtest(i, i8)
lrtest(i, i9) 

#drop bottom from count
i <- hurdle(number ~ year+month+veg | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 

i1 <- hurdle(number ~ month+veg | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
i2 <- hurdle(number ~ year+veg | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
i3 <- hurdle(number ~ year+month | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
i4 <- hurdle(number ~ year+month+veg | month+veg+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
i5 <- hurdle(number ~ year+month+veg | year+veg+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
i6 <- hurdle(number ~ year+month+veg | year+month+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
i7 <- hurdle(number ~ year+month+veg | year+month+veg+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 
i8 <- hurdle(number ~ year+month+veg | year+month+veg+shore, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 

lrtest(i, i1)
lrtest(i, i2)
lrtest(i, i3)
lrtest(i, i4) 
lrtest(i, i5) 
lrtest(i, i6)
lrtest(i, i7)
lrtest(i, i8) 

#all remaining covariates are significant
AIC(i, i1, i2, i3, i4, i5,i6, i7, i8, i9) #AIC says i is best

BestModelIR_zanb <- i

N<-nrow(ir.fl)
EZANB <- resid(BestModelIR_zanb,type="pearson")
Dispersion <- sum(EZANB^2)/(N-71)	#==>be sure to change value for degrees of freedom based on summary output
Dispersion #1.26

sum((resid(BestModelJX_zanb,type="pearson"))^2)/(N-71) > 1+2*(sqrt(2/(N-71)))
sum((resid(BestModelJX_zanb,type="pearson"))^2)/(N-71) > 1+3*(sqrt(2/(N-71)))
#good fit! 

#7. Compare ZINB to ZANB ####
AIC(BestModelAP, BestModelAP_zanb) #ZANB is better
AIC(BestModelCK, BestModelCK_zanb) #ZINB is better
AIC(BestModelTB, BestModelTB_zanb) #ZINB is better
AIC(BestModelCH, BestModelCH_zanb) #ZINB is better
AIC(BestModelJX, BestModelJX_zanb) #ZINB is better
AIC(BestModelIR, BestModelIR_zanb) #ZINB is better

#for every case except for AP, the zero inflated model is better. 

##### MODEL SELECTION BINARY w/ DROP1 command_YOY ######
# # pg 253 Zuur
# # Even though I decided that for this case not to use delta method I will fit the binary model still  
# ##  AP_BIN (Year, Veg = significant)
# summary(Full_ap.bin)
# drop1(Full_ap.bin, test ='Chi')
# # bottom, and shore are not significant so drop them one at a time
# 
# #drop bottom
# M1_ap.bin <- glm(number ~ year+month+veg+shore, data=ap.bin, family=binomial)
# drop1(M1_ap.bin, test ="Chi")
# 
# #drop shore
# M2_ap.bin <- glm(number ~ year+month+veg, data=ap.bin, family=binomial)
# drop1(M2_ap.bin, test ="Chi")
# 
# lrtest(M1_ap.bin, M2_ap.bin) #equal so opt for simpler
# AIC(Full_ap.bin,M1_ap.bin, M2_ap.bin) #AIC agrees
# #M3_ap.bin is best
# 
# ## CK_BIN (Year, Month, Veg, Shore = significant)
# summary(Full_ck.bin)
# drop1(Full_ck.bin, test ='Chi')
# # bottom  is not signficiant
# 
# #drop bottom 
# M1_ck.bin <- glm(number ~ year+month+veg+shore, data=ck.bin, family=binomial, na.action=na.exclude)
# drop1(M1_ck.bin, test ="Chi")
# lrtest(Full_ck.bin, M1_ck.bin)
# 
# ## TB_BIN (Year, Month, Veg, Shore = significant)
# summary(Full_tb.bin)
# drop1(Full_tb.bin, test ='Chi')
# # bottom is not signficiant
# 
# #drop bottom 
# M1_tb.bin <- glm(number ~ year+month+veg+shore, data=tb.bin, family=binomial)
# drop1(M1_tb.bin, test ="Chi")
# lrtest(Full_tb.bin, M1_tb.bin) #equal so go with simpler
# 
# ##  CH_BIN (Year, Month, Veg, Shore = significant)
# summary(Full_ch.bin)
# drop1(Full_ch.bin, test ='Chi')
# # bottom is not signficiant
# 
# #drop bottom 
# M1_ch.bin <- glm(number ~ year+month+veg+shore, data=ch.bin, family=binomial)
# drop1(M1_ch.bin, test ="Chi")
# lrtest(Full_ch.bin, M1_ch.bin) #equal so go with simpler
# 
# ## JX_BIN (Year, Bottom, Veg, Shore +Month =signficant)
# summary(Full_jx.bin)
# drop1(Full_jx.bin, test ='Chi')
# 
# 
# ## IR_BIN (Year, Month, Veg, Shore +bottom = significant)
# summary(Full_ir.bin)
# drop1(Full_ir.bin, test ='Chi')
# 
# #drop bottom 
# M1_ir.bin <- glm(number ~ year+month+veg+shore, data=ir.bin, family=binomial)
# drop1(M1_ir.bin, test ="Chi")
# 
##### ASSIGN FINAL MODELS_YOY ###### 

#After running mixture model (Posisson -> NB -> ZINB) and an alias hurdle model (zero altered (ZA) Poisson -> ZANB)
# the best models for each are below: 
#The zero inflated models were best for each estuary except for AP. For that one, the alias two-part ZANB was best. 


A <- zeroinfl(number ~ year +month+veg+shore | year+month+bottom+veg, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
C= zeroinfl(number ~ year+ month+veg +bottom | year+month+veg+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
t <- zeroinfl(number ~ year+month+veg+shore | year+month, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
ch <- zeroinfl(number ~ year+veg+shore+bottom | year+month+shore+bottom, dist='negbin', link="logit", data=ch.fl, na.action = na.exclude) 
j <- zeroinfl(number ~ year+month+veg+shore | year+month+veg, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
i <- zeroinfl(number ~ year+month+veg | year+month+veg+shore, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 

BestModelAP_zinb <- A 
BestModelCK_zinb <- C
BestModelTB_zinb <- t
BestModelCH_zinb <- ch
BestModelJX_zinb <- j
BestModelIR_zinb <- i


#alternatives
# BestModelAP
A <- hurdle(number ~ year +month+veg+shore | year+month+veg, dist='negbin', link="logit", data=ap.fl, na.action = na.exclude)
C <- hurdle(number ~ year+ month+veg | year+month+veg+shore, dist='negbin', link="logit", data=ck.fl, na.action = na.exclude)
t <- hurdle(number ~ year+month+veg+shore | year+month+shore, dist='negbin', link="logit", data=tb.fl, na.action = na.exclude) 
ch <- hurdle(number ~ year+month+veg+shore+bottom | year+month+veg+shore, dist='negbin', link="logit", data=ch.fl, na.action = na.exclude) 
j <- hurdle(number ~ veg +month | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.fl, na.action = na.exclude) 
i <- hurdle(number ~ year+month+veg | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.fl, na.action = na.exclude) 

BestModelAP_zanb <- A
BestModelCK_zanb <- C
BestModelTB_zanb <- t
BestModelCH_zanb <- ch
BestModelJX_zanb <- j
BestModelIR_zanb <- i 

##### DETERMINE LEAST SQUARE MEANS_YOY ###########
# DETERMINE LEAST SQUARE MEANS 
# Same thing as covariate adjusted means. Basically, determine the mean value of total positive numbers 
# of catch per year controlling for covariates (in this case it would be veg and shore variables). 
# Use lsmeans CRAN document. 

# Looking at the reference grid gives a good idea of over what levels the mean is being averaged. 
ap.rf.grid <- ref.grid(BestModelAP_zanb)

#Can make predictions using the reference grid. It produces a mean value of numbers based on each scenario combination.  
test_ap = summary(ap.rf.grid)

library(dplyr)
lsm.AP       <- lsmeans(BestModelAP_zanb, "year", data = ap.fl, mode="response")
estimate.AP <- as.data.frame(transform(summary(lsm.AP))) %>% select(year, lsmean, SE) 
#lsmeans with alternative model type (ZINB for AP)
alt.AP <- as.data.frame(transform(summary(lsmeans(BestModelAP_zinb,"year", data = ap.fl, mode="response"))))  %>% select(year, lsmean, SE) %>% mutate(alternative.mean=lsmean) %>% select(-lsmean)

lsm.CK       <- lsmeans(BestModelCK_zinb, "year", data = ck.fl, mode="response")
estimate.CK <- as.data.frame(transform(summary(lsm.CK))) %>% select(year, lsmean, SE) 
#lsmeans with alternative model type (ZANB for CK)
alt.CK <- as.data.frame(transform(summary(lsmeans(BestModelCK_zanb,"year", data = ck.fl, mode="response"))))  %>% select(year, lsmean, SE) %>% mutate(alternative.mean=lsmean) %>% select(-lsmean)

lsm.TB       <- lsmeans(BestModelTB_zinb, "year", data = tb.fl, mode="response")
estimate.TB <- as.data.frame(transform(summary(lsm.TB))) %>% select(year, lsmean, SE) 
#lsmeans with alternative model type (ZANB for TB)
alt.TB <- as.data.frame(transform(summary(lsmeans(BestModelTB_zanb,"year", data = tb.fl, mode="response"))))  %>% select(year, lsmean, SE) %>% mutate(alternative.mean=lsmean) %>% select(-lsmean)

lsm.CH       <- lsmeans(BestModelCH_zinb, "year", data = ch.fl, mode="response")
estimate.CH <- as.data.frame(transform(summary(lsm.CH))) %>% select(year, lsmean, SE) 
#lsmeans with alternative model type (ZANB for CH)
alt.CH <- as.data.frame(transform(summary(lsmeans(BestModelCH_zanb,"year", data = ch.fl, mode="response"))))  %>% select(year, lsmean, SE) %>% mutate(alternative.mean=lsmean) %>% select(-lsmean)

lsm.JX       <- lsmeans(BestModelJX_zinb, "year", data = jx.fl, mode="response")
estimate.JX <- as.data.frame(transform(summary(lsm.JX))) %>% select(year, lsmean, SE) 
#lsmeans with alternative model type (ZANB for JX)
alt.JX <- as.data.frame(transform(summary(lsmeans(BestModelJX_zanb,"year", data = jx.fl, mode="response"))))  %>% select(year, lsmean, SE) %>% mutate(alternative.mean=lsmean) %>% select(-lsmean)

lsm.IR       <- lsmeans(BestModelIR_zinb, "year", data = ir.fl, mode="response")
estimate.IR <- as.data.frame(transform(summary(lsm.IR))) %>% select(year, lsmean, SE) 
#lsmeans with alternative model type (ZANB for IR)
alt.IR <- as.data.frame(transform(summary(lsmeans(BestModelIR_zanb,"year", data = ir.fl, mode="response"))))  %>% select(year, lsmean, SE) %>% mutate(alternative.mean=lsmean) %>% select(-lsmean)


#Make a few thing that we need below for Monte Carlo 
# determining the number of trips per year
sample.sizeAP        <- as.data.frame(table(ap.fl$year))
names(sample.sizeAP) <- c("year","num.trips")
sample.sizeAP$year   <- as.character(sample.sizeAP$year)

sample.sizeCK        <- as.data.frame(table(ck.fl$year))
names(sample.sizeCK) <- c("year","num.trips")
sample.sizeCK$year   <- as.character(sample.sizeCK$year)

sample.sizeTB        <- as.data.frame(table(tb.fl$year))
names(sample.sizeTB) <- c("year","num.trips")
sample.sizeTB$year   <- as.character(sample.sizeTB$year)

sample.sizeCH        <- as.data.frame(table(ch.fl$year))
names(sample.sizeCH) <- c("year","num.trips")
sample.sizeCH$year   <- as.character(sample.sizeCH$year)

sample.sizeJX        <- as.data.frame(table(jx.fl$year))
names(sample.sizeJX) <- c("year","num.trips")
sample.sizeJX$year   <- as.character(sample.sizeJX$year)

sample.sizeIR        <- as.data.frame(table(ir.fl$year))
names(sample.sizeIR) <- c("year","num.trips")
sample.sizeIR$year   <- as.character(sample.sizeIR$year)

# Use lsmeans to determine the least square mean of values. 
#to display lsmeans on the response scale (with a glm model) as oppose to the log/logit scale you can do this:
# LSM_ap.pos <- summary(lsmeans(BestModelAP, 'year', data=ap.fl), type="response")
#However, this above method doesnt apply to every model supported by lsmeans. For example, zeroinfl and hurdle models 
# are returned on the response scale (on the scale of the observed counts) as the default. For more info see below source.
# https://rdrr.io/cran/lsmeans/man/models.html#heading-17

# CREATE INDEX WITH MONTE CARLO IF YOU WANT #####
# creating the index through Monte Carlo simulations
# Technically don't have to if you just want the estimates and the standard errors. 
#Monte Carlo simulations are necessary here because its the only way to incorporate the lsmeans error into the estimate and 
#then determine CV, and 95% conf intervals. Without a distribution you can't get to the confidence intervals which are what people actually care about. 
# This is borrowed from tuning_index_update.R from Hogfish2017 at FWRI. This monte carlo method seems standard among the group. 

# AP Monte Carlo ####
#making an empty matrix to fill with results 
num.yr <- length(estimate.AP$year)
index.dist <- matrix(data=NA,nrow=num.yr,ncol=8) 
num.iter=10000

#build the random distribution (i.e the random deviates)
for (i in 1:num.yr) {
  rand.1    <- qt(runif(num.iter,0.001,0.999),sample.sizeAP$num.trips) 
  
  #create the distribution of the data 
  APdist <- estimate.AP$lsmean[i] + estimate.AP$SE[i] * rand.1 
  
  index.dist[i,1] <- mean(APdist) #take the mean of the dist (aka temp)
  index.dist[i,2] <- sd(APdist) #take the sd of the dist
  index.dist[i,3] <- sd(APdist)/mean(APdist) #create CV
  index.dist[i,4:8] <- quantile(APdist,probs=c(0.025,0.25,0.50,0.75,0.975))
}
Upper        <- index.dist[ ,8] - index.dist[ ,7]
Lower        <- index.dist[ ,5] - index.dist[ ,4]
NominalMean  = as.vector(tapply(ap.fl$number, ap.fl$year, mean))
NominalSD    = as.vector(tapply(ap.fl$number, ap.fl$year, sd))
NominalCV    = NominalSD /NominalMean
APindex <- as.data.frame(cbind(sample.sizeAP, index.dist, Upper, Lower, NominalMean, NominalSD, NominalCV ))
names(APindex) <- c("year","Total.num.trips","Mean","std.dev", "CV","Low.95","Qtr.1","Median","Qtr.3","Up.95","Upper",
                  "Lower", "NominalMean", "NominalSD", "NominalCV" )
APindex$year = as.numeric(APindex$year)
  
      # CK Monte Carlo ####
#making an empty matrix to fill with results 
num.yr <- length(estimate.CK$year)
index.dist <- matrix(data=NA,nrow=num.yr,ncol=8) 
num.iter=10000

#build the random distribution (i.e the random deviates)
for (i in 1:num.yr) {
  rand.1    <- qt(runif(num.iter,0.001,0.999),sample.sizeCK$num.trips) 
  
  #create the distribution of the data 
  CKdist <- estimate.CK$lsmean[i] + estimate.CK$SE[i] * rand.1 
  
  index.dist[i,1] <- mean(CKdist) #take the mean of the dist (aka temp)
  index.dist[i,2] <- sd(CKdist) #take the sd of the dist
  index.dist[i,3] <- sd(CKdist)/mean(CKdist) #create CV
  index.dist[i,4:8] <- quantile(CKdist,probs=c(0.025,0.25,0.50,0.75,0.975))
}
Upper        <- index.dist[ ,8] - index.dist[ ,7]
Lower        <- index.dist[ ,5] - index.dist[ ,4]
NominalMean  = as.vector(tapply(ck.fl$number, ck.fl$year, mean))
NominalSD    = as.vector(tapply(ck.fl$number, ck.fl$year, sd))
NominalCV    = NominalSD /NominalMean
CKindex <- as.data.frame(cbind(sample.sizeCK, index.dist, Upper, Lower, NominalMean, NominalSD, NominalCV ))
names(CKindex) <- c("year","Total.num.trips","Mean","std.dev", "CV","Low.95","Qtr.1","Median","Qtr.3","Up.95","Upper",
                    "Lower", "NominalMean", "NominalSD", "NominalCV" )
CKindex$year = as.numeric(CKindex$year)

# TB Monte Carlo ####
num.yr <- length(estimate.TB$year)
index.dist <- matrix(data=NA,nrow=num.yr,ncol=8) 
num.iter=10000

#build the random distribution (i.e the random deviates)
for (i in 1:num.yr) {
  rand.1    <- qt(runif(num.iter,0.001,0.999),sample.sizeTB$num.trips) 
  
  #create the distribution of the data 
  TBdist <- estimate.TB$lsmean[i] + estimate.TB$SE[i] * rand.1 
  
  index.dist[i,1] <- mean(TBdist) #take the mean of the dist (aka temp)
  index.dist[i,2] <- sd(TBdist) #take the sd of the dist
  index.dist[i,3] <- sd(TBdist)/mean(TBdist) #create CV
  index.dist[i,4:8] <- quantile(TBdist,probs=c(0.025,0.25,0.50,0.75,0.975))
}
Upper        <- index.dist[ ,8] - index.dist[ ,7]
Lower        <- index.dist[ ,5] - index.dist[ ,4]
NominalMean  = as.vector(tapply(tb.fl$number, tb.fl$year, mean))
NominalSD    = as.vector(tapply(tb.fl$number, tb.fl$year, sd))
NominalCV    = NominalSD /NominalMean
TBindex <- as.data.frame(cbind(sample.sizeTB, index.dist, Upper, Lower, NominalMean, NominalSD, NominalCV ))
names(TBindex) <- c("year","Total.num.trips","Mean","std.dev", "CV","Low.95","Qtr.1","Median","Qtr.3","Up.95","Upper",
                    "Lower", "NominalMean", "NominalSD", "NominalCV" )
TBindex$year = as.numeric(TBindex$year)

# CH Monte Carlo ####
num.yr <- length(estimate.CH$year)
index.dist <- matrix(data=NA,nrow=num.yr,ncol=8) 
num.iter=10000

#build the random distribution (i.e the random deviates)
for (i in 1:num.yr) {
  rand.1    <- qt(runif(num.iter,0.001,0.999),sample.sizeCH$num.trips) 
  
  #create the distribution of the data 
  CHdist <- estimate.CH$lsmean[i] + estimate.CH$SE[i] * rand.1 
  
  index.dist[i,1] <- mean(CHdist) #take the mean of the dist (aka temp)
  index.dist[i,2] <- sd(CHdist) #take the sd of the dist
  index.dist[i,3] <- sd(CHdist)/mean(CHdist) #create CV
  index.dist[i,4:8] <- quantile(CHdist,probs=c(0.025,0.25,0.50,0.75,0.975))
}
Upper        <- index.dist[ ,8] - index.dist[ ,7]
Lower        <- index.dist[ ,5] - index.dist[ ,4]
NominalMean  = as.vector(tapply(ch.fl$number, ch.fl$year, mean))
NominalSD    = as.vector(tapply(ch.fl$number, ch.fl$year, sd))
NominalCV    = NominalSD /NominalMean
CHindex <- as.data.frame(cbind(sample.sizeCH, index.dist, Upper, Lower, NominalMean, NominalSD, NominalCV ))
names(CHindex) <- c("year","Total.num.trips","Mean","std.dev", "CV","Low.95","Qtr.1","Median","Qtr.3","Up.95","Upper",
                    "Lower", "NominalMean", "NominalSD", "NominalCV" )
CHindex$year = as.numeric(CHindex$year)

# JX Monte Carlo ####
num.yr <- length(estimate.JX$year)
index.dist <- matrix(data=NA,nrow=num.yr,ncol=8) 
num.iter=10000

#build the random distribution (i.e the random deviates)
for (i in 1:num.yr) {
  rand.1    <- qt(runif(num.iter,0.001,0.999),sample.sizeJX$num.trips) 
  
  #create the distribution of the data 
  JXdist <- estimate.JX$lsmean[i] + estimate.JX$SE[i] * rand.1 
  
  index.dist[i,1] <- mean(JXdist) #take the mean of the dist (aka temp)
  index.dist[i,2] <- sd(JXdist) #take the sd of the dist
  index.dist[i,3] <- sd(JXdist)/mean(JXdist) #create CV
  index.dist[i,4:8] <- quantile(JXdist,probs=c(0.025,0.25,0.50,0.75,0.975))
}
Upper        <- index.dist[ ,8] - index.dist[ ,7]
Lower        <- index.dist[ ,5] - index.dist[ ,4]
NominalMean  = as.vector(tapply(jx.fl$number, jx.fl$year, mean))
NominalSD    = as.vector(tapply(jx.fl$number, jx.fl$year, sd))
NominalCV    = NominalSD /NominalMean
JXindex <- as.data.frame(cbind(sample.sizeJX, index.dist, Upper, Lower, NominalMean, NominalSD, NominalCV ))
names(JXindex) <- c("year","Total.num.trips","Mean","std.dev", "CV","Low.95","Qtr.1","Median","Qtr.3","Up.95","Upper",
                    "Lower", "NominalMean", "NominalSD", "NominalCV" )
JXindex$year = as.numeric(JXindex$year)

# IR Monte Carlo ####

num.yr <- length(estimate.IR$year)
index.dist <- matrix(data=NA,nrow=num.yr,ncol=8) 
num.iter=10000

#build the random distribution (i.e the random deviates)
for (i in 1:num.yr) {
  rand.1    <- qt(runif(num.iter,0.001,0.999),sample.sizeIR$num.trips) 
  
  #create the distribution of the data 
  IRdist <- estimate.IR$lsmean[i] + estimate.IR$SE[i] * rand.1 
  
  index.dist[i,1] <- mean(IRdist) #take the mean of the dist (aka temp)
  index.dist[i,2] <- sd(IRdist) #take the sd of the dist
  index.dist[i,3] <- sd(IRdist)/mean(IRdist) #create CV
  index.dist[i,4:8] <- quantile(IRdist,probs=c(0.025,0.25,0.50,0.75,0.975))
}
Upper        <- index.dist[ ,8] - index.dist[ ,7]
Lower        <- index.dist[ ,5] - index.dist[ ,4]
NominalMean  = as.vector(tapply(ir.fl$number, ir.fl$year, mean))
NominalSD    = as.vector(tapply(ir.fl$number, ir.fl$year, sd))
NominalCV    = NominalSD /NominalMean
IRindex <- as.data.frame(cbind(sample.sizeIR, index.dist, Upper, Lower, NominalMean, NominalSD, NominalCV ))
names(IRindex) <- c("year","Total.num.trips","Mean","std.dev", "CV","Low.95","Qtr.1","Median","Qtr.3","Up.95","Upper",
                    "Lower", "NominalMean", "NominalSD", "NominalCV" )
IRindex$year = as.numeric(IRindex$year)

##### EXPORT PREDICTED INDEX_YOY ####
#export to csv _PERSONAL COMPUTER
#write.csv(Mean_AP, "~/Desktop/PhD project/Projects/Seatrout/Data/Indices/UpdatedIndices/AP_yoy_index.csv")
#write.csv(Mean_IR, "~/Desktop/PhD project/Projects/Seatrout/Data/Indices/UpdatedIndices/IR_yoy_index.csv")
#write.csv(Mean_JX, "~/Desktop/PhD project/Projects/Seatrout/Data/Indices/UpdatedIndices/JX_yoy_index.csv")
#write.csv(Mean_CH, "~/Desktop/PhD project/Projects/Seatrout/Data/Indices/UpdatedIndices/CH_yoy_index.csv")
#write.csv(Mean_TB, "~/Desktop/PhD project/Projects/Seatrout/Data/Indices/UpdatedIndices/TB_yoy_index.csv")
#write.csv(Mean_CK, "~/Desktop/PhD project/Projects/Seatrout/Data/Indices/UpdatedIndices/CK_yoy_index.csv")

#export to csv _WORK COMPUTER
write.csv(APindex, "U:/PhD_projectfiles/Exported_R_Datafiles/Indices/UpdatedIndices/AP_yoy_index.csv")
write.csv(CKindex, "U:/PhD_projectfiles/Exported_R_Datafiles/Indices/UpdatedIndices/CK_yoy_index.csv")
write.csv(TBindex, "U:/PhD_projectfiles/Exported_R_Datafiles/Indices/UpdatedIndices/TB_yoy_index.csv")
write.csv(CHindex, "U:/PhD_projectfiles/Exported_R_Datafiles/Indices/UpdatedIndices/CH_yoy_index.csv")
write.csv(JXindex, "U:/PhD_projectfiles/Exported_R_Datafiles/Indices/UpdatedIndices/JX_yoy_index.csv")
write.csv(IRindex, "U:/PhD_projectfiles/Exported_R_Datafiles/Indices/UpdatedIndices/IR_yoy_index.csv")

# PLOT INDICES ####

library(reshape2)
library(ggplot2)

#add on the alternative index 
APred <- APindex[,c(1,3,13)] %>% cbind(alt.AP$alternative.mean) %>% melt(id=c("year"))
CKred <- CKindex[,c(1,3,13)] %>% cbind(alt.CK$alternative.mean) %>% melt(id=c("year"))
TBred <- TBindex[,c(1,3,13)] %>% cbind(alt.TB$alternative.mean) %>% melt(id=c("year"))
CHred <- CHindex[,c(1,3,13)] %>% cbind(alt.CH$alternative.mean) %>% melt(id=c("year"))
JXred <- JXindex[,c(1,3,13)] %>% cbind(alt.JX$alternative.mean) %>% melt(id=c("year"))
IRred <- IRindex[,c(1,3,13)] %>% cbind(alt.IR$alternative.mean) %>% melt(id=c("year"))


library(ggplot2)

ggplot(APred, aes(x=year, y=value, color=variable))+
geom_line()  

ggplot(CKred, aes(x=year, y=value, color=variable))+
  geom_line()  

ggplot(TBred, aes(x=year, y=value, color=variable))+
  geom_line()  

ggplot(CHred, aes(x=year, y=value, color=variable))+
  geom_line()  

ggplot(JXred, aes(x=year, y=value, color=variable))+
  geom_line()  

ggplot(IRred, aes(x=year, y=value, color=variable))+
  geom_line()  

##### ADULT SECTION ####
# Now that the YOY index is done... I need to load the adult data because I need this to predict
# the SR curve and the expected number of survivals above it (residuals)
# After determining numbers, I will then apply an age schedule to the numbers and then apply weight schedule to
# convert total numbers to SSB necessary to create the stock recruitment curves

##### SET WORKING DIRECTORY_ADULT ####
setwd("U:/PhD_projectfiles/Raw_Data/Seatrout_FIM_Data/FIMData")
setwd("~/Desktop/PhD project/Projects/Seatrout/FWRI SCRATCH FOLDER/Elizabeth Herdter/SAS data sets/FIMData")

##### IMPORT DATASETS_ADULT ####
#
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

##### SELECT CATEGORICAL HABITAT VALUES_ ADULT ########
# to be used in the model to predict positive numbers

#Based on FWRI code the three variables were bottom type (bStr, bsan, bmud), bottom vegetation (bveg), and shoreline (Shore)
#There are three different bottom type variables each of them coded in a binary form.
#I want to take bStr, bsan, and bmud and put them into 1 variable so I will make a new variable entirely = 'bottom'
#I also want to turn bveg into a new variable = 'veg' based on the entries. If alg or Sav then turn to SAV because there are only 9 entries for Alg. 
#Same thing for the shore variable = 'shore'. Decided to have only emergent, structure, terrestrial, and mangrove. 
#Removed old variables (bStr, bSan, bMud, bveg, Shore)
#Removed rows when there was no shoreline variable. 

full_ad <-  select(full_ad, c(number,year, month, bio_reference, bStr, bSan, bMud, bveg, Shore, bay, Zone)) %>% 
  mutate(season = ifelse(month<5, "spring", ifelse(month>4 & month <9, "summer", ifelse(month>=9, "fall", "Non"))), bottom = ifelse(full_ad$bStr ==1, "structure", ifelse(full_ad$bSan>0 | full_ad$bMud>0, "mudsand", "unknown")), 
         veg= ifelse(full_ad$bveg == "SAVAlg", "SAV", ifelse(full_ad$bveg == "Alg", "SAV", ifelse(full_ad$bveg =="SAV", "SAV", "Noveg"))),
         shore = ifelse(substr(full_ad$Shore,1,3)=="Eme", "Emerge", ifelse(substr(full_ad$Shore,1,3) =="Man", "Mangrove", ifelse(substr(full_ad$Shore,1,3)=="Str", "Structure", 
            ifelse(substr(full_ad$Shore, 1,3)=="Ter", "Terrestrial", "Non"))))) %>% select(-c(bStr, bSan, bMud, bveg, Shore)) %>% subset(!shore=="Non") 

#Turn habitat variables into factors so they can be treated as categorical
full_ad[,c(2,5:10)] <- lapply(full_ad[,c(2,5:10)], factor)

ap.ad <- droplevels(full_ad %>% subset(bay =='AP'))
ck.ad <- droplevels(full_ad %>% subset(bay =='CK'))
tb.ad <- droplevels(full_ad %>% subset(bay =='TB'))
ch.ad <- droplevels(full_ad %>% subset(bay =='CH'))
jx.ad <- droplevels(full_ad %>% subset(bay =='JX'))
ir.ad <- droplevels(full_ad %>% subset(bay =='IR'))

##### AGGREGATE_ADULT ##########

#AP
#month
xtabs(number ~month, ap.ad)
ap.ad$month <- as.factor(as.character(ap.ad$month))
#veg
xtabs(number ~veg, ap.ad)
#bottom
xtabs(number ~bottom, ap.ad)
#shore
xtabs(number ~shore, ap.ad)
with(ap.ad,tapply(number, list(year,shore),sum))
ap.ad <- subset(ap.ad, shore !=  "Mangrove") %>% droplevels(ap.ad$shore)
ap.ad <- na.omit(ap.ad)

#CK
#month
xtabs(number ~ year + month, ck.ad)
ck.ad$month <- as.factor(as.character(ck.ad$month))
#veg
xtabs(number ~year +veg, ck.ad)
#bottom
xtabs(number ~year + bottom,ck.ad)
ck.ad <- droplevels(subset(ck.ad, bottom != "unknown"))
#shore 
xtabs(number ~year + shore, ck.ad)
ck.ad <- droplevels(subset(ck.ad, shore != "Terrestrial"))
ck.ad <- na.omit(ck.ad)

#TB
xtabs(number ~ year + month, tb.ad)
tb.ad$month <- as.factor(as.character(tb.ad$month))


with(tb.ad,tapply(number, list(year,veg),sum))
with(tb.ad,tapply(number, list(year,bottom),sum))
tb.ad <- droplevels(subset(tb.ad, bottom != "unknown"))
with(tb.ad,tapply(number, list(year,shore),sum))
tb.ad <- na.omit(tb.ad)

#CH - no aggregation needed
with(ch.ad,tapply(number, list(year,month),sum))
with(ch.ad,tapply(number, list(year,veg),sum))
with(ch.ad,tapply(number, list(year,bottom),sum))
with(ch.ad,tapply(number, list(year,shore),sum))
ch.ad <- droplevels(subset(ch.ad, shore != "Emerge"))
ch.ad$month <- as.factor(as.character(ch.ad$month))
ch.ad <- na.omit(ch.ad)

#JX 
with(jx.ad,tapply(number, list(year,month),sum))
with(jx.ad,tapply(number, list(year,veg),sum))
with(jx.ad,tapply(number, list(year,bottom),sum))
jx.ad <- droplevels(subset(jx.ad, bottom != "unknown"))
with(jx.ad,tapply(number, list(year,shore),sum))
jx.ad <- na.omit(jx.ad) #randomly an NA present
jx.ad$month <- as.factor(as.character(jx.ad$month))

#IR - no aggregation needed 
with(ir.ad,tapply(number, list(year,month),sum))
with(ir.ad,tapply(number, list(year,veg),sum))
with(ir.ad,tapply(number, list(year,bottom),sum))
ir.ad <- droplevels(subset(ir.ad, bottom != "unknown"))
with(ir.ad,tapply(number, list(year,shore),sum))
ir.ad <- na.omit(ir.ad)
ir.ad$month <- as.factor(as.character(ir.ad$month))


##### BUILD MODELS & VALIDATION_ADULT #########
# Build full Poisson, NB, ZINB, and ZANB

#1. Poisson and NB ####
apP.ad <- glm(number ~ year +month+veg+bottom+shore, data=ap.ad, family=poisson, na.action = na.exclude)
ckP.ad <- glm(number ~ year +month+veg+bottom+shore, data=ck.ad, family=poisson,  na.action = na.exclude)
tbP.ad <- glm(number ~ year +month+veg+bottom+shore, data=tb.ad, family=poisson,  na.action = na.exclude)
chP.ad <- glm(number ~ year +month+veg+bottom+shore, data=ch.ad, family=poisson, na.action = na.exclude)
jxP.ad <- glm(number ~ year +month+veg+bottom+shore, data=jx.ad, family=poisson, na.action = na.exclude)
irP.ad <- glm(number ~ year +month+veg+bottom+shore, data=ir.ad, family=poisson, na.action = na.exclude)

apNB.ad <- glm.nb(number ~ year +month+veg+bottom+shore, data=ap.ad, na.action = na.exclude)
ckNB.ad <- glm.nb(number ~ year +month+veg+bottom+shore, data=ck.ad,  na.action = na.exclude)
tbNB.ad <- glm.nb(number ~ year +month+veg+bottom+shore, data=tb.ad,  na.action = na.exclude)
chNB.ad <- glm.nb(number ~ year +month+veg+bottom+shore, data=ch.ad,  na.action = na.exclude)
jxNB.ad <- glm.nb(number ~ year +month+veg+bottom+shore, data=jx.ad,  na.action = na.exclude)
irNB.ad <- glm.nb(number ~ year +month+veg+bottom+shore, data=ir.ad,  na.action = na.exclude)


#odds test shows negative binomial better than poisson
odTest(apNB.ad)
odTest(ckNB.ad)
odTest(tbNB.ad)
odTest(chNB.ad)
odTest(jxNB.ad)
odTest(irNB.ad)

#model validation with NB
drop1(apNB.ad, test="Chi") #drop veg
apNB.ad1 <- glm.nb(number ~ year +month+bottom+shore, data=ap.ad, na.action = na.exclude)
drop1(apNB.ad1, test="Chi")
lrtest(apNB.ad, apNB.ad1) #equal so go with simpler (apNB.ad1)

sum((resid(apNB.ad1,type="pearson"))^2)/apNB.ad1$df.resid #1.95
sum((resid(apNB.ad1,type="pearson"))^2)/apNB.ad1$df.resid> 1+2*(sqrt(2/(apNB.ad1$df.residual)))
sum((resid(apNB.ad1,type="pearson"))^2)/apNB.ad1$df.resid> 1+3*(sqrt(2/(apNB.ad1$df.residual)))
#model selection did not improve fit of NB model to the ap.ad data 

drop1(ckNB.ad, test="Chi") #drop bottom 
ckNB.ad1 <- glm.nb(number ~ year +month+veg +shore, data=ck.ad, na.action = na.exclude)
drop1(ckNB.ad1, test="Chi")
lrtest(ckNB.ad, ckNB.ad1) #equal so go with simpler (ckNB.ad1)

sum((resid(ckNB.ad1,type="pearson"))^2)/ckNB.ad1$df.resid #1.26
sum((resid(ckNB.ad1,type="pearson"))^2)/ckNB.ad1$df.resid> 1+2*(sqrt(2/(ckNB.ad1$df.residual)))
sum((resid(ckNB.ad1,type="pearson"))^2)/ckNB.ad1$df.resid> 1+3*(sqrt(2/(ckNB.ad1$df.residual)))
#model selection did not improve fit of NB model to the ck.ad data 

drop1(tbNB.ad, test="Chi") #drop bottom and shore
tbNB.ad1 <- glm.nb(number ~ year +month+veg, data=tb.ad, na.action = na.exclude)
drop1(tbNB.ad1, test="Chi")
lrtest(tbNB.ad, tbNB.ad1) #equal so go with simpler (tbNB.ad1)

sum((resid(tbNB.ad1,type="pearson"))^2)/tbNB.ad1$df.resid #1.24
sum((resid(tbNB.ad1,type="pearson"))^2)/tbNB.ad1$df.resid> 1+2*(sqrt(2/(tbNB.ad1$df.residual)))
sum((resid(tbNB.ad1,type="pearson"))^2)/tbNB.ad1$df.resid> 1+3*(sqrt(2/(tbNB.ad1$df.residual)))
#model selection did not improve fit of NB model to the tb.ad data 

drop1(chNB.ad, test="Chi") #drop bottom and veg
chNB.ad1 <- glm.nb(number ~ year +month+shore, data=ch.ad, na.action = na.exclude)
drop1(chNB.ad1, test="Chi")
lrtest(chNB.ad, chNB.ad1) #equal so go with simpler (chNB.ad1)

sum((resid(chNB.ad1,type="pearson"))^2)/chNB.ad1$df.resid #1.29
sum((resid(chNB.ad1,type="pearson"))^2)/chNB.ad1$df.resid> 1+2*(sqrt(2/(chNB.ad1$df.residual)))
sum((resid(chNB.ad1,type="pearson"))^2)/chNB.ad1$df.resid> 1+3*(sqrt(2/(chNB.ad1$df.residual)))
#model selection did not improve fit of NB model to the ch.ad data 

drop1(jxNB.ad, test="Chi") #drop bottom 
jxNB.ad1 <- glm.nb(number ~ year +month+shore+veg, data=jx.ad, na.action = na.exclude)
drop1(jxNB.ad1, test="Chi")
lrtest(jxNB.ad, jxNB.ad1) #equal so go with simpler (jxNB.ad1)

sum((resid(jxNB.ad1,type="pearson"))^2)/jxNB.ad1$df.resid #1.55
sum((resid(jxNB.ad1,type="pearson"))^2)/jxNB.ad1$df.resid> 1+2*(sqrt(2/(jxNB.ad1$df.residual)))
sum((resid(jxNB.ad1,type="pearson"))^2)/jxNB.ad1$df.resid> 1+3*(sqrt(2/(jxNB.ad1$df.residual)))
#model selection did not improve fit of NB model to the jx.ad data 

drop1(irNB.ad, test="Chi") #drop bottom and shore
irNB.ad1 <- glm.nb(number ~ year +month+veg, data=ir.ad, na.action = na.exclude)
drop1(irNB.ad1, test="Chi")
lrtest(irNB.ad, irNB.ad1) #equal so go with simpler (irNB.ad1)

sum((resid(irNB.ad1,type="pearson"))^2)/irNB.ad1$df.resid #1.63
sum((resid(irNB.ad1,type="pearson"))^2)/irNB.ad1$df.resid> 1+2*(sqrt(2/(irNB.ad1$df.residual)))
sum((resid(irNB.ad1,type="pearson"))^2)/irNB.ad1$df.resid> 1+3*(sqrt(2/(irNB.ad1$df.residual)))
#model selection did not improve fit of NB model to the ir.ad data 

#NB is still overdispersed. Lets try ZINB then

#2. ZIP and ZINB and ZAP and ZANB ####

ap.ad_zip= zeroinfl(number ~ year+month+veg+bottom+shore | year+month+veg+bottom+shore,  link="logit", data=ap.ad, na.action = na.exclude)
ck.ad_zip= zeroinfl(number ~ year+month+veg+bottom+shore | year+month+veg+bottom+shore,  link="logit", data=ck.ad, na.action = na.exclude)
tb.ad_zip= zeroinfl(number ~ year+month+veg+bottom+shore | year+month+veg+bottom+shore,  link="logit", data=tb.ad, na.action = na.exclude)
ch.ad_zip= zeroinfl(number ~ year+month+veg+bottom+shore | year+month+veg+bottom+shore,  link="logit", data=ch.ad, na.action = na.exclude)
jx.ad_zip= zeroinfl(number ~ year+month+veg+bottom+shore | year+month+veg+bottom+shore,  link="logit", data=jx.ad, na.action = na.exclude)
ir.ad_zip= zeroinfl(number ~ year+month+veg+bottom+shore | year+month+veg+bottom+shore,  link="logit", data=ir.ad, na.action = na.exclude)


ap.ad_zinb= zeroinfl(number ~ year+month+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ap.ad, na.action = na.exclude)
ck.ad_zinb= zeroinfl(number ~ year+month+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
tb.ad_zinb= zeroinfl(number ~ year+month+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=tb.ad, na.action = na.exclude)
ch.ad_zinb= zeroinfl(number ~ year+month+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude)
jx.ad_zinb= zeroinfl(number ~ year+month+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude)
ir.ad_zinb= zeroinfl(number ~ year+month+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ir.ad, na.action = na.exclude)

lrtest(ap.ad_zip, ap.ad_zinb) 
lrtest(ck.ad_zip, ck.ad_zinb)
lrtest(tb.ad_zip, tb.ad_zinb)
lrtest(ch.ad_zip, ch.ad_zinb)
lrtest(jx.ad_zip, jx.ad_zinb)
lrtest(ir.ad_zip, ir.ad_zinb)
#zip is nested in zinb. More complicated model is best. Chose ZINB over ZIP

ap.ad_zap= hurdle(number ~ year+month+veg+bottom+shore | year+month+veg+bottom+shore, link="logit", data=ap.ad, na.action = na.exclude)
ck.ad_zap= hurdle(number ~ year+month+veg+bottom+shore | year+month+veg+bottom+shore, link="logit", data=ck.ad, na.action = na.exclude)
tb.ad_zap= hurdle(number ~ year+month+veg+bottom+shore | year+month+veg+bottom+shore,  link="logit", data=tb.ad, na.action = na.exclude)
ch.ad_zap= hurdle(number ~ year+month+veg+bottom+shore | year+month+veg+bottom+shore,  link="logit", data=ch.ad, na.action = na.exclude)
jx.ad_zap= hurdle(number ~ year+month+veg+bottom+shore | year+month+veg+bottom+shore,  link="logit", data=jx.ad, na.action = na.exclude)
ir.ad_zap= hurdle(number ~ year+month+veg+bottom+shore | year+month+veg+bottom+shore,  link="logit", data=ir.ad, na.action = na.exclude)

ap.ad_zanb= hurdle(number ~ year+month+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ap.ad, na.action = na.exclude)
ck.ad_zanb= hurdle(number ~ year+month+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
tb.ad_zanb= hurdle(number ~ year+month+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=tb.ad, na.action = na.exclude)
ch.ad_zanb= hurdle(number ~ year+month+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude)
jx.ad_zanb= hurdle(number ~ year+month+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude)
ir.ad_zanb= hurdle(number ~ year+month+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ir.ad, na.action = na.exclude)

lrtest(ap.ad_zap, ap.ad_zanb) 
lrtest(ck.ad_zap, ck.ad_zanb)
lrtest(tb.ad_zap, tb.ad_zanb)
lrtest(ch.ad_zap, ch.ad_zanb)
lrtest(jx.ad_zap, jx.ad_zanb)
lrtest(ir.ad_zap, ir.ad_zanb)
#zap is nested in zanb. More complicated model is best. Chose ZANB over ZAP

#ZINB vs ZANB?? 
AIC(ap.ad_zinb, ap.ad_zanb) #ZANB is likely best
AIC(ck.ad_zinb, ck.ad_zanb) #ZANB is likely best
AIC(tb.ad_zinb, tb.ad_zanb) #ZANB is likely best
#AIC(ch.ad_zinb, ch.ad_zanb) #ZANB is likely best
#AIC(jx.ad_zinb, jx.ad_zanb) #ZANB is likely best
AIC(ir.ad_zinb, ir.ad_zanb) #ZANB is likely best

# AP adult model selection ZINB ####
A= zeroinfl(number ~ year+month+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ap.ad, na.action = na.exclude)

#drop year from count
A1= zeroinfl(number ~ month+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ap.ad, na.action = na.exclude)
#drop month from count
A2= zeroinfl(number ~ year+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ap.ad, na.action = na.exclude)
#drop veg from count
A3= zeroinfl(number ~ year+month+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ap.ad, na.action = na.exclude)
#drop bottom from count
A4= zeroinfl(number ~ year+month+veg+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ap.ad, na.action = na.exclude)
#drop shore from count
A5= zeroinfl(number ~ year+month+veg+bottom | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ap.ad, na.action = na.exclude)

#drop year from binary
A6= zeroinfl(number ~ year+month+veg+bottom+shore |month+veg+bottom+shore , dist='negbin', link="logit", data=ap.ad, na.action = na.exclude)
#drop month from binary
A7= zeroinfl(number ~ year+month+veg+bottom+shore | year+veg+bottom+shore, dist='negbin', link="logit", data=ap.ad, na.action = na.exclude)
#drop veg from binary
A8= zeroinfl(number ~ year+month+veg+bottom+shore | year+month+bottom+shore, dist='negbin', link="logit", data=ap.ad, na.action = na.exclude)
#drop bottom from binary
A9= zeroinfl(number ~ year+month+veg+bottom+shore | year+month+veg+shore, dist='negbin', link="logit", data=ap.ad, na.action = na.exclude)
#drop shore from binary
A10= zeroinfl(number ~ year+month+veg+bottom+shore | year+month+veg+bottom, dist='negbin', link="logit", data=ap.ad, na.action = na.exclude)

lrtest(A, A1)
lrtest(A, A2)
lrtest(A, A3) #same so drop veg from count
lrtest(A, A4)
lrtest(A, A5)
lrtest(A, A6)
lrtest(A, A7)
lrtest(A, A8)
lrtest(A, A9) #same so drop bottom from binary
lrtest(A, A10)

#drop veg from count
A= zeroinfl(number ~ year+month+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ap.ad, na.action = na.exclude)

A1= zeroinfl(number ~ month+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ap.ad, na.action = na.exclude)
A2= zeroinfl(number ~ year+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ap.ad, na.action = na.exclude)
A3= zeroinfl(number ~ year+month+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ap.ad, na.action = na.exclude)
A4= zeroinfl(number ~ year+month+bottom | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ap.ad, na.action = na.exclude)
A5= zeroinfl(number ~ year+month+bottom+shore | month+veg+bottom+shore, dist='negbin', link="logit", data=ap.ad, na.action = na.exclude)
A6= zeroinfl(number ~ year+month+bottom+shore | year+veg+bottom+shore, dist='negbin', link="logit", data=ap.ad, na.action = na.exclude)
A7= zeroinfl(number ~ year+month+bottom+shore | year+month+bottom+shore, dist='negbin', link="logit", data=ap.ad, na.action = na.exclude)
A8= zeroinfl(number ~ year+month+bottom+shore | year+month+veg+shore, dist='negbin', link="logit", data=ap.ad, na.action = na.exclude)
A9= zeroinfl(number ~ year+month+bottom+shore | year+month+veg+bottom, dist='negbin', link="logit", data=ap.ad, na.action = na.exclude)

lrtest(A, A1)
lrtest(A, A2)
lrtest(A, A3) 
lrtest(A, A4)
lrtest(A, A5)
lrtest(A, A6)
lrtest(A, A7)
lrtest(A, A8) #equal so drop bottom from binary
lrtest(A, A9)

#drop bottom from binary
A= zeroinfl(number ~ year+month+bottom+shore | year+month+veg+shore, dist='negbin', link="logit", data=ap.ad, na.action = na.exclude)

A1= zeroinfl(number ~ month+bottom+shore | year+month+veg+shore, dist='negbin', link="logit", data=ap.ad, na.action = na.exclude)
A2= zeroinfl(number ~ year+bottom+shore | year+month+veg+shore, dist='negbin', link="logit", data=ap.ad, na.action = na.exclude)
A3= zeroinfl(number ~ year+month+shore | year+month+veg+shore, dist='negbin', link="logit", data=ap.ad, na.action = na.exclude)
A4= zeroinfl(number ~ year+month+bottom | year+month+veg+shore, dist='negbin', link="logit", data=ap.ad, na.action = na.exclude)
A5= zeroinfl(number ~ year+month+bottom+shore | month+veg+shore, dist='negbin', link="logit", data=ap.ad, na.action = na.exclude)
A6= zeroinfl(number ~ year+month+bottom+shore | year+veg+shore, dist='negbin', link="logit", data=ap.ad, na.action = na.exclude)
A7= zeroinfl(number ~ year+month+bottom+shore | year+month+shore, dist='negbin', link="logit", data=ap.ad, na.action = na.exclude)
A8= zeroinfl(number ~ year+month+bottom+shore | year+month+veg, dist='negbin', link="logit", data=ap.ad, na.action = na.exclude)

lrtest(A, A1)
lrtest(A, A2)
lrtest(A, A3) 
lrtest(A, A4)
#lrtest(A, A5)
#lrtest(A, A6)
lrtest(A, A7)
lrtest(A, A8)

AIC(A, A1, A2, A3, A4, A7, A8)

#all remaining covariates significant; A best model; can be confirmed using AIC

summary(A)	#==>examine standard errors of covariates; unusually large values indicate problem with model

N<-nrow(ap.ad)
EZIP <- resid(A,type="pearson")
Dispersion <- sum(EZIP^2)/(N-65)	#==>be sure to change value for degrees of freedom based on summary output
Dispersion #2.32

#==>Dispersion is 1.12 for model which is an improvement from  1.39 in the NB 
sum((resid(A,type="pearson"))^2)/(N-65) #1.12 
sum((resid(A,type="pearson"))^2)/(N-65) > 1+3*(sqrt(2/(N-65)))
BestModelAP.ad_zinb <- A

# CK adult model selection ZINB ####
C= zeroinfl(number ~ year+month+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)

C1= zeroinfl(number ~ month+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C2= zeroinfl(number ~ year+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
#C3= zeroinfl(number ~ year+month+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
#C4= zeroinfl(number ~ year+month+veg+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C5= zeroinfl(number ~ year+month+veg+bottom| year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C6= zeroinfl(number ~ year+month+veg+bottom+shore | month+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
#C7= zeroinfl(number ~ year+month+veg+bottom+shore | year+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C8= zeroinfl(number ~ year+month+veg+bottom+shore | year+month+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C9= zeroinfl(number ~ year+month+veg+bottom+shore | year+month+veg+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C10= zeroinfl(number ~ year+month+veg+bottom+shore | year+month+veg+bottom, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)

lrtest(C, C1)
lrtest(C, C2)
lrtest(C, C5) #drop shore from count
lrtest(C, C6)
lrtest(C, C8) #drop veg from binary
lrtest(C, C9) #drop bottom from binary
lrtest(C, C10) #drop shore from binary

#drop shore from count
C= zeroinfl(number ~ year+month+veg+bottom| year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)

C1= zeroinfl(number ~ month+veg+bottom| year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C2= zeroinfl(number ~ year+veg+bottom| year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
#C3= zeroinfl(number ~ year+month+bottom| year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C4= zeroinfl(number ~ year+month+veg| year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C5= zeroinfl(number ~ year+month+veg+bottom| month+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
#C6= zeroinfl(number ~ year+month+veg+bottom| year+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C7= zeroinfl(number ~ year+month+veg+bottom| year+month+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C8= zeroinfl(number ~ year+month+veg+bottom| year+month+veg+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C9= zeroinfl(number ~ year+month+veg+bottom| year+month+veg+bottom, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)

lrtest(C, C1)
lrtest(C, C2)
lrtest(C, C4) #drop bottom from count
lrtest(C, C5)
lrtest(C, C7) #drop veg from binary
lrtest(C, C8) #drop bottom 
lrtest(C, C9) #drop shore

#drop bottom from count
C= zeroinfl(number ~ year+month+veg | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)

C1= zeroinfl(number ~ month+veg | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C2= zeroinfl(number ~ year+veg | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
#C3= zeroinfl(number ~ year+month | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
#C4= zeroinfl(number ~ year+month+veg | month+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
#C5= zeroinfl(number ~ year+month+veg | year+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
#C6= zeroinfl(number ~ year+month+veg | year+month+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C7= zeroinfl(number ~ year+month+veg | year+month+veg+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C8= zeroinfl(number ~ year+month+veg | year+month+veg+bottom, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)

lrtest(C, C1)
lrtest(C, C2)
lrtest(C, C7) #drop bottom
lrtest(C, C8) #drop shore

#drop bottom from binary
C= zeroinfl(number ~ year+month+veg | year+month+veg+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)

C1= zeroinfl(number ~ month+veg | year+month+veg+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C2= zeroinfl(number ~ year+veg | year+month+veg+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
#C3= zeroinfl(number ~ year+month| year+month+veg+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
#C4= zeroinfl(number ~ year+month+veg | month+veg+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
#C5= zeroinfl(number ~ year+month+veg | year+veg+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C6= zeroinfl(number ~ year+month+veg | year+month+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C7= zeroinfl(number ~ year+month+veg | year+month+veg, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)

lrtest(C, C1)
lrtest(C, C2)
lrtest(C, C6) #drop veg from binary
lrtest(C, C7) #drop shore from binary

#drop veg from binary
C= zeroinfl(number ~ year+month+veg | year+month+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)

#C1= zeroinfl(number ~ month+veg | year+month+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C2= zeroinfl(number ~ year+veg | year+month+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C3= zeroinfl(number ~ year+month | year+month+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C4= zeroinfl(number ~ year+month+veg | month+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C5= zeroinfl(number ~ year+month+veg | year+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
#C6= zeroinfl(number ~ year+month+veg | year+month, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)

lrtest(C, C2)
lrtest(C, C3) 
lrtest(C, C5)








 #C= zeroinfl(number ~ year+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)

#C1= zeroinfl(number ~ veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C2= zeroinfl(number ~ year+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C3= zeroinfl(number ~ year+veg+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C4= zeroinfl(number ~ year+veg+bottom | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
#C5= zeroinfl(number ~ year+veg+bottom+shore | month+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C6= zeroinfl(number ~ year+veg+bottom+shore | year+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C7= zeroinfl(number ~ year+veg+bottom+shore | year+month+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C8= zeroinfl(number ~ year+veg+bottom+shore | year+month+veg+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C9= zeroinfl(number ~ year+veg+bottom+shore | year+month+veg+bottom, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)

#lrtest(C, C1)
lrtest(C, C2)
lrtest(C, C3)
lrtest(C, C4) 
#lrtest(C, C5) 
lrtest(C, C6) #equal so drop month from binary
lrtest(C, C7) 
lrtest(C, C8) 
lrtest(C, C9) 

#drop month from binary
C= zeroinfl(number ~ year+veg+bottom+shore | year+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)

#C1= zeroinfl(number ~ veg+bottom+shore | year+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
#C2= zeroinfl(number ~ year+bottom+shore | year+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
#C3= zeroinfl(number ~ year+veg+shore | year+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
#C4= zeroinfl(number ~ year+veg+bottom | year+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C5= zeroinfl(number ~ year+veg+bottom+shore | veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
#C6= zeroinfl(number ~ year+veg+bottom+shore | year+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
#C7= zeroinfl(number ~ year+veg+bottom+shore | year+veg+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
#C8= zeroinfl(number ~ year+veg+bottom+shore | year+veg+bottom, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)

lrtest(C, C5) #equal so drop year from binary

#drop year from binary
C= zeroinfl(number ~ year+veg+bottom+shore | veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)

C1= zeroinfl(number ~ veg+bottom+shore | veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C2= zeroinfl(number ~ year+bottom+shore | veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C3= zeroinfl(number ~ year+veg+shore | veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C4= zeroinfl(number ~ year+veg+bottom | veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C5= zeroinfl(number ~ year+veg+bottom+shore | bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C6= zeroinfl(number ~ year+veg+bottom+shore | veg+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C7= zeroinfl(number ~ year+veg+bottom+shore | veg+bottom, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)

lrtest(C, C1)
lrtest(C, C2)
lrtest(C, C3) #equal sp drop bottom from count
lrtest(C, C4) #equal so drop shore from count
lrtest(C, C5) #equal so drop veg from binary
lrtest(C, C6)
lrtest(C, C7) #equal so drop shore form binary

#drop bottom from count
C= zeroinfl(number ~ year+veg+shore | veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)

C1= zeroinfl(number ~ veg+shore | veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C2= zeroinfl(number ~ year+shore | veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C3= zeroinfl(number ~ year+veg| veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C4= zeroinfl(number ~ year+veg+shore |bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C5= zeroinfl(number ~ year+veg+shore | veg+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C6= zeroinfl(number ~ year+veg+shore | veg+bottom, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)

lrtest(C, C1)
lrtest(C, C2)
lrtest(C, C3) #drop shore from count
lrtest(C, C4) #drop veg from binary
lrtest(C, C5)
lrtest(C, C6) #drop shore from binary

#drop shore from count
C= zeroinfl(number ~ year+veg | veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)

C1= zeroinfl(number ~ veg | veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C2= zeroinfl(number ~ year| veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C3= zeroinfl(number ~ year+veg | bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C4= zeroinfl(number ~ year+veg | veg+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C5= zeroinfl(number ~ year+veg | veg+bottom, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)

lrtest(C, C1)
lrtest(C, C2)
lrtest(C, C3) #drop veg from binary
lrtest(C, C4)
lrtest(C, C5) #drop shore from binary

#drop veg from binary
C= zeroinfl(number ~ year+veg | bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)

C1= zeroinfl(number ~ veg | bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C2= zeroinfl(number ~ year| bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C3= zeroinfl(number ~ year+veg | shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C4= zeroinfl(number ~ year+veg | bottom, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)

lrtest(C, C1)
lrtest(C, C2)
lrtest(C, C3)
lrtest(C, C4) #drop shore from bottom

#drop veg from binary
C= zeroinfl(number ~ year+veg | bottom, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)




ck.ad <- as.data.frame(ck.ad %>% mutate(bin = ifelse(ck.ad$number>= 1, 1, ck.ad$number)))






lrtest(C, C4) #equal so drop year from binary 
#lrtest(C, C5) 
lrtest(C, C6) #drop veg from binary
lrtest(C, C7) #equal so drop bottom from binary
lrtest(C, C8) #equal so drop shore from binary

#drop veg from binary 
C= zeroinfl(number ~ year+month+veg | year+month+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)

#C1= zeroinfl(number ~ month+veg | year+month+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C2= zeroinfl(number ~ year+veg | year+month+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C3= zeroinfl(number ~ year+month| year+month+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C4= zeroinfl(number ~ year+month+veg | month+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
#C5= zeroinfl(number ~ year+month+veg | yearbottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C6= zeroinfl(number ~ year+month+veg | year+month+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C7= zeroinfl(number ~ year+month+veg | year+month+bottom, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)

#lrtest(C, C1) 
lrtest(C, C2) 
lrtest(C, C3) 
lrtest(C, C4) #drop year from binary
lrtest(C, C6) #drop bottom from binary
lrtest(C, C7) #equal so drop shore from binary

#drop year from binary 
C= zeroinfl(number ~ year+month+veg | month+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)

C1= zeroinfl(number ~ month+veg | month+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
#C2= zeroinfl(number ~ year+veg | month+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C3= zeroinfl(number ~ year+month | month+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C4= zeroinfl(number ~ year+month+veg | bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C5= zeroinfl(number ~ year+month+veg | month+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C6= zeroinfl(number ~ year+month+veg | month+bottom, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)

lrtest(C, C1) 
#lrtest(C, C2) 
lrtest(C, C3) 
lrtest(C, C4) #drop month from binary
lrtest(C, C5) 
lrtest(C, C6) #equal so drop shore from binary

#drop month from binary 
C= zeroinfl(number ~ year+month+veg | bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)

C1= zeroinfl(number ~ month+veg | bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C2= zeroinfl(number ~ year+veg | bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C3= zeroinfl(number ~ year+month | bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C4= zeroinfl(number ~ year+month+veg | shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
C5= zeroinfl(number ~ year+month+veg | bottom, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)

lrtest(C, C1) 
lrtest(C, C2) 
lrtest(C, C3) 
lrtest(C, C4) 
lrtest(C, C5) #drop shore from binary

#drop shore from binary 

AIC(C, C1, C2,C3, C4, C5) 

# model C appears best

N<-nrow(ck.ad)
EZIP <- resid(C5,type="pearson")
Dispersion <- sum(EZIP^2)/(N-34)	#==>be sure to change value for degrees of freedom based on summary output
Dispersion #1.27

BestModelCK.ad_zinb <- C5



# TB adult model validation ZINB ####
#full
t <- zeroinfl(number ~ year+month+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=tb.ad, na.action = na.exclude) 
#drop year
t1 <- zeroinfl(number ~month+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=tb.ad, na.action = na.exclude) 
#drop month
t2 <- zeroinfl(number ~ year+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=tb.ad, na.action = na.exclude) 
#drop veg
#t3 <- zeroinfl(number ~ year+month+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=tb.ad, na.action = na.exclude) 
#drop bottom
#t4 <- zeroinfl(number ~ year+month+veg+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=tb.ad, na.action = na.exclude) 
#drop shore
t5 <- zeroinfl(number ~ year+month+veg+bottom | year+month+veg+bottom+shore, dist='negbin', link="logit", data=tb.ad, na.action = na.exclude) 
#drop year binary
t6 <- zeroinfl(number ~ year+month+veg+bottom+shore | month+veg+bottom+shore, dist='negbin', link="logit", data=tb.ad, na.action = na.exclude) 
#drop month binary
#t7 <- zeroinfl(number ~ year+month+veg+bottom+shore | year+veg+bottom+shore, dist='negbin', link="logit", data=tb.ad, na.action = na.exclude) 
#drop veg binary
#t8 <- zeroinfl(number ~ year+month+veg+bottom+shore | year+month+bottom+shore, dist='negbin', link="logit", data=tb.ad, na.action = na.exclude) 
#drop bottom binary
#t9 <- zeroinfl(number ~ year+month+veg+bottom+shore | year+month+veg+shore, dist='negbin', link="logit", data=tb.ad, na.action = na.exclude) 
#drop shore binary
t10 <- zeroinfl(number ~ year+month+veg+bottom+shore | year+month+veg+bottom, dist='negbin', link="logit", data=tb.ad, na.action = na.exclude) 

lrtest(t, t1)
lrtest(t, t2)
#lrtest(t, t4) 
lrtest(t, t5)
lrtest(t, t6)
lrtest(t, t10) 

AIC(t, t1, t2, t5, t6, t10)
# All remaining covariates are significant

 N<-nrow(tb.ad)
 EZIP <- resid(t,type="pearson")
 Dispersion <- sum(EZIP^2)/(N-73)	#==>be sure to change value for degrees of freedom based on summary output
 Dispersion #1.15.. improvement from 1.67 with negative binomial 

#BestModelTB.ad_zinb <- t

# CH adult model validation ZINB ####
#ch <- zeroinfl(number ~ year+month+veg+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 
#drop year
#ch1 <- zeroinfl(number ~ month+veg+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 
#drop month
ch2 <- zeroinfl(number ~ year+veg+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 
#drop veg
#ch3 <- zeroinfl(number ~ year+month+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 
#drop shore
#ch4 <- zeroinfl(number ~ year+month+veg+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 
#drop bottom
#ch5 <- zeroinfl(number ~ year+month+veg+shore | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 
#drop year binary
ch6 <- zeroinfl(number ~ year+month+veg+shore+bottom | month+veg+shore+bottom, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 
#drop month binary
ch7 <- zeroinfl(number ~ year+month+veg+shore+bottom | year+veg+shore+bottom, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 
#drop veg binary
#ch8 <- zeroinfl(number ~ year+month+veg+shore+bottom | year+month+shore+bottom, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 
#drop shore binary
#ch9 <- zeroinfl(number ~ year+month+veg+shore+bottom | year+month+veg+bottom, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 
#drop bottom binary
#ch10 <- zeroinfl(number ~ year+month+veg+shore+bottom | year+month+veg+shore, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 

#few of the above will solve so lets start with this one. 
#drop month from count
ch <- zeroinfl(number ~ year+veg+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 

#ch1 <- zeroinfl(number ~ veg+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 
ch2 <- zeroinfl(number ~ year+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 
ch3 <- zeroinfl(number ~ year+veg+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 
ch4 <- zeroinfl(number ~ year+veg+shore | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 
ch5 <- zeroinfl(number ~ year+veg+shore+bottom |month+veg+shore+bottom, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 
ch6 <- zeroinfl(number ~ year+veg+shore+bottom | year+veg+shore+bottom, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 
ch7 <- zeroinfl(number ~ year+veg+shore+bottom | year+month+shore+bottom, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 
ch8 <- zeroinfl(number ~ year+veg+shore+bottom | year+month+veg+bottom, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 
ch9 <- zeroinfl(number ~ year+veg+shore+bottom | year+month+veg+shore, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 

lrtest(ch, ch2) #drop veg from count
lrtest(ch, ch3)
lrtest(ch, ch4)
lrtest(ch, ch5)
lrtest(ch, ch6)
lrtest(ch, ch7) #drop veg from binary
lrtest(ch, ch8) #drop shore from binary
lrtest(ch, ch9)

#drop veg from count
ch <- zeroinfl(number ~ year+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 

ch1 <- zeroinfl(number ~ shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 
ch2 <- zeroinfl(number ~ year+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 
#ch3 <- zeroinfl(number ~ year+shore | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 
ch4 <- zeroinfl(number ~ year+shore+bottom | month+veg+shore+bottom, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 
ch5 <- zeroinfl(number ~ year+shore+bottom | year+veg+shore+bottom, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 
ch6 <- zeroinfl(number ~ year+shore+bottom | year+month+shore+bottom, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 
ch7 <- zeroinfl(number ~ year+shore+bottom | year+month+veg+bottom, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 
ch8 <- zeroinfl(number ~ year+shore+bottom | year+month+veg+shore, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 

lrtest(ch, ch1)
lrtest(ch, ch2)
lrtest(ch, ch4)
lrtest(ch, ch5)
lrtest(ch, ch6) # drop veg form binary 
lrtest(ch, ch7) #drop shore from binary
lrtest(ch, ch8)

#drop veg from binary
ch <- zeroinfl(number ~ year+shore+bottom | year+month+shore+bottom, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 

ch1 <- zeroinfl(number ~ shore+bottom | year+month+shore+bottom, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 
ch2 <- zeroinfl(number ~ year+bottom | year+month+shore+bottom, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 
ch3 <- zeroinfl(number ~ year+shore | year+month+shore+bottom, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 
ch4 <- zeroinfl(number ~ year+shore+bottom | month+shore+bottom, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 
ch5 <- zeroinfl(number ~ year+shore+bottom | year+shore+bottom, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 
ch6 <- zeroinfl(number ~ year+shore+bottom | year+month+bottom, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 
ch7 <- zeroinfl(number ~ year+shore+bottom | year+month+shore, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 

lrtest(ch, ch1)
lrtest(ch, ch2)
lrtest(ch, ch3)
lrtest(ch, ch4)
lrtest(ch, ch5)
lrtest(ch, ch6) # drop shore from binary
lrtest(ch, ch7) 

#drop shore from binary
ch <- zeroinfl(number ~ year+shore+bottom | year+month+bottom, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 

ch1 <- zeroinfl(number ~ shore+bottom | year+month+bottom, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 
ch2 <- zeroinfl(number ~ year+bottom | year+month+bottom, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 
#ch3 <- zeroinfl(number ~ year+shore | year+month+bottom, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 
ch4 <- zeroinfl(number ~ year+shore+bottom | month+bottom, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 
ch5 <- zeroinfl(number ~ year+shore+bottom | year+bottom, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 
ch6 <- zeroinfl(number ~ year+shore+bottom | year+month, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 

lrtest(ch, ch1)
lrtest(ch, ch2)
#lrtest(ch, ch3)
lrtest(ch, ch4)
lrtest(ch, ch5)
lrtest(ch, ch6) 

AIC(ch, ch1, ch2, ch4, ch5, ch6)
#model ch seems to be the best
#All covariates are significant
#Huge standard errors for a few years in the binary. Weird. 

N<-nrow(ch.ad)
EZIP <- resid(ch,type="pearson")
Dispersion <- sum(EZIP^2)/(N-56)	#==>be sure to change value for degrees of freedom based on summary output
Dispersion #1.18.. improvement from 1.52 with negative binomial 

BestModelCH.ad_zinb <- ch

# JX adult model validation ZINB ####
#j <- zeroinfl(number ~ year+month+veg+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 

#drop year
#j1 <- zeroinfl(number ~ month+veg+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 
#drop month
j2 <- zeroinfl(number ~ year+veg+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 
#drop veg
#j3 <- zeroinfl(number ~ year+month+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 
#drop shore
j4 <- zeroinfl(number ~ year+month+veg+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 
#drop bottom
#j5 <- zeroinfl(number ~ year+month+veg+shore | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 
#drop year binary
j6 <- zeroinfl(number ~ year+month+veg+shore+bottom | month+veg+shore+bottom, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 
#drop month binary
j7 <- zeroinfl(number ~ year+month+veg+shore+bottom | year+veg+shore+bottom, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 
#drop veg binary
#j8 <- zeroinfl(number ~ year+month+veg+shore+bottom | year+month+shore+bottom, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 
#drop shore binary
#j9 <- zeroinfl(number ~ year+month+veg+shore+bottom | year+month+veg+bottom, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 
#drop bottom binary
#j10 <- zeroinfl(number ~ year+month+veg+shore+bottom | year+month+veg+shore, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 


#many wont solve above so just drop month
j <- zeroinfl(number ~ year+veg+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 

j1 <- zeroinfl(number ~ veg+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 
j2 <- zeroinfl(number ~ year+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 
j3 <- zeroinfl(number ~ year+veg+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 
j4 <- zeroinfl(number ~ year+veg+shore | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 
j5 <- zeroinfl(number ~ year+veg+shore+bottom | month+veg+shore+bottom, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 
j6 <- zeroinfl(number ~ year+veg+shore+bottom | year+veg+shore+bottom, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 
j7 <- zeroinfl(number ~ year+veg+shore+bottom | year+month+shore+bottom, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 
#j8 <- zeroinfl(number ~ year+veg+shore+bottom | year+month+veg+bottom, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 
j9 <- zeroinfl(number ~ year+veg+shore+bottom | year+month+veg+shore, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 

lrtest(j, j1)
lrtest(j, j2) #drop veg from count 
lrtest(j, j3) #drop shore from count
lrtest(j, j4) #drop bottom form count
lrtest(j, j5)
lrtest(j, j6) #drop month from binary
lrtest(j, j7)
#lrtest(j, j8)
lrtest(j, j9) #drop bottom from binary 

#drop veg from count
j <- zeroinfl(number ~ year+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 

j1 <- zeroinfl(number ~ shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 
j2 <- zeroinfl(number ~ year+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 
j3 <- zeroinfl(number ~ year+shore| year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 
j4 <- zeroinfl(number ~ year+shore+bottom | month+veg+shore+bottom, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 
j5 <- zeroinfl(number ~ year+shore+bottom | year+veg+shore+bottom, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 
j6 <- zeroinfl(number ~ year+shore+bottom | year+month+shore+bottom, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 
#j7 <- zeroinfl(number ~ year+shore+bottom | year+month+veg+bottom, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 
j8 <- zeroinfl(number ~ year+shore+bottom | year+month+veg+shore, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 

lrtest(j, j1)
lrtest(j, j2) #drop shore from count 
lrtest(j, j3) #drop bottom from count
lrtest(j, j4)
lrtest(j, j5) #drop month from binary
lrtest(j, j6) 
#lrtest(j, j7)
lrtest(j, j8) #drop bottom from binary

#drop shore from count
j <- zeroinfl(number ~ year+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 

j1 <- zeroinfl(number ~ bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 
j2 <- zeroinfl(number ~ year | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 
j3 <- zeroinfl(number ~ year+bottom | month+veg+shore+bottom, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 
j4 <- zeroinfl(number ~ year+bottom | year+veg+shore+bottom, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 
j5 <- zeroinfl(number ~ year+bottom | year+month+shore+bottom, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 
#j6 <- zeroinfl(number ~ year+bottom | year+month+veg+bottom, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 
j7 <- zeroinfl(number ~ year+bottom | year+month+veg+shore, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 

lrtest(j, j1)
lrtest(j, j2) #drop bottom
lrtest(j, j3) 
lrtest(j, j4) #drop momth from binary
lrtest(j, j5) 
#lrtest(j, j6) 
lrtest(j, j7) #drop bottom from binary

#drop bottom from count
j <- zeroinfl(number ~ year | year+month+veg+shore+bottom, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 

j1 <- zeroinfl(number ~ year | month+veg+shore+bottom, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 
j2 <- zeroinfl(number ~ year | year+veg+shore+bottom, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 
j3 <- zeroinfl(number ~ year | year+month+shore+bottom, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 
#j4 <- zeroinfl(number ~ year | year+month+veg+bottom, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 
j5 <- zeroinfl(number ~ year | year+month+veg+shore, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 


lrtest(j, j1)
lrtest(j, j2) 
lrtest(j, j3) 
#lrtest(j, j4) 
lrtest(j, j5) #drop bottom from binary

#drop bottom from count
j <- zeroinfl(number ~ year | year+month+veg+shore, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 

j1 <- zeroinfl(number ~ year | month+veg+shore, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 
j2 <- zeroinfl(number ~ year | year+veg+shore, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 
j3 <- zeroinfl(number ~ year | year+month+shore, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 
#j4 <- zeroinfl(number ~ year | year+month+veg, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 

lrtest(j, j1)
lrtest(j, j2) 
lrtest(j, j3) 

AIC(j, j1, j2, j3) #undecided. says J2 might be better
#best model is j

N<-nrow(jx.ad)
EZIP <- resid(j2,type="pearson")
Dispersion <- sum(EZIP^2)/(N-34)	#==>be sure to change value for degrees of freedom based on summary output
Dispersion #1.08.. improvement from 1.34 with negative binomial 

BestModelJX.ad_zinb <- j2

# IR adult model validation ZINB ####
i <- zeroinfl(number ~ year+month+veg+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action = na.exclude) 

i1 <- zeroinfl(number ~ month+veg+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action = na.exclude) 
i2 <- zeroinfl(number ~ year+veg+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action = na.exclude) 
i3 <- zeroinfl(number ~ year+month+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action = na.exclude) 
#i4 <- zeroinfl(number ~ year+month+veg+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action = na.exclude) 
i5 <- zeroinfl(number ~ year+month+veg+shore | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action = na.exclude) 
i6 <- zeroinfl(number ~ year+month+veg+shore+bottom | month+veg+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action = na.exclude) 
i7 <- zeroinfl(number ~ year+month+veg+shore+bottom | year+veg+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action = na.exclude) 
#i8 <- zeroinfl(number ~ year+month+veg+shore+bottom | year+month+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action = na.exclude) 
#i9 <- zeroinfl(number ~ year+month+veg+shore+bottom | year+month+veg+bottom, dist='negbin', link="logit", data=ir.ad, na.action = na.exclude) 
i10 <- zeroinfl(number ~ year+month+veg+shore+bottom | year+month+veg+shore, dist='negbin', link="logit", data=ir.ad, na.action = na.exclude) 

lrtest(i, i1)
lrtest(i, i2)
lrtest(i, i3)
#lrtest(i, i4)
lrtest(i, i5) #equal so drop bottom from count
lrtest(i, i6)
lrtest(i, i7)
#lrtest(i, i8)
#lrtest(i, i9)
lrtest(i, i10) #equal so drop bottom from binary

#drop bottom from count
i <- zeroinfl(number ~ year+month+veg+shore | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action = na.exclude) 

i1 <- zeroinfl(number ~ month+veg+shore | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action = na.exclude) 
i2 <- zeroinfl(number ~ year+veg+shore | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action = na.exclude) 
i3 <- zeroinfl(number ~ year+month+shore | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action = na.exclude) 
i4 <- zeroinfl(number ~ year+month+veg | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action = na.exclude) 
i5 <- zeroinfl(number ~ year+month+veg+shore | month+veg+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action = na.exclude) 
i6 <- zeroinfl(number ~ year+month+veg+shore | year+veg+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action = na.exclude) 
i7 <- zeroinfl(number ~ year+month+veg+shore | year+month+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action = na.exclude) 
i8 <- zeroinfl(number ~ year+month+veg+shore | year+month+veg+bottom, dist='negbin', link="logit", data=ir.ad, na.action = na.exclude) 
i9 <- zeroinfl(number ~ year+month+veg+shore | year+month+veg+shore, dist='negbin', link="logit", data=ir.ad, na.action = na.exclude) 

lrtest(i, i1)
lrtest(i, i2)
lrtest(i, i3)
lrtest(i, i4) #equal so drop shore from count
lrtest(i, i5) 
lrtest(i, i6)
lrtest(i, i7)
lrtest(i, i8) #drop shore from binary
lrtest(i, i9) #drop bottom from binary

#drop shore from count
i <- zeroinfl(number ~ year+month+veg | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action = na.exclude) 

i1 <- zeroinfl(number ~ month+veg | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action = na.exclude) 
i2 <- zeroinfl(number ~ year+veg | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action = na.exclude) 
i3 <- zeroinfl(number ~ year+month| year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action = na.exclude) 
i4 <- zeroinfl(number ~ year+month+veg |month+veg+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action = na.exclude) 
i5 <- zeroinfl(number ~ year+month+veg | year+veg+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action = na.exclude) 
i6 <- zeroinfl(number ~ year+month+veg | year+month+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action = na.exclude) 
i7 <- zeroinfl(number ~ year+month+veg | year+month+veg+bottom, dist='negbin', link="logit", data=ir.ad, na.action = na.exclude) 
i8 <- zeroinfl(number ~ year+month+veg | year+month+veg+shore, dist='negbin', link="logit", data=ir.ad, na.action = na.exclude) 

lrtest(i, i1)
lrtest(i, i2)
lrtest(i, i3)
lrtest(i, i4) 
lrtest(i, i5) 
lrtest(i, i6)
#lrtest(i, i7)
#lrtest(i, i8)

#all reamining covariates are significant

AIC (i, i1, i2, i3, i4, i5, i6)
#model i is best

N<-nrow(ir.ad)
EZIP <- resid(i,type="pearson")
Dispersion <- sum(EZIP^2)/(N-67)	#==>be sure to change value for degrees of freedom based on summary output
Dispersion #1.59

BestModelIR.ad_zinb <- i


#3. ZANB ####
#AP adult model validation ZANB ####
A = hurdle(number ~ year+month+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ap.ad, na.action=na.exclude)

A1 = hurdle(number ~ month+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ap.ad, na.action=na.exclude)
A2 = hurdle(number ~ year+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ap.ad, na.action=na.exclude)
A3 = hurdle(number ~ year+month+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ap.ad, na.action=na.exclude)
A4 = hurdle(number ~ year+month+veg+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ap.ad, na.action=na.exclude)
A5 = hurdle(number ~ year+month+veg+bottom | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ap.ad, na.action=na.exclude)
A6 = hurdle(number ~ year+month+veg+bottom+shore | month+veg+bottom+shore, dist='negbin', link="logit", data=ap.ad, na.action=na.exclude)
A7 = hurdle(number ~ year+month+veg+bottom+shore | year+veg+bottom+shore, dist='negbin', link="logit", data=ap.ad, na.action=na.exclude)
A8 = hurdle(number ~ year+month+veg+bottom+shore | year+month+bottom+shore, dist='negbin', link="logit", data=ap.ad, na.action=na.exclude)
A9 = hurdle(number ~ year+month+veg+bottom+shore | year+month+veg+shore, dist='negbin', link="logit", data=ap.ad, na.action=na.exclude)
A10 = hurdle(number ~ year+month+veg+bottom+shore | year+month+veg+bottom, dist='negbin', link="logit", data=ap.ad, na.action=na.exclude)

lrtest(A, A1)
lrtest(A, A2)
lrtest(A, A3) #equal so drop veg from count
lrtest(A, A4)
lrtest(A, A5) #equal so drop shore from count
lrtest(A, A6)
lrtest(A, A7)
lrtest(A, A8)
lrtest(A, A9) #equal so frop bottom from binary
lrtest(A, A10)

#equal so drop veg from count
A = hurdle(number ~ year+month+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ap.ad, na.action=na.exclude)

A1 = hurdle(number ~ month+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ap.ad, na.action=na.exclude)
A2 = hurdle(number ~ year+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ap.ad, na.action=na.exclude)
A3 = hurdle(number ~ year+month+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ap.ad, na.action=na.exclude)
A4 = hurdle(number ~ year+month+bottom | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ap.ad, na.action=na.exclude)
A5 = hurdle(number ~ year+month+bottom+shore | month+veg+bottom+shore, dist='negbin', link="logit", data=ap.ad, na.action=na.exclude)
A6 = hurdle(number ~ year+month+bottom+shore | year+veg+bottom+shore, dist='negbin', link="logit", data=ap.ad, na.action=na.exclude)
A7 = hurdle(number ~ year+month+bottom+shore | year+month+bottom+shore, dist='negbin', link="logit", data=ap.ad, na.action=na.exclude)
A8 = hurdle(number ~ year+month+bottom+shore | year+month+veg+shore, dist='negbin', link="logit", data=ap.ad, na.action=na.exclude)
A9 = hurdle(number ~ year+month+bottom+shore | year+month+veg+bottom, dist='negbin', link="logit", data=ap.ad, na.action=na.exclude)

lrtest(A, A1)
lrtest(A, A2)
lrtest(A, A3) 
lrtest(A, A4) #equal so drop shore from count
lrtest(A, A5) 
lrtest(A, A6)
lrtest(A, A7)
lrtest(A, A8)#equal so frop bottom from binary
lrtest(A, A9) 

#equal so drop shore from count
A = hurdle(number ~ year+month+bottom | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ap.ad, na.action=na.exclude)

A1 = hurdle(number ~ month+bottom | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ap.ad, na.action=na.exclude)
A2 = hurdle(number ~ year+bottom | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ap.ad, na.action=na.exclude)
A3 = hurdle(number ~ year+month| year+month+veg+bottom+shore, dist='negbin', link="logit", data=ap.ad, na.action=na.exclude)
A4 = hurdle(number ~ year+month+bottom | month+veg+bottom+shore, dist='negbin', link="logit", data=ap.ad, na.action=na.exclude)
A5 = hurdle(number ~ year+month+bottom | year+veg+bottom+shore, dist='negbin', link="logit", data=ap.ad, na.action=na.exclude)
A6 = hurdle(number ~ year+month+bottom | year+month+bottom+shore, dist='negbin', link="logit", data=ap.ad, na.action=na.exclude)
A7 = hurdle(number ~ year+month+bottom | year+month+veg+shore, dist='negbin', link="logit", data=ap.ad, na.action=na.exclude)
A8 = hurdle(number ~ year+month+bottom | year+month+veg+bottom, dist='negbin', link="logit", data=ap.ad, na.action=na.exclude)

lrtest(A, A1)
lrtest(A, A2)
lrtest(A, A3) 
lrtest(A, A4) 
lrtest(A, A5) 
lrtest(A, A6)
lrtest(A, A7)#equal so frop bottom from binary
lrtest(A, A8)

#equal so drop bottom from binary
A = hurdle(number ~ year+month+bottom | year+month+veg+shore, dist='negbin', link="logit", data=ap.ad, na.action=na.exclude)

A1 = hurdle(number ~ month+bottom | year+month+veg+shore, dist='negbin', link="logit", data=ap.ad, na.action=na.exclude)
A2 = hurdle(number ~ year+bottom | year+month+veg+shore, dist='negbin', link="logit", data=ap.ad, na.action=na.exclude)
A3 = hurdle(number ~ year+month | year+month+veg+shore, dist='negbin', link="logit", data=ap.ad, na.action=na.exclude)
A4 = hurdle(number ~ year+month+bottom | month+veg+shore, dist='negbin', link="logit", data=ap.ad, na.action=na.exclude)
A5 = hurdle(number ~ year+month+bottom | year+veg+shore, dist='negbin', link="logit", data=ap.ad, na.action=na.exclude)
A6 = hurdle(number ~ year+month+bottom | year+month+shore, dist='negbin', link="logit", data=ap.ad, na.action=na.exclude)
A7 = hurdle(number ~ year+month+bottom | year+month+veg, dist='negbin', link="logit", data=ap.ad, na.action=na.exclude)

lrtest(A, A1)
lrtest(A, A2)
lrtest(A, A3) 
lrtest(A, A4) 
lrtest(A, A5) 
lrtest(A, A6)
lrtest(A, A7)

#all covariates are significant
#model A is best

AIC (A, A1, A2, A3, A4, A5, A6)
#model A is best

N<-nrow(ap.ad)
EZIP <- resid(A,type="pearson")
Dispersion <- sum(EZIP^2)/(N-63)	#==>be sure to change value for degrees of freedom based on summary output
Dispersion #1.45

BestModelAP.ad_zinb <- A


#CK adult model validation ZANB ####
C = hurdle(number ~ year+month+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action=na.exclude)

C1 = hurdle(number ~ month+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action=na.exclude)
C2 = hurdle(number ~ year+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action=na.exclude)
C3 = hurdle(number ~ year+month+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action=na.exclude)
C4 = hurdle(number ~ year+month+veg+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action=na.exclude)
C5 = hurdle(number ~ year+month+veg+bottom | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action=na.exclude)
C6 = hurdle(number ~ year+month+veg+bottom+shore | month+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action=na.exclude)
C7 = hurdle(number ~ year+month+veg+bottom+shore | year+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action=na.exclude)
C8 = hurdle(number ~ year+month+veg+bottom+shore | year+month+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action=na.exclude)
C9 = hurdle(number ~ year+month+veg+bottom+shore | year+month+veg+shore, dist='negbin', link="logit", data=ck.ad, na.action=na.exclude)
C10 = hurdle(number ~ year+month+veg+bottom+shore | year+month+veg+bottom, dist='negbin', link="logit", data=ck.ad, na.action=na.exclude)

lrtest(C, C1)
lrtest(C, C2)
lrtest(C, C3)
lrtest(C, C4) #equal so drop bottom from count
lrtest(C, C5) #equal so drop shore from count
lrtest(C, C6)
lrtest(C, C7)
lrtest(C, C8)
lrtest(C, C9)
lrtest(C, C10)

#drop bottom from count
C = hurdle(number ~ year+month+veg+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action=na.exclude)

C1 = hurdle(number ~ month+veg+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action=na.exclude)
C2 = hurdle(number ~ year+veg+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action=na.exclude)
C3 = hurdle(number ~ year+month+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action=na.exclude)
C4 = hurdle(number ~ year+month+veg | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action=na.exclude)
C5 = hurdle(number ~ year+month+veg+shore | month+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action=na.exclude)
C6 = hurdle(number ~ year+month+veg+shore | year+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action=na.exclude)
C7 = hurdle(number ~ year+month+veg+shore | year+month+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action=na.exclude)
C8 = hurdle(number ~ year+month+veg+shore | year+month+veg+shore, dist='negbin', link="logit", data=ck.ad, na.action=na.exclude)
C9 = hurdle(number ~ year+month+veg+shore | year+month+veg+bottom, dist='negbin', link="logit", data=ck.ad, na.action=na.exclude)

lrtest(C, C1)
lrtest(C, C2)
lrtest(C, C3)
lrtest(C, C4) #drop shore from count
lrtest(C, C5) 
lrtest(C, C6)
lrtest(C, C7)
lrtest(C, C8)
lrtest(C, C9)

#drop shore from count
C = hurdle(number ~ year+month+veg | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action=na.exclude)

C1 = hurdle(number ~ month+veg | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action=na.exclude)
C2 = hurdle(number ~ year+veg | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action=na.exclude)
C3 = hurdle(number ~ year+month| year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action=na.exclude)
C4 = hurdle(number ~ year+month+veg | month+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action=na.exclude)
C5 = hurdle(number ~ year+month+veg | year+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action=na.exclude)
C6 = hurdle(number ~ year+month+veg | year+month+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action=na.exclude)
C7 = hurdle(number ~ year+month+veg | year+month+veg+shore, dist='negbin', link="logit", data=ck.ad, na.action=na.exclude)
C8 = hurdle(number ~ year+month+veg | year+month+veg+bottom, dist='negbin', link="logit", data=ck.ad, na.action=na.exclude)

lrtest(C, C1)
lrtest(C, C2)
lrtest(C, C3)
lrtest(C, C4) 
lrtest(C, C5) 
lrtest(C, C6)
lrtest(C, C7)
lrtest(C, C8)

#all covariates are significant
#model C is best

AIC (C, C1, C2, C3, C4, C5, C6)
#model C is best

N<-nrow(ck.ad)
EZIP <- resid(C,type="pearson")
Dispersion <- sum(EZIP^2)/(N-67)	#==>be sure to change value for degrees of freedom based on summary output
Dispersion #1.10

BestModelCK.ad_zinb <- C

#TB adult model validation ZANB ####
t = hurdle(number ~ year+month+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)

t1 = hurdle(number ~ month+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)
t2 = hurdle(number ~ year+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)
t3 = hurdle(number ~ year+month+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)
t4 = hurdle(number ~ year+month+veg+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)
t5 = hurdle(number ~ year+month+veg+bottom | year+month+veg+bottom+shore, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)
t6 = hurdle(number ~ year+month+veg+bottom+shore | month+veg+bottom+shore, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)
t7 = hurdle(number ~ year+month+veg+bottom+shore | year+veg+bottom+shore, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)
t8 = hurdle(number ~ year+month+veg+bottom+shore | year+month+bottom+shore, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)
t9 = hurdle(number ~ year+month+veg+bottom+shore | year+month+veg+shore, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)
t10 = hurdle(number ~ year+month+veg+bottom+shore | year+month+veg+bottom, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)

lrtest(t, t1)
lrtest(t, t2)
lrtest(t, t3)
lrtest(t, t4) #equal but drop bottom from count
lrtest(t, t5)
lrtest(t, t6)
lrtest(t, t7)
lrtest(t, t8) #equal so drop veg from binary
lrtest(t, t9) # equal so drop bottom from binary
lrtest(t, t10) #equa so drop shore from binary

#drop bottom from count
t = hurdle(number ~ year+month+veg+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)

t1 = hurdle(number ~ month+veg+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)
t2 = hurdle(number ~ year+veg+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)
t3 = hurdle(number ~ year+month+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)
t4 = hurdle(number ~ year+month+veg | year+month+veg+bottom+shore, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)
t5 = hurdle(number ~ year+month+veg+shore | month+veg+bottom+shore, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)
t6 = hurdle(number ~ year+month+veg+shore | year+veg+bottom+shore, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)
t7 = hurdle(number ~ year+month+veg+shore | year+month+bottom+shore, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)
t8 = hurdle(number ~ year+month+veg+shore | year+month+veg+shore, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)
t9 = hurdle(number ~ year+month+veg+shore | year+month+veg+bottom, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)

lrtest(t, t1)
lrtest(t, t2)
lrtest(t, t3)
lrtest(t, t4) 
lrtest(t, t5)
lrtest(t, t6)
lrtest(t, t7)#equal so drop veg from binary
lrtest(t, t8) # equal so drop bottom from binary
lrtest(t, t9) #equal so drop shore from binary

#drop veg from binary
t = hurdle(number ~ year+month+veg+shore | year+month+bottom+shore, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)

t1 = hurdle(number ~ month+veg+shore | year+month+bottom+shore, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)
t2 = hurdle(number ~ year+veg+shore | year+month+bottom+shore, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)
t3 = hurdle(number ~ year+month+shore | year+month+bottom+shore, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)
t4 = hurdle(number ~ year+month+veg| year+month+bottom+shore, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)
t5 = hurdle(number ~ year+month+veg+shore | month+bottom+shore, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)
t6 = hurdle(number ~ year+month+veg+shore | year+bottom+shore, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)
t7 = hurdle(number ~ year+month+veg+shore | year+month+shore, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)
t8 = hurdle(number ~ year+month+veg+shore | year+month+bottom, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)

lrtest(t, t1)
lrtest(t, t2)
lrtest(t, t3)
lrtest(t, t4) 
lrtest(t, t5)
lrtest(t, t6)
lrtest(t, t7)
lrtest(t, t8) #equal so drop shore from binary 

#drop shore from binary
t = hurdle(number ~ year+month+veg+shore | year+month+bottom, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)

t1 = hurdle(number ~ month+veg+shore | year+month+bottom, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)
t2 = hurdle(number ~ year+veg+shore | year+month+bottom, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)
t3 = hurdle(number ~ year+month+shore | year+month+bottom, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)
t4 = hurdle(number ~ year+month+veg | year+month+bottom, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)
t5 = hurdle(number ~ year+month+veg+shore | month+bottom, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)
t6 = hurdle(number ~ year+month+veg+shore | year+bottom, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)
t7 = hurdle(number ~ year+month+veg+shore | year+month, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)

lrtest(t, t1)
lrtest(t, t2)
lrtest(t, t3)
lrtest(t, t4) 
lrtest(t, t5)
lrtest(t, t6)
lrtest(t, t7) #drop bottom from binary

#drop bottom from binary
t = hurdle(number ~ year+month+veg+shore | year+month, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)

t1 = hurdle(number ~ month+veg+shore | year+month, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)
t2 = hurdle(number ~ year+veg+shore | year+month, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)
t3 = hurdle(number ~ year+month+shore | year+month, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)
t4 = hurdle(number ~ year+month+veg | year+month, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)
t5 = hurdle(number ~ year+month+veg+shore | month, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)
t6 = hurdle(number ~ year+month+veg+shore | year, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)

lrtest(t, t1)
lrtest(t, t2)
lrtest(t, t3)
lrtest(t, t4) 
lrtest(t, t5)
lrtest(t, t6)

#all covariates are significant
#model t is best

AIC (t, t1, t2, t3, t4, t5, t6)
#model C is best

N<-nrow(tb.ad)
EZIP <- resid(t,type="pearson")
Dispersion <- sum(EZIP^2)/(N-68)	#==>be sure to change value for degrees of freedom based on summary output
Dispersion #1.13

BestModelTB.ad_zinb <- t


#CH adult model validation ZANB ####
#ch = hurdle(number ~ year+month+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ch.ad, na.action=na.exclude)

#ch1 = hurdle(number ~ month+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ch.ad, na.action=na.exclude)
#ch2 = hurdle(number ~ year+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ch.ad, na.action=na.exclude)
#ch3 = hurdle(number ~ year+month+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ch.ad, na.action=na.exclude)
ch4 = hurdle(number ~ year+month+veg+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ch.ad, na.action=na.exclude)
ch5 = hurdle(number ~ year+month+veg+bottom | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ch.ad, na.action=na.exclude)
#ch6 = hurdle(number ~ year+month+veg+bottom+shore | month+veg+bottom+shore, dist='negbin', link="logit", data=ch.ad, na.action=na.exclude)
#ch7 = hurdle(number ~ year+month+veg+bottom+shore | year+veg+bottom+shore, dist='negbin', link="logit", data=ch.ad, na.action=na.exclude)
#ch8 = hurdle(number ~ year+month+veg+bottom+shore | year+month+bottom+shore, dist='negbin', link="logit", data=ch.ad, na.action=na.exclude)
#ch9 = hurdle(number ~ year+month+veg+bottom+shore | year+month+veg+shore, dist='negbin', link="logit", data=ch.ad, na.action=na.exclude)
#ch10 = hurdle(number ~ year+month+veg+bottom+shore | year+month+veg+bottom, dist='negbin', link="logit", data=ch.ad, na.action=na.exclude)

#many of above did not solve so using ch4 as starting which did solve
ch = hurdle(number ~ year+month+veg+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ch.ad, na.action=na.exclude)

#ch1 = hurdle(number ~ month+veg+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ch.ad, na.action=na.exclude)
#ch2 = hurdle(number ~ year+veg+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ch.ad, na.action=na.exclude)
#ch3 = hurdle(number ~ year+month+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ch.ad, na.action=na.exclude)
ch4 = hurdle(number ~ year+month+veg | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ch.ad, na.action=na.exclude)
ch5 = hurdle(number ~ year+month+veg+shore | month+veg+bottom+shore, dist='negbin', link="logit", data=ch.ad, na.action=na.exclude)
#ch6 = hurdle(number ~ year+month+veg+shore | yearveg+bottom+shore, dist='negbin', link="logit", data=ch.ad, na.action=na.exclude)
ch7 = hurdle(number ~ year+month+veg+shore | year+month+bottom+shore, dist='negbin', link="logit", data=ch.ad, na.action=na.exclude)
ch8 = hurdle(number ~ year+month+veg+shore | year+month+veg+shore, dist='negbin', link="logit", data=ch.ad, na.action=na.exclude)
ch9 = hurdle(number ~ year+month+veg+shore | year+month+veg+bottom, dist='negbin', link="logit", data=ch.ad, na.action=na.exclude)

lrtest(ch, ch4)
lrtest(ch, ch5) #equal so drop year from binary
lrtest(ch, ch7)
lrtest(ch, ch8) #equal sp drop bottom frim binary
lrtest(ch, ch9)

#drop year from binary
ch = hurdle(number ~ year+month+veg+shore | month+veg+bottom+shore, dist='negbin', link="logit", data=ch.ad, na.action=na.exclude)

#ch1 = hurdle(number ~ month+veg+shore | month+veg+bottom+shore, dist='negbin', link="logit", data=ch.ad, na.action=na.exclude)
#ch2 = hurdle(number ~ year+veg+shore | month+veg+bottom+shore, dist='negbin', link="logit", data=ch.ad, na.action=na.exclude)
#ch3 = hurdle(number ~ year+month+shore | month+veg+bottom+shore, dist='negbin', link="logit", data=ch.ad, na.action=na.exclude)
ch4 = hurdle(number ~ year+month+veg | month+veg+bottom+shore, dist='negbin', link="logit", data=ch.ad, na.action=na.exclude)
ch5 = hurdle(number ~ year+month+veg+shore | veg+bottom+shore, dist='negbin', link="logit", data=ch.ad, na.action=na.exclude)
ch6 = hurdle(number ~ year+month+veg+shore | month+bottom+shore, dist='negbin', link="logit", data=ch.ad, na.action=na.exclude)
ch7 = hurdle(number ~ year+month+veg+shore | month+veg+shore, dist='negbin', link="logit", data=ch.ad, na.action=na.exclude)
ch8 = hurdle(number ~ year+month+veg+shore | month+veg+bottom, dist='negbin', link="logit", data=ch.ad, na.action=na.exclude)

lrtest(ch, ch4)
lrtest(ch, ch5) 
lrtest(ch, ch6)
lrtest(ch, ch7) #equal sp drop bottom frim binary
lrtest(ch, ch8) 
lrtest(ch, ch9) #equal so drop shore form bottom 

#drop bottom from binary
ch = hurdle(number ~ year+month+veg+shore | month+veg+shore, dist='negbin', link="logit", data=ch.ad, na.action=na.exclude)

#ch1 = hurdle(number ~ month+veg+shore | month+veg+shore, dist='negbin', link="logit", data=ch.ad, na.action=na.exclude)
#ch2 = hurdle(number ~ year+veg+shore | month+veg+shore, dist='negbin', link="logit", data=ch.ad, na.action=na.exclude)
#ch3 = hurdle(number ~ year+month+shore | month+veg+shore, dist='negbin', link="logit", data=ch.ad, na.action=na.exclude)
ch4 = hurdle(number ~ year+month+veg | month+veg+shore, dist='negbin', link="logit", data=ch.ad, na.action=na.exclude)
ch5 = hurdle(number ~ year+month+veg+shore | veg+shore, dist='negbin', link="logit", data=ch.ad, na.action=na.exclude)
ch6 = hurdle(number ~ year+month+veg+shore | month+shore, dist='negbin', link="logit", data=ch.ad, na.action=na.exclude)
ch7 = hurdle(number ~ year+month+veg+shore | month+veg, dist='negbin', link="logit", data=ch.ad, na.action=na.exclude)

lrtest(ch, ch4)
lrtest(ch, ch5) 
lrtest(ch, ch6)
lrtest(ch, ch7) 

#all covariates are significant
#model ch is best

AIC (ch, ch4,ch5,ch6,ch7)
#model ch is best

N<-nrow(ch.ad)
EZIP <- resid(ch,type="pearson")
Dispersion <- sum(EZIP^2)/(N-59)	#==>be sure to change value for degrees of freedom based on summary output
Dispersion #1.12

BestModelCH.ad_zinb <- ch

#JX adult model validation ZANB ####
jx = hurdle(number ~ year+month+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)

jx1 = hurdle(number ~ month+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)
jx2 = hurdle(number ~ year+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)
jx3 = hurdle(number ~ year+month+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)
jx4 = hurdle(number ~ year+month+veg+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)
jx5 = hurdle(number ~ year+month+veg+bottom | year+month+veg+bottom+shore, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)
jx6 = hurdle(number ~ year+month+veg+bottom+shore | month+veg+bottom+shore, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)
jx7 = hurdle(number ~ year+month+veg+bottom+shore | year+veg+bottom+shore, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)
jx8 = hurdle(number ~ year+month+veg+bottom+shore | year+month+bottom+shore, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)
jx9 = hurdle(number ~ year+month+veg+bottom+shore | year+month+veg+shore, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)
jx10 = hurdle(number ~ year+month+veg+bottom+shore | year+month+veg+bottom, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)

lrtest(jx, jx1)
lrtest(jx, jx2)
lrtest(jx, jx3) #equal so drop veg from count
lrtest(jx, jx4) #equal so drop bottom from count
lrtest(jx, jx5) #equal so drop shore from count
lrtest(jx, jx6) #equal sp drop year from count
lrtest(jx, jx7)
lrtest(jx, jx8)
lrtest(jx, jx9)
lrtest(jx, jx10)

#equal so drop veg from count
jx = hurdle(number ~ year+month+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)

jx1 = hurdle(number ~ month+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)
jx2 = hurdle(number ~ year+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)
jx3 = hurdle(number ~ year+month+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)
jx4 = hurdle(number ~ year+month+bottom | year+month+veg+bottom+shore, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)
jx5 = hurdle(number ~ year+month+bottom+shore | month+veg+bottom+shore, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)
jx6 = hurdle(number ~ year+month+bottom+shore | year+veg+bottom+shore, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)
jx7 = hurdle(number ~ year+month+bottom+shore | year+month+bottom+shore, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)
jx8 = hurdle(number ~ year+month+bottom+shore | year+month+veg+shore, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)
jx9 = hurdle(number ~ year+month+bottom+shore | year+month+veg+bottom, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)

lrtest(jx, jx1)
lrtest(jx, jx2)
lrtest(jx, jx3) #equal so drop bottom from count
lrtest(jx, jx4) #equal so drop shore from count
lrtest(jx, jx5) #equal so drop year from count
lrtest(jx, jx6)
lrtest(jx, jx7)
lrtest(jx, jx8)
lrtest(jx, jx9)

#equal so drop bottom from count
jx = hurdle(number ~ year+month+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)

jx1 = hurdle(number ~ month+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)
jx2 = hurdle(number ~ year+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)
jx3 = hurdle(number ~ year+month| year+month+veg+bottom+shore, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)
jx4 = hurdle(number ~ year+month+shore | month+veg+bottom+shore, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)
jx5 = hurdle(number ~ year+month+shore | year+veg+bottom+shore, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)
jx6 = hurdle(number ~ year+month+shore | year+month+bottom+shore, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)
jx7 = hurdle(number ~ year+month+shore | year+month+veg+shore, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)
jx8 = hurdle(number ~ year+month+shore | year+month+veg+bottom, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)

lrtest(jx, jx1)
lrtest(jx, jx2)
lrtest(jx, jx3) #equal so drop shore from count
lrtest(jx, jx4) #equal so drop year from binary
lrtest(jx, jx5) 
lrtest(jx, jx6)
lrtest(jx, jx7)
lrtest(jx, jx8)

#equal so drop shore from count
jx = hurdle(number ~ year+month | year+month+veg+bottom+shore, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)

jx1 = hurdle(number ~ month | year+month+veg+bottom+shore, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)
jx2 = hurdle(number ~ year | year+month+veg+bottom+shore, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)
jx3 = hurdle(number ~ year+month |month+veg+bottom+shore, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)
jx4 = hurdle(number ~ year+month | year+veg+bottom+shore, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)
jx5 = hurdle(number ~ year+month | year+month+bottom+shore, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)
jx6 = hurdle(number ~ year+month | year+month+veg+shore, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)
jx7 = hurdle(number ~ year+month | year+month+veg+bottom, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)

lrtest(jx, jx1)
lrtest(jx, jx2)
lrtest(jx, jx3) #equal so drop year from binary
lrtest(jx, jx4) 
lrtest(jx, jx5) 
lrtest(jx, jx6)
lrtest(jx, jx7)

#equal so drop year from binary
jx = hurdle(number ~ year+month | month+veg+bottom+shore, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)

jx1 = hurdle(number ~ month | month+veg+bottom+shore, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)
jx2 = hurdle(number ~ year | month+veg+bottom+shore, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)
jx3 = hurdle(number ~ year+month | veg+bottom+shore, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)
jx4 = hurdle(number ~ year+month | month+bottom+shore, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)
jx5 = hurdle(number ~ year+month | month+veg+shore, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)
jx6 = hurdle(number ~ year+month | month+veg+bottom, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)

lrtest(jx, jx1)
lrtest(jx, jx2)
lrtest(jx, jx3) 
lrtest(jx, jx4) 
lrtest(jx, jx5) 
lrtest(jx, jx6)

#all covariates are significant
#model jx is best

AIC (jx, jx4,jx5,jx6,jx7)
#model jx is best

N<-nrow(jx.ad)
EZIP <- resid(jx,type="pearson")
Dispersion <- sum(EZIP^2)/(N-43)	#==>be sure to change value for degrees of freedom based on summary output
Dispersion #1.21

BestModelJX.ad_zinb <- jx


#IR adult model validation ZANB ####
ir = hurdle(number ~ year+month+veg+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)

ir1 = hurdle(number ~ month+veg+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)
ir2 = hurdle(number ~ year+veg+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)
ir3 = hurdle(number ~ year+month+shore+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)
ir4 = hurdle(number ~ year+month+veg+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)
ir5 = hurdle(number ~ year+month+veg+shore| year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)
ir6 = hurdle(number ~ year+month+veg+shore+bottom | month+veg+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)
ir7 = hurdle(number ~ year+month+veg+shore+bottom | year+veg+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)
ir8 = hurdle(number ~ year+month+veg+shore+bottom | year+month+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)
ir9 = hurdle(number ~ year+month+veg+shore+bottom | year+month+veg+bottom, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)
ir10 = hurdle(number ~ year+month+veg+shore+bottom | year+month+veg+shore, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)

lrtest(ir, ir1)
lrtest(ir, ir2)
lrtest(ir, ir3)
lrtest(ir, ir4) #equal so drop shore from count
lrtest(ir, ir5) #equal so drop bottom from count
lrtest(ir, ir6)
lrtest(ir, ir7)
lrtest(ir, ir8) #equal so drop veg from count
lrtest(ir, ir9) #equal so drop shore from count
lrtest(ir, ir10) #equal so drop bottom from count

#drop shore from count
ir = hurdle(number ~ year+month+veg+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)

ir1 = hurdle(number ~ month+veg+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)
ir2 = hurdle(number ~ year+veg+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)
ir3 = hurdle(number ~ year+month+bottom | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)
ir4 = hurdle(number ~ year+month+veg| year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)
ir5 = hurdle(number ~ year+month+veg+bottom | month+veg+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)
ir6 = hurdle(number ~ year+month+veg+bottom | year+veg+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)
ir7 = hurdle(number ~ year+month+veg+bottom | year+month+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)
ir8 = hurdle(number ~ year+month+veg+bottom | year+month+veg+bottom, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)
ir9 = hurdle(number ~ year+month+veg+bottom | year+month+veg+shore, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)

lrtest(ir, ir1)
lrtest(ir, ir2)
lrtest(ir, ir3)
lrtest(ir, ir4) #equal so drop bottom from count
lrtest(ir, ir5) 
lrtest(ir, ir6)
lrtest(ir, ir7) #equal sp drp veg from binary
lrtest(ir, ir8) #equal so drop shore from binary
lrtest(ir, ir9) #equal so drop bottom from binary 

#drop bottom from count
ir = hurdle(number ~ year+month+veg | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)

ir1 = hurdle(number ~ month+veg | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)
ir2 = hurdle(number ~ year+veg | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)
ir3 = hurdle(number ~ year+month | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)
ir4 = hurdle(number ~ year+month+veg | month+veg+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)
ir5 = hurdle(number ~ year+month+veg | year+veg+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)
ir6 = hurdle(number ~ year+month+veg | year+month+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)
ir7 = hurdle(number ~ year+month+veg | year+month+veg+bottom, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)
ir8 = hurdle(number ~ year+month+veg | year+month+veg+shore, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)

lrtest(ir, ir1)
lrtest(ir, ir2)
lrtest(ir, ir3)
lrtest(ir, ir4) 
lrtest(ir, ir5) 
lrtest(ir, ir6) #drop veg from binary
lrtest(ir, ir7) #equal sp drp shore from binary
lrtest(ir, ir8) #equal so drop bottom from binary

#drop veg from binary
ir = hurdle(number ~ year+month+veg | year+month+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)

ir1 = hurdle(number ~ month+veg | year+month+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)
ir2 = hurdle(number ~ year+veg | year+month+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)
ir3 = hurdle(number ~ year+month | year+month+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)
ir4 = hurdle(number ~ year+month+veg | month+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)
ir5 = hurdle(number ~ year+month+veg | year+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)
ir6 = hurdle(number ~ year+month+veg | year+month+bottom, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)
ir7 = hurdle(number ~ year+month+veg | year+month+shore, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)

lrtest(ir, ir1)
lrtest(ir, ir2)
lrtest(ir, ir3)
lrtest(ir, ir4) 
lrtest(ir, ir5) 
lrtest(ir, ir6) #drop shore from binary
lrtest(ir, ir7) #equal sp drp bottom from binary

#drop shore from binary
ir = hurdle(number ~ year+month+veg | year+month+bottom, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)

ir1 = hurdle(number ~ month+veg | year+month+bottom, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)
ir2 = hurdle(number ~ year+veg | year+month+bottom, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)
ir3 = hurdle(number ~ year+month| year+month+bottom, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)
ir4 = hurdle(number ~ year+month+veg | month+bottom, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)
ir5 = hurdle(number ~ year+month+veg | year+bottom, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)
ir6 = hurdle(number ~ year+month+veg | year+month, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)

lrtest(ir, ir1)
lrtest(ir, ir2)
lrtest(ir, ir3)
lrtest(ir, ir4) 
lrtest(ir, ir5) 
lrtest(ir, ir6) #drop bottom from binary

#drop bottom from binary
ir = hurdle(number ~ year+month+veg | year+month, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)

ir1 = hurdle(number ~ month+veg | year+month, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)
ir2 = hurdle(number ~ year+veg | year+month, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)
ir3 = hurdle(number ~ year+month| year+month, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)
ir4 = hurdle(number ~ year+month+veg | month, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)
ir5 = hurdle(number ~ year+month+veg | year, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)

lrtest(ir, ir1)
lrtest(ir, ir2)
lrtest(ir, ir3)
lrtest(ir, ir4) 
lrtest(ir, ir5) 

#all covariates are significant
#model ir is best

AIC (ir,ir1, ir2, ir3, ir4, ir5)
#model ir is best

N<-nrow(ir.ad)
EZIP <- resid(ir,type="pearson")
Dispersion <- sum(EZIP^2)/(N-62)	#==>be sure to change value for degrees of freedom based on summary output
Dispersion #1.39

BestModelIR.ad_zinb <- ir

##### ASSIGN FINAL MODELS_ADULT ###### 

# BOTH the zinb and the zanb produce standard error estimates that are huge! 
#zinb
A= zeroinfl(number ~ year+bottom+month+shore | year+month+veg+shore, dist='negbin', link="logit", data=ap.ad, na.action = na.exclude)
C= zeroinfl(number ~ year+month+veg | bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action = na.exclude)
t <- zeroinfl(number ~ year+month+veg+bottom+shore | year+month+veg+bottom+shore, dist='negbin', link="logit", data=tb.ad, na.action = na.exclude) 
ch <- zeroinfl(number ~ year+shore+bottom | year+month+bottom, dist='negbin', link="logit", data=ch.ad, na.action = na.exclude) 
j2 <- zeroinfl(number ~ year | year+veg+shore, dist='negbin', link="logit", data=jx.ad, na.action = na.exclude) 
i <- zeroinfl(number ~ year+month+veg | year+month+veg+shore+bottom, dist='negbin', link="logit", data=ir.ad, na.action = na.exclude) 

BestModelAP.ad_zinb <- A
BestModelCK.ad_zinb <- C
BestModelTB.ad_zinb <- t
BestModelCH.ad_zinb <- ch
BestModelJX.ad_zinb <- j2
BestModelIR.ad_zinb <- i

#zanb
A = hurdle(number ~ year+month+bottom | year+month+veg+shore, dist='negbin', link="logit", data=ap.ad, na.action=na.exclude)
C = hurdle(number ~ year+month+veg | year+month+veg+bottom+shore, dist='negbin', link="logit", data=ck.ad, na.action=na.exclude)
t = hurdle(number ~ year+month+veg+shore | year+month, dist='negbin', link="logit", data=tb.ad, na.action=na.exclude)
ch = hurdle(number ~ year+month+veg+shore | month+veg+shore, dist='negbin', link="logit", data=ch.ad, na.action=na.exclude)
jx = hurdle(number ~ year+month | month+veg+bottom+shore, dist='negbin', link="logit", data=jx.ad, na.action=na.exclude)
ir = hurdle(number ~ year+month+veg | year+month, dist='negbin', link="logit", data=ir.ad, na.action=na.exclude)

BestModelAP.ad_zanb <- A
BestModelCK.ad_zanb <- C
BestModelTB.ad_zanb <- t
BestModelCH.ad_zanb <- ch
BestModelJX.ad_zanb <- jx
BestModelIR.ad_zanb <- ir

#Compare ZANB to ZINB ####
AIC(BestModelAP.ad_zanb, BestModelAP.ad_zinb, apP.ad, apNB.ad1) #zanb is better
AIC(BestModelCK.ad_zanb, BestModelCK.ad_zinb, ckP.ad, ckNB.ad1) #zanb is better
AIC(BestModelTB.ad_zanb, BestModelTB.ad_zinb, tbP.ad, tbNB.ad1) #zanb is better
AIC(BestModelCH.ad_zanb, BestModelCH.ad_zinb, chNB.ad1) #zanb is better
AIC(BestModelJX.ad_zanb, BestModelJX.ad_zinb, jxNB.ad1) #zanb is better
AIC(BestModelIR.ad_zanb, BestModelIR.ad_zinb, irNB.ad1) #zanb is better

##### DETERMINE LEAST SQUARE MEANS_ADULT######
# Same thing as covariate adjusted means. Basically, determine the mean value of total positive numbers 
# of catch per year controlling for covariates (in this case it would be veg and shore variables). 
# Use lsmeans CRAN document. 


library(dplyr)
lsm.APad       <- lsmeans(BestModelAP.ad_zanb, "year", data = ap.ad, mode="response")
estimate.APad <- as.data.frame(transform(summary(lsm.APad))) %>% select(year, lsmean, SE) 
#lsmeans with alternative model type (ZINB for AP)
alt.APad <- as.data.frame(transform(summary(lsmeans(BestModelAP.ad_zinb,"year", data = ap.ad, mode="response"))))  %>% select(year, lsmean, SE) %>% mutate(alternative.mean=lsmean) %>% select(-lsmean)
estimate.APad$year <- as.numeric(as.character(estimate.APad$year))

lsm.CKad       <- lsmeans(BestModelCK.ad_zanb, "year", data = ck.ad, mode="response")
estimate.CKad <- as.data.frame(transform(summary(lsm.CKad))) %>% select(year, lsmean, SE) 
#lsmeans with alternative model type (ZANB for CK)
alt.CKad <- as.data.frame(transform(summary(lsmeans(BestModelCK.ad_zinb,"year", data = ck.ad, mode="response"))))  %>% select(year, lsmean, SE) %>% mutate(alternative.mean=lsmean) %>% select(-lsmean)
estimate.CKad$year <- as.numeric(as.character(estimate.CKad$year))

lsm.TBad       <- lsmeans(BestModelTB.ad_zanb, "year", data = tb.ad, mode="response")
estimate.TBad <- as.data.frame(transform(summary(lsm.TBad))) %>% select(year, lsmean, SE) 
#lsmeans with alternative model type (ZANB for TB)
alt.TBad <- as.data.frame(transform(summary(lsmeans(BestModelTB.ad_zinb,"year", data = tb.ad, mode="response"))))  %>% select(year, lsmean, SE) %>% mutate(alternative.mean=lsmean) %>% select(-lsmean)
estimate.TBad$year <- as.numeric(as.character(estimate.TBad$year))

lsm.CHad       <- lsmeans(BestModelCH.ad_zanb, "year", data = ch.ad, mode="response")
estimate.CHad <- as.data.frame(transform(summary(lsm.CHad))) %>% select(year, lsmean, SE) 
#lsmeans with alternative model type (ZANB for CH)
alt.CHad <- as.data.frame(transform(summary(lsmeans(BestModelCH.ad_zinb,"year", data = ch.ad, mode="response"))))  %>% select(year, lsmean, SE) %>% mutate(alternative.mean=lsmean) %>% select(-lsmean)
estimate.CHad$year <- as.numeric(as.character(estimate.CHad$year))

lsm.JXad       <- lsmeans(BestModelJX.ad_zanb, "year", data = jx.ad, mode="response")
estimate.JXad <- as.data.frame(transform(summary(lsm.JXad))) %>% select(year, lsmean, SE) 
#lsmeans with alternative model type (ZANB for JX)
alt.JXad <- as.data.frame(transform(summary(lsmeans(BestModelJX.ad_zinb,"year", data = jx.ad, mode="response"))))  %>% select(year, lsmean, SE) %>% mutate(alternative.mean=lsmean) %>% select(-lsmean)
estimate.JXad$year <- as.numeric(as.character(estimate.JXad$year))

lsm.IRad       <- lsmeans(BestModelIR.ad_zanb, "year", data = ir.ad, mode="response")
estimate.IRad <- as.data.frame(transform(summary(lsm.IRad))) %>% select(year, lsmean, SE) 
#lsmeans with alternative model type (ZANB for IR)
alt.IRad <- as.data.frame(transform(summary(lsmeans(BestModelIR.ad_zinb,"year", data = ir.ad, mode="response"))))  %>% select(year, lsmean, SE) %>% mutate(alternative.mean=lsmean) %>% select(-lsmean)
estimate.IRad$year <- as.numeric(as.character(estimate.IRad$year))


NominalMeanAP  = as.vector(tapply(ap.ad$number, ap.ad$year, mean))
NominalMeanCK  = as.vector(tapply(ck.ad$number, ck.ad$year, mean))
NominalMeanTB  = as.vector(tapply(tb.ad$number, tb.ad$year, mean))
NominalMeanCH  = as.vector(tapply(ch.ad$number, ch.ad$year, mean))
NominalMeanJX  = as.vector(tapply(jx.ad$number, jx.ad$year, mean))
NominalMeanIR  = as.vector(tapply(ir.ad$number, ir.ad$year, mean))

# PLOT INDICES_adult ####

library(reshape2)
library(ggplot2)

#add on the alternative index 
APadred <- estimate.APad[,c(1,2)] %>% cbind(alt.APad$alternative.mean, NominalMeanAP) %>% melt(id=c("year"))
CKadred <- estimate.CKad[,c(1,2)] %>% cbind(alt.CKad$alternative.mean, NominalMeanCK) %>% melt(id=c("year"))
TBadred <- estimate.TBad[,c(1,2)] %>% cbind(alt.TBad$alternative.mean, NominalMeanTB) %>% melt(id=c("year"))
CHadred <- estimate.CHad[,c(1,2)] %>% cbind(alt.CHad$alternative.mean, NominalMeanCH) %>% melt(id=c("year"))
JXadred <- estimate.JXad[,c(1,2)] %>% cbind(alt.JXad$alternative.mean, NominalMeanJX) %>% melt(id=c("year"))
IRadred <- estimate.IRad[,c(1,2)] %>% cbind(alt.IRad$alternative.mean, NominalMeanIR) %>% melt(id=c("year"))


library(ggplot2)

ggplot(APadred, aes(x=year, y=value, color=variable))+
  geom_line()  

ggplot(CKadred, aes(x=year, y=value, color=variable))+
  geom_line()  

ggplot(TBadred, aes(x=year, y=value, color=variable))+
  geom_line()  

ggplot(CHadred, aes(x=year, y=value, color=variable))+
  geom_line()  

ggplot(JXadred, aes(x=year, y=value, color=variable))+
  geom_line()  

ggplot(IRadred, aes(x=year, y=value, color=variable))+
  geom_line()  



##### EXPORT PREDICTED INDEX (NUMBERS)_ADULT ######
#export to csv _PERSONAL COMPUTER
write.csv(Mean_AP_ad, "~/Desktop/PhD project/Projects/Seatrout/Data/Indices/DeltaMethod Indices/AP_adult_index.csv")
write.csv(Mean_IR_ad, "~/Desktop/PhD project/Projects/Seatrout/Data/Indices/DeltaMethod Indices/CK_adult_index.csv")
write.csv(Mean_JX_ad, "~/Desktop/PhD project/Projects/Seatrout/Data/Indices/DeltaMethod Indices/TB_adult_index.csv")
write.csv(Mean_CH_ad, "~/Desktop/PhD project/Projects/Seatrout/Data/Indices/DeltaMethod Indices/CH_adult_index.csv")
write.csv(Mean_TB_ad, "~/Desktop/PhD project/Projects/Seatrout/Data/Indices/DeltaMethod Indices/JX_adult_index.csv")
write.csv(Mean_CK_ad, "~/Desktop/PhD project/Projects/Seatrout/Data/Indices/DeltaMethod Indices/IR_adult_index.csv")

#export to csv _WORK

write.csv(estimate.APad, "U:/PhD_projectfiles/Exported_R_Datafiles/Indices/UpdatedIndices/AP_adult_index.csv")
write.csv(estimate.CKad, "U:/PhD_projectfiles/Exported_R_Datafiles/Indices/UpdatedIndices/CK_adult_index.csv")
write.csv(estimate.TBad, "U:/PhD_projectfiles/Exported_R_Datafiles/Indices/UpdatedIndices/TB_adult_index.csv")
write.csv(estimate.CHad, "U:/PhD_projectfiles/Exported_R_Datafiles/Indices/UpdatedIndices/CH_adult_index.csv")
write.csv(estimate.JXad, "U:/PhD_projectfiles/Exported_R_Datafiles/Indices/UpdatedIndices/JX_adult_index.csv")
write.csv(estimate.IRad, "U:/PhD_projectfiles/Exported_R_Datafiles/Indices/UpdatedIndices/IR_adult_index.csv")


##### APPLY AGE & WEIGHT SCHEDULE TO OBTAIN PREDICTED SSB OF ADULTS #######
#obtain proportion at age schedule from the output files of ALK_analysis.R

#setwd to get the age proportion data on my work computer
setwd("U:/PhD_projectfiles/Exported_R_Datafiles")
setwd("~/Desktop/PhD project/Projects/Seatrout/Data/Exported R Dataframes")

#Import the datafiles
prop_AP <- read.csv("PropAtAge_APadult_FIMdata.csv", header=T)
prop_CK <- read.csv("PropAtAge_CKadult_FIMdata.csv", header=T)
prop_TB <- read.csv("PropAtAge_TBadult_FIMdata.csv", header=T)
prop_CH <- read.csv("PropAtAge_CHadult_FIMdata.csv", header=T)
prop_JX <- read.csv("PropAtAge_JXadult_FIMdata.csv", header=T)
prop_IR <- read.csv("PropAtAge_IRadult_FIMdata.csv", header=T)

#Apply proportions to the adult indices created above using a for loop
#create empty dataframe for results to be stored
#cycle through rows of the adult index data frame and multiply by the proportion at age

#AP
 num.yr = length(estimate.APad$year)
 num.age =length(prop_AP$Freq)
 Pred_Numbers_at_age_AP_adult <- data.frame(matrix(data=NA, nrow=num.yr, ncol=num.age)) #make a dataframe for the loop to store results in
 for (i in 1:num.yr) {
   Pred_Numbers_at_age_AP_adult[i,] <- estimate.APad$lsmean[i]*prop_AP$Freq
 }
colnames(Pred_Numbers_at_age_AP_adult) <- c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10")

#CK
num.yr = length(estimate.CKad$year)
num.age =length(prop_CK$Freq)
Pred_Numbers_at_age_CK_adult <- data.frame(matrix(data=NA, nrow=num.yr, ncol=num.age)) #make a dataframe for the loop to store results in
for (i in 1:num.yr) {
  Pred_Numbers_at_age_CK_adult[i,] <- estimate.CKad$lsmean[i]*prop_CK$Freq
}
colnames(Pred_Numbers_at_age_CK_adult) <- c("P1", "P2", "P3", "P4", "P5", "P6", "P7") 

#TB
num.yr = length(estimate.TBad$year)
num.age =length(prop_TB$Freq)
Pred_Numbers_at_age_TB_adult <- data.frame(matrix(data=NA, nrow=num.yr, ncol=num.age)) #make a dataframe for the loop to store results in
for (i in 1:num.yr) {
  Pred_Numbers_at_age_TB_adult[i,] <- estimate.TBad$lsmean[i]*prop_TB$Freq
}
colnames(Pred_Numbers_at_age_TB_adult) <- c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9") 

#CH
num.yr = length(estimate.CHad$year)
num.age =length(prop_CH$Freq)
Pred_Numbers_at_age_CH_adult <- data.frame(matrix(data=NA, nrow=num.yr, ncol=num.age)) #make a dataframe for the loop to store results in
for (i in 1:num.yr) {
  Pred_Numbers_at_age_CH_adult[i,] <- estimate.CHad$lsmean[i]*prop_CH$Freq
}
colnames(Pred_Numbers_at_age_CH_adult) <- c("P1", "P2", "P3", "P4", "P5", "P6", "P7")

#JX
num.yr = length(estimate.JXad$year)
num.age =length(prop_JX$Freq)
Pred_Numbers_at_age_JX_adult <- data.frame(matrix(data=NA, nrow=num.yr, ncol=num.age)) #make a dataframe for the loop to store results in
for (i in 1:num.yr) {
  Pred_Numbers_at_age_JX_adult[i,] <- estimate.JXad$lsmean[i]*prop_JX$Freq
}
colnames(Pred_Numbers_at_age_JX_adult) <- c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8") 

#IR
num.yr = length(estimate.IRad$year) 
num.age =length(prop_IR$Freq)
Pred_Numbers_at_age_IR_adult <- data.frame(matrix(data=NA, nrow=num.yr, ncol=num.age)) #make a dataframe for the loop to store results in
for (i in 1:num.yr) {
  Pred_Numbers_at_age_IR_adult[i,] <- estimate.IRad$lsmean[i]*prop_IR$Freq
}

colnames(Pred_Numbers_at_age_IR_adult) <- c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9")

#Now apply weight (per individual)-at-age schedules to the index of numbers of adult-at-age to get index of SSB-at-age () and then total biomass (sumbiomass)

AP_weight <- read.csv("Weight_at_Age_AP_ad_FIMdata.csv", header=TRUE)
CK_weight <- read.csv("Weight_at_Age_CK_ad_FIMdata.csv", header=TRUE)
TB_weight <- read.csv("Weight_at_Age_TB_ad_FIMdata.csv", header=TRUE)
CH_weight <- read.csv("Weight_at_Age_CH_ad_FIMdata.csv", header=TRUE)
IR_weight <- read.csv("Weight_at_Age_IR_ad_FIMdata.csv", header=TRUE)
JX_weight <- read.csv("Weight_at_Age_JX_ad_FIMdata.csv", header=TRUE)

apmxwt = 2.1812083
AP_weight$mean_wt[8:10] = apmxwt

ckmxwt=1.89
CK_weight$mean_wt[7] =ckmxwt

mxwt=1.45
TB_weight$mean_wt[8:9] =mxwt

mxwt=1.21
CH_weight$mean_wt[7] =mxwt

mxwt=2.54
JX_weight$mean_wt[7:8] =mxwt

mxwt=2.36
IR_weight$mean_wt[9] =mxwt


#AP
num.yr = length(estimate.APad$year) 
num.age =length(prop_AP$Freq)
Pred_Biomass_AP_adult <- data.frame(matrix(data=NA, nrow=num.yr, ncol=num.age)) #make a dataframe for the loop to store results in
for (i in 1:num.yr) {
  Pred_Biomass_AP_adult[i,] <- Pred_Numbers_at_age_AP_adult[i,]*AP_weight$mean_wt
}

Pred_Biomass_AP_adult <- Pred_Biomass_AP_adult %>% mutate(sumbiomass = rowSums(Pred_Biomass_AP_adult))



#CK
num.yr = length(estimate.CKad$year) 
num.age =length(prop_CK$Freq)
Pred_Biomass_CK_adult <- data.frame(matrix(data=NA, nrow=num.yr, ncol=num.age)) #make a dataframe for the loop to store results in
for (i in 1:num.yr) {
  Pred_Biomass_CK_adult[i,] <- Pred_Numbers_at_age_CK_adult[i,]*CK_weight$mean_wt
}

Pred_Biomass_CK_adult <- Pred_Biomass_CK_adult %>% mutate(sumbiomass = rowSums(Pred_Biomass_CK_adult))

#TB
num.yr = length(estimate.TBad$year) 
num.age =length(prop_TB$Freq)
Pred_Biomass_TB_adult <- data.frame(matrix(data=NA, nrow=num.yr, ncol=num.age)) #make a dataframe for the loop to store results in
for (i in 1:num.yr) {
  Pred_Biomass_TB_adult[i,] <- Pred_Numbers_at_age_TB_adult[i,]*TB_weight$mean_wt
}

Pred_Biomass_TB_adult <- Pred_Biomass_TB_adult %>% mutate(sumbiomass = rowSums(Pred_Biomass_TB_adult))

#CH
num.yr = length(estimate.CHad$year) 
num.age =length(prop_CH$Freq)
Pred_Biomass_CH_adult <- data.frame(matrix(data=NA, nrow=num.yr, ncol=num.age)) #make a dataframe for the loop to store results in
for (i in 1:num.yr) {
  Pred_Biomass_CH_adult[i,] <- Pred_Numbers_at_age_CH_adult[i,]*CH_weight$mean_wt
}

Pred_Biomass_CH_adult <- Pred_Biomass_CH_adult %>% mutate(sumbiomass = rowSums(Pred_Biomass_CH_adult))

#JX
num.yr = length(estimate.JXad$year) 
num.age =length(prop_JX$Freq)
Pred_Biomass_JX_adult <- data.frame(matrix(data=NA, nrow=num.yr, ncol=num.age)) #make a dataframe for the loop to store results in
for (i in 1:num.yr) {
  Pred_Biomass_JX_adult[i,] <- Pred_Numbers_at_age_JX_adult[i,]*JX_weight$mean_wt
}

Pred_Biomass_JX_adult <- Pred_Biomass_JX_adult %>% mutate(sumbiomass = rowSums(Pred_Biomass_JX_adult))


#IR
num.yr = length(estimate.IRad$year) 
num.age =length(prop_IR$Freq)
Pred_Biomass_IR_adult <- data.frame(matrix(data=NA, nrow=num.yr, ncol=num.age)) #make a dataframe for the loop to store results in
for (i in 1:num.yr) {
  Pred_Biomass_IR_adult[i,] <- Pred_Numbers_at_age_IR_adult[i,]*IR_weight$mean_wt
}

Pred_Biomass_IR_adult <- Pred_Biomass_IR_adult %>% mutate(sumbiomass = rowSums(Pred_Biomass_IR_adult))

###### FIT SR CURVES TO DETERMINE RESIDUALS #####
#Fit SR Curves to yoy:adult index
# First, make combined indices for YOY (Mean_AP) and Adults (Mean_AP_ad), and Adult biomass (Pred_Biomass_AP_adult)
#I will trim the yoy data because the adult timeseries are shorter- only for CK, TB, CH, and IR

#load you index so I dont have to re run all of the above

APindex <- read.csv(file="U:\\PhD_projectfiles\\Exported_R_Datafiles\\Indices\\UpdatedIndices\\AP_yoy_index.csv", header=T)
CKindex <- read.csv(file="U:\\PhD_projectfiles\\Exported_R_Datafiles\\Indices\\UpdatedIndices\\CK_yoy_index.csv", header=T)
TBindex <- read.csv(file="U:\\PhD_projectfiles\\Exported_R_Datafiles\\Indices\\UpdatedIndices\\TB_yoy_index.csv", header=T)
CHindex <- read.csv(file="U:\\PhD_projectfiles\\Exported_R_Datafiles\\Indices\\UpdatedIndices\\CH_yoy_index.csv", header=T)
JXindex <- read.csv(file="U:\\PhD_projectfiles\\Exported_R_Datafiles\\Indices\\UpdatedIndices\\JX_yoy_index.csv", header=T)
IRindex <- read.csv(file="U:\\PhD_projectfiles\\Exported_R_Datafiles\\Indices\\UpdatedIndices\\IR_yoy_index.csv", header=T)



AP_ind <- data.frame(cbind(APindex$Mean, estimate.APad$lsmean, Pred_Biomass_AP_adult$sumbiomass)) %>% mutate(logyoy=log(X1), logadult=log(X2), logadultbio=log(X3))
names(AP_ind) <- c("yoy", "adult", "adult_bio", "logyoy", "logadult", "logadult_bio")
plot(yoy~adult, data=AP_ind)

CK_ind <- data.frame(cbind(CKindex$Mean[2:20], estimate.CKad$lsmean, Pred_Biomass_CK_adult$sumbiomass)) %>% mutate(logyoy=log(X1), logadult=log(X2), logadultbio=log(X3))
names(CK_ind) <- c("yoy", "adult", "adult_bio", "logyoy", "logadult", "logadult_bio")
plot(yoy~adult, data=CK_ind)

TB_ind <- data.frame(cbind(TBindex$Mean[8:27], estimate.TBad$lsmean, Pred_Biomass_TB_adult$sumbiomass)) %>% mutate(logyoy=log(X1), logadult=log(X2), logadultbio=log(X3))
names(TB_ind) <- c("yoy", "adult", "adult_bio", "logyoy", "logadult", "logadult_bio")
plot(yoy~adult, data=TB_ind)

CH_ind <- data.frame(cbind(CHindex$Mean[8:27], estimate.CHad$lsmean, Pred_Biomass_CH_adult$sumbiomass)) %>% mutate(logyoy=log(X1), logadult=log(X2), logadultbio=log(X3))
names(CH_ind) <- c("yoy", "adult", "adult_bio", "logyoy", "logadult", "logadult_bio")
plot(yoy~adult, data=CH_ind)

JX_ind <- data.frame(cbind(JXindex$Mean, estimate.JXad$lsmean, Pred_Biomass_JX_adult$sumbiomass)) %>% mutate(logyoy=log(X1), logadult=log(X2), logadultbio=log(X3))
names(JX_ind) <- c("yoy", "adult", "adult_bio", "logyoy", "logadult", "logadult_bio")
plot(yoy~adult, data=JX_ind)

IR_ind <- data.frame(cbind(IRindex$Mean[8:26], estimate.IRad$lsmean, Pred_Biomass_IR_adult$sumbiomass)) %>% mutate(logyoy=log(X1), logadult=log(X2), logadultbio=log(X3))
names(IR_ind) <- c("yoy", "adult", "adult_bio", "logyoy", "logadult", "logadult_bio")
plot(yoy~adult, data=IR_ind)

library(FSA)
library(nlstools)
#BEVERTON-HOLT, then density idependent, then immediately followed by a Ricker Model


### AP STOCK RECRUITMENT Fitting #####
####Ricker
srStarts(yoy ~ adult, data=AP_ind, type="Ricker") #determine starting values 
svR_ap <- list(a=3, b=0.5) #putting starting values into a named list for later use
RK <- srFuns("Ricker") #define stock recruit function that ill be using 
srRK_ap <- nls(logyoy~log(RK(adult,a,b)), data=AP_ind, start=svR_ap) #stock recruitment function with multiplicative errors is fit with nls 
overview(srRK_ap) #produces parameter estimates, confidence intervals, Residual sums squares

#visualize the model fit
x=seq(0,5, length.out=999) #many S for predictions
pR <- RK(x, a=coef(srRK_ap)) #predicted mean R for the values of S using the coefficients of the fitted model above (srRK_ap)
xlmts=range(c(x,AP_ind$adult)) #make xlmts on the plot
plot(yoy~adult, data=AP_ind, xlim=xlmts) #plot 
lines(pR~x, lwd=2) #add the line of the predicted mean R

##### Beverton Holt 
srStarts(yoy~adult, data=AP_ind, type="BevertonHolt", param=1)  #determine starting values
svR_ap <- list(a=-52, b=-65)
BH <- srFuns("BevertonHolt")
srBH_ap <- nls(logyoy~log(BH(adult,a,b)), data=AP_ind, start=svR_ap)
overview(srBH_ap)

#visualize the model fit
x=seq(0,5, length.out=999) #many S for predictions
pR <- BH(x, a=coef(srBH_ap)) #predicted mean R for the values of S using the coefficients of the fitted model above (srRK_ap)
xlmts=range(c(x,AP_ind$adult)) #make xlmts on the plot
plot(yoy~adult, data=AP_ind, xlim=xlmts) #plot 
lines(pR~x, lwd=2) #add the line of the predicted mean R

##### Density Independent
ind <- srFuns("independence")
svI <- srStarts(yoy~adult, data=AP_ind, type= "independence")
srI <- nls(logyoy~log(ind(adult,a)), data=AP_ind, start=svI)

#test whether density independent are better than either 
extraSS(srI, com=srRK_ap)
extraSS(srI, com=srBH_ap) #better fit than Bh
              
#test whether models are better
AIC(srBH_ap, srRK_ap)
#says BH is better fit bit the plot looks like absolute crap. the b ricker was not significant

#do not export 
#write.csv(data.frame(residuals(srBH_ap)) %>% mutate(year = c(1998:2015)), "~/Desktop/PhD project/Projects/Seatrout/Data/Exported R Dataframes/ap_resid_BH.csv")
#write.csv(data.frame(residuals(srRK_ap)) %>% mutate(year = c(1998:2015)), "~/Desktop/PhD project/Projects/Seatrout/Data/Exported R Dataframes/ap_resid_RK.csv")

#CK STOCK RECRUITMENT FITTING #####
####Ricker
srStarts(yoy ~ adult, data=CK_ind, type="Ricker") #determine starting values 
svR_ck <- list(a=0.9, b=-0.41) #putting starting values into a named list for later use
RK <- srFuns("Ricker") #define stock recruit function that ill be using 
srRK_ck <- nls(logyoy~log(RK(adult,a,b)), data=CK_ind, start=svR_ck) #stock recruitment function with multiplicative errors is fit with nls 
overview(srRK_ck) #produces parameter estimates, confidence intervals, Residual sums squares
#b parameter is not significant

#visualize the model fit
x=seq(0,5, length.out=999) #many S for predictions
pR <- RK(x, a=coef(srRK_ck)) #predicted mean R for the values of S using the coefficients of the fitted model above (srRK_ap)
xlmts=range(c(x,CK_ind$adult)) #make xlmts on the plot
plot(yoy~adult, data=CK_ind, xlim=xlmts) #plot 
lines(pR~x, lwd=2) #add the line of the predicted mean R

##### Beverton Holt 
srStarts(yoy~adult, data=CK_ind, type="BevertonHolt")  #determine starting values
svR_ck <- list(a=0.5, b=0.73)
BH <- srFuns("BevertonHolt")
srBH_ck <- nls(logyoy~log(BH(adult,a,b)), data=CK_ind, start=svR_ck)
overview(srBH_ck)
#not significant here - wont even fit it  

#visualize the model fit
x=seq(0,5, length.out=999) #many S for predictions
pR <- BH(x, a=coef(srBH_ck)) #predicted mean R for the values of S using the coefficients of the fitted model above (srRK_ap)
xlmts=range(c(x,CK_ind$adult)) #make xlmts on the plot
plot(yoy~adult, data=CK_ind, xlim=xlmts) #plot 
lines(pR~x, lwd=2) #add the line of the predicted mean R
#weird

##### Density Independent
ind <- srFuns("independence")
svI <- srStarts(yoy~adult, data=CK_ind, type= "independence")
srI <- nls(logyoy~log(ind(adult,a)), data=CK_ind, start=svI)

#test whether density independent are better than either 
extraSS(srI, com=srRK_ck) #a density independent is even better 
extraSS(srI, com=srBH_ck)

#test whether models are better
AIC(srBH_ck, srRK_ck)  #both are crap

# #both weren't even significant. both were crap. do not export
# write.csv(data.frame(residuals(srBH_ch)) %>% mutate(year = c(1997:2015)), "~/Desktop/PhD project/Projects/Seatrout/Data/Exported R Dataframes/ap_resid_BH.csv")
# write.csv(data.frame(residuals(srRK_ch)) %>% mutate(year = c(1997:2015)), "~/Desktop/PhD project/Projects/Seatrout/Data/Exported R Dataframes/ap_resid_RK.csv")

# TB STOCK RECRUITMENT FITTING #####
####Ricker
srStarts(yoy ~ adult, data=TB_ind, type="Ricker") #determine starting values 
svR_tb <- list(a=7.7, b=2.25) #putting starting values into a named list for later use
RK <- srFuns("Ricker") #define stock recruit function that ill be using 
srRK_tb <- nls(logyoy~log(RK(adult,a,b)), data=TB_ind, start=svR_tb) #stock recruitment function with multiplicative errors is fit with nls 
overview(srRK_tb) #produces parameter estimates, confidence intervals, Residual sums squares

#visualize the model fit
x=seq(0,max(TB_ind$adult), length.out=999) #many S for predictions
pR <- RK(x, a=coef(srRK_tb)) #predicted mean R for the values of S using the coefficients of the fitted model above (srRK_ap)
xlmts=range(c(x,TB_ind$adult)) #make xlmts on the plot
plot(yoy~adult, data=TB_ind, xlim=xlmts) #plot 
lines(pR~x, lwd=2) #add the line of the predicted mean R

##### Beverton Holt 
srStarts(yoy~adult, data=TB_ind, type="BevertonHolt")  #determine starting values
svR_tb <- list(a=-52, b=-39)
BH <- srFuns("BevertonHolt")
srBH_tb <- nls(logyoy~log(BH(adult,a,b)), data=TB_ind, start=svR_tb)
overview(srBH_tb)
#not significant here 

#visualize the model fit
x=seq(0,max(TB_ind$adult), length.out=999) #many S for predictions
pR <- BH(x, a=coef(srBH_tb)) #predicted mean R for the values of S using the coefficients of the fitted model above (srRK_ap)
xlmts=range(c(x,TB_ind$adult)) #make xlmts on the plot
plot(yoy~adult, data=TB_ind, xlim=xlmts) #plot 
lines(pR~x, lwd=2) #add the line of the predicted mean R


##### Density Independent
ind <- srFuns("independence")
svI <- srStarts(yoy~adult, data=TB_ind, type= "independence")
srI <- nls(logyoy~log(ind(adult,a)), data=TB_ind, start=svI)

#test whether density independent are better than either 
extraSS(srI, com=srRK_tb) #ricker is better
extraSS(srI, com=srBH_tb) #bh is better but its not significant 

#test whether models are better
AIC(srBH_tb, srRK_tb)  #Ricker is better. But still isnt good. 

#not exporting. ricker is partly significant. BH is shit
#write.csv(data.frame(residuals(srBH_tb)) %>% mutate(year = c(1996:2015)), "~/Desktop/PhD project/Projects/Seatrout/Data/Exported R Dataframes/ap_resid_BH.csv")
#write.csv(data.frame(residuals(srRK_tb)) %>% mutate(year = c(1996:2015)), "~/Desktop/PhD project/Projects/Seatrout/Data/Exported R Dataframes/ap_resid_RK.csv")

#CH STOCK RECRUITMENT FITTING #####
####Ricker
srStarts(yoy ~ adult, data=CH_ind, type="Ricker") #determine starting values 
svR_ch <- list(a=15, b=5) #putting starting values into a named list for later use
RK <- srFuns("Ricker") #define stock recruit function that ill be using 
srRK_ch <- nls(logyoy~log(RK(adult,a,b)), data=CH_ind, start=svR_ch) #stock recruitment function with multiplicative errors is fit with nls 
overview(srRK_ch) #produces parameter estimates, confidence intervals, Residual sums squares

#visualize the model fit
x=seq(0,max(CH_ind$adult), length.out=999) #many S for predictions
pR <- RK(x, a=coef(srRK_ch)) #predicted mean R for the values of S using the coefficients of the fitted model above (srRK_ch)
xlmts=range(c(x,CH_ind$adult)) #make xlmts on the plot
plot(yoy~adult, data=CH_ind, xlim=xlmts) #plot 
lines(pR~x, lwd=2) #add the line of the predicted mean R

##### Beverton Holt 
srStarts(yoy~adult, data=CH_ind, type="BevertonHolt", param=1)  #determine starting values
svR_ch <- list(a=-13, b=-19)
BH <- srFuns("BevertonHolt")
srBH_ch <- nls(logyoy~log(BH(adult,a,b)), data=CH_ind, start=svR_ch)
overview(srBH_ch)
#not significant here 

#visualize the model fit
x=seq(0,max(CH_ind$adult), length.out=999) #many S for predictions
pR <- BH(x, a=coef(srBH_ch)) #predicted mean R for the values of S using the coefficients of the fitted model above (srRK_ap)
xlmts=range(c(x,CH_ind$adult)) #make xlmts on the plot
plot(yoy~adult, data=CH_ind, xlim=xlmts) #plot 
lines(pR~x, lwd=2) #add the line of the predicted mean R
#weird

##### Density Independent
ind <- srFuns("independence")
svI <- srStarts(yoy~adult, data=CH_ind, type= "independence")
srI <- nls(logyoy~log(ind(adult,a)), data=CH_ind, start=svI)

#test whether density independent are better than either 
extraSS(srI, com=srRK_ch)
extraSS(srI, com=srBH_ch)

#test whether models are better
AIC(srBH_ch, srRK_ch)  #dont believe it

#both are shit. Ricker is a bit better. 
#write.csv(data.frame(residuals(srBH_ch)) %>% mutate(year = c(1996:2015)), "~/Desktop/PhD project/Projects/Seatrout/Data/Exported R Dataframes/ap_resid_BH.csv")
#write.csv(data.frame(residuals(srRK_ch)) %>% mutate(year = c(1996:2015)), "~/Desktop/PhD project/Projects/Seatrout/Data/Exported R Dataframes/ap_resid_RK.csv")

#JX STOCK RECRUITMENT CURVE FITTING ####
####Ricker
srStarts(yoy ~ adult, data=JX_ind, type="Ricker") #determine starting values 
svR_jx <- list(a=1.2, b=103) #putting starting values into a named list for later use
RK <- srFuns("Ricker") #define stock recruit function that ill be using 
srRK_jx <- nls(logyoy~log(RK(adult,a,b)), data=JX_ind, start=svR_jx) #stock recruitment function with multiplicative errors is fit with nls 
overview(srRK_jx) #produces parameter estimates, confidence intervals, Residual sums squares
#not significant here 

#visualize the model fit
x=seq(0,max(JX_ind$adult), length.out=999) #many S for predictions
pR <- RK(x, a=coef(srRK_jx)) #predicted mean R for the values of S using the coefficients of the fitted model above (srRK_ap)
xlmts=range(c(x,JX_ind$adult)) #make xlmts on the plot
plot(yoy~adult, data=JX_ind, xlim=xlmts) #plot 
lines(pR~x, lwd=2) #add the line of the predicted mean R

##### Beverton Holt 
srStarts(yoy~adult,data=JX_ind, type="BevertonHolt")  #determine starting values
svR_jx <- list(a=0.63, b=40)
BH <- srFuns("BevertonHolt")
srBH_jx <- nls(logyoy~log(BH(adult,a,b)), data=JX_ind, start=svR_jx)
overview(srBH_jx)
#not significant here 

#visualize the model fit
x=seq(0,5, length.out=999) #many S for predictions
pR <- BH(x, a=coef(srBH_jx)) #predicted mean R for the values of S using the coefficients of the fitted model above (srRK_ap)
xlmts=range(c(x,JX_ind$adult)) #make xlmts on the plot
plot(yoy~adult, data=JX_ind, xlim=xlmts) #plot 
lines(pR~x, lwd=2) #add the line of the predicted mean R

##### Density Independent
ind <- srFuns("independence")
svI <- srStarts(yoy~adult, data=JX_ind, type= "independence")
srI <- nls(logyoy~log(ind(adult,a)), data=JX_ind, start=svI)

#test whether density independent are better than either 
extraSS(srI, com=srRK_jx) #independent is better
extraSS(srI, com=srBH_jx) #independent is better

#test whether models are better
AIC(srBH_jx, srRK_jx)  #Ricker is better. But still isnt good. 

#both are shit
#write.csv(data.frame(residuals(srBH_jx)) %>% mutate(year = c(2001:2015)), "~/Desktop/PhD project/Projects/Seatrout/Data/Exported R Dataframes/ap_resid_BH.csv")
#write.csv(data.frame(residuals(srRK_jx)) %>% mutate(year = c(2001:2015)), "~/Desktop/PhD project/Projects/Seatrout/Data/Exported R Dataframes/ap_resid_RK.csv")

# IR STOCK RECRUITMENT CURVE FITTING #####
####Ricker
srStarts(yoy ~ adult, data=IR_ind, type="Ricker") #determine starting values 
svR_ir <- list(a=2.8, b=1.01) #putting starting values into a named list for later use
RK <- srFuns("Ricker") #define stock recruit function that ill be using 
srRK_ir <- nls(logyoy~log(RK(adult,a,b)), data=IR_ind, start=svR_ir) #stock recruitment function with multiplicative errors is fit with nls 
overview(srRK_ir) #produces parameter estimates, confidence intervals, Residual sums squares

#visualize the model fit
x=seq(0,5, length.out=999) #many S for predictions
pR <- RK(x, a=coef(srRK_ir)) #predicted mean R for the values of S using the coefficients of the fitted model above (srRK_ap)
xlmts=range(c(x,IR_ind$adult)) #make xlmts on the plot
plot(yoy~adult, data=IR_ind, xlim=xlmts) #plot 
lines(pR~x, lwd=2) #add the line of the predicted mean R

##### Beverton Holt 
srStarts(yoy~adult, data=IR_ind, type="BevertonHolt")  #determine starting values
svR_ir <- list(a=4.9, b=4.3)
BH <- srFuns("BevertonHolt")
srBH_ir <- nls(logyoy~log(BH(adult,a,b)), data=IR_ind, start=svR_ir)
overview(srBH_ir)
#not significant here 

#visualize the model fit
x=seq(0,5, length.out=999) #many S for predictions
pR <- BH(x, a=coef(srBH_ir)) #predicted mean R for the values of S using the coefficients of the fitted model above (srRK_ap)
xlmts=range(c(x,IR_ind$adult)) #make xlmts on the plot
plot(yoy~adult, data=IR_ind, xlim=xlmts) #plot 
lines(pR~x, lwd=2) #add the line of the predicted mean R

##### Density Independent
ind <- srFuns("independence")
svI <- srStarts(yoy~adult, data=IR_ind, type= "independence")
srI <- nls(logyoy~log(ind(adult,a)), data=IR_ind, start=svI)

#test whether density independent are better than either 
extraSS(srI, com=srRK_ir) #ricker is better
extraSS(srI, com=srBH_ir) #bh is better but its not significant 

#test whether models are better
AIC(srBH_tb, srRK_tb)  #Ricker is better. But still isnt good. 

#Ricker is ok but continuity not going to export because I didn't export the others.... 
#write.csv(data.frame(residuals(srBH_ir)) %>% mutate(year = c(1997:2015)), "~/Desktop/PhD project/Projects/Seatrout/Data/Exported R Dataframes/ap_resid_BH.csv")
#write.csv(data.frame(residuals(srRK_ir)) %>% mutate(year = c(1997:2015)), "~/Desktop/PhD project/Projects/Seatrout/Data/Exported R Dataframes/ap_resid_RK.csv")


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



# Produce Adjusted Indices for Adult


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


# Make combined biomass indices for YOY and Adults 

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

#plot bh and density independent
plot(yoy~adult,data=AP_ind)
curve((coef(srBH_ap)[1]*x)/(1+coef(srBH_ap)[2]*x),from=0,to=120,col="red",lwd=2,add=TRUE)
curve(coef(bh0nls)[1]*x,from=0,to=120,col="blue",lwd=2,add=TRUE)
legend("topleft",legend=c("density independent","density dependent"),col=c("blue","red"),lwd=2,cex=0.6)

#plot bh
x=seq(0,5, length.out=999)
pR <- bh(x, a=coef(srBH_ap))
xlmts=range(c(x,AP_bio$adult))
plot(yoy~adult, data=AP_bio, xlim=xlmts)
lines(pR~x, lwd=2)

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

##### ERROR PROPAGATION TO FIND FINAL VALUE (pos * prop.pos)_YOY 
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
  x = c(LSM_ap.pos$rate[i], LSM_ap.pos$SE[i])
  y= c(LSM_ap.bin$prob[i], LSM_ap.bin$SE[i])
  EXPR <- expression(x*y)
  DF <- cbind(x,y)
  RES <- propagate(expr=EXPR, data=DF, type='stat', do.sim=TRUE, verbose=TRUE)
  df[i,] <- t(matrix(RES$sim))
}
Mean_AP <- df %>% cbind(LSM_ap.pos$year)
colnames(Mean_AP) <- c("Mean", "SD", "Median", "MAD", "2.5%", "97.5%", "Year")

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

# OLD _ CREATE INDEX WITH MONTE CARLO_ADULT #####
#DO NOT DO THIS. THe standard errors of the lsmeans are fucked (can't figure out how to get them to a reasonable range)
# so I'm just going to use the lsmeans estimates and sort of ignore the standard error for the adult
#I only need the adult series to fit the curves so I'm giving up on figuring out the problems


# creating the index through Monte Carlo simulations
# Technically don't have to if you just want the estimates and the standard errors. 
#Monte Carlo simulations are necessary here because its the only way to incorporate the lsmeans error into the estimate and 
#then determine CV, and 95% conf intervals. Without a distribution you can't get to the confidence intervals which are what people actually care about. 
# This is borrowed from tuning_index_update.R from Hogfish2017 at FWRI. This monte carlo method seems standard among the group. 

# # AP Monte Carlo 
# #making an empty matrix to fill with results 
# num.yr <- length(estimate.APad$year)
# index.dist <- matrix(data=NA,nrow=num.yr,ncol=8) 
# num.iter=10000
# 
# #build the random distribution (i.e the random deviates)
# for (i in 1:num.yr) {
#   rand.1    <- qt(runif(num.iter,0.001,0.999),sample.sizeAPad$num.trips) 
#   
#   #create the distribution of the data 
#   APdist <- estimate.APad$lsmean[i] + estimate.APad$SE[i] * rand.1 
#   
#   index.dist[i,1] <- mean(APdist) #take the mean of the dist (aka temp)
#   index.dist[i,2] <- sd(APdist) #take the sd of the dist
#   index.dist[i,3] <- sd(APdist)/mean(APdist) #create CV
#   index.dist[i,4:8] <- quantile(APdist,probs=c(0.025,0.25,0.50,0.75,0.975))
# }
# Upper        <- index.dist[ ,8] - index.dist[ ,7]
# Lower        <- index.dist[ ,5] - index.dist[ ,4]
# NominalMean  = as.vector(tapply(ap.ad$number, ap.ad$year, mean))
# NominalSD    = as.vector(tapply(ap.ad$number, ap.ad$year, sd))
# NominalCV    = NominalSD /NominalMean
# APadindex <- as.data.frame(cbind(sample.sizeAPad, index.dist, Upper, Lower, NominalMean, NominalSD, NominalCV ))
# names(APadindex) <- c("year","Total.num.trips","Mean","std.dev", "CV","Low.95","Qtr.1","Median","Qtr.3","Up.95","Upper",
#                     "Lower", "NominalMean", "NominalSD", "NominalCV" )
# APadindex$year = as.numeric(APadindex$year)
# 
# # CK Monte Carlo
# #making an empty matrix to fill with results 
# num.yr <- length(estimate.CKad$year)
# index.dist <- matrix(data=NA,nrow=num.yr,ncol=8) 
# num.iter=10000
# 
# #build the random distribution (i.e the random deviates)
# for (i in 1:num.yr) {
#   rand.1    <- qt(runif(num.iter,0.001,0.999),sample.sizeCKad$num.trips) 
#   
#   #create the distribution of the data 
#   CKdist <- estimate.CKad$lsmean[i] + estimate.CKad$SE[i] * rand.1 
#   
#   index.dist[i,1] <- mean(CKdist) #take the mean of the dist (aka temp)
#   index.dist[i,2] <- sd(CKdist) #take the sd of the dist
#   index.dist[i,3] <- sd(CKdist)/mean(CKdist) #create CV
#   index.dist[i,4:8] <- quantile(CKdist,probs=c(0.025,0.25,0.50,0.75,0.975))
# }
# Upper        <- index.dist[ ,8] - index.dist[ ,7]
# Lower        <- index.dist[ ,5] - index.dist[ ,4]
# NominalMean  = as.vector(tapply(ck.ad$number, ck.ad$year, mean))
# NominalSD    = as.vector(tapply(ck.ad$number, ck.ad$year, sd))
# NominalCV    = NominalSD /NominalMean
# CKadindex <- as.data.frame(cbind(sample.sizeCKad, index.dist, Upper, Lower, NominalMean, NominalSD, NominalCV ))
# names(CKadindex) <- c("year","Total.num.trips","Mean","std.dev", "CV","Low.95","Qtr.1","Median","Qtr.3","Up.95","Upper",
#                     "Lower", "NominalMean", "NominalSD", "NominalCV" )
# CKadindex$year = as.numeric(CKadindex$year)
# 
# # TB Monte Carlo
# num.yr <- length(estimate.TBad$year)
# index.dist <- matrix(data=NA,nrow=num.yr,ncol=8) 
# num.iter=10000
# 
# #build the random distribution (i.e the random deviates)
# for (i in 1:num.yr) {
#   rand.1    <- qt(runif(num.iter,0.001,0.999),sample.sizeTBad$num.trips) 
#   
#   #create the distribution of the data 
#   TBdist <- estimate.TBad$lsmean[i] + estimate.TBad$SE[i] * rand.1 
#   
#   index.dist[i,1] <- mean(TBdist) #take the mean of the dist (aka temp)
#   index.dist[i,2] <- sd(TBdist) #take the sd of the dist
#   index.dist[i,3] <- sd(TBdist)/mean(TBdist) #create CV
#   index.dist[i,4:8] <- quantile(TBdist,probs=c(0.025,0.25,0.50,0.75,0.975))
# }
# Upper        <- index.dist[ ,8] - index.dist[ ,7]
# Lower        <- index.dist[ ,5] - index.dist[ ,4]
# NominalMean  = as.vector(tapply(tb.ad$number, tb.ad$year, mean))
# NominalSD    = as.vector(tapply(tb.ad$number, tb.ad$year, sd))
# NominalCV    = NominalSD /NominalMean
# TBadindex <- as.data.frame(cbind(sample.sizeTBad, index.dist, Upper, Lower, NominalMean, NominalSD, NominalCV ))
# names(TBadindex) <- c("year","Total.num.trips","Mean","std.dev", "CV","Low.95","Qtr.1","Median","Qtr.3","Up.95","Upper",
#                     "Lower", "NominalMean", "NominalSD", "NominalCV" )
# TBadindex$year = as.numeric(TBadindex$year)
# 
# # CH Monte Carlo
# num.yr <- length(estimate.CHad$year)
# index.dist <- matrix(data=NA,nrow=num.yr,ncol=8) 
# num.iter=10000
# 
# #build the random distribution (i.e the random deviates)
# for (i in 1:num.yr) {
#   rand.1    <- qt(runif(num.iter,0.001,0.999),sample.sizeCHad$num.trips) 
#   
#   #create the distribution of the data 
#   CHdist <- estimate.CHad$lsmean[i] + estimate.CHad$SE[i] * rand.1 
#   
#   index.dist[i,1] <- mean(CHdist) #take the mean of the dist (aka temp)
#   index.dist[i,2] <- sd(CHdist) #take the sd of the dist
#   index.dist[i,3] <- sd(CHdist)/mean(CHdist) #create CV
#   index.dist[i,4:8] <- quantile(CHdist,probs=c(0.025,0.25,0.50,0.75,0.975))
# }
# Upper        <- index.dist[ ,8] - index.dist[ ,7]
# Lower        <- index.dist[ ,5] - index.dist[ ,4]
# NominalMean  = as.vector(tapply(ch.ad$number, ch.ad$year, mean))
# NominalSD    = as.vector(tapply(ch.ad$number, ch.ad$year, sd))
# NominalCV    = NominalSD /NominalMean
# CHadindex <- as.data.frame(cbind(sample.sizeCHad, index.dist, Upper, Lower, NominalMean, NominalSD, NominalCV ))
# names(CHadindex) <- c("year","Total.num.trips","Mean","std.dev", "CV","Low.95","Qtr.1","Median","Qtr.3","Up.95","Upper",
#                     "Lower", "NominalMean", "NominalSD", "NominalCV" )
# CHadindex$year = as.numeric(CHadindex$year)
# 
# # JX Monte Carlo 
# num.yr <- length(estimate.JXad$year)
# index.dist <- matrix(data=NA,nrow=num.yr,ncol=8) 
# num.iter=10000
# 
# #build the random distribution (i.e the random deviates)
# for (i in 1:num.yr) {
#   rand.1    <- qt(runif(num.iter,0.001,0.999),sample.sizeJXad$num.trips) 
#   
#   #create the distribution of the data 
#   JXdist <- estimate.JXad$lsmean[i] + estimate.JXad$SE[i] * rand.1 
#   
#   index.dist[i,1] <- mean(JXdist) #take the mean of the dist (aka temp)
#   index.dist[i,2] <- sd(JXdist) #take the sd of the dist
#   index.dist[i,3] <- sd(JXdist)/mean(JXdist) #create CV
#   index.dist[i,4:8] <- quantile(JXdist,probs=c(0.025,0.25,0.50,0.75,0.975))
# }
# Upper        <- index.dist[ ,8] - index.dist[ ,7]
# Lower        <- index.dist[ ,5] - index.dist[ ,4]
# NominalMean  = as.vector(tapply(jx.fl$number, jx.fl$year, mean))
# NominalSD    = as.vector(tapply(jx.fl$number, jx.fl$year, sd))
# NominalCV    = NominalSD /NominalMean
# JXadindex <- as.data.frame(cbind(sample.sizeJXad, index.dist, Upper, Lower, NominalMean, NominalSD, NominalCV ))
# names(JXadindex) <- c("year","Total.num.trips","Mean","std.dev", "CV","Low.95","Qtr.1","Median","Qtr.3","Up.95","Upper",
#                     "Lower", "NominalMean", "NominalSD", "NominalCV" )
# JXadindex$year = as.numeric(JXadindex$year)
# 
# # IR Monte Carlo
# 
# num.yr <- length(estimate.IRad$year)
# index.dist <- matrix(data=NA,nrow=num.yr,ncol=8) 
# num.iter=10000
# 
# #build the random distribution (i.e the random deviates)
# for (i in 1:num.yr) {
#   rand.1    <- qt(runif(num.iter,0.001,0.999),sample.sizeIRad$num.trips) 
#   
#   #create the distribution of the data 
#   IRdist <- estimate.IRad$lsmean[i] + estimate.IRad$SE[i] * rand.1 
#   
#   index.dist[i,1] <- mean(IRdist) #take the mean of the dist (aka temp)
#   index.dist[i,2] <- sd(IRdist) #take the sd of the dist
#   index.dist[i,3] <- sd(IRdist)/mean(IRdist) #create CV
#   index.dist[i,4:8] <- quantile(IRdist,probs=c(0.025,0.25,0.50,0.75,0.975))
# }
# Upper        <- index.dist[ ,8] - index.dist[ ,7]
# Lower        <- index.dist[ ,5] - index.dist[ ,4]
# NominalMean  = as.vector(tapply(ir.ad$number, ir.ad$year, mean))
# NominalSD    = as.vector(tapply(ir.ad$number, ir.ad$year, sd))
# NominalCV    = NominalSD /NominalMean
# IRadindex <- as.data.frame(cbind(sample.sizeIRad, index.dist, Upper, Lower, NominalMean, NominalSD, NominalCV ))
# names(IRadindex) <- c("year","Total.num.trips","Mean","std.dev", "CV","Low.95","Qtr.1","Median","Qtr.3","Up.95","Upper",
#                     "Lower", "NominalMean", "NominalSD", "NominalCV" )
# IRadindex$year = as.numeric(IRadindex$year)

