# 10/17/2016 This script imports recreational (MRFSS, MRIP) length data 
# Main Objectives of this script: 
# 1. Imports mrfss and mrip length data.
# 4. One way ANOVA to determine whether there is a significant difference in total length among estuaries
# 2. Plots of ECDF
# 3. Statistical comparison of bay-bay length frequencies. 
#     - K-S test
#     - Chi-Square Test
###############################################################

# ******Note about sample sizes affecting results of hypothesis testing:
# http://stats.stackexchange.com/questions/47498/mann-whitney-u-test-and-k-s-test-with-unequal-sample-sizes
# "With such large sample sizes both tests will have high power to detect minor differences. 
#    The 2 distributions could be almost identical with a small difference in shape location that is not of practical importance 
#    and the tests would reject (because they are different). If all you really care about is a statistically significant difference 
#    then you can be happy with the results of the KS test (and others, even a t-test will be meaningful with non-normal data of 
#    those sample sizes due to the Central Limit Theorem)."


library(FSA)
library(magrittr)
library(dplyr)
library(plotrix)
library(haven)
library(Matching)
library(dplyr)
library(ggplot2)
library(scales)

# set working directory
setwd("~/Desktop/PhD project/Projects/Seatrout/Data")
setwd("U:/PhD_projectfiles/Raw_Data/Age_Length_Data")

# **************Length frequency distributions of MRFSS and MRIP

#import the sas data **************************
mrfss <- read_sas("mrfss_lens_8115.sas7bdat")
mrip <- read_sas("mrip_lens_20042015.sas7bdat")

# subset each estuary using the Country of Encounter codes
# From the definitions defined by FWRI
#FIPS county codes
# County_Bay_assignments.xlsx in GitHub>Seatrout>Data

#load data
#select the counties that most closely align with the bays in question
# convert fork length to tl using the equation TLmm = 1.00467 * FL + 0.04850 from the SAS mrfss leng_freq
# chose data from 2000 on because of major regulatory changes 

mrfss_AP <- subset(read_sas("mrfss_lens_8115.sas7bdat"), CNTY %in% c(5,33,37,45,65,77,91,113,123,129,131,133) & YEAR>=2000 & YEAR <= 2003) %>% mutate(tl = (1.00467*flmm+0.04850)/10) %>% filter(!is.na(tl)) %>% mutate(bay=rep("AP",2234)) 
  mrfss_AP <- droplevels(subset(mrfss_AP, select=c(tl, bay)))
mrfss_TB <- subset(read_sas("mrfss_lens_8115.sas7bdat"), CNTY %in% c(53,101,57,103,81) & YEAR>=2000 & YEAR<= 2003) %>% mutate(tl = (1.00467*flmm+0.04850)/10) %>% filter(!is.na(tl)) %>% mutate(bay=rep("TB",1741))
  mrfss_TB <- droplevels(subset(mrfss_TB, select=c(tl, bay)))
mrfss_CH <- subset(read_sas("mrfss_lens_8115.sas7bdat"), CNTY %in% c(27,115,7,15,71,21,87,51) & YEAR>=2000 & YEAR<= 2003) %>% mutate(tl = (1.00467*flmm+0.04850)/10) %>% filter(!is.na(tl))%>% mutate(bay=rep("CH",1085))
  mrfss_CH <- droplevels(subset(mrfss_CH, select=c(tl, bay)))
mrfss_CK <- subset(read_sas("mrfss_lens_8115.sas7bdat"), CNTY %in% c(75,29,17,1) & YEAR>=2000 & YEAR<= 2003) %>% mutate(tl = (1.00467*flmm+0.04850)/10) %>% filter(!is.na(tl))%>% mutate(bay=rep("CK",1225))
  mrfss_CK <- droplevels(subset(mrfss_CK, select=c(tl, bay)))
mrfss_JX <- subset(read_sas("mrfss_lens_8115.sas7bdat"), CNTY %in% c(19,31,35,89,107,109) & YEAR>=2000 & YEAR<= 2003) %>% mutate(tl = (1.00467*flmm+0.04850)/10) %>% filter(!is.na(tl))%>% mutate(bay=rep("JX",291))
  mrfss_JX <- droplevels(subset(mrfss_JX, select=c(tl, bay)))
mrfss_IR <- subset(read_sas("mrfss_lens_8115.sas7bdat"), CNTY %in% c(9,11,25,61,85,99,111,127) & YEAR>=2000 & YEAR <=2003) %>% mutate(tl = (1.00467*flmm+0.04850)/10) %>% filter(!is.na(tl))%>% mutate(bay=rep("IR",812))
  mrfss_IR <- droplevels(subset(mrfss_IR, select=c(tl, bay)))
  
mrip_AP <- subset(read_sas("mrip_lens_20042015.sas7bdat"), CNTY %in% c(5,33,37,45,65,77,91,113,123,129,131,133)) %>% mutate(tl = (1.00467*LNGTH+0.04850)/10) %>% filter(!is.na(tl)) %>% mutate(bay=rep("AP",9265))
  mrip_AP <- droplevels(subset(mrip_AP, select=c(tl, bay)))
mrip_TB <- subset(read_sas("mrip_lens_20042015.sas7bdat"), CNTY %in% c(53,101,57,103,81)) %>% mutate(tl = (1.00467*LNGTH+0.04850)/10) %>% filter(!is.na(tl)) %>% mutate(bay=rep("TB",8220))
  mrip_TB <- droplevels(subset(mrip_TB, select=c(tl, bay)))
mrip_CH <- subset(read_sas("mrip_lens_20042015.sas7bdat"), CNTY %in% c(27,115,7,15,71,21,87,51)) %>% mutate(tl = (1.00467*LNGTH+0.04850)/10) %>% filter(!is.na(tl))%>% mutate(bay=rep("CH",4832))
  mrip_CH <- droplevels(subset(mrip_CH, select=c(tl, bay)))
mrip_CK <- subset(read_sas("mrip_lens_20042015.sas7bdat"), CNTY %in% c(75,29,17,1)) %>% mutate(tl = (1.00467*LNGTH+0.04850)/10) %>% filter(!is.na(tl))%>% mutate(bay=rep("CK",5613))
  mrip_CK <- droplevels(subset(mrip_CK, select=c(tl, bay)))
mrip_JX <- subset(read_sas("mrip_lens_20042015.sas7bdat"), CNTY %in% c(19,31,35,89,107,109)) %>% mutate(tl = (1.00467*LNGTH+0.04850)/10) %>% filter(!is.na(tl))%>% mutate(bay=rep("JX",1510))
  mrip_JX <- droplevels(subset(mrip_JX, select=c(tl, bay)))
mrip_IR <- subset(read_sas("mrip_lens_20042015.sas7bdat"), CNTY %in% c(9,11,25,61,85,99,111,127)) %>% mutate(tl = (1.00467*LNGTH+0.04850)/10) %>% filter(!is.na(tl))%>% mutate(bay=rep("IR",2595))
  mrip_IR <- droplevels(subset(mrip_IR, select=c(tl, bay)))

# COMBINE MRFSS And MRIP DATA

AP <- rbind(mrfss_AP, mrip_AP)
TB <- rbind(mrfss_TB, mrip_TB)
CK <- rbind(mrfss_CK, mrip_CK)
CH <- rbind(mrfss_CH, mrip_CH)
JX <- rbind(mrfss_JX, mrip_JX)
IR <- rbind(mrfss_IR, mrip_IR)

All= rbind(AP, TB, CK, CH, JX, IR)
All$bay <- as.factor(All$bay)
All <- na.omit(All)

All_sum <- summarise(group_by(All, bay), N = length(tl), mean_tl = mean(tl), sd_tl= sd(tl), se_tl= sd_tl/(sqrt(length(tl))))


N = nrow(AP) +nrow(TB) +nrow(CK) +nrow(CH) +nrow(JX) +nrow(IR)


######################
#PLOT ECDF (empirical cumulative distribution function) - used for K-S test
######################

#MRFSS
clr<- c("black", "gray50", "red", "blue", "green", "pink")
plot(ecdf(mrfss_TB$tl), xlab="Total Length (cm)", do.points=FALSE, verticals=TRUE, main="", col.01line=NULL)
plot(ecdf(mrfss_AP$tl), add=TRUE, do.points=FALSE, verticals=TRUE, col=clr[2], col.01line=NULL)
plot(ecdf(mrfss_CK$tl), add=TRUE, do.points=FALSE, verticals=TRUE, col=clr[3], col.01line=NULL)
plot(ecdf(mrfss_CH$tl), add=TRUE, do.points=FALSE, verticals=TRUE, col=clr[4], col.01line=NULL)
plot(ecdf(mrfss_IR$tl), add=TRUE, do.points=FALSE, verticals=TRUE, col=clr[5], col.01line=NULL)
plot(ecdf(mrfss_JX$tl), add=TRUE, do.points=FALSE, verticals=TRUE, col=clr[6], col.01line=NULL)

#MRIP
clr<- c("black", "gray50", "red", "blue", "green", "pink")
plot(ecdf(mrip_TB$tl), xlab="Total Length (cm)", do.points=FALSE, verticals=TRUE, main="", col.01line=NULL)
plot(ecdf(mrip_AP$tl), add=TRUE, do.points=FALSE, verticals=TRUE, col=clr[2], col.01line=NULL)
plot(ecdf(mrip_CK$tl), add=TRUE, do.points=FALSE, verticals=TRUE, col=clr[3], col.01line=NULL)
plot(ecdf(mrip_CH$tl), add=TRUE, do.points=FALSE, verticals=TRUE, col=clr[4], col.01line=NULL)
plot(ecdf(mrip_IR$tl), add=TRUE, do.points=FALSE, verticals=TRUE, col=clr[5], col.01line=NULL)
plot(ecdf(mrip_JX$tl), add=TRUE, do.points=FALSE, verticals=TRUE, col=clr[6], col.01line=NULL)

#ALL
clr<- c("black", "gray50", "red", "blue", "green", "pink")
plot(ecdf(TB$tl), xlab="Total Length (cm)", do.points=FALSE, verticals=TRUE, main="", col.01line=NULL)
plot(ecdf(AP$tl), add=TRUE, do.points=FALSE, verticals=TRUE, col=clr[2], col.01line=NULL)
plot(ecdf(CK$tl), add=TRUE, do.points=FALSE, verticals=TRUE, col=clr[3], col.01line=NULL)
plot(ecdf(CH$tl), add=TRUE, do.points=FALSE, verticals=TRUE, col=clr[4], col.01line=NULL)
plot(ecdf(IR$tl), add=TRUE, do.points=FALSE, verticals=TRUE, col=clr[5], col.01line=NULL)
plot(ecdf(JX$tl), add=TRUE, do.points=FALSE, verticals=TRUE, col=clr[6], col.01line=NULL)

test <- ecdf(TB$tl)
#########################################                                                                  
#AMONG GROUP STATISTICAL COMPARISONS
#########################################

####################
#First with K-S Test
####################
# used to determine whether ECDF are the same between two groups and can detect differences in location, dispersion and shape of the ditrubtions
#bootstrapped version of the K-S test that is insensitive to ties with noncontinuous data is implemented in ks.boot
ks<- c(ks.boot(mrfss_TB$tl, mrfss_CK$tl, nboots=1000)$p.value, ks.boot(mrfss_TB$tl, mrfss_CH$tl, nboots=5000)$p.value)
c <- ks.boot(mrfss_TB$tl, mrfss_TB$tl)
#if comparing more than two groups then each pair must be compared separately and the p value from each comparison must
# be adjusted for an increasing experimentwise error rate due to multiple comparisons. 

#MRFSS
ks_ALL_mrfss <-c(ks.boot(mrfss_TB$tl, mrfss_CK$tl, nboots=1000)$ks.boot.pvalue, 
  ks.boot(mrfss_TB$tl, mrfss_CH$tl, nboots=1000)$ks.boot.pvalue,
  ks.boot(mrfss_TB$tl, mrfss_AP$tl, nboots=1000)$ks.boot.pvalue,
  ks.boot(mrfss_TB$tl, mrfss_JX$tl, nboots=1000)$ks.boot.pvalue,
  ks.boot(mrfss_TB$tl, mrfss_IR$tl, nboots=1000)$ks.boot.pvalue,
  ks.boot(mrfss_CK$tl, mrfss_CH$tl, nboots=1000)$ks.boot.pvalue,
  ks.boot(mrfss_CK$tl, mrfss_AP$tl, nboots=1000)$ks.boot.pvalue,
  ks.boot(mrfss_CK$tl, mrfss_JX$tl, nboots=1000)$ks.boot.pvalue,
  ks.boot(mrfss_CK$tl, mrfss_IR$tl, nboots=1000)$ks.boot.pvalue,
  ks.boot(mrfss_CH$tl, mrfss_AP$tl, nboots=1000)$ks.boot.pvalue,
  ks.boot(mrfss_CH$tl, mrfss_JX$tl, nboots=1000)$ks.boot.pvalue,
  ks.boot(mrfss_CH$tl, mrfss_IR$tl, nboots=1000)$ks.boot.pvalue,
  ks.boot(mrfss_AP$tl, mrfss_JX$tl, nboots=1000)$ks.boot.pvalue,
  ks.boot(mrfss_AP$tl, mrfss_IR$tl, nboots=1000)$ks.boot.pvalue,
  ks.boot(mrfss_JX$tl, mrfss_IR$tl, nboots=1000)$ks.boot.pvalue)

p.adjust(ks_ALL_mrfss) 

#MRIP
ks_ALL_mrip <-c(ks.boot(mrip_TB$tl, mrip_CK$tl, nboots=1000)$ks.boot.pvalue, 
                 ks.boot(mrip_TB$tl, mrip_CH$tl, nboots=1000)$ks.boot.pvalue,
                 ks.boot(mrip_TB$tl, mrip_AP$tl, nboots=1000)$ks.boot.pvalue,
                 ks.boot(mrip_TB$tl, mrip_JX$tl, nboots=1000)$ks.boot.pvalue,
                 ks.boot(mrip_TB$tl, mrip_IR$tl, nboots=1000)$ks.boot.pvalue,
                 ks.boot(mrip_CK$tl, mrip_CH$tl, nboots=1000)$ks.boot.pvalue,
                 ks.boot(mrip_CK$tl, mrip_AP$tl, nboots=1000)$ks.boot.pvalue,
                 ks.boot(mrip_CK$tl, mrip_JX$tl, nboots=1000)$ks.boot.pvalue,
                 ks.boot(mrip_CK$tl, mrip_IR$tl, nboots=1000)$ks.boot.pvalue,
                 ks.boot(mrip_CH$tl, mrip_AP$tl, nboots=1000)$ks.boot.pvalue,
                 ks.boot(mrip_CH$tl, mrip_JX$tl, nboots=1000)$ks.boot.pvalue,
                 ks.boot(mrip_CH$tl, mrip_IR$tl, nboots=1000)$ks.boot.pvalue,
                 ks.boot(mrip_AP$tl, mrip_JX$tl, nboots=1000)$ks.boot.pvalue,
                 ks.boot(mrip_AP$tl, mrip_IR$tl, nboots=1000)$ks.boot.pvalue,
                 ks.boot(mrip_JX$tl, mrip_IR$tl, nboots=1000)$ks.boot.pvalue)

p.adjust(ks_ALL_mrip)

#ALL
ks_ALL <-c(ks.boot(TB$tl, CK$tl, nboots=1000)$ks.boot.pvalue, 
                ks.boot(TB$tl, CH$tl, nboots=1000)$ks.boot.pvalue,
                ks.boot(TB$tl, AP$tl, nboots=1000)$ks.boot.pvalue,
                ks.boot(TB$tl, JX$tl, nboots=1000)$ks.boot.pvalue,
                ks.boot(TB$tl, IR$tl, nboots=1000)$ks.boot.pvalue,
                ks.boot(CK$tl, CH$tl, nboots=1000)$ks.boot.pvalue,
                ks.boot(CK$tl, AP$tl, nboots=1000)$ks.boot.pvalue,
                ks.boot(CK$tl, JX$tl, nboots=1000)$ks.boot.pvalue,
                ks.boot(CK$tl, IR$tl, nboots=1000)$ks.boot.pvalue,
                ks.boot(CH$tl, AP$tl, nboots=1000)$ks.boot.pvalue,
                ks.boot(CH$tl, JX$tl, nboots=1000)$ks.boot.pvalue,
                ks.boot(CH$tl, IR$tl, nboots=1000)$ks.boot.pvalue,
                ks.boot(AP$tl, JX$tl, nboots=1000)$ks.boot.pvalue,
                ks.boot(AP$tl, IR$tl, nboots=1000)$ks.boot.pvalue,
                ks.boot(JX$tl, IR$tl, nboots=1000)$ks.boot.pvalue)

p.adjust(ks_ALL)



# Chi-Square Test #####
###########################
# used to test for differences in length frequencies among different groups (among, times, locations or gears)
# chi square will generally fail where the contingency table contains zeroes or very small frequencies 
#  so best to make the length categories to ensure appropriate sample size 
# exclude some low and high values so that I can make more appropriate length intervals 

#length frequencies of the recreational catch will have changed by year based upon regulatory changes so that
# it is invalid to do yearly comparisons of length frequencies because they were affected by slot limit changes
#just compare by bays so combine the data and then make length frequencies
# current is 38 cm to 50 cm  (slot of 15-20)
# started in 2000 but specified above already

#MRFSS
mrfss_all <- subset(rbind(mrfss_AP, mrfss_TB, mrfss_CH, mrfss_CK, mrfss_JX, mrfss_IR), tl>=30 & tl<= 65) %>% mutate(lcat5=lencat(tl, w=5)) 
mrfss_all$bay <- as.factor(mrfss_all$bay)

(mrfss_lf <- xtabs(~bay+lcat5, data=mrfss_all))

#multiple comparisons are conducted by isolating pairs of groups, accumulating p-values for each pair and then adjusting the p values for multiple comparisons
ps_mrfss<- c(chisq.test(mrfss_lf[1:2,])$p.value,
chisq.test(mrfss_lf[1:3,])$p.value,
chisq.test(mrfss_lf[1:4,])$p.value,
chisq.test(mrfss_lf[1:5,])$p.value,
chisq.test(mrfss_lf[1:6,])$p.value,
chisq.test(mrfss_lf[2:3,])$p.value,
chisq.test(mrfss_lf[2:4,])$p.value,
chisq.test(mrfss_lf[2:5,])$p.value,
chisq.test(mrfss_lf[2:6,])$p.value,
chisq.test(mrfss_lf[3:4,])$p.value,
chisq.test(mrfss_lf[3:5,])$p.value,
chisq.test(mrfss_lf[3:6,])$p.value,
chisq.test(mrfss_lf[4:5,])$p.value,
chisq.test(mrfss_lf[4:6,])$p.value,
chisq.test(mrfss_lf[5:6,])$p.value)

p.adjust(ps_mrfss)

#MRIP
mrip_all <- subset(rbind(mrip_AP, mrip_TB, mrip_CH, mrip_CK, mrip_JX, mrip_IR), tl>=30 & tl<= 65) %>% mutate(lcat5=lencat(tl, w=5)) 
mrip_all$bay <- as.factor(mrip_all$bay)

(mrip_lf <- xtabs(~bay+lcat5, data=mrip_all))

#multiple comparisons are conducted by isolating pairs of groups, accumulating p-values for each pair and then adjusting the p values for multiple comparisons
ps_mrip<- c(chisq.test(mrip_lf[1:2,])$p.value,
       chisq.test(mrip_lf[1:3,])$p.value,
       chisq.test(mrip_lf[1:4,])$p.value,
       chisq.test(mrip_lf[1:5,])$p.value,
       chisq.test(mrip_lf[1:6,])$p.value,
       chisq.test(mrip_lf[2:3,])$p.value,
       chisq.test(mrip_lf[2:4,])$p.value,
       chisq.test(mrip_lf[2:5,])$p.value,
       chisq.test(mrip_lf[2:6,])$p.value,
       chisq.test(mrip_lf[3:4,])$p.value,
       chisq.test(mrip_lf[3:5,])$p.value,
       chisq.test(mrip_lf[3:6,])$p.value,
       chisq.test(mrip_lf[4:5,])$p.value,
       chisq.test(mrip_lf[4:6,])$p.value,
       chisq.test(mrip_lf[5:6,])$p.value)

p.adjust(ps_mrip)

#ALL
ALL <- subset(rbind(AP, TB, CH, CK, JX, IR), tl>=30 & tl<= 65) %>% mutate(lcat5=lencat(tl, w=5)) 
ALL$bay <- as.factor(ALL$bay)

(ALL_lf <- xtabs(~bay+lcat5, data=ALL))

#multiple comparisons are conducted by isolating pairs of groups, accumulating p-values for each pair and then adjusting the p values for multiple comparisons
ps_ALL<- c(chisq.test(ALL_lf[1:2,])$p.value, #Ap to CH X-squared = 100.06, df = 6, p-value < 0.00000000000000022
            chisq.test(ALL_lf[1:3,])$p.value, #AP to CK X-squared = 152.07, df = 12, p-value < 0.00000000000000022
            chisq.test(ALL_lf[1:4,])$p.value, #AP to IR X-squared = 702.63, df = 18, p-value < 0.00000000000000022
            chisq.test(ALL_lf[1:5,])$p.value, #AP to JX X-squared = 773.19, df = 24, p-value < 0.00000000000000022
            chisq.test(ALL_lf[1:6,])$p.value, #AP to TB X-squared = 954.13, df = 30, p-value < 0.00000000000000022
            chisq.test(ALL_lf[2:3,])$p.value, #CH to CK X-squared = 75.102, df = 6, p-value = 0.00000000000003656
            chisq.test(ALL_lf[2:4,])$p.value, #CH to IR X-squared = 555.33, df = 12, p-value < 0.00000000000000022
            chisq.test(ALL_lf[2:5,])$p.value, #CH to JX X-squared = 655.97, df = 18, p-value < 0.00000000000000022
            chisq.test(ALL_lf[2:6,])$p.value, #CH to TB X-squared = 769.09, df = 24, p-value < 0.00000000000000022
            chisq.test(ALL_lf[3:4,])$p.value, #CK to IR X-squared = 481.84, df = 6, p-value < 0.00000000000000022
            chisq.test(ALL_lf[3:5,])$p.value, #CK to JX X-squared = 566.06, df = 12, p-value < 0.00000000000000022
            chisq.test(ALL_lf[3:6,])$p.value, #CK to TB X-squared = 723.83, df = 18, p-value < 0.00000000000000022
            chisq.test(ALL_lf[4:5,])$p.value, #IR to JX X-squared = 174.63, df = 6, p-value < 0.00000000000000022
            chisq.test(ALL_lf[4:6,])$p.value, #IR to TB X-squared = 397.14, df = 12, p-value < 0.00000000000000022
            chisq.test(ALL_lf[5:6,])$p.value) #JX to TB X-squared = 145.54, df = 6, p-value < 0.00000000000000022

p.adjust(ps_ALL)

#MULTIPLOT LENGTH HISTOGRAMS #####
length_AP <- ggplot(AP, aes(x=tl))+ 
  geom_histogram(aes(y=(..count..)/sum(..count..)), binwidth=1, boundary=-0.5, fill="white", color="black")+    # plotting in percent frequency
  scale_y_continuous(labels=percent_format(), name="Frequency (%)", limits=c(0,.15), breaks=seq(0,.15, .03))+  ## plotting in percent frequency
  scale_x_continuous(name="Length (cm)", limits=c(30,65), breaks=seq(30,65,5))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_rect(colour="black",fill="white"), 
        axis.text.x=element_text(colour="black", size=24), #changing  colour of x axisaxis.text.y=element_text(colour="black"), #changing colour of y acis
        axis.text.y=element_text(colour="black", size=24), #changing colour of y acis
        axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24),
        plot.title=element_text(size=14))+
  annotate("text", x=60, y= .13, label="AP", size=10)

length_CK <- ggplot(CK, aes(x=tl))+ 
  geom_histogram(aes(y=(..count..)/sum(..count..)), binwidth=1, boundary=-0.5, fill="white", color="black")+    # plotting in percent frequency
  scale_y_continuous(labels=percent_format(), name="Frequency (%)", limits=c(0,.15), breaks=seq(0,.15, .03))+  ## plotting in percent frequency
  scale_x_continuous(name="Length (cm)", limits=c(30,65), breaks=seq(30,65,5))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_rect(colour="black",fill="white"), 
        axis.text.x=element_text(colour="black", size=24), #changing  colour of x axisaxis.text.y=element_text(colour="black"), #changing colour of y acis
        axis.text.y=element_text(colour="black", size=24), #changing colour of y acis
        axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24),
        plot.title=element_text(size=14))+
  annotate("text", x=60, y= .13, label="CK", size=10)

length_TB <- ggplot(TB, aes(x=tl))+ 
  geom_histogram(aes(y=(..count..)/sum(..count..)), binwidth=1, boundary=-0.5, fill="white", color="black")+    # plotting in percent frequency
  scale_y_continuous(labels=percent_format(), name="Frequency (%)", limits=c(0,.15), breaks=seq(0,.15, .03))+  ## plotting in percent frequency
  scale_x_continuous(name="Length (cm)", limits=c(30,65), breaks=seq(30,65,5))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_rect(colour="black",fill="white"), 
        axis.text.x=element_text(colour="black", size=24), #changing  colour of x axisaxis.text.y=element_text(colour="black"), #changing colour of y acis
        axis.text.y=element_text(colour="black", size=24), #changing colour of y acis
        axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24),
        plot.title=element_text(size=14))+
  annotate("text", x=60, y= .13, label="TB", size=10)

length_CH <- ggplot(CH, aes(x=tl))+ 
  geom_histogram(aes(y=(..count..)/sum(..count..)), binwidth=1, boundary=-0.5, fill="white", color="black")+    # plotting in percent frequency
  scale_y_continuous(labels=percent_format(), name="Frequency (%)", limits=c(0,.15), breaks=seq(0,.15, .03))+  ## plotting in percent frequency
  scale_x_continuous(name="Length (cm)", limits=c(30,65), breaks=seq(30,65,5))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_rect(colour="black",fill="white"), 
        axis.text.x=element_text(colour="black", size=24), #changing  colour of x axisaxis.text.y=element_text(colour="black"), #changing colour of y acis
        axis.text.y=element_text(colour="black", size=24), #changing colour of y acis
        axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24),
        plot.title=element_text(size=14))+
  annotate("text", x=60, y= .13, label="CH", size=10)

length_JX <- ggplot(JX, aes(x=tl))+ 
  geom_histogram(aes(y=(..count..)/sum(..count..)), binwidth=1, boundary=-0.5, fill="white", color="black")+    # plotting in percent frequency
  scale_y_continuous(labels=percent_format(), name="Frequency (%)", limits=c(0,.15), breaks=seq(0,.15, .03))+  ## plotting in percent frequency
  scale_x_continuous(name="Length (cm)", limits=c(30,65), breaks=seq(30,65,5))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_rect(colour="black",fill="white"), 
        axis.text.x=element_text(colour="black", size=24), #changing  colour of x axisaxis.text.y=element_text(colour="black"), #changing colour of y acis
        axis.text.y=element_text(colour="black", size=24), #changing colour of y acis
        axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24),
        plot.title=element_text(size=14))+
  annotate("text", x=60, y= .13, label="JX", size=10)

length_IR <- ggplot(IR, aes(x=tl))+ 
  geom_histogram(aes(y=(..count..)/sum(..count..)), binwidth=1, boundary=-0.5, fill="white", color="black")+    # plotting in percent frequency
  scale_y_continuous(labels=percent_format(), name="Frequency (%)", limits=c(0,0.15), breaks=seq(0,0.15, .03))+  ## plotting in percent frequency
  scale_x_continuous(name="Length (cm)", limits=c(30,65), breaks=seq(30,65,5))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_rect(colour="black",fill="white"), 
        axis.text.x=element_text(colour="black", size=24), #changing  colour of x axisaxis.text.y=element_text(colour="black"), #changing colour of y acis
        axis.text.y=element_text(colour="black", size=24), #changing colour of y acis
        axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24),
        plot.title=element_text(size=14))+
  annotate("text", x=60, y= .13, label="IR", size=10)




