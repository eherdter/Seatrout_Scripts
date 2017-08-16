# Script for figuring out the medoids in all zones so they can be used as distance references
# 6/1/2016 - Edited script for 
# 1. Determining medoids of each estuary so they can be used as great circle (GC) distance references
# 2. ESTIMATE SPATIAL DECAY WITH NLS MODEL for both recruitment indices and RESIDUALS
# 3. Determines lats and longs for adult data as well. 
# 4. Produces a CSV file with all of the lats and longs


setwd("~/Desktop/Github Repo/Seatrout/FWRI SCRATCH FOLDER/Elizabeth Herdter/SAS data sets/FIMData/NEWNov7")
library(haven) #for loading SAS data
library(cluster) #for finding medoid
library(geosphere) # for calculating great circle distances between each medoid
library(ggplot2) #for plotting
library(nlstools) #for fitting nls models
library(dplyr)

################################
# SELECT DATA for YOY 
# same way as in Delta Method for Producing Nominal Indices.R
########################################
ap = subset(read_sas("ap_yoy_cn_c.sas7bdat"), month %in% c(6,7,8,9,10,11))
apl = read_sas("ap_yoy_cn_l.sas7bdat")

ck = subset(read_sas("ck_yoy_cn_c.sas7bdat"),  month %in% c(5,6,7,8,9,10,11)) 
ckl = read_sas("ck_yoy_cn_l.sas7bdat")

ch = subset(read_sas("ch_yoy_cn_c.sas7bdat"), month %in% c(4,5,6,7,8,9,10)) 
chl = read_sas("ch_yoy_cn_l.sas7bdat")

tb = subset(read_sas("tb_yoy_cn_c.sas7bdat"), month %in% c(4,5,6,7,8,9,10)) 
tbl = read_sas("tb_yoy_cn_l.sas7bdat")

ir = subset(read_sas("ir_yoy_cn_c.sas7bdat"), month %in% c(5,6,7,8,9,10,11)) 
irl = read_sas("ir_yoy_cn_l.sas7bdat")

jx = subset(read_sas("jx_yoy_cn_c.sas7bdat") , month %in% c(5,6,7,8,9,10,11))
jxl = read_sas("jx_yoy_cn_l.sas7bdat")


ap_all = ap
tb_all=tb
ch_all=ch
ck_all=ck
ir_all=ir
jx_all=jx


########################################
# MEDOID CALCULATION
#########################################


### Figuring out whether to use the mediod or the mean. The mediod is an actual point in the data set whereas the mean is a number not included. 
AP_LL <- subset(ap_all, select=c("Longitude", "Latitude")) 
  APMed <- pam(AP_LL,1)$medoids
    #Lat <- APmed[,2]
    #Long <- TBAmed[,1]
  #plot(TBA_LL$Longitude ~ TBA_LL$Latitude)
  #  points(Lat, Long, pch=16, col="red")

#   meanLL <- data.frame(colMeans(LatLong))
#   test <- cbind(meanLL[1,], meanLL[2,])

#   plot(TB_BAY_AUn$Longitude ~ TB_BAY_AUn$Latitude)
#   cords <- xy.coords(-82.63272,27.9578)
#   points(test, pch=16, col="red")
CH_LL <- subset(ch_all, select=c("Longitude", "Latitude")) 
  CHMed <- pam(CH_LL,1)$medoids

CK_LL <- subset(ck_all, select=c("Longitude", "Latitude")) 
  CKMed <- pam(CK_LL,1)$medoids

IR_LL <-subset(ir_all, select=c("Longitude", "Latitude"))
  IRMed <- pam(IR_LL,1)$medoids

TB_LL <-na.omit(subset(tb_all, select=c("Longitude", "Latitude")))
  TBMed <- pam(TB_LL,1)$medoids

JX_LL <-subset(jx_all, select=c("Longitude", "Latitude"))
  JXMed <- pam(JX_LL,1)$medoids

########################################################
## CALCULATING GREAT CIRCLE DISTANCES
#########################################################

AP_compars = rbind(distGeo(APMed, CHMed, a=6378137, f=1/298.257223563),
                        distGeo(APMed, CKMed, a=6378137, f=1/298.257223563),
                        distGeo(APMed, IRMed, a=6378137, f=1/298.257223563),
                        distGeo(APMed, TBMed, a=6378137, f=1/298.257223563),
                        distGeo(APMed, JXMed, a=6378137, f=1/298.257223563))

CH_compars = rbind(distGeo(CHMed, CKMed, a=6378137, f=1/298.257223563),
                   distGeo(CHMed, IRMed, a=6378137, f=1/298.257223563),
                   distGeo(CHMed, TBMed, a=6378137, f=1/298.257223563),
                   distGeo(CHMed, JXMed, a=6378137, f=1/298.257223563))


CK_compars = rbind(distGeo(CKMed, IRMed, a=6378137, f=1/298.257223563),
                   distGeo(CKMed, TBMed, a=6378137, f=1/298.257223563),
                   distGeo(CKMed, JXMed, a=6378137, f=1/298.257223563))


IR_compars = rbind(distGeo(IRMed, TBMed, a=6378137, f=1/298.257223563),
                    distGeo(IRMed, JXMed, a=6378137, f=1/298.257223563))

TB_compars = distGeo(TBMed, JXMed, a=6378137, f=1/298.257223563)



#join distances together. make sure they are in the correct order so that they can bind together with the csv imported below. 
distances <- rbind(AP_compars, CH_compars, CK_compars, IR_compars, TB_compars)
rownames(distances) <- c("AP_CH", "AP_CK", "AP_IR", "AP_TB", "AP_JX", 
                         "CH_CK", "CH_IR", "CH_TB", "CH_JX",
                         "CK_IR", "CK_TB", "CK_JX",
                         "IR_TB", "IR_JX",
                         "TB_JX")

#bring in the rho_P_vector and combinind with the distances
rho_P_vector <- read.csv('~/Desktop/Github Repo/Seatrout/Data/Exported R Dataframes/rho_P_vector.csv')
rho_P_vector_residuals_RK <- read.csv('~/Desktop/Github Repo/Seatrout/Data/Exported R Dataframes/rho_P_vector_residuals_RK.csv')
rho_P_vector_residuals_BH <- read.csv('~/Desktop/Github Repo/Seatrout/Data/Exported R Dataframes/rho_P_vector_residuals_BH.csv')


#combine rho and distances vectors and then turn distance into kilometers, name x and y for plotting convention below
t<-cbind(rho_P_vector, distances) %>% mutate(distancesKM=distances/1000) %>% arrange(distancesKM) %>% rename(y=rho, x=distancesKM)
r<-cbind(rho_P_vector_residuals_RK, distances) %>% mutate(distancesKM=distances/1000) %>% arrange(distancesKM) %>% rename(y=rho, x=distancesKM)
b<-cbind(rho_P_vector_residuals_BH, distances) %>% mutate(distancesKM=distances/1000) %>% arrange(distancesKM) %>% rename(y=rho, x=distancesKM)


write.csv(t,'~/Desktop/Github Repo/Seatrout/Data/Exported R Dataframes/rho_vs_distance.csv' )
write.csv(r,'~/Desktop/Github Repo/Seatrout/Data/Exported R Dataframes/rho_resid_RK_vs_distance.csv' )
write.csv(b,'~/Desktop/Github Repo/Seatrout/Data/Exported R Dataframes/rho_resid_BH_vs_distance.csv' )
r_min = r[-13,]  #testing out what happens if we remove the outlier that is CH-JX


############################################
#ESTIMATE SPATIAL DECAY WITH NLS MODEL
#1. Fix p0= 1, just exclude it from the equation
#2. Estimate p0 
#Pyper et al. 2001
#Peterman et al. 1998
#P(d)= p0e(-d/v)- constrain p0 to 1
#P(d)=p0e(-d/v) - estimate P0
#v (e-folding scale) where e-folding scale tells distance 
#############################################
setwd('~/Desktop/Github Repo/Seatrout/Data/Exported R Dataframes')
t <- read.csv('~/Desktop/Github Repo/Seatrout/Data/Exported R Dataframes/rho_vs_distance.csv', header=T)

# WITH RECRUITMENT INDICES
#must find acceptable starting values (can play around with these in excel-rho_vs_distance_param_start_estimation.csv)
p0 = 1
v = 125
#p0 constrained to 1, weighted by the total sample number which is shared between each comparison
m1 = nls(y ~exp(-(x/v)), start=list(v=v), data=t, weights=(N))
overview(m1)
#estimated v 174.04(SE 38.29)

# p0 estimated
p0=0.9
v=125
m2 = nls(y ~p0*exp(-(x/v)), start=list(p0=p0,v=v), data=t, weights=(N))
summary(m2)
yfitted <- predict(m2)
boot <- nlsBoot(m2, niter=2000) #bootstrapping 
bootCI <- boot$bootCI #gives confidence intervals and median values for p0 and v (e-folding scale) where e-folding scale tells distance 

#evalute fit with anova and F test. if they are significantly different than F test will indicated significant p value
anova(m1, m2)
# Interpreting F table: https://sakai.duke.edu/access/content/group/25e08a3d-9fc4-41b0-a7e9-815732c1c4ba/New%20folder/Stat%20Topic%20Files/Non-Linear%20Regression/FTestTutorial.pdf
# If p-value is significant, then it indicates the more complex model fits the data significantly better than the simpler
# In our case the p value is not significant so we don't need to estimate p0- we can constrain it to 1


#WITH RESIDUAL TIMESERIES
#RK
p0 = 1
v = 125
#p0 constrained to 1, weighted by the total sample number which is shared between each comparison
m3 = nls(y ~exp(-(x/v)), start=list(v=v), data=r, weights=(N))
overview(m3)
#estimated v=146.65 SE(25.14)


# p0 estimated
p0=0.9
v=125
m4 = nls(y ~p0*exp(-(x/v)), start=list(p0=p0,v=v), data=r, weights=(N))
summary(m4)
yfitted <- predict(m4)
boot <- nlsBoot(m4, niter=2000) #bootstrapping 
bootCI <- boot$bootCI #gives confidence intervals and median values for p0 and v (e-folding scale) where e-folding scale tells distance 

#evalute fit with anova and F test. if they are significantly different than F test will indicated significant p value
anova(m3, m4)

#BH
p0 = 1
v = 125
#p0 constrained to 1, weighted by the total sample number which is shared between each comparison
m5 = nls(y ~exp(-(x/v)), start=list(v=v), data=b, weights=(N))
overview(m5)
#estimated v=148.54 SE(25.64)


# p0 estimated
p0=0.9
v=125
m6 = nls(y ~p0*exp(-(x/v)), start=list(p0=p0,v=v), data=b, weights=(N))
overview(m6)
yfitted <- predict(m4)
boot <- nlsBoot(m4, niter=2000) #bootstrapping 
bootCI <- boot$bootCI #gives confidence intervals and median values for p0 and v (e-folding scale) where e-folding scale tells distance 

#evalute fit with anova and F test. if they are significantly different than F test will indicated significant p value
anova(m5, m6)


# Interpreting F table: https://sakai.duke.edu/access/content/group/25e08a3d-9fc4-41b0-a7e9-815732c1c4ba/New%20folder/Stat%20Topic%20Files/Non-Linear%20Regression/FTestTutorial.pdf
# If p-value is significant, then it indicates the more complex model fits the data significantly better than the simpler
# In our case the p value is not significant so we don't need to estimate p0- we can constrain it to 1


#################################
#PLOT CORRELATION BY DISTANCE
#################################
#Plot correlation by distance and add the fitted curve to the plot(here p0 is estimated. see below for explanation)
#http://stackoverflow.com/questions/25030653/fitting-with-ggplot2-geom-smooth-and-nls
# NOTE: within geom_smooth  it does not recognize variable names; they must be named x and y; weight must be in an aesthetic in ggplot

rhobydis <- ggplot(data=t, aes(x=x, y=y))+geom_point()+ 
  geom_smooth(method="nls",formula=y ~exp(-(x/v)), method.args=list(start=c(v=150)), aes(weight=N), se=FALSE, color="black", size=0.5)+                                           
  ylab("Correlation") +
  xlab("Distance (km)")+ 
  geom_vline(xintercept = 174, linetype="dotted")+
  scale_x_continuous(limits=c(75,500))+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(),  
        panel.background=element_rect(fill='white', colour='black'),                                                    
        axis.text.x=element_text(colour="black"), #changing  colour of x axis
        axis.text.y=element_text(colour="black"), #changing colour of y axis
        plot.title=element_text(size=14), # changing size of plot title)+
        legend.text=element_text(size=10))+
  annotate("text", x=90, y=-0.35, label="(A)", size=5, family="Times New Roman")




#check to make sure that ggplot is producing correct plot
t_fitted <- cbind(t, yfitted)
plot(y ~x, data=t)


#### PLOT CORRELATION BY DISTANCE WITH RESIDUAL RHOs


resid_RK_rhobydis <- ggplot(data=r, aes(x=x, y=y))+geom_point()+ 
  geom_smooth(method="nls",formula=y ~exp(-(x/v)), method.args=list(start=c(v=150)), aes(weight=N), se=FALSE, color="black", size=0.5)+                                           
  ylab("Correlation") +
  xlab("Distance (km)")+ 
  geom_vline(xintercept = 146, linetype="dotted")+
  scale_x_continuous(limits=c(75,500))+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(),  
        panel.background=element_rect(fill='white', colour='black'),                                                    
        axis.text.x=element_text(colour="black"), #changing  colour of x axis
        axis.text.y=element_text(colour="black"), #changing colour of y axis
        plot.title=element_text(size=14), # changing size of plot title)+
        legend.text=element_text(size=10))+
  annotate("text", x=90, y=-0.35, label="(B)", size=5, family="Times New Roman")


resid_BH_rhobydis <- ggplot(data=b, aes(x=x, y=y))+geom_point()+ 
  geom_smooth(method="nls",formula=y ~exp(-(x/v)), method.args=list(start=c(v=150)), aes(weight=N), se=FALSE, color="black", size=0.5)+                                           
  ylab("Correlation") +
  xlab("Distance (km)")+ 
  geom_vline(xintercept = 148, linetype="dotted")+
  scale_x_continuous(limits=c(75,500))+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(),  
        panel.background=element_rect(fill='white', colour='black'),                                                    
        axis.text.x=element_text(colour="black"), #changing  colour of x axis
        axis.text.y=element_text(colour="black"), #changing colour of y axis
        plot.title=element_text(size=14), # changing size of plot title)+
        legend.text=element_text(size=10))+
  annotate("text", x=90, y=-0.35, label="(C)", size=5, family="Times New Roman")


multi <- multiplot(rhobydis, resid_RK_rhobydis, resid_BH_rhobydis, cols=1)
#############################################
# PRODUCE LAT AND LONG DATA FOR YOY AND ADULT
#############################################

#yoy
yoy_LL=rbind(AP_LL, CH_LL, CK_LL, TB_LL, IR_LL, JX_LL) %>% mutate(group=rep("yoy", 48179) )

#adult 
setwd("~/Desktop/Github Repo/Seatrout/FWRI SCRATCH FOLDER/Elizabeth Herdter/SAS data sets/FIMData")
ap_ad = subset(read_sas("ap_adult_cn_c.sas7bdat"))
APad_LL <- subset(ap_ad, select=c("Longitude", "Latitude")) 

ch_ad = subset(read_sas("ch_adult_cn_c.sas7bdat")) # *******
CHad_LL <- subset(ch_ad, select=c("Longitude", "Latitude")) 

ck_ad = subset(read_sas("ck_adult_cn_c.sas7bdat"))
CKad_LL <- subset(ck_ad, select=c("Longitude", "Latitude")) 

tb_ad = subset(read_sas("tb_adult_cn_c.sas7bdat"))
TBad_LL <- subset(tb_ad, select=c("Longitude", "Latitude")) 

jx_ad = subset(read_sas("jx_adult_cn_c.sas7bdat"))
JXad_LL <- subset(jx_ad, select=c("Longitude", "Latitude")) 

ir_ad = subset(read_sas("ir_adult_cn_c.sas7bdat"))
IRad_LL <- subset(ir_ad, select=c("Longitude", "Latitude")) 

adult_LL=rbind(APad_LL, CHad_LL, CKad_LL, TBad_LL, IRad_LL, JXad_LL) %>% mutate(group=rep("adult", 26957) )

combined_LL <- rbind(yoy_LL, adult_LL)
write.csv(combined_LL,'~/Desktop/Github Repo/Seatrout/Data/Exported R Dataframes/FIM_lats_longs.csv' )

