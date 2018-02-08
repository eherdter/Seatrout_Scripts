##### ABOUT ######
#5/24/2016 This file imports the indices of abundance produced from the DeltaLogNormal....R file.
# Updated 10/21/16
# ***** Formerly this was called  "Autocorrelation and correlation with DeltaLogNorm Indices.R"
# Main Objectives of this script: 
# 1. Imports Indices
# 2. Plots the indices
# 3. Tests interseries correlation (autocorrelation). If interseries correlation exists 
#   then significance tests must be adjusted based on degrees of freedom (Pyper et al. 2001)
# 4. Pearson product-moment correlation (Mueter et al. 2002, Field and Ralston 2005, Pyper et al. 2001, Peterman et al. 1998, Myers at al. 1997)
# 5. Cluster analysis
# 6. DFA
# 7. PCA
# 8. Tests indices for linear trends

##### LOAD PACKAGES #####
# library(dplyr) NOTE: loading dplyr conflicts with the lag1.plot in astsa. DO NOT LOAD before using ASTSA.
# or load it and then unload it
library(cluster)
library(Hmisc) #for correlation
library(ggplot2)
library(astsa)
library(MARSS) #for DFA
library(dplyr)

setwd("~/Desktop/PhD project/Projects/Seatrout/Data/Indices/DeltaMethod Indices")
setwd("U:/PhD_projectfiles/Exported_R_Datafiles/Indices/UpdatedIndices")

##### 1. IMPORT INDICES #################################
# IMPORT INDICES
# indices produced with Delta Method for Producing Nominal Indices.R
# and select the year, assign bay name,
# scale each index (z score, 0 mean 1 unit variance)

AP<- subset(read.csv("AP_yoy_index.csv", header=T), year>1998, select=c(year, Mean)) %>% mutate(Est= rep("AP",17)) #1997 seems like a crazy outlier so I removed it becuse it was affecting the scaled mean values 
  names(AP) <- c("year", "Mean", "Est") #assign names, est as in estuary not estimate
  AP$Est <- as.factor(AP$Est) #turn bay into a factor
 AP$Mean_scaled<-  as.numeric(scale(AP$Mean, scale=TRUE)) #z-score (mean= 0, var=1), as.numeric is important for full_join command below. doesnt like it when it isnt a true numeric
  # ap does not have a riv component because of so many 1 or 0 positive trips in the FIM data. The csv is available but its mostly NaNs

 CK<- subset(read.csv("CK_yoy_index.csv", header=T), select=c(year, Mean)) %>% mutate(Est= rep("CK",20)) 
 names(CK) <- c("year", "Mean", "Est") #assign names, est as in estuary not estimate
 CK$Est <- as.factor(CK$Est) #turn bay into a factor
 CK$Mean_scaled<-  as.numeric(scale(CK$Mean, scale=TRUE)) #z-score (mean= 0, var=1), as.numeric is important for full_join command below. doesnt like it when it isnt a true numeric
 
 TB<- subset(read.csv("TB_yoy_index.csv", header=T), select=c(year, Mean)) %>% mutate(Est= rep("TB",27)) 
 names(TB) <- c("year", "Mean", "Est") #assign names, est as in estuary not estimate
 TB$Est <- as.factor(TB$Est) #turn bay into a factor
 TB$Mean_scaled<-  as.numeric(scale(TB$Mean, scale=TRUE)) #z-score (mean= 0, var=1), as.numeric is important for full_join command below. doesnt like it when it isnt a true numeric
 

 CH<- subset(read.csv("CH_yoy_index.csv", header=T), select=c(year, Mean)) %>% mutate(Est= rep("CH",27)) 
 names(CH) <- c("year", "Mean", "Est") #assign names, est as in estuary not estimate
 CH$Est <- as.factor(CH$Est) #turn bay into a factor
 CH$Mean_scaled<-  as.numeric(scale(CH$Mean, scale=TRUE)) #z-score (mean= 0, var=1), as.numeric is important for full_join command below. doesnt like it when it isnt a true numeric
 
 IR<- subset(read.csv("IR_yoy_index.csv", header=T), select=c(year, Mean)) %>% mutate(Est= rep("IR",26)) 
 names(IR) <- c("year", "Mean", "Est") #assign names, est as in estuary not estimate
 IR$Est <- as.factor(IR$Est) #turn bay into a factor
 IR$Mean_scaled<-  as.numeric(scale(IR$Mean, scale=TRUE)) #z-score (mean= 0, var=1), as.numeric is important for full_join command below. doesnt like it when it isnt a true numeric
 

 JX<- subset(read.csv("JX_yoy_index.csv", header=T), select=c(year, Mean)) %>% mutate(Est= rep("JX",15)) 
 names(JX) <- c("year", "Mean", "Est") #assign names, est as in estuary not estimate
 JX$Est <- as.factor(JX$Est) #turn bay into a factor
 JX$Mean_scaled<-  as.numeric(scale(JX$Mean, scale=TRUE)) #z-score (mean= 0, var=1), as.numeric is important for full_join command below. doesnt like it when it isnt a true numeric
 

All <- rbind(AP,CH,CK,IR,JX,TB)
All_NW <- rbind(AP, CK)
All_SW <- rbind(TB, CH)
All_East <- rbind(IR, JX)

##### 2. PLOT TIMES SERIES DATA ######
# PLOT TIME SERIES DATA 

# plot either the standardized or non standardized time series on the same plot
# might do this with ggplot
library(ggplot2)


All_mean <- ggplot(All, aes(x=year, y=Mean, group=Est)) + #color=est
  #geom_line(aes(linetype=Est), size =.5)+ # make line types based on the different labels- this will be our workaround because in a few stps we will specify the first label (obserseved) be a blank line (therefore a scatter plot)
  geom_line(aes(color=Est), size=0.5)+
  geom_point(aes(color=Est), size=2) + #, color=bay))+ # groups the points together correctly and then gives them a unique shape them differently based on the line type 
  ylab("LSMean estimate- YOY per haul") +
  xlab("Year")+ 
  scale_y_continuous(limits=c(-0.5,8.5))+
  scale_x_continuous(limits=c(1989, 2015), breaks=seq(1989, 2015, 2))+
  scale_colour_discrete(name="Location")+
  theme(panel.grid.minor=element_blank(), 
        panel.grid.major=element_blank(), 
        panel.background=element_rect(colour="black", fill="white"),
        axis.title.x =element_text(colour="black"),
        axis.text.x = element_text(colour="black"),
        axis.title.y =element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title=element_text(size=14))


All_scaled <- ggplot(All, aes(x=year, y=Mean_scaled, group=Est)) + #color=est
        #geom_line(aes(linetype=Est), size =.5)+ # make line types based on the different labels- this will be our workaround because in a few stps we will specify the first label (obserseved) be a blank line (therefore a scatter plot)
        geom_line(aes(color=Est), size=0.5)+
      geom_point(aes(color=Est), size=2) + #, color=bay))+ # groups the points together correctly and then gives them a unique shape them differently based on the line type 
      ylab("LSMean estimate- YOY per haul- Scaled") +
      xlab("Year")+ 
  scale_y_continuous(limits=c(-2,4.5))+
  scale_x_continuous(limits=c(1989, 2015), breaks=seq(1989, 2015, 2))+
  scale_colour_discrete(name="Location")+
  theme(panel.grid.minor=element_blank(), 
        panel.grid.major=element_blank(), 
        panel.background=element_rect(colour="black", fill="white"),
        axis.title.x =element_text(colour="black"),
        axis.text.x = element_text(colour="black"),
        axis.title.y =element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title=element_text(size=14))

#multiplot function
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

File <- ("U:/PhD_projectfiles/Figures/abundance_timeseries_multiplot.tiff")
if (file.exists(File)) stop(File, " already exists")
dir.create(dirname(File), showWarnings = FALSE)

tiff(File, units="in", width=7, height=7, res=300)

multiplot(All_mean, All_scaled, cols=1)

dev.off()













ggplot(All_NW, aes(x=year, y=Mean, group=Est)) + #color=est
  #geom_line(aes(linetype=Est), size =.5)+ # make line types based on the different labels- this will be our workaround because in a few stps we will specify the first label (obserseved) be a blank line (therefore a scatter plot)
  geom_line(aes(color=Est), size=1)+
  geom_point(aes(color=Est), size=2.5) + #, color=bay))+ # groups the points together correctly and then gives them a unique shape them differently based on the line type 
  ylab("LSMean  #/haul") +
  xlab("Year")+ 
  scale_y_continuous(limits=c(-0.5,5), breaks=seq(-0.5,5, 0.5))+
  scale_x_continuous(limits=c(1995, 2015), breaks=seq(1996, 2015, 3))+
  scale_colour_discrete(name="Location")+
  theme(panel.grid.minor=element_blank(), 
        panel.grid.major=element_blank(), 
        panel.background=element_rect(colour="black", fill="white"),
        axis.title.x =element_text(colour="black", size=20),
        axis.text.x = element_text(colour="black", size=18),
        axis.title.y =element_text(colour="black", size=20),
        axis.text.y = element_text(colour="black", size=18),
        plot.title=element_text(size=14),
        legend.title=element_text(size=14),
        legend.text=element_text(size=14))

ggplot(All_SW, aes(x=Year, y=Mean, group=Est)) + #color=est
  #geom_line(aes(linetype=Est), size =.5)+ # make line types based on the different labels- this will be our workaround because in a few stps we will specify the first label (obserseved) be a blank line (therefore a scatter plot)
  geom_line(aes(color=Est), size=1)+
  geom_point(aes(color=Est), size=2.5) + #, color=bay))+ # groups the points together correctly and then gives them a unique shape them differently based on the line type 
  ylab("LSMean  #/haul") +
  xlab("Year")+ 
  scale_y_continuous(limits=c(0,3.5), breaks=seq(0,3.5, 0.5))+
  scale_x_continuous(limits=c(1989, 2015), breaks=seq(1989, 2015, 4))+
  scale_colour_discrete(name="Location")+
  theme(panel.grid.minor=element_blank(), 
        panel.grid.major=element_blank(), 
        panel.background=element_rect(colour="black", fill="white"),
        axis.title.x =element_text(colour="black", size=20),
        axis.text.x = element_text(colour="black", size=18),
        axis.title.y =element_text(colour="black", size=20),
        axis.text.y = element_text(colour="black", size=18),
        plot.title=element_text(size=14),
        legend.title=element_text(size=14),
        legend.text=element_text(size=14))

ggplot(IR, aes(x=Year, y=Mean, group=Est)) + #color=est
  #geom_line(aes(linetype=Est), size =.5)+ # make line types based on the different labels- this will be our workaround because in a few stps we will specify the first label (obserseved) be a blank line (therefore a scatter plot)
  geom_line(aes(color=Est), size=1)+
  geom_point(aes(color=Est), size=2.5) + #, color=bay))+ # groups the points together correctly and then gives them a unique shape them differently based on the line type 
  ylab("LSMean #/haul") +
  xlab("Year")+ 
  scale_y_continuous(limits=c(0,3), breaks=seq(0,3, 0.5))+
  scale_x_continuous(limits=c(1989, 2015), breaks=seq(1989, 2018, 4))+
  scale_colour_discrete(name="Location")+
  theme(panel.grid.minor=element_blank(), 
        panel.grid.major=element_blank(), 
        panel.background=element_rect(colour="black", fill="white"),
        axis.title.x =element_text(colour="black", size=20),
        axis.text.x = element_text(colour="black", size=18),
        axis.title.y =element_text(colour="black", size=20),
        axis.text.y = element_text(colour="black", size=18),
        plot.title=element_text(size=14),
        legend.title=element_text(size=14),
        legend.text=element_text(size=14))


ggplot(JX, aes(x=Year, y=Mean, group=Est)) + #color=est
  #geom_line(aes(linetype=Est), size =.5)+ # make line types based on the different labels- this will be our workaround because in a few stps we will specify the first label (obserseved) be a blank line (therefore a scatter plot)
  geom_line(aes(color=Est), size=1)+
  geom_point(aes(color=Est), size=2.5) + #, color=bay))+ # groups the points together correctly and then gives them a unique shape them differently based on the line type 
  ylab("LSMean  #/haul") +
  xlab("Year")+ 
  scale_y_continuous(limits=c(-0.01,0.02), breaks=seq(-0.01,0.02, 0.005))+
  scale_x_continuous(limits=c(2000, 2015), breaks=seq(2000, 2015, 3))+
  scale_colour_discrete(name="Location")+
  theme(panel.grid.minor=element_blank(), 
        panel.grid.major=element_blank(), 
        panel.background=element_rect(colour="black", fill="white"),
        axis.title.x =element_text(colour="black", size=20),
        axis.text.x = element_text(colour="black", size=18),
        axis.title.y =element_text(colour="black", size=20),
        axis.text.y = element_text(colour="black", size=18),
        plot.title=element_text(size=14),
        legend.title=element_text(size=14),
        legend.text=element_text(size=14))





#plot each index separately but have to turn year into numeric
#also plotted in Delta_Method R script
AP$Year <- as.numeric(AP$Year)


plot(Mean~year, data=AP, type="l")
plot(Mean_scaled~year, data=CK, type='l')
plot(Mean~year, data=TB, type='l')
plot(Mean~year, data=CH, type='l')
plot(Mean~year, data=JX, type='l')
plot(Mean~year, data=IR, type='l')

#### 8. TEST INDICES FOR LINEAR TRENDS #####
#some linear trends are present for only two estuaries. This may not matter in duture steps because the other steps are done with the scaled mean values
summary(lm(Mean ~ year, data=AP))
summary(lm(Mean ~ year, data=CK))
summary(lm(Mean ~ year, data=TB))
summary(lm(Mean ~ year, data=CH))
summary(lm(Mean ~ year, data=JX))
summary(lm(Mean ~ year, data=IR))
#some have slight linear trends but the adjusted R squared values are really small <0.25 

#test scaled indices for linear trends
summary(lm(Mean_scaled ~ year, data=AP))
summary(lm(Mean_scaled ~ year, data=CK))
summary(lm(Mean_scaled ~ year, data=TB))
summary(lm(Mean_scaled ~ year, data=CH))
summary(lm(Mean_scaled ~ year, data=JX))
summary(lm(Mean_scaled ~ year, data=IR))





##### 3. LAGGED SCATTERPLOTS OF EACH TIMESERIES#######
#LAGGED SCATTERPLOTS of each time series
# From Notes_3, GEOS 585A, Spring 2015 handout printed from arizona (website above)
# " An attribute of the lagged scatterplot is that it can display autocorrelation regardless of the form of the
# dependence on the past values. An assumption of linear dependence is not necessary."
# Such nonlinear dependence might not be effectively summarized by other methods (such as the acf function)

#1A.1 Determine Critical level of correlation for 95% significance (alpha = 0.5) r= 0+- 2/sqrt(N)
#AP = +-2/sqrt(17)= 0.48
#CH = +-2/sqrt(27)= 0.38 1989 ->
#CK = +-2/sqrt(20)= 0.44 1996- > 
#IR = +- 2/sqrt(26) =0.39
#TB = +- 2/sqrt(27) = 0.38 ***
#JX= +- 2/sqrt(14) = 0.53

#unload dplyr
detach("package:dplyr", unload=TRUE)
library(astsa)
lag1.plot(AP, 3, corr=TRUE) #autocorrelated up to third lag
lag1.plot(CH[8:27,], 4, corr=TRUE) #autocorrelated up to 3rd lag
lag1.plot(CK, 4, corr=TRUE) #autocorrelated to 3rd lag
lag1.plot(IR, 5, corr=TRUE) #autocorrelated to 6th lag
lag1.plot(TB, 5, corr=TRUE) #autocorrelation at 1st and 2nd lag
lag1.plot(JX, 3, corr=TRUE) #autocorrelated to 2nd lag

#1B. ACF and correlograms
acf2(AP$Mean)
acf2(CH$Mean)
acf2(CK$Mean)
acf2(IR$Mean)
acf2(TB$Mean)
acf2(JX$Mean)


#MOdified Chelton ####
#adjust degrees of freedom to number of independent pairs (determine at what lag each series becomse uncorrelated)
# use this table and the corr value to reevaulate significance
# http://users.sussex.ac.uk/~grahamh/RM1web/Pearsonstable.pdf

#moderate evidence for auto correlation so maybe cant use normal significance tests for pairwise correlations??

## 4.PEARSON CORRELATION TESTS######
# PEARSON CORRELATION

# To not end up with a bunch of "est" columns I will just make some new dfs that do not contain the est columns

ap_min <- subset(AP, select=(-Est)) 
ch_min <- subset(CH, select=(-Est))
ck_min <- subset(CK, select=(-Est))
ir_min <- subset(IR, select=(-Est))
tb_min <- subset(TB, select=(-Est))
jx_min <- subset(JX, select=(-Est))

# Make dataframe with just indices by year

# For some reason full_join stops working after 4 dataframes at which point it starts duplicating things so 
# I just made two and joined them together

ind1 <- full_join(ap_min, ch_min, by='year') %>% 
            full_join(.,ck_min, by='year' ) 
names(ind1) <- c("year", "AP","AP_Scaled" , "CH", "CH_Scaled", "CK", "CK_Scaled")

ind2 <-  full_join(tb_min,ir_min, by='year' ) %>%
            full_join(.,jx_min, by='year' ) 
names(ind2) <- c("year","TB","TB_Scaled", "IR", "IR_Scaled", "JX", "JX_Scaled" )

ind_Year <- full_join(ind1, ind2, by='year')

mat<- as.matrix(arrange(ind_Year, year)) #arrange the data,frame by year and then turn into a matrix
mattest=mat[,-1] #remove the first column which is the year to just have a matrix of scaled indices


# Now pearson correlation with rcorr from hmisc
#rcorr- missing values are deleted in pairs rather than deleting all rows of x having any missing variables

corr_mat_ALL <- rcorr(mattest, type="pearson")
rho_mat_ALL <- as.data.frame(corr_mat_ALL$r) #rho values
n_mat_ALL <- as.data.frame(corr_mat_ALL$n) #number of samples used to compute correlation
P_mat_ALL <- as.data.frame(corr_mat_ALL$P) #P values

#P_mat_adjust_ALL <- as.data.frame(p.adjust(P_mat_ALL)) #adjust the p-values for multiple comparisons but is this the right thing to do???

#select the non scaled ones
rho_mat <- (as.data.frame(rho_mat_ALL %>% select(1,3,5,7,9,11)))[c(1,3,5,7,9,11),] #select columns and rows 
n_mat <- (as.data.frame(n_mat_ALL %>% select(1,3,5,7,9,11)))[c(1,3,5,7,9,11),]
P_mat <- (P_mat_ALL %>% select(1,3,5,7,9,11))[c(1,3,5,7,9,11),]

#select scaled ones
rho_mat_scl <- (as.data.frame(rho_mat_ALL %>% select(2,4,6,8,10,12)))[c(2,4,6,8,10,12),] #select columns and rows 
n_mat_scl <- (as.data.frame(n_mat_ALL %>% select(2,4,6,8,10,12)))[c(2,4,6,8,10,12),]



#unwrap each matrix
library(gdata)
rho_vec <- as.data.frame(lowerTriangle(rho_mat, diag=FALSE, byrow=FALSE))
p_vec   <- lowerTriangle(P_mat, diag=FALSE, byrow=FALSE)
p.adjust(p_vec)

#produce sample number (impt for fitting the nls function and weighting by sample number in Medoid_GC_NLSfit_plot.R)
n_mat <- as.data.frame(lowerTriangle(n_mat, diag=FALSE, byrow=FALSE)) 
rho_P_vec <- cbind(rho_vec, p_vec, n_mat)
ronam <- c( "AP_CH","AP_CK",  "AP_TB",  "AP_IR", "AP_JX", 
              "CH_CK", "CH_TB", "CH_IR", "CH_JX",
              "CK_TB", "CK_IR", "CK_JX",
              "TB_IR", "TB_JX", "IR_JX")

row.names(rho_P_vec)<- ronam
data.frame(rho_P_vec)
colnames(rho_P_vec)<-  c("rho", "P", "N")
## export the rho_dataset to be used by Medoid_GC_NLSfit_plot.R when comparing great circle distances to rho
write.csv(rho_P_vec, "U:/PhD_projectfiles/Exported_R_Datafiles/Rho_P_vector.csv")


corr()








### 5. CLUSTER ANALYSIS####
# CLUSTER ANALYSIS

#using compelete linkage method and the correlation matrix as the distance matrix

library(cluster)
ac <- agnes(rho_mat)
plot(ac, ask=TRUE)

##### 7. PRINCIPAL COMPONENTS######
# Principal Components 

library(vegan)

pcamat <- mat[11:25,2:7]
rownames(pcamat) <- c("2001", "2002", "2003", "2004", "2005", "2006", "2007", "2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015")

pcaind <- rda(pcamat)


plot.new()
plot(pcaind, type="n", xlab="", ylab="", scaling="species") #when type="n" no points. it just sets the frame of the plot
#remove the generic plot labels and then add them at the end- see below
text(data.rda, dis="si") #dis is short for display- when dis="sp" it will plot response variables(species); when its ="si" it will plot the years(sites); when its = "cn" it will plot the predictor vectors
#points(data.rda, pch=21, col="red", bg="yellow", cex=1.2) - this is if we want the years(sites) plotted as unidentified points
text(pcaind, "species",  col="blue", cex=0.8) #plots the response data (species)
text(pcaind, "sites", col="red", cex=0.8)
title(xlab="PCA 1 (33.45%)", ylab= "PCA 2 (28.8%)")  #labeled axes with %variance explained. Unfortunately it doesnt look like this is a default in the vegan package so I had to hand calculcate below. 

##### 6. DFA######
# DFA

#page 120 in Holmes-Analysis of multivariate using MARSS
# first fit two trends
# form argument for a MARSS call to specify that we want to it a
# DFA model. This will set up the Z matrix and the other parameters for you.
# Specify how many trends you want by passing in model=list(m=x). You can
# also pass in different forms for the R matrix in the usual way.

indices <- mat[,2:7] #select just the indices- remove year 
indices<- t(indices) #transpose the matrix (not sure why but says so in the instructions)
N.ts = dim(indices)[1]

# # Hand build the Z matrix for a model with 6 observed time series and 2 hidden trends
# Z.vals = list("z11", 0,
#               "z21", "z22",
#               "z31", "z32",
#               "z41", "z42",
#               "z51", "z52",
#               "z61", "z62")
# Z=matrix(Z.vals, nrow = N.ts, ncol=3, byrow=TRUE)




model.list=list(m=2, R="diagonal and equal") #define the form of the model- two trends and the R matrix to be different variances but same covariance
cntl.list = list(minit=200, maxit=5000, allow.degen=FALSE) #set control parameters
cntl.list = list(maxit=50)

mod.2=MARSS(indices, model=model.list, z.score=TRUE, form='dfa', control=cntl.list)

#now fit three trends
model.list=list(m=3, R="diagonal and equal") #define the form of the model- two trends and the R matrix to be different variances but same covariance
mod.3=MARSS(indices, model=model.list, z.score=TRUE, form='dfa', control=cntl.list)


#now fit 1 trend
model.list=list(m=1, R="diagonal and equal") #define the form of the model- two trends and the R matrix to be different variances but same covariance
mod.1=MARSS(indices, model=model.list, z.score=TRUE, form='dfa', control=cntl.list)


# now compare the AIC value from the two trend model to the three trend model and the one trend model
print(cbind(model=c("3 trends", "2 trends", "1 trend"),
            AICc=round(c(mod.3$AICc, mod.2$AICc, mod.1$AICc))),
      quote=FALSE)


# page 121 Holmes- Analysis of multivariate using MARSS
# build a few models 
# test from 1 to 3 underlying trends

# Following code builds the model matrices
# different structures for the R matrix in the DFA analysis equation
# 1. same variances & no covariance (``diagonal and equal'');
# 2. different variances & no covariance (``diagonal and unequal'');
# 3. same variances & same covariance (``equalvarcov''); and
# 4. different variances & covariances (``unconstrained'').

cntl.list = list(minit=200, maxit=5000, allow.degen=FALSE) #set control parameters
levels.R = c("diagonal and equal",
             "diagonal and unequal",
             "equalvarcov",
             "unconstrained")   

model.data = data.frame()
# fit lots of models & store results
# NOTE: this will take a long time to run!
for(R in levels.R) {
  for(m in 1:(N.ts-1)) {
    dfa.model = list(A="zero", R=R, m=m)
    kemz = MARSS(indices, model=dfa.model, control=cntl.list,
                 form="dfa", z.score=TRUE)
    model.data = rbind(model.data,
                       data.frame(R=R,
                                  m=m,
                                  logLik=kemz$logLik,
                                  K=kemz$num.params,
                                  AICc=kemz$AICc,
                                  stringsAsFactors=FALSE))
    assign(paste("kemz", m, R, sep="."), kemz)
  } # end m loop
} # end R loop

#get the "best" model

best.model <- mod.3[1,]
fitname = paste("mod.3",best.model$m, best.model$R,sep=".")
best.fit = get(fitname)

# get the inverse of the rotation matrix
H.inv = varimax(coef(best.fit, type="matrix")$Z)$rotmat
H.inv = varimax(coef(mod.3, type="matrix")$Z)$rotmat

# rotate factor loadings
Z.rot = coef(mod.3, type="matrix")$Z %*% H.inv

# rotate trends
trends.rot = solve(H.inv) %*% mod.3$states


##### Plotting Example- to produce plot on page 118 of Holmes-MARSS #####
#https://github.com/cran/MARSS/blob/master/inst/doc/Chapter_DFA.R

fit = kemz.3
spp = rownames(dat.z)
par(mfcol=c(3,2), mar=c(3,4,1.5,0.5), oma=c(0.4,1,1,1))
for(i in 1:length(spp)){
  plot(dat.z[i,],xlab="",ylab="abundance index",bty="L", xaxt="n", ylim=c(-4,3), pch=16, col="blue")
  axis(1,12*(0:dim(dat.z)[2])+1,1980+0:dim(dat.z)[2])
  par.mat=coef(fit,type="matrix")
  lines(as.vector(par.mat$Z[i,,drop=FALSE]%*%fit$states+par.mat$A[i,]), lwd=2)
  title(spp[i])
}


#Mine# 
fit=mod.1 
spp=rownames(indices)
par(mfcol=c(3,3), mar=c(3,4,1.5,0.5), oma=c(0.4,1,1,1))
for(i in 1:length(spp)){
  plot(indices[i,],xlab="",ylab="abundance index",bty="L", xaxt="n", ylim=c(-4,3), pch=16, col="blue")
  axis(1,12*(0:dim(indices)[2])+1,1990+0:dim(indices)[2])
  par.mat=coef(fit,type="matrix")
  lines(as.vector(par.mat$Z[i,,drop=FALSE]%*%fit$states+par.mat$A[i,]), lwd=2)
  title(spp[i])
}




### code chunk number 29: Cs23_plotfacloadings

spp = rownames(dat.z)
minZ = 0.05
ylims = c(-1.1*max(abs(Z.rot)), 1.1*max(abs(Z.rot)))
par(mfrow=c(ceiling(dim(trends.rot)[1]/2),2), mar=c(3,4,1.5,0.5), oma=c(0.4,1,1,1))
for(i in 1:best.model$m) {
  plot(c(1:N.ts)[abs(Z.rot[,i])>minZ], as.vector(Z.rot[abs(Z.rot[,i])>minZ,i]),
       type="h", lwd=2, xlab="", ylab="", xaxt="n", ylim=ylims, xlim=c(0,N.ts+1))
  for(j in 1:N.ts) {
    if(Z.rot[j,i] > minZ) {text(j, -0.05, spp[j], srt=90, adj=1, cex=0.9)}
    if(Z.rot[j,i] < -minZ) {text(j, 0.05, spp[j], srt=90, adj=0, cex=0.9)}
    abline(h=0, lwd=1, col="gray")
  } # end j loop
  mtext(paste("Factor loadings on trend",i,sep=" "),side=3,line=.5)
} # end i loop



### code chunk number 30: Cs24_plottrends

# get ts of trends
ts.trends = t(trends.rot)
par(mfrow=c(ceiling(dim(ts.trends)[2]/2),2), mar=c(3,4,1.5,0.5), oma=c(0.4,1,1,1))
# loop over each trend
for(i in 1:dim(ts.trends)[2]) {
  # set up plot area
  plot(ts.trends[,i],
       ylim=c(-1.1,1.1)*max(abs(ts.trends)), 
       type="n", lwd=2, bty="L", 
       xlab="", ylab="", xaxt="n", yaxt="n")
  # draw zero-line
  abline(h=0, col="gray")
  # plot trend line
  par(new=TRUE)
  plot(ts.trends[,i],
       ylim=c(-1.1,1.1)*max(abs(ts.trends)), 
       type="l", lwd=2, bty="L", 
       xlab="", ylab="", xaxt="n")
  # add panel labels
  mtext(paste("Trend",i,sep=" "), side=3, line=0.5)
  axis(1,12*(0:dim(dat.spp.1980)[2])+1,1980+0:dim(dat.spp.1980)[2])
} # end i loop (trends)



### code chunk number 31: Cs25_plotbestfits

par.mat=coef(best.fit, type="matrix")
fit.b = par.mat$Z %*% best.fit$states + matrix(par.mat$A, nrow=N.ts, ncol=TT)
spp = rownames(dat.z)
par(mfcol=c(3,2), mar=c(3,4,1.5,0.5), oma=c(0.4,1,1,1))
for(i in 1:length(spp)){
  plot(dat.z[i,],xlab="",ylab="abundance index",bty="L", xaxt="n", ylim=c(-4,3), pch=16, col="blue")
  axis(1,12*(0:dim(dat.z)[2])+1,1980+0:dim(dat.z)[2])
  lines(fit.b[i,], lwd=2)
  title(spp[i])
}


##### Example from the Text: HOLMES #####


# data(lakeWAplankton)
# # we want lakeWAplanktonTrans, which has been transformed
# # so the 0s are replaced with NAs and the data z-scored
# dat = lakeWAplanktonTrans
# # use only the 10 years from 1980-1989
# plankdat = dat[dat[,"Year"]>=1980 & dat[,"Year"]<1990,]
# # create vector of phytoplankton group names
# phytoplankton = c("Cryptomonas", "Diatoms", "Greens",
#                   "Unicells", "Other.algae")
# # get only the phytoplankton
# dat.spp.1980 = plankdat[,phytoplankton]
# 
# # transpose data so time goes across columns
# dat.spp.1980 = t(dat.spp.1980)
# # get number of time series
# N.ts = dim(dat.spp.1980)[1]
# # get length of time series
# TT = dim(dat.spp.1980)[2]
# 
# model.list = list(m=3, R="diagonal and unequal")
# kemz.3 = MARSS(dat.spp.1980, model=model.list,
#                z.score=TRUE, form="dfa", control=cntl.list)
# 
# best.model <- kemz.3[1,]
# fitname = paste("kemz.2",best.model$m, best.model$R,sep=".")
# best.fit = get(fitname)
# 
# # get the inverse of the rotation matrix
# H.inv = varimax(coef(best.fit, type="matrix")$Z)$rotmat
# H.inv = varimax(coef(mod.3, type="matrix")$Z)$rotmat
# 
# # rotate factor loadings
# Z.rot = coef(mod.3, type="matrix")$Z %*% H.inv
# 
# # rotate trends
# trends.rot = solve(H.inv) %*% mod.3$states
# 
