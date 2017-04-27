# 5/24/2016 This file imports the indices of abundance produced from the DeltaLogNormal....R file.
# Updated 10/21/16
# ***** Formerly this was called  "Autocorrelation and correlation with DeltaLogNorm Indices.R"
# Main Objectives of this script: 
# 1. Plots the indices
# 2. Tests interseries correlation (autocorrelation). If interseries correlation exists 
#   then significance tests must be adjusted based on degrees of freedom (Pyper et al. 2001)
# 3. Pearson product-moment correlation (Mueter et al. 2002, Field and Ralston 2005, Pyper et al. 2001, Peterman et al. 1998, Myers at al. 1997)
# 4. Cluster analysis
# 5. DFA
# 6. PCA

#Load packages
# library(dplyr) NOTE: loading dplyr conflicts with the lag1.plot in astsa. DO NOT LOAD before using ASTSA.
# or load it and then unload it
library(cluster)
library(dplyr)
library(Hmisc) #for correlation
library(ggplot2)
library(astsa)
library(MARSS) #for DFA
setwd("~/Desktop/Github Repo/Seatrout/Data/Indices/DeltaMethod Indices")

######################################
# IMPORT INDICES
# indices produced with Delta Method for Producing Nominal Indices.R
# and select the year, assign bay name,
# scale each index (z score, 0 mean 1 unit variance)
##############################
AP<- subset(read.csv("AP_yoy_index.csv", header=T), select=c(year, index)) %>% mutate(bay= rep("AP",18)) 
  names(AP) <- c("year", "index", "est") #assign names
  AP$est <- as.factor(AP$est) #turn bay into a factor
 AP$index<-  as.numeric(scale(AP$index, scale=TRUE)) #z-score (mean= 0, var=1), as.numeric is important for full_join command below. doesnt like it when it isnt a true numeric
  # ap does not have a riv component because of so many 1 or 0 positive trips in the FIM data. The csv is available but its mostly NaNs

CH<- subset(read.csv("CH_yoy_index.csv", header=T), select=c(year, index)) %>% mutate(bay= rep("CH",27)) 
  names(CH) <- c("year", "index", "est")
  CH$est <- as.factor(CH$est) 
  CH$index <- as.numeric(scale(CH$index, scale=TRUE)) 
  
# #ch_riv<- subset(read.csv("CH_yoy_index.csv", header=T), select=c(year, index)) %>% mutate(bay= rep("CHR",27))
#   names(ch_riv) <- c("year", "index", "est")
#   ch_riv$est <- as.factor(ch_riv$est)
#   ch_riv$index <- as.numeric(scale(ch_riv$index, scale=TRUE))

CK<- subset(read.csv("CK_yoy_index.csv", header=T), select=c(year, index)) %>% mutate(bay= rep("CK",20))
  names(CK) <- c("year", "index", "est")
  CK$est <- as.factor(CK$est)
  CK$index <- as.numeric(scale(CK$index , scale=TRUE))
  # Riv is same as AP above
  
IR<- subset(read.csv("IR_yoy_index.csv", header=T), select=c(year, index)) %>% mutate(bay= rep("IR",26))
  names(IR) <- c("year", "index", "est")
  IR$est <- as.factor(IR$est)
  IR$index <- as.numeric(scale(IR$index , scale=TRUE))

TB<- subset(read.csv("TB_yoy_index.csv", header=T), select=c(year, index))%>% mutate(bay= rep("TB",27))
  names(TB) <- c("year", "index", "est")
  TB$est <- as.factor(TB$est)
  TB$index <- as.numeric(scale(TB$index , scale=TRUE))
  
# tb_riv<- subset(read.csv("index_tbriv_yoy.csv", header=T), select=c(year, Mean))%>% mutate(bay= rep("TBR",27))
#   names(tb_riv) <- c("year", "mean", "est")
#   tb_riv$est <- as.factor(tb_riv$est)
#   tb_riv$mean <- as.numeric(scale(tb_riv$mean , scale=TRUE))
  
  
JX<- subset(read.csv("JX_yoy_index.csv", header=T), select=c(year, index))%>% mutate(bay= rep("JX",15))
  names(JX) <- c("year", "index", "est")
  JX$est <- as.factor(JX$est)
  JX$index <- as.numeric(scale(JX$index , scale=TRUE))
  
All <- rbind(AP,CH,CK,IR,JX, TB)


######################################################
# PLOT TIME SERIES DATA
######################################################
# plot either the standardized or non standardized time series on the same plot
# might do this with ggplot
library(ggplot2)

recruitment <- ggplot(All, aes(x=year, y=index, group=est)) + #color=est
        geom_line(aes(linetype=est), size =.5)+ # make line types based on the different labels- this will be our workaround because in a few stps we will specify the first label (obserseved) be a blank line (therefore a scatter plot)
        #geom_line(aes(color=bay))+
      geom_point(aes(shape=est), size=2) + #, color=bay))+ # groups the points together correctly and then gives them a unique shape them differently based on the line type 
      ylab("Scaled YOY abundance index") +
      xlab("Year")+ 
  scale_y_continuous(limits=c(-2,3))+
  scale_colour_discrete(name="Location")+
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(), 
        panel.background=element_rect(colour="black", fill="white"),
        axis.title.x =element_text(colour="black"),
        axis.text.x = element_text(colour="black"),
        axis.title.y =element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        axis.text.x=element_text(colour="black"), #changing  colour of x axis
        axis.text.y=element_text(colour="black"), #changing colour of y acis
        plot.title=element_text(size=14))

############################################
#1A. LAGGED SCATTERPLOTS of each time series
#################################################
# From Notes_3, GEOS 585A, Spring 2015 handout printed from arizona (website above)
# " An attribute of the lagged scatterplot is that it can display autocorrelation regardless of the form of the
# dependence on the past values. An assumption of linear dependence is not necessary."
# Such nonlinear dependence might not be effectively summarized by other methods (such as the acf function)

#1A.1 Determine Critical level of correlation for 95% significance (alpha = 0.5) r= 0+- 2/sqrt(N)
#AP = +-2/sqrt(18)= 0.47
#CH = +-2/sqrt(27)= 0.38
#CK = +-2/sqrt(20)= 0.44
#IR = +- 2/sqrt(26) =0.39
#TB = +- 2/sqrt(27) = 0.38 ***
#Jx= +- 2/sqrt(14) = 0.53

#unload dplyr
detach("package:dplyr", unload=TRUE)
library(astsa)
lag1.plot(AP, 4, corr=TRUE)
lag1.plot(CH, 4, corr=TRUE)
lag1.plot(CK, 4, corr=TRUE)
lag1.plot(IR, 4, corr=TRUE)
lag1.plot(TB, 4, corr=TRUE) #autocorrelation at 1st and 2nd lag
lag1.plot(JX, 4, corr=TRUE)

#1B. ACF and correlograms
acf2(AP$index)
acf2(CH$index)
acf2(CK$index)
acf2(IR$index)
acf2(TB$index)
acf2(JX$index)


#little evidence of interseries correlation so can use standard significance test for pairwise correlations

#########################
# PEARSON CORRELATION
#########################
# To not end up with a bunch of "est" columns I will just make some new dfs that do not contain the est columns

ap_min <- subset(AP, select=(-est)) 
ch_min <- subset(CH, select=(-est))
ck_min <- subset(CK, select=(-est))
ir_min <- subset(IR, select=(-est))
tb_min <- subset(TB, select=(-est))
jx_min <- subset(JX, select=(-est))

# Make dataframe with just indices by year

# For some reason full_join stops working after 4 dataframes at which point it starts duplicating things so 
# I just made two and joined them together

ind1 <- full_join(ap_min, ch_min, by='year') %>% 
            full_join(.,ck_min, by='year' ) 
names(ind1) <- c("year", "AP", "CH", "CK")

ind2 <-  full_join(tb_min,ir_min, by='year' ) %>%
            full_join(.,jx_min, by='year' ) 
names(ind2) <- c("year","IR", "TB", "JX" )

ind_year <- full_join(ind1, ind2, by='year')

mat<- as.matrix(arrange(ind_year, year)) #arrange the data,frame by year and then turn into a matrix
mattest=mat[,-1] #remove the first column which is the year to just have a matrix of scaled indices


# Now pearson correlation with rcorr from hmisc
#rcorr- missing values are deleted in pairs rather than deleting all rows of x having any missing variables

corr_mat <- rcorr(mattest, type="pearson")
rho_mat <- corr_mat$r #rho values
n_mat <- corr_mat$n #number of samples used to compute correlation
P_mat <- corr_mat$P #P values

P_mat_adjust <- as.data.frame(p.adjust(P_mat)) #adjust the p-values for multiple comparisons

#unwrap each matrix
library(gdata)
rho_vec <- as.data.frame(lowerTriangle(rho_mat, diag=FALSE, byrow=FALSE))
p_vec   <- as.data.frame(lowerTriangle(P_mat, diag=FALSE, byrow=FALSE))

#produce sample number (impt for fitting the nls function and weighting by sample number in Medoid_GC_NLSfit_plot.R)
n_mat <- as.data.frame(lowerTriangle(n_mat, diag=FALSE, byrow=FALSE)) 
rho_P_vec <- cbind(rho_vec, p_vec, n_mat)
rownames <- c("AP_CH", "AP_CK", "AP_IR", "AP_TB", "AP_JX", 
              "CH_CK", "CH_IR", "CH_TB", "CH_JX",
              "CK_IR", "CK_TB", "CK_JX",
              "IR_TB", "IR_JX",
              "TB_JX")

row.names(rho_P_vec)<- rownames
data.frame(rho_P_vec)
colnames(rho_P_vec)<-  c("rho", "P", "N")
## export the rho_dataset to be used by Medoid_GC_NLSfit_plot.R when comparing great circle distances to rho
write.csv(rho_P_vec, "~/Desktop/Github Repo/Seatrout/Data/Exported R Dataframes/Rho_P_vector.csv")


#############################
# CLUSTER ANALYSIS
#############################

#using compelete linkage method and the correlation matrix as the distance matrix

library(cluster)
ac <- agnes(rho_mat)
plot(ac, ask=TRUE)

################################
# Principal Components 
#################################
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

#############################
# DFA
#############################

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


#Plotting Example- to produce plot on page 118
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



###################################################
### code chunk number 29: Cs23_plotfacloadings
###################################################
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


###################################################
### code chunk number 30: Cs24_plottrends
###################################################
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


###################################################
### code chunk number 31: Cs25_plotbestfits
###################################################
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











#Example from the Text: 


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
