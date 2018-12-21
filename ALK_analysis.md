R Notebook for Age Length Key Analysis
================

### **Main Objectives of this script**

1.  Imports otolith data. Summarizes data, mean age, mean length, sd and se and total sample number.
2.  Makes bay-specific observed ALK, calculates some summary statistics.
3.  Makes bay specific smoothed (modeled) ALK with multinomial modeling methods in Ogle (87-).
4.  Likelihood ratio testing to do among group statistical comparisons - Ogle (102-103)
5.  Plots observed and smoothed ALK for each bay.
6.  Calculates proportional age distribution
7.  Determines mean length -at-age (otolith database) and produces plots
8.  One way anova to determine if the mean lengths (and mean ages) of males and females are significantly different among estuaries
9.  T test to determine whether there is a significant difference between male and female age for each estuary
10. T test to determine whether there is a significant difference in male and female length for each estuary. 11.Calculates proportional age distribution for ADULTS to be used in Delta\_Method script (different from \#6)
11. Calculates weight at age for ADULTS to be used in Delta\_Method script

### *Start Here*

Load packages and set my working directory.

``` r
library(FSA)
library(magrittr)
library(nnet)
library(plotrix)
library(haven)
library(ggplot2)
library(scales)
library(fishmethods) #masking select from dplyr 
library(dplyr)
library(gridExtra)
```

#### 1. LOAD DATA & DO SOME WRANGLING

1.  load the csv file
2.  subset by which bay I want
3.  subset the FIM program
4.  make sure I have just the "aged sample"
5.  turn mm to cm
6.  select just a few variables to make it more manageable
7.  then drop the remaining bay levels still sticking around (droplevels)
8.  turn tl from mm to cm
9.  create length categories with FSA package

*Note: Because each dataset has its peculiarities I did this for each dataset and not within a loop or function*

``` r
setwd("U:/PhD_projectfiles/Raw_Data/Age_Length_Data")
Agelength_TB<- droplevels(subset(as.data.frame(read.csv("ALK_Bay_and_weight.csv", header=T)), bay=="TB" & tl>14 & final_age >0 & program== 'FIM', select=c(sex,SpecimenNumber, bay, tl,sl, final_age, wt_total, date, program))) %>% mutate(tl=tl/10, sl=sl/10, lcat2 =lencat(tl, w=2)) %>% filter(sex %in% c("m", "M", "F")) #, as.fact=TRUE))- can include this when determing ALK below but the smoothed ALK needs to be the nonfactored version of the length categorization variable. 
Agelength_TB$sex[which(Agelength_TB$sex == "m")] = "M"
Agelength_TB$sex <- droplevels(Agelength_TB$sex)
```

Change date format into a factor so that I can do summary statistics by year later on

``` r
Agelength_TB$date=as.character(Agelength_TB$date)
Agelength_TB$DateNew = as.POSIXct(strptime(Agelength_TB$date, format="%d-%B-%y", tz=""))  #B is the selection for when month is spelled out
Agelength_TB = mutate(Agelength_TB, year = strftime(DateNew, format="%Y")) %>% dplyr::select(c(-date, -DateNew))
Agelength_TB$year = as.factor(Agelength_TB$year)
```

Now do the same thing for other estuaries. R code will not be displayed because it was fairly similar.

#### 2. Basic data summarization

``` r
#total sample number of FIM data
All= rbind(Agelength_AP, Agelength_CK, Agelength_TB, Agelength_CH, Agelength_JX, Agelength_IR) #omits na
All <- na.omit(All)
```

    ## [1] 10.6

    ## [1] 75.7

``` r
#summarize the entire set and group by sex and then bay
Combined_sum <- dplyr::summarize(group_by(All, sex,bay), mean_tl=mean(tl), sd_tl=sd(tl), se_tl= sd_tl/(sqrt(length(final_age))),min_tl=min(tl), max_tl=max(tl), mean_age=mean(final_age), median_age=median(final_age),sd_age= sd(final_age), se_age=sd_age/sqrt(length(final_age)), min_age= min(final_age), max_age=max(final_age))

Combined_sum
```

    ## # A tibble: 12 x 13
    ## # Groups:   sex [?]
    ##    sex   bay   mean_tl sd_tl se_tl min_tl max_tl mean_age median_age
    ##    <fct> <fct>   <dbl> <dbl> <dbl>  <dbl>  <dbl>    <dbl>      <dbl>
    ##  1 F     AP       41.6  9.16 0.328   15.7   71.0     2.58       2.00
    ##  2 F     CK       39.0  7.80 0.302   16.4   64.5     2.04       2.00
    ##  3 F     TB       40.0  7.85 0.232   16.3   70.4     2.58       2.00
    ##  4 F     CH       40.5  7.08 0.276   20.8   64.7     2.59       2.00
    ##  5 F     JX       35.7  7.64 0.343   15.2   69.8     1.83       2.00
    ##  6 F     IR       39.7 10.6  0.236   10.6   75.7     2.50       2.00
    ##  7 M     AP       35.1  5.26 0.324   15.2   50.5     2.66       2.00
    ##  8 M     CK       34.4  5.55 0.397   20.5   52.4     2.07       2.00
    ##  9 M     TB       35.8  5.25 0.201   16.0   61.9     2.86       3.00
    ## 10 M     CH       35.3  5.09 0.255   22.1   62.9     2.66       3.00
    ## 11 M     JX       32.6  4.61 0.239   18.0   51.9     2.01       2.00
    ## 12 M     IR       34.1  7.50 0.290   15.1   60.3     2.86       3.00
    ## # ... with 4 more variables: sd_age <dbl>, se_age <dbl>, min_age <dbl>,
    ## #   max_age <dbl>

Summarize each to get sex ratios etc.

``` r
#an example with Tampa Bay data
TB_sum <- summarise(Agelength_TB %>% filter(sex %in% c("M", "F")),mean_age=mean(final_age), median_age=median(final_age),sd_age= sd(final_age), se_age=sd_age/sqrt(length(final_age)), mean_tl = mean(tl), median_tl = median(tl), min_tl=min(tl), max_tl=max(tl), sd_tl= sd(tl), se_tl= sd_tl/(sqrt(length(final_age))), prop_F=nrow(subset(Agelength_TB, sex=="F"))/length(final_age),prop_M=nrow(subset(Agelength_TB, sex=="M"))/length(final_age), N_F=nrow(subset(Agelength_TB, sex=="F")), N_M=nrow(subset(Agelength_TB, sex=="M")), totN =N_F+N_M)
```

Now for all other estuaries but code not displayed.

Add them together and make sure rownames are added

``` r
All_sum <- rbind(AP_sum, CK_sum, TB_sum, CH_sum, JX_sum, IR_sum) 
rownames(All_sum) <- c("AP", "CK", "TB", "CH", "JX", "IR")
```

#### TWO WAY ANOVAS

To test for differences in age and length by bay and sex.

``` r
All_MF <- droplevels(subset(All, sex %in% c("M", "F")))

#differences in final age and length by estuary and sex
aov <- aov(final_age ~ sex + bay, data= All_MF)
summary(aov)
```

    ##               Df Sum Sq Mean Sq F value   Pr(>F)    
    ## sex            1     70   69.72   42.54 7.32e-11 ***
    ## bay            5    602  120.38   73.46  < 2e-16 ***
    ## Residuals   8325  13642    1.64                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Tukey's test for pairwise comparison
TukeyHSD(aov)
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = final_age ~ sex + bay, data = All_MF)
    ## 
    ## $sex
    ##        diff       lwr       upr p adj
    ## M-F 0.19776 0.1383262 0.2571938     0
    ## 
    ## $bay
    ##                diff         lwr          upr     p adj
    ## CK-AP -0.5491323628 -0.71694930 -0.381315428 0.0000000
    ## TB-AP  0.0626960378 -0.07875504  0.204147118 0.8049849
    ## CH-AP -0.0095706787 -0.16882131  0.149679951 0.9999797
    ## JX-AP -0.7274866204 -0.89509184 -0.559881396 0.0000000
    ## IR-AP -0.0087666919 -0.14199094  0.124457557 0.9999681
    ## TB-CK  0.6118284006  0.46127858  0.762378226 0.0000000
    ## CH-CK  0.5395616841  0.37217708  0.706946286 0.0000000
    ## JX-CK -0.1783542576 -0.35370634 -0.003002175 0.0435139
    ## IR-CK  0.5403656709  0.39751770  0.683213640 0.0000000
    ## CH-TB -0.0722667165 -0.21320461  0.068671176 0.6889091
    ## JX-TB -0.7901826582 -0.94049646 -0.639868861 0.0000000
    ## IR-TB -0.0714627297 -0.18214889  0.039223434 0.4397216
    ## JX-CH -0.7179159417 -0.88508829 -0.550743598 0.0000000
    ## IR-CH  0.0008039868 -0.13187526  0.133483231 1.0000000
    ## IR-JX  0.7187199285  0.57612074  0.861319122 0.0000000

``` r
aov2 <- aov(tl ~ bay + sex, data= All_MF)
summary(aov2)
```

    ##               Df Sum Sq Mean Sq F value Pr(>F)    
    ## bay            5  16439    3288   49.74 <2e-16 ***
    ## sex            1  41446   41446  627.04 <2e-16 ***
    ## Residuals   8325 550266      66                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
TukeyHSD(aov2)
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = tl ~ bay + sex, data = All_MF)
    ## 
    ## $bay
    ##             diff        lwr        upr     p adj
    ## CK-AP -2.0057352 -3.0715401 -0.9399302 0.0000012
    ## TB-AP -1.5358818 -2.4342373 -0.6375263 0.0000165
    ## CH-AP -1.4232857 -2.4346860 -0.4118853 0.0008630
    ## JX-AP -5.6318568 -6.6963172 -4.5673965 0.0000000
    ## IR-AP -1.7030310 -2.5491379 -0.8569241 0.0000001
    ## TB-CK  0.4698534 -0.4862882  1.4259949 0.7266964
    ## CH-CK  0.5824495 -0.4806097  1.6455087 0.6239606
    ## JX-CK -3.6261217 -4.7397823 -2.5124610 0.0000000
    ## IR-CK  0.3027042 -0.6045229  1.2099313 0.9330243
    ## CH-TB  0.1125961 -0.7825001  1.0076923 0.9992268
    ## JX-TB -4.0959750 -5.0506176 -3.1413325 0.0000000
    ## IR-TB -0.1671492 -0.8701167  0.5358184 0.9844179
    ## JX-CH -4.2085712 -5.2702823 -3.1468600 0.0000000
    ## IR-CH -0.2797453 -1.1223909  0.5629002 0.9343789
    ## IR-JX  3.9288258  3.0231787  4.8344730 0.0000000
    ## 
    ## $sex
    ##          diff       lwr       upr p adj
    ## M-F -4.763083 -5.140547 -4.385619     0

``` r
#differences in final age and length by sex within each estuary
summary(aov(final_age ~  sex, data= AL_TB))
```

    ##               Df Sum Sq Mean Sq F value   Pr(>F)    
    ## sex            1     32   32.15   16.62 4.75e-05 ***
    ## Residuals   1858   3593    1.93                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(aov(final_age ~  sex, data= AL_AP))
```

    ##               Df Sum Sq Mean Sq F value Pr(>F)  
    ## sex            1      5   5.037   2.835 0.0925 .
    ## Residuals   1245   2212   1.777                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(aov(final_age ~  sex, data= AL_CK))
```

    ##              Df Sum Sq Mean Sq F value Pr(>F)
    ## sex           1    0.1  0.1025   0.081  0.777
    ## Residuals   818 1040.9  1.2725

``` r
summary(aov(final_age ~  sex, data= AL_CH))
```

    ##               Df Sum Sq Mean Sq F value Pr(>F)
    ## sex            1    0.7  0.7301   0.479  0.489
    ## Residuals   1144 1744.4  1.5248

``` r
summary(aov(final_age ~  sex, data= AL_JX))
```

    ##              Df Sum Sq Mean Sq F value  Pr(>F)   
    ## sex           1    7.9   7.878   7.952 0.00492 **
    ## Residuals   819  811.5   0.991                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(aov(final_age ~  sex, data= AL_IR))
```

    ##               Df Sum Sq Mean Sq F value  Pr(>F)    
    ## sex            1     95   95.27   58.78 2.5e-14 ***
    ## Residuals   2494   4042    1.62                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(aov(tl ~  sex, data= AL_TB))
```

    ##               Df Sum Sq Mean Sq F value Pr(>F)    
    ## sex            1   7986    7986   168.3 <2e-16 ***
    ## Residuals   1858  88159      47                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(aov(tl  ~  sex, data= AL_AP))
```

    ##               Df Sum Sq Mean Sq F value Pr(>F)    
    ## sex            1   9085    9085   153.5 <2e-16 ***
    ## Residuals   1245  73669      59                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(aov(tl  ~  sex, data= AL_CK))
```

    ##              Df Sum Sq Mean Sq F value Pr(>F)    
    ## sex           1   3114  3114.1   72.98 <2e-16 ***
    ## Residuals   818  34905    42.7                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(aov(tl  ~  sex, data= AL_CH))
```

    ##               Df Sum Sq Mean Sq F value Pr(>F)    
    ## sex            1   8584    8584   213.7 <2e-16 ***
    ## Residuals   1144  45954      40                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(aov(tl  ~  sex, data= AL_JX))
```

    ##              Df Sum Sq Mean Sq F value   Pr(>F)    
    ## sex           1   2235  2234.6   64.41 3.49e-15 ***
    ## Residuals   819  28413    34.7                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(aov(tl  ~  sex, data= AL_IR))
```

    ##               Df Sum Sq Mean Sq F value Pr(>F)    
    ## sex            1  12253   12253     153 <2e-16 ***
    ## Residuals   2494 199671      80                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

CHI SQ TEST To test differencs in age distribution between bays
===============================================================

non-parameteric like anova is below so it doesnt assume normality of the age distribution
can do a test of variance to determine whether can use anova
anova assumes homogenous variances
test for variances
H0= ratio of variances is equal to 1
Ha = ratio of variances is not equal to 1

``` r
var.test(Agelength_TB$final_age, Agelength_AP$final_age) #equal variances
```

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  Agelength_TB$final_age and Agelength_AP$final_age
    ## F = 1.0974, num df = 1873, denom df = 1254, p-value = 0.07313
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##  0.991309 1.213639
    ## sample estimates:
    ## ratio of variances 
    ##           1.097419

``` r
var.test(Agelength_TB$final_age, Agelength_CK$final_age)
```

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  Agelength_TB$final_age and Agelength_CK$final_age
    ## F = 1.5524, num df = 1873, denom df = 863, p-value = 2.34e-13
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##  1.383334 1.737926
    ## sample estimates:
    ## ratio of variances 
    ##           1.552416

``` r
var.test(Agelength_TB$final_age, Agelength_CH$final_age)
```

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  Agelength_TB$final_age and Agelength_CH$final_age
    ## F = 1.2721, num df = 1873, denom df = 1168, p-value = 6.399e-06
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##  1.146447 1.409787
    ## sample estimates:
    ## ratio of variances 
    ##           1.272115

``` r
var.test(Agelength_TB$final_age, Agelength_IR$final_age) 
```

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  Agelength_TB$final_age and Agelength_IR$final_age
    ## F = 1.1415, num df = 1873, denom df = 2765, p-value = 0.001685
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##  1.050929 1.240627
    ## sample estimates:
    ## ratio of variances 
    ##           1.141462

``` r
var.test(Agelength_TB$final_age, Agelength_JX$final_age)
```

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  Agelength_TB$final_age and Agelength_JX$final_age
    ## F = 1.9811, num df = 1873, denom df = 867, p-value < 2.2e-16
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##  1.765692 2.217492
    ## sample estimates:
    ## ratio of variances 
    ##           1.981129

``` r
var.test(Agelength_AP$final_age, Agelength_CK$final_age)
```

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  Agelength_AP$final_age and Agelength_CK$final_age
    ## F = 1.4146, num df = 1254, denom df = 863, p-value = 4.533e-08
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##  1.250425 1.598094
    ## sample estimates:
    ## ratio of variances 
    ##           1.414607

``` r
var.test(Agelength_AP$final_age, Agelength_CH$final_age) 
```

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  Agelength_AP$final_age and Agelength_CH$final_age
    ## F = 1.1592, num df = 1254, denom df = 1168, p-value = 0.01035
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##  1.035452 1.297415
    ## sample estimates:
    ## ratio of variances 
    ##           1.159189

``` r
var.test(Agelength_AP$final_age, Agelength_IR$final_age)
```

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  Agelength_AP$final_age and Agelength_IR$final_age
    ## F = 1.0401, num df = 1254, denom df = 2765, p-value = 0.4093
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##  0.9472387 1.1440803
    ## sample estimates:
    ## ratio of variances 
    ##           1.040134

``` r
var.test(Agelength_AP$final_age, Agelength_JX$final_age) #equal variances
```

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  Agelength_AP$final_age and Agelength_JX$final_age
    ## F = 1.8053, num df = 1254, denom df = 867, p-value < 2.2e-16
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##  1.596025 2.039100
    ## sample estimates:
    ## ratio of variances 
    ##           1.805263

``` r
var.test(Agelength_CK$final_age, Agelength_CH$final_age)
```

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  Agelength_CK$final_age and Agelength_CH$final_age
    ## F = 0.81944, num df = 863, denom df = 1168, p-value = 0.00184
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##  0.7239490 0.9286267
    ## sample estimates:
    ## ratio of variances 
    ##          0.8194424

``` r
var.test(Agelength_CK$final_age, Agelength_IR$final_age)
```

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  Agelength_CK$final_age and Agelength_IR$final_age
    ## F = 0.73528, num df = 863, denom df = 2765, p-value = 6.204e-08
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##  0.6609450 0.8205232
    ## sample estimates:
    ## ratio of variances 
    ##          0.7352814

``` r
var.test(Agelength_CK$final_age, Agelength_JX$final_age)
```

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  Agelength_CK$final_age and Agelength_JX$final_age
    ## F = 1.2762, num df = 863, denom df = 867, p-value = 0.0003424
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##  1.116831 1.458247
    ## sample estimates:
    ## ratio of variances 
    ##           1.276159

``` r
var.test(Agelength_CH$final_age, Agelength_IR$final_age)
```

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  Agelength_CH$final_age and Agelength_IR$final_age
    ## F = 0.89729, num df = 1168, denom df = 2765, p-value = 0.02982
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##  0.8153206 0.9894154
    ## sample estimates:
    ## ratio of variances 
    ##          0.8972947

``` r
var.test(Agelength_CH$final_age, Agelength_JX$final_age)
```

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  Agelength_CH$final_age and Agelength_JX$final_age
    ## F = 1.5574, num df = 1168, denom df = 867, p-value = 5.956e-12
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##  1.374485 1.762501
    ## sample estimates:
    ## ratio of variances 
    ##            1.55735

``` r
var.test(Agelength_IR$final_age, Agelength_JX$final_age)
```

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  Agelength_IR$final_age and Agelength_JX$final_age
    ## F = 1.7356, num df = 2765, denom df = 867, p-value < 2.2e-16
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##  1.555611 1.930462
    ## sample estimates:
    ## ratio of variances 
    ##           1.735606

``` r
#mostly all variances are unequal so use a chi squared test

#truncate to less than 8 to satisfy requirements of the test
ALL_to7 <- All %>% subset(final_age<8)
(age_freq <- xtabs(~bay+final_age, data=ALL_to7))
```

    ##     final_age
    ## bay    1   2   3   4   5   6   7
    ##   AP 243 327 225 157  58  19  12
    ##   CK 347 264 154  71  23   3   2
    ##   TB 389 540 473 235 119  51  23
    ##   CH 207 329 282 163  55  16   4
    ##   JX 338 360 113  37  12   6   1
    ##   IR 629 737 702 388 136  52  16

``` r
chisq.test(age_freq[1:2,]) 
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  age_freq[1:2, ]
    ## X-squared = 89.013, df = 6, p-value < 2.2e-16

``` r
chisq.test(age_freq[1:3,]) 
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  age_freq[1:3, ]
    ## X-squared = 162.27, df = 12, p-value < 2.2e-16

``` r
chisq.test(age_freq[1:4,]) 
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  age_freq[1:4, ]
    ## X-squared = 198.05, df = 18, p-value < 2.2e-16

``` r
chisq.test(age_freq[1:5,]) 
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  age_freq[1:5, ]
    ## X-squared = 404.56, df = 24, p-value < 2.2e-16

``` r
chisq.test(age_freq[1:6,]) 
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  age_freq[1:6, ]
    ## X-squared = 449.2, df = 30, p-value < 2.2e-16

``` r
chisq.test(age_freq[2:3,]) 
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  age_freq[2:3, ]
    ## X-squared = 144.78, df = 6, p-value < 2.2e-16

``` r
chisq.test(age_freq[2:4,]) 
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  age_freq[2:4, ]
    ## X-squared = 188.52, df = 12, p-value < 2.2e-16

``` r
chisq.test(age_freq[2:5,]) 
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  age_freq[2:5, ]
    ## X-squared = 384.93, df = 18, p-value < 2.2e-16

``` r
chisq.test(age_freq[2:6,]) 
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  age_freq[2:6, ]
    ## X-squared = 435.87, df = 24, p-value < 2.2e-16

``` r
chisq.test(age_freq[3:4,])
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  age_freq[3:4, ]
    ## X-squared = 16.99, df = 6, p-value = 0.00932

``` r
chisq.test(age_freq[3:5,]) 
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  age_freq[3:5, ]
    ## X-squared = 282.68, df = 12, p-value < 2.2e-16

``` r
chisq.test(age_freq[3:6,]) 
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  age_freq[3:6, ]
    ## X-squared = 311.04, df = 18, p-value < 2.2e-16

``` r
chisq.test(age_freq[4:5,]) 
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  age_freq[4:5, ]
    ## X-squared = 201.89, df = 6, p-value < 2.2e-16

``` r
chisq.test(age_freq[4:6,]) 
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  age_freq[4:6, ]
    ## X-squared = 266.71, df = 12, p-value < 2.2e-16

``` r
chisq.test(age_freq[5:6,]) 
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  age_freq[5:6, ]
    ## X-squared = 235.72, df = 6, p-value < 2.2e-16

``` r
#Adjusting p-values for the multople comparisons
ps_age<- c(chisq.test(age_freq[1:2,])$p.value,
           chisq.test(age_freq[1:3,])$p.value,
           chisq.test(age_freq[1:4,])$p.value,
           chisq.test(age_freq[1:5,])$p.value,
           chisq.test(age_freq[1:6,])$p.value,
           chisq.test(age_freq[2:3,])$p.value,
           chisq.test(age_freq[2:4,])$p.value,
           chisq.test(age_freq[2:5,])$p.value,
           chisq.test(age_freq[2:6,])$p.value,
           chisq.test(age_freq[3:4,])$p.value,
           chisq.test(age_freq[3:5,])$p.value,
           chisq.test(age_freq[3:6,])$p.value,
           chisq.test(age_freq[4:5,])$p.value,
           chisq.test(age_freq[4:6,])$p.value,
           chisq.test(age_freq[5:6,])$p.value)

pdf <- as.data.frame(p.adjust(ps_age))

#all significantly different excpet for TB to AP and CH and IR
```

### CHI SQ TEST To test differencs in length distribution between bays

test for equal variances

``` r
var.test(Agelength_TB$tl, Agelength_AP$tl)
```

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  Agelength_TB$tl and Agelength_AP$tl
    ## F = 0.77476, num df = 1873, denom df = 1254, p-value = 6.243e-07
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##  0.6998504 0.8568126
    ## sample estimates:
    ## ratio of variances 
    ##          0.7747624

``` r
var.test(Agelength_TB$tl, Agelength_CK$tl)
```

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  Agelength_TB$tl and Agelength_CK$tl
    ## F = 0.93581, num df = 1873, denom df = 863, p-value = 0.2494
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##  0.8338865 1.0476374
    ## sample estimates:
    ## ratio of variances 
    ##          0.9358103

``` r
var.test(Agelength_TB$tl, Agelength_CH$tl)
```

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  Agelength_TB$tl and Agelength_CH$tl
    ## F = 1.07, num df = 1873, denom df = 1168, p-value = 0.2023
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##  0.9642759 1.1857711
    ## sample estimates:
    ## ratio of variances 
    ##           1.069975

``` r
var.test(Agelength_TB$tl, Agelength_IR$tl) #equal variances
```

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  Agelength_TB$tl and Agelength_IR$tl
    ## F = 0.52453, num df = 1873, denom df = 2765, p-value < 2.2e-16
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##  0.4829248 0.5700948
    ## sample estimates:
    ## ratio of variances 
    ##          0.5245267

``` r
var.test(Agelength_TB$tl, Agelength_JX$tl)
```

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  Agelength_TB$tl and Agelength_JX$tl
    ## F = 1.2036, num df = 1873, denom df = 867, p-value = 0.001665
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##  1.072683 1.347158
    ## sample estimates:
    ## ratio of variances 
    ##           1.203564

``` r
var.test(Agelength_AP$tl, Agelength_CK$tl)
```

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  Agelength_AP$tl and Agelength_CK$tl
    ## F = 1.2079, num df = 1254, denom df = 863, p-value = 0.002743
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##  1.067680 1.364538
    ## sample estimates:
    ## ratio of variances 
    ##           1.207867

``` r
var.test(Agelength_AP$tl, Agelength_CH$tl) #equal variances
```

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  Agelength_AP$tl and Agelength_CH$tl
    ## F = 1.381, num df = 1254, denom df = 1168, p-value = 2.282e-08
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##  1.233619 1.545717
    ## sample estimates:
    ## ratio of variances 
    ##           1.381037

``` r
var.test(Agelength_AP$tl, Agelength_IR$tl)
```

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  Agelength_AP$tl and Agelength_IR$tl
    ## F = 0.67702, num df = 1254, denom df = 2765, p-value = 3.248e-15
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##  0.6165513 0.7446742
    ## sample estimates:
    ## ratio of variances 
    ##          0.6770162

``` r
var.test(Agelength_AP$tl, Agelength_JX$tl)
```

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  Agelength_AP$tl and Agelength_JX$tl
    ## F = 1.5535, num df = 1254, denom df = 867, p-value = 4.289e-12
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##  1.373410 1.754684
    ## sample estimates:
    ## ratio of variances 
    ##           1.553463

``` r
var.test(Agelength_CK$tl, Agelength_CH$tl)
```

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  Agelength_CK$tl and Agelength_CH$tl
    ## F = 1.1434, num df = 863, denom df = 1168, p-value = 0.03405
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##  1.010126 1.295713
    ## sample estimates:
    ## ratio of variances 
    ##           1.143368

``` r
var.test(Agelength_CK$tl, Agelength_IR$tl)
```

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  Agelength_CK$tl and Agelength_IR$tl
    ## F = 0.56051, num df = 863, denom df = 2765, p-value < 2.2e-16
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##  0.5038387 0.6254853
    ## sample estimates:
    ## ratio of variances 
    ##          0.5605054

``` r
var.test(Agelength_CK$tl, Agelength_JX$tl)
```

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  Agelength_CK$tl and Agelength_JX$tl
    ## F = 1.2861, num df = 863, denom df = 867, p-value = 0.0002201
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##  1.125549 1.469629
    ## sample estimates:
    ## ratio of variances 
    ##            1.28612

``` r
var.test(Agelength_CH$tl, Agelength_IR$tl)
```

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  Agelength_CH$tl and Agelength_IR$tl
    ## F = 0.49022, num df = 1168, denom df = 2765, p-value < 2.2e-16
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##  0.4454378 0.5405518
    ## sample estimates:
    ## ratio of variances 
    ##          0.4902231

``` r
var.test(Agelength_CH$tl, Agelength_JX$tl)
```

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  Agelength_CH$tl and Agelength_JX$tl
    ## F = 1.1249, num df = 1168, denom df = 867, p-value = 0.06484
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##  0.9927713 1.2730302
    ## sample estimates:
    ## ratio of variances 
    ##           1.124852

``` r
var.test(Agelength_IR$tl, Agelength_JX$tl)
```

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  Agelength_IR$tl and Agelength_JX$tl
    ## F = 2.2946, num df = 2765, denom df = 867, p-value < 2.2e-16
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##  2.056609 2.552183
    ## sample estimates:
    ## ratio of variances 
    ##           2.294572

``` r
#unequal variances so will use chi.squared again

#truncate lengths to satisfy requirements
All_len_freq_30to58 <- All %>% subset(lcat2>=32 & lcat2<= 58)
(len_freq <- xtabs(~bay+lcat2, data=All_len_freq_30to58))
```

    ##     lcat2
    ## bay   32  34  36  38  40  42  44  46  48  50  52  54  56  58
    ##   AP  99 111 120  84  95  75  62  44  37  34  22  36  16  12
    ##   CK  84  94 122  95  75  59  44  33  31  22  15   7   4   4
    ##   TB 241 279 242 157 138 142 103  61  55  35  24  24  18  16
    ##   CH 111 137 172 118  85  63  59  41  41  27  16  11   9   6
    ##   JX 144 142  79  58  28  27  21  12   8   6   7   3   3   4
    ##   IR 201 215 265 245 198 129 117  91  82  80  65  46  48  47

``` r
chisq.test(len_freq[1:2,]) 
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  len_freq[1:2, ]
    ## X-squared = 31.487, df = 13, p-value = 0.00286

``` r
chisq.test(len_freq[1:3,]) 
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  len_freq[1:3, ]
    ## X-squared = 76.934, df = 26, p-value = 6.171e-07

``` r
chisq.test(len_freq[1:4,]) 
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  len_freq[1:4, ]
    ## X-squared = 97.572, df = 39, p-value = 6.343e-07

``` r
chisq.test(len_freq[1:5,]) 
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  len_freq[1:5, ]
    ## X-squared = 241.7, df = 52, p-value < 2.2e-16

``` r
chisq.test(len_freq[1:6,]) 
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  len_freq[1:6, ]
    ## X-squared = 366.66, df = 65, p-value < 2.2e-16

``` r
chisq.test(len_freq[2:3,]) 
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  len_freq[2:3, ]
    ## X-squared = 26.538, df = 13, p-value = 0.01438

``` r
chisq.test(len_freq[2:4,]) 
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  len_freq[2:4, ]
    ## X-squared = 39.529, df = 26, p-value = 0.04337

``` r
chisq.test(len_freq[2:5,]) 
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  len_freq[2:5, ]
    ## X-squared = 158.14, df = 39, p-value = 2.781e-16

``` r
chisq.test(len_freq[2:6,]) 
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  len_freq[2:6, ]
    ## X-squared = 322.9, df = 52, p-value < 2.2e-16

``` r
chisq.test(len_freq[3:4,]) 
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  len_freq[3:4, ]
    ## X-squared = 23.467, df = 13, p-value = 0.0364

``` r
chisq.test(len_freq[3:5,]) 
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  len_freq[3:5, ]
    ## X-squared = 125.33, df = 26, p-value = 5.735e-15

``` r
chisq.test(len_freq[3:6,]) 
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  len_freq[3:6, ]
    ## X-squared = 296.02, df = 39, p-value < 2.2e-16

``` r
chisq.test(len_freq[4:5,]) 
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  len_freq[4:5, ]
    ## X-squared = 102.49, df = 13, p-value = 5.455e-16

``` r
chisq.test(len_freq[4:6,]) 
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  len_freq[4:6, ]
    ## X-squared = 251, df = 26, p-value < 2.2e-16

``` r
chisq.test(len_freq[5:6,]) 
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  len_freq[5:6, ]
    ## X-squared = 210.68, df = 13, p-value < 2.2e-16

``` r
ps_length <- c(chisq.test(len_freq[1:2,])$p.value, 
chisq.test(len_freq[1:3,])$p.value, 
chisq.test(len_freq[1:4,])$p.value, 
chisq.test(len_freq[1:5,])$p.value, 
chisq.test(len_freq[1:6,])$p.value, 
chisq.test(len_freq[2:3,])$p.value, 
chisq.test(len_freq[2:4,])$p.value, 
chisq.test(len_freq[2:5,])$p.value, 
chisq.test(len_freq[2:6,])$p.value, 
chisq.test(len_freq[3:4,])$p.value, 
chisq.test(len_freq[3:5,])$p.value, 
chisq.test(len_freq[3:6,])$p.value, 
chisq.test(len_freq[4:5,])$p.value, 
chisq.test(len_freq[4:6,])$p.value, 
chisq.test(len_freq[5:6,])$p.value) 

plength <- as.data.frame(p.adjust(ps_length))
```

### MAKE ALKS

Make table with observed total numbers at length by age

``` r
(rawfreq_TB <- xtabs(~lcat2+final_age, data=Agelength_TB)) 
```

    ##      final_age
    ## lcat2   1   2   3   4   5   6   7   8   9
    ##    14   1   0   0   0   0   0   0   0   0
    ##    16   4   0   0   0   0   0   0   0   0
    ##    18   1   0   0   0   0   0   0   0   0
    ##    20   2   0   0   0   0   0   0   0   0
    ##    22   6   0   0   0   0   0   0   0   0
    ##    24  12   1   0   0   0   0   0   0   0
    ##    26  17   1   1   0   0   0   0   0   0
    ##    28  41  21   9   0   0   0   0   0   0
    ##    30  62  67  33   4   1   0   0   0   0
    ##    32  80  78  63  18   2   0   0   0   0
    ##    34  76 105  72  29   3   1   0   0   0
    ##    36  61  83  68  23   7   2   1   0   0
    ##    38  18  61  51  22  11   2   0   0   0
    ##    40   6  38  46  33   9   5   3   0   0
    ##    42   4  47  48  31  10   2   1   1   0
    ##    44   1  23  33  25  13   9   1   0   0
    ##    46   0  16  22  12   8   3   1   1   0
    ##    48   0   4  15  16  15   4   1   0   0
    ##    50   1   0  11   9  10   7   0   0   0
    ##    52   0   1   5   7   7   2   2   0   0
    ##    54   0   0   5   8   8   1   2   0   0
    ##    56   0   0   2   2   6   4   4   0   0
    ##    58   0   0   2   4   7   3   0   1   0
    ##    60   0   0   0   0   2   6   3   0   0
    ##    62   0   0   0   2   3   1   2   1   0
    ##    66   0   0   0   0   0   1   1   0   1
    ##    68   0   0   0   0   1   0   0   0   0
    ##    70   0   0   0   0   0   0   1   0   0

``` r
#rawfreq_TB_test_df <- as.data.frame(as.matrix(xtabs(~lcat2+final_age, data=TB_test)))
# there appears to be a fish that was assigned an age of 3 but is in the 0-2 length category. Going to remove this because its probably a typo. Specified in above subsetting step as tl>20mm =(2cm).   
colSums(rawfreq_TB) #number of lengths obs per age
```

    ##   1   2   3   4   5   6   7   8   9 
    ## 393 546 486 245 123  53  23   4   1

``` r
rowSums(rawfreq_TB) #number age obs per length
```

    ##  14  16  18  20  22  24  26  28  30  32  34  36  38  40  42  44  46  48 
    ##   1   4   1   2   6  13  19  71 167 241 286 245 165 140 144 105  63  55 
    ##  50  52  54  56  58  60  62  66  68  70 
    ##  38  24  24  18  17  11   9   3   1   1

Now do for other estuaries

    ##      final_age
    ## lcat2  1  2  3  4  5  6  7  8  9 10
    ##    14  3  0  0  0  0  0  0  0  0  0
    ##    16  4  0  0  0  0  0  0  0  0  0
    ##    18  1  0  0  0  0  0  0  0  0  0
    ##    20  4  0  0  0  0  0  0  0  0  0
    ##    22  3  0  0  0  0  0  0  0  0  0
    ##    24 10  0  0  0  0  0  0  0  0  0
    ##    26 12  0  0  1  0  0  0  0  0  0
    ##    28 41 10  2  1  0  0  0  0  0  0
    ##    30 47 24  4  0  0  0  0  0  0  0
    ##    32 44 41 10  3  1  0  0  0  0  0
    ##    34 30 53 26  5  0  1  0  0  0  0
    ##    36 38 68 31  9  2  0  0  0  0  0
    ##    38 22 70 28  8  3  0  0  0  0  0
    ##    40  5 69 31 18  7  1  0  0  1  0
    ##    42  3 53 36 15  7  0  0  0  0  0
    ##    44  1 30 28 15  1  0  0  0  0  0
    ##    46  0 13 20 15  4  3  0  0  0  0
    ##    48  0  9 19 17  3  0  0  0  0  0
    ##    50  0  5 14 14  7  2  1  0  0  0
    ##    52  0  0  9 11  5  0  2  0  0  0
    ##    54  0  0 10 15  8  5  0  0  0  0
    ##    56  0  1  3  8  5  1  0  0  0  0
    ##    58  0  0  0  5  2  4  1  2  0  0
    ##    60  0  0  0  4  5  3  3  0  0  0
    ##    62  0  0  0  3  1  1  0  0  0  0
    ##    64  0  0  0  2  1  1  4  0  0  0
    ##    66  0  0  0  0  0  0  1  0  0  0
    ##    70  0  0  0  0  0  0  1  0  0  1

    ##  14  16  18  20  22  24  26  28  30  32  34  36  38  40  42  44  46  48 
    ##   3   4   1   4   3  10  13  54  75  99 115 148 131 132 114  75  55  48 
    ##  50  52  54  56  58  60  62  64  66  70 
    ##  43  27  38  18  14  15   5   8   1   2

    ##   1   2   3   4   5   6   7   8   9  10 
    ## 268 446 271 169  62  22  13   2   1   1

    ##      final_age
    ## lcat2  1  2  3  4  5  6  7
    ##    16  4  0  0  0  0  0  0
    ##    18  4  0  0  0  0  0  0
    ##    20  7  2  0  0  0  0  0
    ##    22 12  0  0  0  0  0  0
    ##    24 14  1  0  0  0  0  0
    ##    26 11  1  0  0  0  0  0
    ##    28 31  6  0  0  0  0  0
    ##    30 55 17  1  0  0  0  0
    ##    32 56 23  5  0  0  0  0
    ##    34 51 31 12  0  0  0  0
    ##    36 51 53 15  2  1  0  0
    ##    38 28 42 20  5  0  0  0
    ##    40 14 37 18  6  0  0  0
    ##    42  7 24 16 11  1  0  0
    ##    44  2 10 22  9  1  0  0
    ##    46  0  6 18  5  3  0  1
    ##    48  0  6 12  9  4  0  0
    ##    50  0  4  9  5  2  1  1
    ##    52  0  0  2  8  5  0  0
    ##    54  0  1  3  2  1  0  0
    ##    56  0  0  0  3  1  0  0
    ##    58  0  0  0  2  1  1  0
    ##    60  0  0  1  3  2  0  0
    ##    62  0  0  0  1  1  0  0
    ##    64  0  0  0  0  0  1  0

    ##  16  18  20  22  24  26  28  30  32  34  36  38  40  42  44  46  48  50 
    ##   4   4   9  12  15  12  37  73  84  94 122  95  75  59  44  33  31  22 
    ##  52  54  56  58  60  62  64 
    ##  15   7   4   4   6   2   1

    ##   1   2   3   4   5   6   7 
    ## 347 264 154  71  23   3   2

    ##      final_age
    ## lcat2  1  2  3  4  5  6  7
    ##    20  0  1  0  0  0  0  0
    ##    22  5  0  0  0  0  0  0
    ##    24 14  1  0  0  0  0  0
    ##    26 29  3  0  0  0  0  0
    ##    28 36 23  5  0  0  0  0
    ##    30 40 40 14  5  0  0  0
    ##    32 45 32 33  6  1  1  0
    ##    34 37 55 37 13  2  0  1
    ##    36 43 65 52 17  2  0  0
    ##    38  9 42 40 25  5  0  0
    ##    40  5 34 30 18  5  0  0
    ##    42  3 23 24 10  5  0  1
    ##    44  1  7 27 23  7  1  0
    ##    46  0 10 13 14  3  3  1
    ##    48  1  5 13 16  6  2  0
    ##    50  0  3  9 10  3  2  0
    ##    52  0  1  3  7  3  3  0
    ##    54  0  0  1  2  7  1  0
    ##    56  0  0  1  5  3  1  0
    ##    58  0  0  0  3  3  0  1
    ##    60  0  0  0  2  2  1  0
    ##    62  0  0  0  0  1  0  0
    ##    64  0  0  0  0  0  1  0

    ##  20  22  24  26  28  30  32  34  36  38  40  42  44  46  48  50  52  54 
    ##   1   5  15  32  64  99 118 145 179 121  92  66  66  44  43  27  17  11 
    ##  56  58  60  62  64 
    ##  10   7   5   1   1

    ##   1   2   3   4   5   6   7 
    ## 268 345 302 176  58  16   4

    ##      final_age
    ## lcat2   1   2   3   4   5   6   7   8   9
    ##    10   2   0   0   0   0   0   0   0   0
    ##    12   3   0   0   0   0   0   0   0   0
    ##    14   5   0   0   0   0   0   0   0   0
    ##    16  14   1   0   0   0   0   0   0   0
    ##    18  11   2   0   0   0   0   0   0   0
    ##    20  17   5   1   0   0   0   0   0   0
    ##    22  45  17   2   0   0   0   0   0   0
    ##    24 105  26   5   1   0   0   0   0   0
    ##    26  87  40  20   4   0   0   0   0   0
    ##    28  74  48  23   3   0   0   1   0   0
    ##    30  78  69  44   8   1   0   0   0   0
    ##    32  58  76  55  14   1   0   0   0   0
    ##    34  39  84  77  21   2   1   0   0   0
    ##    36  43 110  89  35   4   1   0   0   0
    ##    38  28  91  90  40   7   1   0   0   0
    ##    40  16  72  64  39  14   1   2   0   0
    ##    42  11  39  43  29   7   4   1   0   0
    ##    44   1  35  49  23  12   1   0   0   0
    ##    46   2  18  43  26   4   2   0   0   0
    ##    48   1  18  30  26   7   3   0   0   0
    ##    50   1  13  24  30   9   5   2   0   1
    ##    52   0   7  25  21  13   2   3   0   0
    ##    54   0   0  18  15  10   4   0   0   0
    ##    56   0   0  14  16  16   3   2   0   0
    ##    58   0   1   5  21  12   7   1   1   0
    ##    60   0   0   6  16   7   5   1   0   0
    ##    62   0   0   1   5   8   2   2   2   0
    ##    64   1   0   0   5   4   5   0   0   0
    ##    66   0   1   0   4   4   0   0   1   0
    ##    68   0   0   0   1   1   2   0   0   0
    ##    70   0   0   0   2   0   0   0   0   0
    ##    72   0   0   0   0   0   2   0   0   0
    ##    74   0   0   0   0   1   2   1   0   0

    ##  10  12  14  16  18  20  22  24  26  28  30  32  34  36  38  40  42  44 
    ##   2   3   5  15  13  23  64 137 151 149 200 204 224 282 257 208 134 121 
    ##  46  48  50  52  54  56  58  60  62  64  66  68  70  72  74 
    ##  95  85  85  71  47  51  48  35  20  15  10   4   2   2   4

    ##   1   2   3   4   5   6   7   8   9 
    ## 642 773 728 405 144  53  16   4   1

    ##      final_age
    ## lcat2  1  2  3  4  5  6  7  8
    ##    14  3  0  0  0  0  0  0  0
    ##    16  2  0  0  0  0  0  0  0
    ##    18  7  0  0  0  0  0  0  0
    ##    20  9  0  0  0  0  0  0  0
    ##    22  9  1  0  0  0  0  0  0
    ##    24 13  2  0  0  0  0  0  0
    ##    26 25  6  0  0  0  0  0  0
    ##    28 58 40  0  0  0  0  0  0
    ##    30 74 63  9  1  0  0  0  0
    ##    32 53 77 12  2  0  0  0  0
    ##    34 44 68 25  4  1  0  0  0
    ##    36 19 44 13  2  1  0  0  0
    ##    38 12 26 17  2  1  0  0  0
    ##    40  7 11  6  3  1  0  0  0
    ##    42  2  7 10  6  2  0  0  0
    ##    44  1  9  6  3  1  1  0  0
    ##    46  0  3  6  2  0  0  1  0
    ##    48  0  3  2  3  0  0  0  0
    ##    50  0  0  2  3  0  1  0  0
    ##    52  0  0  3  2  2  0  0  0
    ##    54  0  0  1  2  0  0  0  0
    ##    56  0  0  0  1  1  1  0  0
    ##    58  0  0  1  1  1  1  0  0
    ##    62  0  0  0  0  0  1  0  0
    ##    64  0  0  0  0  1  0  0  0
    ##    66  0  0  0  0  0  0  0  1
    ##    68  0  0  0  0  0  1  0  0

    ##  14  16  18  20  22  24  26  28  30  32  34  36  38  40  42  44  46  48 
    ##   3   2   7   9  10  15  31  98 147 144 142  79  58  28  27  21  12   8 
    ##  50  52  54  56  58  62  64  66  68 
    ##   6   7   3   3   4   1   1   1   1

    ##   1   2   3   4   5   6   7   8 
    ## 338 360 113  37  12   6   1   1

### MAKE OBSERVED ALK FROM ABOVE TABLES

The conditional proportions that form the ALK are calculated by dividing ecah cell of the frequency table by the sum of the corresponding row. These row proportions are constructed by submitting the xtabs() object to prop.table() and including margin=1 to indicate that the proportions are computed by row (page 92).

The alkPlot command used for plotting the observed ALK is unable to extend the x axis to the bounds of c(0,80) because xlim is not working. Therefore, in order to produce a plot with an x axis that can span from 0-80 (like what is happening with the length frequency and the smoothed ALK) I need to add in "observed" proportions for length categories that were not sampled. I could have added them to the original data frame but I was concerned that in the process of proportion calculations the extra entries would affect the proportions or result in proportions that I didn't want. The smoothing process can estimate proportions outside of the range but I wanted to keep the observed plot with just the proportions from the observed data. Therefore, I added in the zero proportion data by writing and editing a csv file which I then read below.

``` r
as.data.frame.matrix((prop.table(rawfreq_TB, margin=1))) %>% write.csv("U:/PhD_projectfiles/Exported_R_Datafiles/alk_TB.csv")
alk_TB <- read.csv("U:/PhD_projectfiles/Exported_R_Datafiles/alk_TB_edit.csv", row.names=1)
names(alk_TB) <- c(1,2,3,4,5,6,7,8,9)
round(alk_TB,3)
```

    ##        1     2     3     4     5     6     7     8     9
    ## 0  0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 2  0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 4  0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 6  0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 8  0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 10 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 12 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 14 0.500 0.500 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 16 1.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 18 0.500 0.500 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 20 1.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 22 1.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 24 0.923 0.077 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 26 0.900 0.050 0.050 0.000 0.000 0.000 0.000 0.000 0.000
    ## 28 0.569 0.306 0.125 0.000 0.000 0.000 0.000 0.000 0.000
    ## 30 0.371 0.401 0.198 0.024 0.006 0.000 0.000 0.000 0.000
    ## 32 0.332 0.324 0.261 0.075 0.008 0.000 0.000 0.000 0.000
    ## 34 0.267 0.368 0.250 0.101 0.010 0.003 0.000 0.000 0.000
    ## 36 0.249 0.339 0.278 0.094 0.029 0.008 0.004 0.000 0.000
    ## 38 0.109 0.370 0.309 0.133 0.067 0.012 0.000 0.000 0.000
    ## 40 0.050 0.270 0.326 0.234 0.064 0.035 0.021 0.000 0.000
    ## 42 0.028 0.324 0.331 0.221 0.069 0.014 0.007 0.007 0.000
    ## 44 0.010 0.219 0.314 0.238 0.124 0.086 0.010 0.000 0.000
    ## 46 0.000 0.254 0.349 0.190 0.127 0.048 0.016 0.016 0.000
    ## 48 0.000 0.073 0.273 0.291 0.273 0.073 0.018 0.000 0.000
    ## 50 0.026 0.000 0.289 0.237 0.263 0.184 0.000 0.000 0.000
    ## 52 0.000 0.042 0.208 0.292 0.292 0.083 0.083 0.000 0.000
    ## 54 0.000 0.000 0.200 0.360 0.320 0.040 0.080 0.000 0.000
    ## 56 0.000 0.000 0.111 0.111 0.333 0.222 0.222 0.000 0.000
    ## 58 0.000 0.000 0.118 0.235 0.412 0.176 0.000 0.059 0.000
    ## 60 0.000 0.000 0.000 0.000 0.182 0.545 0.273 0.000 0.000
    ## 62 0.000 0.000 0.000 0.222 0.333 0.111 0.222 0.111 0.000
    ## 64 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 66 0.000 0.000 0.000 0.000 0.000 0.333 0.333 0.000 0.333
    ## 68 0.000 0.000 0.000 0.000 1.000 0.000 0.000 0.000 0.000
    ## 70 0.000 0.000 0.000 0.000 0.000 0.000 1.000 0.000 0.000
    ## 72 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 74 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 76 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 78 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 80 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000

Now do the same with all of the others

    ##        1     2     3     4     5     6     7
    ## 0  0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 2  0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 4  0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 6  0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 8  0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 10 1.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 12 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 14 1.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 16 0.000 0.000 1.000 0.000 0.000 0.000 0.000
    ## 18 1.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 20 0.750 0.250 0.000 0.000 0.000 0.000 0.000
    ## 22 1.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 24 0.933 0.067 0.000 0.000 0.000 0.000 0.000
    ## 26 0.879 0.121 0.000 0.000 0.000 0.000 0.000
    ## 28 0.562 0.359 0.078 0.000 0.000 0.000 0.000
    ## 30 0.404 0.404 0.141 0.051 0.000 0.000 0.000
    ## 32 0.381 0.271 0.280 0.051 0.008 0.008 0.000
    ## 34 0.255 0.379 0.255 0.090 0.014 0.000 0.007
    ## 36 0.240 0.363 0.291 0.095 0.011 0.000 0.000
    ## 38 0.074 0.347 0.331 0.207 0.041 0.000 0.000
    ## 40 0.054 0.370 0.326 0.196 0.054 0.000 0.000
    ## 42 0.045 0.348 0.364 0.152 0.076 0.000 0.015
    ## 44 0.015 0.106 0.409 0.348 0.106 0.015 0.000
    ## 46 0.000 0.227 0.295 0.318 0.068 0.068 0.023
    ## 48 0.023 0.116 0.302 0.372 0.140 0.047 0.000
    ## 50 0.000 0.111 0.333 0.370 0.111 0.074 0.000
    ## 52 0.000 0.059 0.176 0.412 0.176 0.176 0.000
    ## 54 0.000 0.000 0.091 0.182 0.636 0.091 0.000
    ## 56 0.000 0.000 0.100 0.500 0.300 0.100 0.000
    ## 58 0.000 0.000 0.000 0.429 0.429 0.000 0.143
    ## 60 0.000 0.000 0.000 0.400 0.400 0.200 0.000
    ## 62 0.000 0.000 0.000 0.000 1.000 0.000 0.000
    ## 64 0.000 0.000 0.000 0.000 0.000 1.000 0.000
    ## 66 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 68 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 70 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 72 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 74 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 76 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 78 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 80 0.000 0.000 0.000 0.000 0.000 0.000 0.000

    ##        1     2     3     4     5     6     7
    ## 0  0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 2  0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 4  0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 6  0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 8  0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 10 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 12 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 14 1.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 16 0.900 0.000 0.100 0.000 0.000 0.000 0.000
    ## 18 1.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 20 0.875 0.125 0.000 0.000 0.000 0.000 0.000
    ## 22 1.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 24 0.938 0.062 0.000 0.000 0.000 0.000 0.000
    ## 26 0.917 0.083 0.000 0.000 0.000 0.000 0.000
    ## 28 0.838 0.162 0.000 0.000 0.000 0.000 0.000
    ## 30 0.753 0.233 0.014 0.000 0.000 0.000 0.000
    ## 32 0.667 0.274 0.060 0.000 0.000 0.000 0.000
    ## 34 0.543 0.330 0.128 0.000 0.000 0.000 0.000
    ## 36 0.418 0.434 0.123 0.016 0.008 0.000 0.000
    ## 38 0.295 0.442 0.211 0.053 0.000 0.000 0.000
    ## 40 0.184 0.500 0.237 0.079 0.000 0.000 0.000
    ## 42 0.119 0.407 0.271 0.186 0.017 0.000 0.000
    ## 44 0.045 0.227 0.500 0.205 0.023 0.000 0.000
    ## 46 0.000 0.182 0.545 0.152 0.091 0.000 0.030
    ## 48 0.000 0.194 0.387 0.290 0.129 0.000 0.000
    ## 50 0.000 0.182 0.409 0.227 0.091 0.045 0.045
    ## 52 0.000 0.000 0.133 0.533 0.333 0.000 0.000
    ## 54 0.000 0.143 0.429 0.286 0.143 0.000 0.000
    ## 56 0.000 0.000 0.000 0.750 0.250 0.000 0.000
    ## 58 0.000 0.000 0.000 0.500 0.250 0.250 0.000
    ## 60 0.000 0.000 0.167 0.500 0.333 0.000 0.000
    ## 62 0.000 0.000 0.000 0.500 0.500 0.000 0.000
    ## 64 0.000 0.000 0.000 0.000 0.000 1.000 0.000
    ## 66 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 68 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 70 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 72 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 74 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 76 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 78 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 80 0.000 0.000 0.000 0.000 0.000 0.000 0.000

    ##        1     2     3     4     5     6     7     8     9  10
    ## 0  0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.0
    ## 2  0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.0
    ## 4  0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.0
    ## 6  1.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.0
    ## 8  0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.0
    ## 10 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.0
    ## 12 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.0
    ## 14 1.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.0
    ## 16 1.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.0
    ## 18 1.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.0
    ## 20 1.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.0
    ## 22 1.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.0
    ## 24 1.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.0
    ## 26 0.867 0.000 0.000 0.067 0.000 0.067 0.000 0.000 0.000 0.0
    ## 28 0.759 0.185 0.037 0.019 0.000 0.000 0.000 0.000 0.000 0.0
    ## 30 0.627 0.320 0.053 0.000 0.000 0.000 0.000 0.000 0.000 0.0
    ## 32 0.450 0.410 0.100 0.030 0.010 0.000 0.000 0.000 0.000 0.0
    ## 34 0.261 0.461 0.226 0.043 0.000 0.009 0.000 0.000 0.000 0.0
    ## 36 0.257 0.459 0.209 0.061 0.014 0.000 0.000 0.000 0.000 0.0
    ## 38 0.168 0.534 0.214 0.061 0.023 0.000 0.000 0.000 0.000 0.0
    ## 40 0.038 0.523 0.235 0.136 0.053 0.008 0.000 0.000 0.008 0.0
    ## 42 0.026 0.465 0.316 0.132 0.061 0.000 0.000 0.000 0.000 0.0
    ## 44 0.013 0.400 0.373 0.200 0.013 0.000 0.000 0.000 0.000 0.0
    ## 46 0.000 0.236 0.364 0.273 0.073 0.055 0.000 0.000 0.000 0.0
    ## 48 0.000 0.188 0.396 0.354 0.062 0.000 0.000 0.000 0.000 0.0
    ## 50 0.000 0.116 0.326 0.326 0.163 0.047 0.023 0.000 0.000 0.0
    ## 52 0.000 0.000 0.333 0.407 0.185 0.000 0.074 0.000 0.000 0.0
    ## 54 0.000 0.000 0.263 0.395 0.211 0.132 0.000 0.000 0.000 0.0
    ## 56 0.000 0.056 0.167 0.444 0.278 0.056 0.000 0.000 0.000 0.0
    ## 58 0.000 0.000 0.000 0.357 0.143 0.286 0.071 0.143 0.000 0.0
    ## 60 0.000 0.000 0.000 0.267 0.333 0.200 0.200 0.000 0.000 0.0
    ## 62 0.000 0.000 0.000 0.600 0.200 0.200 0.000 0.000 0.000 0.0
    ## 64 0.000 0.000 0.000 0.250 0.125 0.125 0.500 0.000 0.000 0.0
    ## 66 0.000 0.000 0.000 0.000 0.000 0.000 1.000 0.000 0.000 0.0
    ## 68 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.0
    ## 70 0.000 0.000 0.000 0.000 0.000 0.000 0.500 0.000 0.000 0.5
    ## 72 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.0
    ## 73 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.0
    ## 76 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.0
    ## 78 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.0
    ## 80 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.0

    ##        1     2     3     4     5     6     7 8
    ## 0  0.000 0.000 0.000 0.000 0.000 0.000 0.000 0
    ## 2  0.000 0.000 0.000 0.000 0.000 0.000 0.000 0
    ## 4  0.000 0.000 0.000 0.000 0.000 0.000 0.000 0
    ## 6  0.000 0.000 0.000 0.000 0.000 0.000 0.000 0
    ## 8  0.000 0.000 0.000 0.000 0.000 0.000 0.000 0
    ## 10 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0
    ## 12 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0
    ## 14 1.000 0.000 0.000 0.000 0.000 0.000 0.000 0
    ## 16 1.000 0.000 0.000 0.000 0.000 0.000 0.000 0
    ## 18 1.000 0.000 0.000 0.000 0.000 0.000 0.000 0
    ## 20 1.000 0.000 0.000 0.000 0.000 0.000 0.000 0
    ## 22 0.900 0.100 0.000 0.000 0.000 0.000 0.000 0
    ## 24 0.867 0.133 0.000 0.000 0.000 0.000 0.000 0
    ## 26 0.806 0.194 0.000 0.000 0.000 0.000 0.000 0
    ## 28 0.592 0.408 0.000 0.000 0.000 0.000 0.000 0
    ## 30 0.500 0.432 0.061 0.007 0.000 0.000 0.000 0
    ## 32 0.368 0.535 0.083 0.014 0.000 0.000 0.000 0
    ## 34 0.310 0.479 0.176 0.028 0.007 0.000 0.000 0
    ## 36 0.241 0.557 0.165 0.025 0.013 0.000 0.000 0
    ## 38 0.207 0.448 0.293 0.034 0.017 0.000 0.000 0
    ## 40 0.250 0.393 0.214 0.107 0.036 0.000 0.000 0
    ## 42 0.074 0.259 0.370 0.222 0.074 0.000 0.000 0
    ## 44 0.048 0.429 0.286 0.143 0.048 0.048 0.000 0
    ## 46 0.000 0.250 0.500 0.167 0.000 0.000 0.083 0
    ## 48 0.000 0.375 0.250 0.375 0.000 0.000 0.000 0
    ## 50 0.000 0.000 0.333 0.500 0.000 0.167 0.000 0
    ## 52 0.000 0.000 0.429 0.286 0.286 0.000 0.000 0
    ## 54 0.000 0.000 0.333 0.667 0.000 0.000 0.000 0
    ## 56 0.000 0.000 0.000 0.333 0.333 0.333 0.000 0
    ## 58 0.000 0.000 0.250 0.250 0.250 0.250 0.000 0
    ## 60 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0
    ## 62 0.000 0.000 0.000 0.000 0.000 1.000 0.000 0
    ## 64 0.000 0.000 0.000 0.000 1.000 0.000 0.000 0
    ## 66 0.000 0.000 0.000 0.000 0.000 0.000 0.000 1
    ## 68 0.000 0.000 0.000 0.000 0.000 1.000 0.000 0
    ## 70 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0
    ## 72 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0
    ## 74 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0
    ## 76 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0
    ## 78 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0
    ## 80 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0

    ##        1     2     3     4     5     6     7     8     9
    ## 0  0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 2  0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 4  0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 6  0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 8  1.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 10 1.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 12 1.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 14 1.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 16 0.938 0.062 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 18 0.846 0.154 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 20 0.739 0.217 0.043 0.000 0.000 0.000 0.000 0.000 0.000
    ## 22 0.703 0.266 0.031 0.000 0.000 0.000 0.000 0.000 0.000
    ## 24 0.766 0.190 0.036 0.007 0.000 0.000 0.000 0.000 0.000
    ## 26 0.576 0.265 0.132 0.026 0.000 0.000 0.000 0.000 0.000
    ## 28 0.497 0.322 0.154 0.020 0.000 0.000 0.007 0.000 0.000
    ## 30 0.390 0.345 0.220 0.040 0.005 0.000 0.000 0.000 0.000
    ## 32 0.284 0.373 0.270 0.069 0.005 0.000 0.000 0.000 0.000
    ## 34 0.174 0.375 0.344 0.094 0.009 0.004 0.000 0.000 0.000
    ## 36 0.152 0.390 0.316 0.124 0.014 0.004 0.000 0.000 0.000
    ## 38 0.109 0.354 0.350 0.156 0.027 0.004 0.000 0.000 0.000
    ## 40 0.077 0.346 0.308 0.188 0.067 0.005 0.010 0.000 0.000
    ## 42 0.082 0.291 0.321 0.216 0.052 0.030 0.007 0.000 0.000
    ## 44 0.008 0.289 0.405 0.190 0.099 0.008 0.000 0.000 0.000
    ## 46 0.021 0.189 0.453 0.274 0.042 0.021 0.000 0.000 0.000
    ## 48 0.012 0.212 0.353 0.306 0.082 0.035 0.000 0.000 0.000
    ## 50 0.012 0.153 0.282 0.353 0.106 0.059 0.024 0.000 0.012
    ## 52 0.000 0.099 0.352 0.296 0.183 0.028 0.042 0.000 0.000
    ## 54 0.000 0.000 0.383 0.319 0.213 0.085 0.000 0.000 0.000
    ## 56 0.000 0.000 0.275 0.314 0.314 0.059 0.039 0.000 0.000
    ## 58 0.000 0.021 0.104 0.438 0.250 0.146 0.021 0.021 0.000
    ## 60 0.000 0.000 0.167 0.444 0.222 0.139 0.028 0.000 0.000
    ## 62 0.000 0.000 0.050 0.250 0.400 0.100 0.100 0.100 0.000
    ## 64 0.067 0.000 0.000 0.333 0.267 0.333 0.000 0.000 0.000
    ## 66 0.000 0.100 0.000 0.400 0.400 0.000 0.000 0.100 0.000
    ## 68 0.000 0.000 0.000 0.250 0.250 0.500 0.000 0.000 0.000
    ## 70 0.000 0.000 0.000 1.000 0.000 0.000 0.000 0.000 0.000
    ## 72 0.000 0.000 0.000 0.000 0.000 1.000 0.000 0.000 0.000
    ## 74 0.000 0.000 0.000 0.000 0.250 0.500 0.250 0.000 0.000
    ## 76 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 78 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    ## 80 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000

### Apply Proportional Age Distribution and Mean Length at Age

1.  Proportional age distribution
2.  Mean length-at-age. (otolith database)

``` r
#Proportional age distribution
  # might need to make sure I add some dummy data for bays that don't have equal number of ages
ad_TB <- xtabs(~final_age, data=Agelength_TB)
round(prop.table(ad_TB), 3)
```

    ## final_age
    ##     1     2     3     4     5     6     7     8     9 
    ## 0.210 0.291 0.259 0.131 0.066 0.028 0.012 0.002 0.001

Now do for other estuaries

    ## final_age
    ##     1     2     3     4     5     6     7 
    ## 0.402 0.306 0.178 0.082 0.027 0.003 0.002

    ## final_age
    ##     1     2     3     4     5     6     7 
    ## 0.229 0.295 0.258 0.151 0.050 0.014 0.003

    ## final_age
    ##     1     2     3     4     5     6     7     8     9    10 
    ## 0.214 0.355 0.216 0.135 0.049 0.018 0.010 0.002 0.001 0.001

    ## final_age
    ##     1     2     3     4     5     6     7     8     9 
    ## 0.232 0.279 0.263 0.146 0.052 0.019 0.006 0.001 0.000

    ## final_age
    ##     1     2     3     4     5     6     7     8 
    ## 0.389 0.415 0.130 0.043 0.014 0.007 0.001 0.001

#### Mean Length-at-age stats for each

``` r
TB_sumlen <- Agelength_TB %>% group_by(final_age) %>% summarize(n=validn(tl), mn=mean(tl, na.rm=TRUE),
                                                    sd=sd(tl, na.rm=TRUE), se=se(tl, na.rm=TRUE)) %>%
                                                    as.data.frame()
```

Now do the rest of the estuaries

### MULTIPLOT AGE HISTOGRAMS

``` r
age_AP <- ggplot(Agelength_AP, aes(x=final_age))+ 
  geom_histogram(aes(y=(..count..)/sum(..count..)), binwidth=1, boundary=-0.5, fill="white", color="black")+    # plotting in percent frequency
  scale_y_continuous(labels=percent_format(), name="Frequency (%)", limits=c(0,.4))+  ## plotting in percent frequency
  scale_x_continuous(name="Age (years)", limits=c(0,11), breaks=seq(1,10,1), labels=c(1,2,3,4,5,6,7,8,9,"10+"))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_rect(colour="black",fill="white"), 
        axis.text.x=element_text(colour="black", size=24, family="Times New Roman"), #changing  colour of x axisaxis.text.y=element_text(colour="black"), #changing colour of y acis
        axis.text.y=element_text(colour="black", size=24, family="Times New Roman"), #changing colour of y acis
        axis.title.x=element_text(size=24, family="Times New Roman"),
        axis.title.y=element_text(size=24, family="Times New Roman"),
        plot.title=element_text(size=14))+
  annotate("text", x=10, y=.35, label="AP", size=10, family="Times New Roman")

age_AP
```

![](ALK_analysis_files/figure-markdown_github/unnamed-chunk-24-1.png) Now do the others

![](ALK_analysis_files/figure-markdown_github/unnamed-chunk-25-1.png)![](ALK_analysis_files/figure-markdown_github/unnamed-chunk-25-2.png)![](ALK_analysis_files/figure-markdown_github/unnamed-chunk-25-3.png)![](ALK_analysis_files/figure-markdown_github/unnamed-chunk-25-4.png)![](ALK_analysis_files/figure-markdown_github/unnamed-chunk-25-5.png)

### MULTIPLOT LENGTH HISTOGRAMS

``` r
length_AP <- ggplot(Agelength_AP, aes(x=tl))+ 
  geom_histogram(aes(y=(..count..)/sum(..count..)), binwidth=1, boundary=-0.5, fill="white", color="black")+    # plotting in percent frequency
  scale_y_continuous(labels=percent_format(), name="Frequency (%)", limits=c(0,.10), breaks=seq(0,.10, .02))+  ## plotting in percent frequency
  scale_x_continuous(name="Length (cm)", limits=c(0,80), breaks=seq(0,80,20))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_rect(colour="black",fill="white"), 
        axis.text.x=element_text(colour="black", size=24, family="Times New Roman"), #changing  colour of x axisaxis.text.y=element_text(colour="black"), #changing colour of y acis
        axis.text.y=element_text(colour="black", size=24, family="Times New Roman"), #changing colour of y acis
        axis.title.x=element_text(size=24, family="Times New Roman"),
        axis.title.y=element_text(size=24, family="Times New Roman"),
        plot.title=element_text(size=14))+
  annotate("text", x=65, y= .08, label="AP", size=10, family="Times New Roman")
length_AP
```

![](ALK_analysis_files/figure-markdown_github/unnamed-chunk-27-1.png) Now make plots for others.

![](ALK_analysis_files/figure-markdown_github/unnamed-chunk-28-1.png)![](ALK_analysis_files/figure-markdown_github/unnamed-chunk-28-2.png)![](ALK_analysis_files/figure-markdown_github/unnamed-chunk-28-3.png)![](ALK_analysis_files/figure-markdown_github/unnamed-chunk-28-4.png)![](ALK_analysis_files/figure-markdown_github/unnamed-chunk-28-5.png)

When you save the notebook, an HTML file containing the code and output will be saved alongside it (click the *Preview* button or press *Ctrl+Shift+K* to preview the HTML file).
