# VOR Data for poke-an-egg experiment (conducted summer 2015)
# Nov 2015
# Julie Jung

ls()
rm(list=ls())
ls()
setwd('/Users/juliejung/Desktop/VOR m.s.') 
getwd()         


###################################################################################################################################################################################################################################################################################
###################################################################################################################################################################################################################################################################################
###################################################################################################################################################################################################################################################################################
###############                    f u n c t i o n s                     ##########################################################################################################################################################################################################
###################################################################################################################################################################################################################################################################################
###################################################################################################################################################################################################################################################################################
###################################################################################################################################################################################################################################################################################

# 
# #library(reshape2)

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

#test differences in variation

globalVariables(c("stats", "aggregate", "pchisq", "rchisq", "rnorm", "sd"))

# Miller test for raw data ----------------------------

#' Asymptotic test for the equality of coefficients of variation from k populations, using measurement data
#'
#' Test for k samples (k sample populations with unequal sized) from
#' Feltz CJ, Miller GE (1996) An asymptotic test for the equality of coefficients
#' of variation from k population. Stat Med 15:647–658
#'
#' @param x a numeric vector containing individual measurement values
#' @param y a vector of any type containing a grouping variable
#'
#' @return a list with the test statistic and p-value
#' @export
#'
#' @examples
#'
#'  y <- unlist(lapply(letters[1:5], function(i) rep(i, 20)))
#'  x <- rnorm(100)
#'
#'  asymptotic_test(x, y)
#'
asymptotic_test <-
  function(x, y){
    
    if (!is.numeric(x) && (!is.numeric(y) | !is.character(y) | !is.integer(y))) {
      warning("x is not numeric or y is not numeric, character or integer: returning NA")
      return(NA_real_)
    }
    
    if (anyNA(x)) {
      warning("x cannot contain any NA values: returning NA")
      return(NA_real_)
    }
    
    if (anyNA(y)) {
      warning("y cannot contain any NA values: returning NA")
      return(NA_real_)
    }
    
    # k is the number of groups
    k <- length(unique(y))
    # j is an index referring to the group number
    # n_j is the sample size of the jth group
    n_j <- data.frame(table(y))$Freq
    # s_j is the sd of the jth population
    s_j <-  aggregate(x, by = list(y), FUN = sd)[2]
    # x_j is the mean of the jth population
    x_j <- aggregate(x, by = list(y), FUN = mean)[2]
    
    m_j <- n_j - 1
    
    D <- (sum(m_j * (s_j/x_j))) / sum(m_j)
    
    D_AD <- (sum(m_j * (s_j/x_j - D)^2 )) / ( D^2 * (0.5 + D^2) )
    
    # D_AD distributes as a central chi-sq random variable with k - 1 degrees of freedom
    p_value <-  pchisq(D_AD,
                       k - 1,
                       lower.tail = FALSE)
    
    return(list(D_AD = D_AD,
                p_value = p_value))
  }


# Miller test for summary data ------------------------

#' Asymptotic test for the equality of coefficients of variation from k populations, using summary statistics when raw measurement data are not available.
#'
#' Test for k samples (k sample populations with unequal sized) from
#' Feltz CJ, Miller GE (1996) An asymptotic test for the equality of coefficients
#' of variation from k population. Stat Med 15:647–658
#'
#' @param k a numeric vector the number of groups
#' @param n a numeric vector the numer of measurements in each group
#' @param s a numeric vector the standard deviation of each group
#' @param x a numeric vector the mean of each group
#'
#' @return a list with the test statistic and p-value
#' @export
#'
#' @examples
#'
#' # Summary stats from Feltz and Miller 1996
#'
#' miller <- data.frame(test = c('ELISA', 'WEHI', '`Viral inhibition`'),
#'                     Mean = c(6.8, 8.5, 6.0),
#'                     CV =   c(0.090, 0.462, 0.340),
#'                     N =    c(5, 5, 5))
#' # compute SD from mean and cv
#' miller$SD <- with(miller, CV * Mean)
#'
#'  asymptotic_test2(k = nrow(miller), n = miller$N, s = miller$SD, x = miller$Mean)
#'
asymptotic_test2 <-
  function(k, n, s, x){
    
    
    if (!is.integer(k) && !is.numeric(k) && !is.integer(n) && !is.numeric(n) && !is.numeric(s) && !is.numeric(x) ){
      warning("input values are not numeric: returning NA")
      return(NA_real_)
    }
    
    
    # k is the number of groups
    k <- k
    # j is an index referring to the group number
    # n_j is the sample size of the jth group
    n_j <- n
    # s_j is the sd of the jth population
    s_j <-  s
    # x_j is the mean of the jth population
    x_j <- x
    
    m_j <- n_j - 1
    
    D <- (sum(m_j * (s_j/x_j))) / sum(m_j)
    
    D_AD <- (sum(m_j * (s_j/x_j - D)^2 )) / ( D^2 * (0.5 + D^2) )
    
    # D_AD distributes as a central chi-sq random variable with k - 1 degrees of freedom
    p_value <-  pchisq(D_AD,
                       k - 1,
                       lower.tail = FALSE)
    
    return(list(D_AD = D_AD,
                p_value = p_value))
  }



# Krishnamoorthy test for raw data --------------------

#' Modified signed-likelihood ratio test (SLRT) for equality of CVs, using measurement data

#' @references \url{http://link.springer.com/article/10.1007/s00180-013-0445-2}
#' Krishnamoorthy, K. & Lee, M. Comput Stat (2014) 29: 215. doi:10.1007/s00180-013-0445-2
#'
#'
#' @param nr numeric vector length one, number of simulation runs, default is 1e3
#' @param x a numeric vector containing individual measurement values
#' @param y a vector of any type containing a grouping variable
#'
#' @return a list with the test statistic and p-value
#' @export
#'
#' @examples
#'
#'  x <- rnorm(100)
#'  y <- unlist(lapply(letters[1:5], function(i) rep(i, 20)))
#'
#'  mslr_test(nr = 1e3, x, y)
#'
#'
mslr_test <-
  function(nr = 1e3, x, y){
    
    if (!is.numeric(x) && !is.numeric(y) && !is.character(y)) {
      warning("x is not numeric or y is not numeric or character: returning NA")
      return(NA_real_)
    }
    
    if (anyNA(x)) {
      warning("x cannot contain any NA values: returning NA")
      return(NA_real_)
    }
    
    if (anyNA(y)) {
      warning("y cannot contain any NA values: returning NA")
      return(NA_real_)
    }
    
    # j is an index referring to the group number
    # n is the sample size of the jth group
    n <- data.frame(table(y))$Freq
    # s is the sd of the jth population
    s <-  aggregate(x, by = list(y), FUN = sd)$x
    # x_j is the mean of the jth population
    x <- aggregate(x, by = list(y), FUN = mean)$x
    
    
    k <- length(x)
    gv <- as.vector(nr)
    df <- n - 1
    xst0 <- LRT_STAT(n, x, s)
    uh0 <- xst0[1:k]
    tauh0 <- xst0[k+1]; stat0 <- xst0[k+2]
    sh0 <- tauh0*uh0
    se0 <- tauh0*uh0/sqrt(n)
    # PB estimates of the mean and SD of the LRT
    for(ii in 1:nr){
      z <- rnorm(k)
      x <- uh0 + z*se0
      ch <- rchisq(k,df)
      s <- sh0*sqrt(ch/df)
      xst <- LRT_STAT(n,x,s)
      gv[ii] <- xst[k+2]
    }
    am <- mean(gv); sd <- sd(gv)
    # end PB estimates
    statm <- sqrt(2.0*(k-1))*(stat0-am)/sd+(k-1)
    pval <- 1.0-pchisq(statm,k-1)
    
    
    return(list(MSLRT = statm,
                p_value = pval))
  }






#' # Modified signed-likelihood ratio test (SLRT) for equality of CVs, using summary statistics when raw measurement data are not available.
#'
#'
# Krishnamoorthy, K. & Lee, M. Comput Stat (2014) 29: 215. doi:10.1007/s00180-013-0445-2
#' @references \url{http://link.springer.com/article/10.1007/s00180-013-0445-2}
#
#'
#' @param nr numeric vector lenght one, number of simulation runs
#' @param n a numeric vector, the number of observations in each group
#' @param x a numeric vector, the mean of each group
#' @param s a numeric vector, the standard deviation of each group
#'
#' @return a list with the test statistic and p-value
#' @export
#'
#' @examples
#'
#' # Summary stats from Feltz and Miller 1996
#'
#' miller <- data.frame(test = c('ELISA', 'WEHI', '`Viral inhibition`'),
#'                     Mean = c(6.8, 8.5, 6.0),
#'                     CV =   c(0.090, 0.462, 0.340),
#'                     N =    c(5, 5, 5))
#' # compute SD from mean and cv
#' miller$SD <- with(miller, CV * Mean)
#'
#'  mslr_test2(nr = 1e3, n = miller$N, s = miller$SD, x = miller$Mean)
#'
#'
mslr_test2 <-
  function(nr, n, x, s){
    
    if (!is.numeric(nr) &&!is.numeric(n) && !is.integer(n) && !is.numeric(s) && !is.numeric(x)) {
      warning("input values are not numeric: returning NA")
      return(NA_real_)
    }
    
    
    k <- length(x)
    gv <- as.vector(nr)
    df <- n - 1
    xst0 <- LRT_STAT(n, x, s)
    uh0 <- xst0[1:k]
    tauh0 <- xst0[k+1]; stat0 <- xst0[k+2]
    sh0 <- tauh0*uh0
    se0 <- tauh0*uh0/sqrt(n)
    # PB estimates of the mean and SD of the LRT
    for(ii in 1:nr){
      z <- rnorm(k)
      x <- uh0 + z*se0
      ch <- rchisq(k,df)
      s <- sh0*sqrt(ch/df)
      xst <- LRT_STAT(n,x,s)
      gv[ii] <- xst[k+2]
    }
    am <- mean(gv); sd <- sd(gv)
    # end PB estimates
    statm <- sqrt(2.0*(k-1))*(stat0-am)/sd+(k-1)
    pval <- 1.0-pchisq(statm,k-1)
    
    
    return(list(MSLRT = statm,
                p_value = pval))
  }


#' LRT_STAT, required by mlrt_test
#'
#' @param n ... as above
#' @param x ...
#' @param s ...
#'
#' @return xx
#'
LRT_STAT <- function(n, x, s){
  
  k <- length(x)
  df <- n-1
  ssq <- s**2
  vsq <- df*ssq/n
  v <- sqrt(vsq)
  sn <- sum(n)
  #MLES
  tau0 <- sum(n*vsq/x**2)/sn
  l <- 1
  repeat{
    uh <- (-x+sqrt(x**2+4.0*tau0*(vsq+x**2)))/2.0/tau0
    tau <- sum(n*(vsq+(x-uh)**2)/uh**2)/sn
    if(abs(tau-tau0) <= 1.0e-7 || l > 30) break
    l <- l + 1
    tau0 <- tau
  }
  tauh <- sqrt(tau)
  #END MLES
  elf <- 0.0; clf <- 0.0
  for(j in 1:k){
    clf <- clf-n[j]*log(tauh*uh[j])-(n[j]*(vsq[j]+(x[j]-uh[j])**2))/(2.0*tauh**2*uh[j]**2)
    elf <- elf-n[j]*log(v[j])-n[j]/2.0
  }
  stat <- 2.0*(elf-clf)
  return(c(uh,tauh,stat))
}


###################################################################################################################################################################################################################################################################################
###################################################################################################################################################################################################################################################################################
###################################################################################################################################################################################################################################################################################
###############                    start c o d e                     ##############################################################################################################################################################################################################
###################################################################################################################################################################################################################################################################################
###################################################################################################################################################################################################################################################################################
###################################################################################################################################################################################################################################################################################







###################################################################################################################################################################################################################################################################################
###################################################################################################################################################################################################################################################################################
###############        ontogeny of V O R                                    #######################################################################################################################################################################################################
###############     d a i l y   developmental series #1                     #######################################################################################################################################################################################################
###################################################################################################################################################################################################################################################################################
###################################################################################################################################################################################################################################################################################
###################################################################################################################################################################################################################################################################################


Daily.df<-read.csv(file="Daily.csv")
error <- summarySE(Daily.df, measurevar="Average", groupvars="Age")

library(clinfun)

jonckheere.test(as.numeric(Daily.df$Average), as.numeric(Daily.df$Age), alternative="increasing")
jonckheere.test(as.numeric(Daily.df$Average), as.numeric(Daily.df$Stage), alternative="increasing")

cor.test(Daily.df$Age, Daily.df$Average, method="pearson")
cor.test(Daily.df$Stage, Daily.df$Average, method="pearson")

hist(Daily.df$Average) #normal enough to use an ANOVA
lm1 <- lm(Average ~ Age, data=Daily.df)
lm2 <- lm(Average ~ Stage, data=Daily.df)
lm3 <- lm(Average ~ Clutch, data=Daily.df)

library(car)
Anova(lm1)
Anova(lm2)
Anova(lm3)
 
c8 <- subset(Daily.df, Clutch == "C8")
c9 <- subset(Daily.df, Clutch == "C9")
c10 <- subset(Daily.df, Clutch == "C10")
c11 <- subset(Daily.df, Clutch == "C11")
c12 <- subset(Daily.df, Clutch == "C12")
c13 <- subset(Daily.df, Clutch == "C13")
c14 <- subset(Daily.df, Clutch == "C14")
c15 <- subset(Daily.df, Clutch == "C15")
c16 <- subset(Daily.df, Clutch == "C16")
c17 <- subset(Daily.df, Clutch == "C17")
c19 <- subset(Daily.df, Clutch == "C19")
c20 <- subset(Daily.df, Clutch == "C20")
c22 <- subset(Daily.df, Clutch == "C22")
c23 <- subset(Daily.df, Clutch == "C23")

sum(Daily.df$Age==3) #6
sum(Daily.df$Age==4)#10
sum(Daily.df$Age==5) #9
sum(Daily.df$Age==6) #9
sum(Daily.df$Age==7) #2

mean(Daily.df$Age) #4.75

col = rainbow(14)

###############   colorful plot with lines   ###############

# ggplot(Daily.df, aes(x=Age, y=Average, colour = Tad)) + 
#   geom_point(data=c8, colour=col[1], size=3) +
#   geom_line(data=c8, colour=col[1])+
#   geom_point(data=c9, colour=col[2], size=3) +
#   geom_line(data=c9, colour=col[2])+
#   geom_point(data=c10, colour=col[3], size=3) +
#   geom_line(data=c10, colour=col[3])+
#   geom_point(data=c11, colour=col[4], size=3) +
#   geom_line(data=c11, colour=col[4])+
#   geom_point(data=c12, colour=col[5], size=3) +
#   geom_line(data=c12, colour=col[5])+
#   geom_point(data=c13, colour=col[6], size=3) +
#   geom_line(data=c13, colour=col[6])+
#   geom_point(data=c14, colour=col[7], size=3) +
#   geom_line(data=c14, colour=col[7])+
#   geom_point(data=c15, colour=col[8], size=3) +
#   geom_line(data=c15, colour=col[8])+
#   geom_point(data=c16, colour=col[9], size=3) +
#   geom_line(data=c16, colour=col[9])+
#   geom_point(data=c17, colour=col[10], size=3) +
#   geom_line(data=c17, colour=col[10])+
#   geom_point(data=c19, colour=col[11], size=3) +
#   geom_line(data=c19, colour=col[11])+
#   geom_point(data=c20, colour=col[12], size=3) +
#   geom_line(data=c20, colour=col[12])+
#   geom_point(data=c22, colour=col[13], size=3) +
#   geom_line(data=c22, colour=col[13])+
#   geom_point(data=c23, colour=col[14], size=3) +
#   geom_line(data=c23, colour=col[14])+
#   theme_bw(20)+
#   ylab("VOR Amplitude (°)\n") +  
#   xlab("\nDevelopmental age (days)")
# #geom_errorbar(data=error, aes(ymin=Average-se, ymax=Average+se), width=.1)
# #scale_fill_continuous(guide=FALSE)
library(ggplot2)
ggplot(Daily.df, aes(x=Age, y=Average, group=Age)) + 
  geom_boxplot(data=Daily.df, size=1) +
  theme_bw(20)+
  geom_point(data=c8, colour=col[1], size=3) +
  geom_point(data=c9, colour=col[2], size=3) +
  geom_point(data=c10, colour=col[3], size=3) +
  geom_point(data=c11, colour=col[4], size=3) +
  geom_point(data=c12, colour=col[5], size=3) +
  geom_point(data=c13, colour=col[6], size=3) +
  geom_point(data=c14, colour=col[7], size=3) +
  geom_point(data=c15, colour=col[8], size=3) +
  geom_point(data=c16, colour=col[9], size=3) +
  geom_point(data=c17, colour=col[10], size=3) +
  geom_point(data=c19, colour=col[11], size=3) +
  geom_point(data=c20, colour=col[12], size=3) +
  geom_point(data=c22, colour=col[13], size=3) +
  geom_point(data=c23, colour=col[14], size=3) +
  theme_bw(20)+
  ylab("VOR amplitude (°)\n") +  
  xlab("\nDevelopmental age (days)")

stage2 <- subset(Daily.df, Stage == 2)
stage3 <- subset(Daily.df, Stage == 3)
stage4 <- subset(Daily.df, Stage == 4)
stage5 <- subset(Daily.df, Stage == 5)
stage6 <- subset(Daily.df, Stage == 6)
stage7 <- subset(Daily.df, Stage == 7)

ggplot(Daily.df, aes(x=Age, y=Average, group=Age)) + 
  geom_boxplot(data=Daily.df, size=1) +
  theme_bw(20)+
  geom_point(data=stage2, colour="red", size=3) +
  geom_point(data=stage3, colour="orange", size=3) +
  geom_point(data=stage4, colour="yellow", size=3) +
  geom_point(data=stage5, colour="green", size=3) +
  geom_point(data=stage6, colour="blue", size=3) +
  geom_point(data=stage7, colour="purple", size=3) +
  theme_bw(20)+
  ylab("VOR amplitude (°)\n") +  
  xlab("\nDevelopmental age (days)")

# Tukey test
library(MASS)
library(multcomp)
Daily.df$Age<-as.factor(Daily.df$Age)
lm1<-lm(Average~Age, data=Daily.df)
lm2<-glht(lm1, linfct=mcp(Age="Tukey"))
summary(lm2)
cld(lm2) 

VOR_asymptotic_test <- 
  with(Daily.df, 
       asymptotic_test(Average, 
                       Age))

VOR_mlrt_test <- 
  with(Daily.df, 
       mslr_test(nr = 1e4, 
                 Average, 
                 Age))

Daily.4v3567<-Daily.df
Daily.4v3567$Age<-as.numeric(Daily.4v3567$Age)
Daily.4v3567$Age[Daily.4v3567$Age!= 4]<-"other"

VOR_asymptotic_test <- 
  with(Daily.4v3567, 
       asymptotic_test(Average, 
                       Age))

VOR_mlrt_test <- 
  with(Daily.4v3567, 
       mslr_test(nr = 1e4, 
                 Average, 
                 Age))

ggplot(Daily.df, aes(x=as.factor(Stage), y=Average)) + 
  geom_boxplot(data=Daily.df, size=1)+
  theme_bw(20)+
  ylab("VOR amplitude (°)\n") +  
  xlab("\nDevelopmental stage")+
  geom_point(data=c8, colour=col[1], size=3) +
  geom_point(data=c9, colour=col[2], size=3) +
  geom_point(data=c10, colour=col[3], size=3) +
  geom_point(data=c11, colour=col[4], size=3) +
  geom_point(data=c12, colour=col[5], size=3) +
  geom_point(data=c13, colour=col[6], size=3) +
  geom_point(data=c14, colour=col[7], size=3) +
  geom_point(data=c15, colour=col[8], size=3) +
  geom_point(data=c16, colour=col[9], size=3) +
  geom_point(data=c17, colour=col[10], size=3) +
  geom_point(data=c19, colour=col[11], size=3) +
  geom_point(data=c20, colour=col[12], size=3) +
  geom_point(data=c22, colour=col[13], size=3) +
  geom_point(data=c23, colour=col[14], size=3)

sum(Daily.df$Stage==2) #3
sum(Daily.df$Stage==3) #4
sum(Daily.df$Stage==4) #0 
sum(Daily.df$Stage==5) #2
sum(Daily.df$Stage==6) #9
sum(Daily.df$Stage==7) #18

mean(Daily.df$Age) #4.75

###################################################################################################################################################################################################################################################################################
###################################################################################################################################################################################################################################################################################
###############        ontogeny of V O R                                    #######################################################################################################################################################################################################
###############     s i b s h i p s   developmental series #2               #######################################################################################################################################################################################################
###################################################################################################################################################################################################################################################################################
###################################################################################################################################################################################################################################################################################
###################################################################################################################################################################################################################################################################################

SibshipsSonia.df<-read.csv(file="SibshipsSonia.csv")
sibshiperrorstats <- summarySE(SibshipsSonia.df, measurevar="AverageAmp", groupvars="Time")
sibshiperrorstats

c101 <- subset(SibshipsSonia.df, Clutch == "101")
c102 <- subset(SibshipsSonia.df, Clutch == "102")
c104 <- subset(SibshipsSonia.df, Clutch == "104")
c105 <- subset(SibshipsSonia.df, Clutch == "105")
c106 <- subset(SibshipsSonia.df, Clutch == "106")

colour = c("red", "orange", "green", "blue", "purple")
ggplot(SibshipsSonia.df, aes(x=Time, y=AverageAmp, colour=Clutch)) + 
  geom_point(data=c101, colour = colour[1], size=5) +
  geom_line(data=c101, colour = colour[1])+
  geom_errorbar(data=c101, aes(ymin=AverageAmp-SE, ymax=AverageAmp+SE), width=0.05, colour = colour[1])+ 
  geom_point(data=c102, colour = colour[2], size=5) +
  geom_line(data=c102, colour = colour[2])+
  geom_errorbar(data=c102, aes(ymin=AverageAmp-SE, ymax=AverageAmp+SE), width=0.05, colour = colour[2])+ 
  geom_point(data=c104, colour = colour[3], size=5) +
  geom_line(data=c104, colour = colour[3])+
  geom_errorbar(data=c104, aes(ymin=AverageAmp-SE, ymax=AverageAmp+SE), width=0.05, colour = colour[3])+ 
  geom_point(data=c105, colour = colour[4], size=5) +
  geom_line(data=c105, colour = colour[4])+
  geom_errorbar(data=c105, aes(ymin=AverageAmp-SE, ymax=AverageAmp+SE), width=0.05, colour = colour[4])+ 
  geom_point(data=c106, colour = colour[5], size=5) +
  geom_line(data=c106, colour = colour[5])+
  geom_errorbar(data=c106, aes(ymin=AverageAmp-SE, ymax=AverageAmp+SE), width=0.05, colour = colour[5])+ 
  theme_bw(20)+
  #theme (panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())+
  #geom_errorbar(data = SibshipsSonia.df, aes(ymin=AverageAmp-SE, ymax=AverageAmp+SE), width = 0.05, colour = Clutch)+  
  ylab("VOR amplitude (°)\n") +  
  xlab("\nDevelopmental age (days)")

SibshipsSonia.df$Time<-as.factor(SibshipsSonia.df$Time)
lm3<-lm(AverageAmp ~ Time, data=SibshipsSonia.df)
lm4<-glht(lm3, linfct=mcp(Time="Tukey"))
summary(lm4)
cld(lm4) 

##COMPARING COEFFS OF VARIATION
# test for differences in variance for ages 4.25 and 4.5 d (vs all others)

A <-subset(SibshipsSonia.df, Time==4.00)
B <-subset(SibshipsSonia.df, Time==4.25)
C <-subset(SibshipsSonia.df, Time==4.50)
D <-subset(SibshipsSonia.df, Time==4.75)
E <-subset(SibshipsSonia.df, Time==5.75)

sibs <- data.frame(Test = c(4.25, 4.50, 4.75, 5.75),
                   Mean = c(mean(B$AverageAmp), mean(C$AverageAmp), mean(D$AverageAmp), mean(E$AverageAmp)),
                   CV =   c(sd(B$AverageAmp)/mean(B$AverageAmp), sd(C$AverageAmp)/mean(C$AverageAmp), sd(D$AverageAmp)/mean(D$AverageAmp),  sd(E$AverageAmp)/mean(E$AverageAmp)),
                   N =    c(length(B$AverageAmp), length(C$AverageAmp), length(D$AverageAmp), length(E$AverageAmp)))
            
sibs$SD <- with(sibs, CV * Mean)
                   
ggplot(sibs,
       aes(Test,
           Mean)) +
  # points to show mean values
  geom_point(size = 4) +
  # lines to show standard deviations
  geom_linerange(aes(ymin = Mean - SD,
                     ymax = Mean + SD)) +
  theme_bw()


Sibships_asymptotic_test2 <- 
  asymptotic_test2(k = nrow(sibs), 
                   n = sibs$N, 
                   s = sibs$SD, 
                   x = sibs$Mean)

Sibships_mlrt_test2 <-
  mslr_test2(nr = 1e4, 
             n = sibs$N, 
             s = sibs$SD, 
             x = sibs$Mean)


####################### compare variation at 4.25d vs. 4.50 and 4.75 d 

sibs <- data.frame(Test = c(4.25, 4.50, 4.75),
                   Mean = c(mean(B$AverageAmp), mean(C$AverageAmp), mean(D$AverageAmp)),
                   CV =   c(sd(B$AverageAmp)/mean(B$AverageAmp), sd(C$AverageAmp)/mean(C$AverageAmp), sd(D$AverageAmp)/mean(D$AverageAmp)),
                   N =    c(length(B$AverageAmp), length(C$AverageAmp), length(D$AverageAmp)))

sibs$SD <- with(sibs, CV * Mean)

ggplot(sibs,
       aes(Test,
           Mean)) +
  # points to show mean values
  geom_point(size = 4) +
  # lines to show standard deviations
  geom_linerange(aes(ymin = Mean - SD,
                     ymax = Mean + SD)) +
  theme_bw()


Sibships_asymptotic_test2 <- 
  asymptotic_test2(k = nrow(sibs), 
                   n = sibs$N, 
                   s = sibs$SD, 
                   x = sibs$Mean)

Sibships_mlrt_test2 <-
  mslr_test2(nr = 1e4, 
             n = sibs$N, 
             s = sibs$SD, 
             x = sibs$Mean)

####################### compare variation at 4.25 vs. 4.50 d  (REPORTED)
sibs <- data.frame(Test = c(4.25, 4.50),
                   Mean = c(mean(B$AverageAmp), mean(C$AverageAmp)),
                   CV =   c(sd(B$AverageAmp)/mean(B$AverageAmp), sd(C$AverageAmp)/mean(C$AverageAmp)),
                   N =    c(length(B$AverageAmp), length(C$AverageAmp)))

sibs$SD <- with(sibs, CV * Mean)

ggplot(sibs,
       aes(Test,
           Mean)) +
  # points to show mean values
  geom_point(size = 4) +
  # lines to show standard deviations
  geom_linerange(aes(ymin = Mean - SD,
                     ymax = Mean + SD)) +
  theme_bw()


Sibships_asymptotic_test2 <- 
  asymptotic_test2(k = nrow(sibs), 
                   n = sibs$N, 
                   s = sibs$SD, 
                   x = sibs$Mean)

Sibships_mlrt_test2 <-
  mslr_test2(nr = 1e4, 
             n = sibs$N, 
             s = sibs$SD, 
             x = sibs$Mean)

####################### compare variation at 4.50 vs. 4.75 d 
sibs <- data.frame(Test = c(4.50, 4.75),
                   Mean = c(mean(C$AverageAmp), mean(D$AverageAmp)),
                   CV =   c(sd(C$AverageAmp)/mean(C$AverageAmp), sd(D$AverageAmp)/mean(D$AverageAmp)),
                   N =    c(length(C$AverageAmp), length(D$AverageAmp)))

sibs$SD <- with(sibs, CV * Mean)

ggplot(sibs,
       aes(Test,
           Mean)) +
  # points to show mean values
  geom_point(size = 4) +
  # lines to show standard deviations
  geom_linerange(aes(ymin = Mean - SD,
                     ymax = Mean + SD)) +
  theme_bw()


Sibships_asymptotic_test2 <- 
  asymptotic_test2(k = nrow(sibs), 
                   n = sibs$N, 
                   s = sibs$SD, 
                   x = sibs$Mean)

Sibships_mlrt_test2 <-
  mslr_test2(nr = 1e4, 
             n = sibs$N, 
             s = sibs$SD, 
             x = sibs$Mean)

###### test for clutch effects


library(MASS)
library(pscl)

#m1<-hurdle(AverageAmp ~ Time, data=SibshipsSonia.df,
#       dist = "poisson", zero.dist = "binomial", link = "logit",
#       model = TRUE, y = TRUE, x = FALSE)

glm1 <- glm(AverageAmp ~ Time, family= poisson(link="log"), data=SibshipsSonia.df)
glm2 <- glm(AverageAmp ~ Clutch, family= poisson(link="log"), data=SibshipsSonia.df)
glm3 <- glm(AverageAmp ~ Time + Clutch, family= poisson(link="log"), data=SibshipsSonia.df)
glm4 <- glm(AverageAmp ~ Time * Clutch, family= poisson(link="log"), data=SibshipsSonia.df)

library("AICcmodavg")
VORmods<-list(glm1, glm2, glm3, glm4)
aictab(VORmods)

library(car)
Anova(glm2)
Anova(glm3) 
Anova(glm4) 

### test for clutch caused variation

hist(SibshipsSonia.df$AverageAmp, breaks=100) #non-normal
hist(log(SibshipsSonia.df$AverageAmp))#still not normal
shapiro.test(SibshipsSonia.df$AverageAmp) #not normal
# P<0.05, therefore not normal so use non-parametric test

kruskal.test(AverageAmp ~ Clutch, data= SibshipsSonia.df)
# no clutch effect for all ages
kruskal.test(AverageAmp ~ Time, data= SibshipsSonia.df)
# yes time/age effect for all ages

cor.test(SibshipsSonia.df$AverageAmp, SibshipsSonia.df$Clutch, method="pearson")
#r=0.09538258 
0.09538258*0.09538258
#The fraction of the variance in VOR that is “explained” by clutch is r2 = 0.009097837.


greatest.variability <- subset(SibshipsSonia.df, SibshipsSonia.df$Time==4.25|SibshipsSonia.df$Time==4.5)
kruskal.test(AverageAmp ~ Clutch, data= greatest.variability)
# no clutch effect for age 4.25 an 4.5
kruskal.test(AverageAmp ~ Time, data= greatest.variability)
# yes time/age effect for ages 4.25 and 4.5

cor.test(greatest.variability$AverageAmp, greatest.variability$Clutch, method="pearson")
#r=0.606184 
0.606184 * 0.606184 
#The fraction of the variance in VOR that is “explained” by clutch is r2 = 0.367459. 


four.twofive <- subset(SibshipsSonia.df, SibshipsSonia.df$Time==4.25)
kruskal.test(AverageAmp ~ Clutch, data= four.twofive)
# no clutch effect for age 4.25 
cor.test(four.twofive$AverageAmp, four.twofive$Clutch, method="pearson")
#r=0.7107471 
0.7107471  * 0.7107471 
#The fraction of the variance in VOR that is “explained” by clutch is r2 = 0.5051614 


four.five<- subset(SibshipsSonia.df, SibshipsSonia.df$Time==4.5)
kruskal.test(AverageAmp ~ Clutch, data= four.five)
# no clutch effect for age 4.5
cor.test(four.five$AverageAmp, four.five$Clutch, method="pearson")
#r=0.9255137 
0.9255137   * 0.9255137  
#The fraction of the variance in VOR that is “explained” by clutch is r2 = 0.8565756.  


ggplot(SibshipsSonia.df, aes(x=Stage, y=AverageAmp, colour=Clutch)) + 
  geom_point(data=c101, colour = colour[1], size=5) +
  geom_line(data=c101, colour = colour[1])+
  geom_errorbar(data=c101, aes(ymin=AverageAmp-SE, ymax=AverageAmp+SE), width=0.05, colour = colour[1])+ 
  geom_point(data=c102, colour = colour[2], size=5) +
  geom_line(data=c102, colour = colour[2])+
  geom_errorbar(data=c102, aes(ymin=AverageAmp-SE, ymax=AverageAmp+SE), width=0.05, colour = colour[2])+ 
  geom_point(data=c104, colour = colour[3], size=5) +
  geom_line(data=c104, colour = colour[3])+
  geom_errorbar(data=c104, aes(ymin=AverageAmp-SE, ymax=AverageAmp+SE), width=0.05, colour = colour[3])+ 
  geom_point(data=c105, colour = colour[4], size=5) +
  geom_line(data=c105, colour = colour[4])+
  geom_errorbar(data=c105, aes(ymin=AverageAmp-SE, ymax=AverageAmp+SE), width=0.05, colour = colour[4])+ 
  geom_point(data=c106, colour = colour[5], size=5) +
  geom_line(data=c106, colour = colour[5])+
  geom_errorbar(data=c106, aes(ymin=AverageAmp-SE, ymax=AverageAmp+SE), width=0.05, colour = colour[5])+ 
  theme_bw(20)+
  #theme (panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())+
  #geom_errorbar(data = SibshipsSonia.df, aes(ymin=AverageAmp-SE, ymax=AverageAmp+SE), width = 0.05, colour = Clutch)+  
  ylab("VOR amplitude (°)\n") +  
  xlab("\nDevelopmental stage")

###################################################################################################################################################################################################################################################################################
###################################################################################################################################################################################################################################################################################
###############                                                ####################################################################################################################################################################################################################
###############        S H A K E    n'    R O L L              ####################################################################################################################################################################################################################
###############                                                ####################################################################################################################################################################################################################
###################################################################################################################################################################################################################################################################################
###################################################################################################################################################################################################################################################################################
###################################################################################################################################################################################################################################################################################

#ShakeN'RollV1-Kwork.xlsx is the final excel file

ShakeRoll.df<-read.csv(file="ShakeRoll.csv")

ShakeRoll.df$HorNH<-as.numeric(ShakeRoll.df$HorNH)

stage1 <- subset(ShakeRoll.df, ShakeRoll.df$SUM.of.trait.values==1)
stage1.NH <- subset(stage1, stage1$HorNH==0)
stage1.H <- subset(stage1, stage1$HorNH==1)

stage2 <- subset(ShakeRoll.df, ShakeRoll.df$SUM.of.trait.values==2)
stage2.NH <- subset(stage2, stage2$HorNH==0)
stage2.H <- subset(stage2, stage2$HorNH==1)

stage3 <- subset(ShakeRoll.df, ShakeRoll.df$SUM.of.trait.values==3)
stage3.NH <- subset(stage3, stage3$HorNH==0)
stage3.H <- subset(stage3, stage3$HorNH==1)

stage4 <- subset(ShakeRoll.df, ShakeRoll.df$SUM.of.trait.values==4)
stage4.NH <- subset(stage4, stage4$HorNH==0)
stage4.H <- subset(stage4, stage4$HorNH==1)

stage5 <- subset(ShakeRoll.df, ShakeRoll.df$SUM.of.trait.values==5)
stage5.NH <- subset(stage5, stage5$HorNH==0)
stage5.H <- subset(stage5, stage5$HorNH==1)

stage6 <- subset(ShakeRoll.df, ShakeRoll.df$SUM.of.trait.values==6)
stage6.NH <- subset(stage6, stage6$HorNH==0)
stage6.H <- subset(stage6, stage6$HorNH==1)

stage7 <- subset(ShakeRoll.df, ShakeRoll.df$SUM.of.trait.values==7)
stage7.NH <- subset(stage7, stage7$HorNH==0)
stage7.H <- subset(stage7, stage7$HorNH==1)

younger <- subset(ShakeRoll.df, ShakeRoll.df$SUM.of.trait.values<5)
length(younger$AverageAmp) #N=18
sum(younger$PropH) #none hatched
mean(younger$AverageAmp) #average VOR was 2.028089
mean(younger$AverageR2) #average R2 was 0.2722102
sd(younger$AverageAmp) 
sd(younger$AverageR2)

ggplot(ShakeRoll.df, aes(x=as.factor(SUM.of.trait.values), y=PropH)) + 
  geom_boxplot(data=ShakeRoll.df, size=1, na.rm=T) +
  ylab("Proportion of clutch hatched\n") +  
  theme_bw(20) +
  xlab("\n Developmental stage")

# Tukey test PropH
ShakeRoll.df$SUM.of.trait.values  <-as.factor(ShakeRoll.df$SUM.of.trait.values  )
lm5<-lm(PropH~SUM.of.trait.values  , data=ShakeRoll.df)
lm6<-glht(lm5, linfct=mcp(SUM.of.trait.values  ="Tukey"))
summary(lm6)
cld(lm6) 

ggplot(ShakeRoll.df, aes(x=as.factor(SUM.of.trait.values), y=AverageAmp)) + 
  geom_boxplot(data=ShakeRoll.df, size=1, na.rm=T) +
  ylab("VOR amplitude (°)\n") + 
  theme_bw(20) +
  xlab("\n Developmental stage")

# Tukey test AverageAmp
lm7<-lm(AverageAmp~SUM.of.trait.values, data=ShakeRoll.df)
lm8<-glht(lm7, linfct=mcp(SUM.of.trait.values  ="Tukey"))
summary(lm8)
cld(lm8) 

#H vs. did not H at different stages
stage5<- subset(ShakeRoll.df, ShakeRoll.df$SUM.of.trait.values==5)
length(stage5)
stage5.NH <-subset(stage5, ShakeRoll.df$HorNH==0)
length(stage5.NH)
stage5.H <-subset(stage5, ShakeRoll.df$HorNH==1)

lm9<-lm(AverageAmp~HorNH, data=stage5)
lm10<-glht(lm9, linfct=mcp(HorNH  ="Tukey"))
summary(lm10)

lm11<-lm(PropH~HorNH, data=stage5)
lm12<-glht(lm11, linfct=mcp(HorNH  ="Tukey"))
summary(lm12)

stage6<- subset(ShakeRoll.df, ShakeRoll.df$SUM.of.trait.values==6)

lm13<-lm(AverageAmp~HorNH, data=stage6)
lm14<-glht(lm13, linfct=mcp(HorNH  ="Tukey"))
summary(lm14)

lm15<-lm(PropH~HorNH, data=stage6)
lm16<-glht(lm15, linfct=mcp(HorNH  ="Tukey"))
summary(lm16)

stage7<- subset(ShakeRoll.df, ShakeRoll.df$SUM.of.trait.values==7)

lm17<-lm(AverageAmp~HorNH, data=stage7)
lm18<-glht(lm17, linfct=mcp(HorNH  ="Tukey"))
summary(lm18)

lm19<-lm(PropH~HorNH, data=stage7)
lm20<-glht(lm19, linfct=mcp(HorNH  ="Tukey"))
summary(lm20)

# 2 boxplots for figure
ShakeRoll.df$HorNH <- as.factor(ShakeRoll.df$HorNH)

ggplot(ShakeRoll.df, aes(x=as.factor(SUM.of.trait.values), y=PropH, color = HorNH)) + 
  geom_boxplot(data=ShakeRoll.df, size=1, na.rm=T) +
  ylab("Proportion of clutch hatched\n") +  
  theme_bw(20) +
  xlab("\n Developmental stage")

plot<-ggplot(ShakeRoll.df, aes(x=as.factor(SUM.of.trait.values), y=AverageAmp, color = HorNH)) + 
  geom_boxplot(data=ShakeRoll.df, size=1) +
  ylab("VOR amplitude (°)\n")+
  theme_bw(20) +
  xlab("\n Developmental stage")

ggplot_build(plot)$data
#F8766D
#00BFC4

#Correlation tests:
library(psych)
ShakeRoll.df$SUM.of.trait.values<-as.numeric(ShakeRoll.df$SUM.of.trait.values)
cor.test(ShakeRoll.df$AverageAmp, ShakeRoll.df$SUM.of.trait.values, method="pearson")
cor.test(ShakeRoll.df$PropH, ShakeRoll.df$SUM.of.trait.values, method="pearson")

## Subset by older embryos (take out Stage 1-4)
older<- subset(ShakeRoll.df, ShakeRoll.df$SUM.of.trait.values>4)
cor.test(older$AverageAmp, older$SUM.of.trait.values, method="pearson")
cor.test(older$PropH, older$SUM.of.trait.values, method="pearson")


## Subset by partially hatched clutches (take out PropH==0 and PropH==1)
parthatch<-subset(ShakeRoll.df, ShakeRoll.df$PropH<1 & ShakeRoll.df$PropH>0)
cor.test(parthatch$AverageAmp, parthatch$SUM.of.trait.values, method="pearson")
cor.test(parthatch$PropH, parthatch$SUM.of.trait.values, method="pearson")

## signif difference in VOR between H and NH??
# in stages 5-7:
lm21<-lm(AverageAmp~HorNH, data=older)
lm22<-glht(lm21, linfct=mcp(HorNH  ="Tukey"))
summary(lm22) ## no significant differences in VOR

# in partially hatched clutches:
lm23<-lm(AverageAmp~HorNH, data=parthatch)
lm24<-glht(lm23, linfct=mcp(HorNH  ="Tukey"))
summary(lm24)


# variation in stage 5 vs variation in other stages. 

VOR_asymptotic_test <- 
  with(older, 
       asymptotic_test(AverageAmp, 
                       SUM.of.trait.values))

VOR_mlrt_test <- 
  with(older, 
       mslr_test(nr = 1e4, 
                 AverageAmp, 
                 SUM.of.trait.values))

propH_asymptotic_test <- 
  with(older, 
       asymptotic_test(PropH, 
                       SUM.of.trait.values))

propH_mlrt_test <- 
  with(older, 
       mslr_test(nr = 1e4, 
                 PropH, 
                 SUM.of.trait.values))


# Does hatching inc with vestibular function?
cor.test(ShakeRoll.df$AverageAmp, ShakeRoll.df$PropH, method="pearson")

## in older embryos (take out Stage 1-4)
cor.test(older$AverageAmp, older$PropH, method="pearson")

## in partially hatched clutches (take out PropH==0 and PropH==1)
cor.test(parthatch$AverageAmp, parthatch$PropH, method="pearson")


## shake and roll, colored according to developmental stage. 

ShakeRoll.df$HorNH <- factor(ShakeRoll.df$HorNH)
str(ShakeRoll.df$HorNH)
colour = c("red", "orange", "yellow", "green", "blue", "purple", "brown")
ggplot(ShakeRoll.df, aes(x=AverageAmp, y=PropH, colour=SUM.of.trait.values)) + 
  geom_point(data=stage1, colour = colour[1], size=5) +
  geom_point(data=stage2, colour = colour[2], size=5) +
  geom_point(data=stage3, colour = colour[3], size=5) +
  geom_point(data=stage4, colour = colour[4], size=5) +
  geom_point(data=stage5, colour = colour[5], size=5) +
  geom_point(data=stage6, colour = colour[6], size=5) +
  geom_point(data=stage7, colour = colour[7], size=5) +
  theme_bw(20)+
  #theme (panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())+
  #geom_errorbar(data = SibshipsSonia.df, aes(ymin=AverageAmp-SE, ymax=AverageAmp+SE), width = 0.05, colour = Clutch)+  
  ylab("Proportion of clutch hatched\n") +  
  xlab("VOR amplitude (°)\n")

ggplot(ShakeRoll.df, aes(x=SUM.of.trait.values, y=AverageAmp, colour=HorNH)) + 
  geom_point(size=3, position = "jitter") +
  theme_bw(20)+
  #theme (panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())+
  #geom_errorbar(data = SibshipsSonia.df, aes(ymin=AverageAmp-SE, ymax=AverageAmp+SE), width = 0.05, colour = Clutch)+  
  ylab("VOR amplitude (°)\n")+
  xlab("\nDevelopmental stage")

Hat<- subset(ShakeRoll.df, HorNH == "1")
NoHat<- subset(ShakeRoll.df, HorNH == "0")

ggplot(ShakeRoll.df, aes(x=AverageAmp, y=PropH, colour=HorNH)) + 
  geom_point(data=Hat, colour = "#00BFC4", size=3) +
  geom_point(data=NoHat, colour = "#F8766D", size=3) +
  theme_bw(20)+
  #theme (panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())+
  #geom_errorbar(data = SibshipsSonia.df, aes(ymin=AverageAmp-SE, ymax=AverageAmp+SE), width = 0.05, colour = Clutch)+  
  ylab("Proportion of clutch hatched\n") +  
  xlab("VOR amplitude (°)\n")

###################################################################################################################################################################################################################################################################################
###################################################################################################################################################################################################################################################################################
###############                                                ####################################################################################################################################################################################################################
###############        J I G G L E   n'   R O L L              ####################################################################################################################################################################################################################
###############                                                ####################################################################################################################################################################################################################
###################################################################################################################################################################################################################################################################################
###################################################################################################################################################################################################################################################################################
###################################################################################################################################################################################################################################################################################

# read.csv(file="my.csv.filename")
VOR.df<-read.csv(file="2cues2015.csv")
tactile <- subset(VOR.df, Stimulus == "T", na.rm=T, select=c(Clutch, Age.Block, Developmental.Stage, Individual, Response, Average.R2.Value, Average.Amplitude, Time.of.First.Hatch, Age.Hatched, Hours.Since.First.Hatch, Latency.to.Hatch.in.Minutes))

hist(log(tactile$Average.Amplitude)) #log distribution is normal-ish

shapiro.test(log(tactile$Average.Amplitude)) #not normal
# P<0.05, therefore not normal so use non-parametric test (like JT test)

jonckheere.test(as.numeric(tactile$Average.Amplitude), as.numeric(tactile$Age.Block), nperm=1000, alternative="increasing")
jonckheere.test(as.numeric(tactile$Average.Amplitude), as.numeric(tactile$Developmental.Stage), nperm=1000, alternative="increasing")

cor.test(as.numeric(tactile$Average.Amplitude), as.numeric(tactile$Age.Block), method="pearson")
cor.test(as.numeric(tactile$Average.Amplitude), as.numeric(tactile$Developmental.Stage), method="pearson")

glm1 <- glm(Average.Amplitude ~ Age.Block, data=tactile)
glm2 <- glm(Average.Amplitude ~ Developmental.Stage, data=tactile)
glm3 <- glm(Average.Amplitude ~ Clutch, data=tactile)
glm4 <- glm(Average.Amplitude ~ Developmental.Stage + Age.Block, data=tactile)
glm5 <- glm(Average.Amplitude ~ Developmental.Stage * Age.Block, data=tactile)

glms<-list(glm1, glm2, glm3, glm4, glm5)
aictab(glms)

library(car)
Anova(glm1)
Anova(glm2)
Anova(glm3)
Anova(glm4)
Anova(glm5)

# PROP H of EACH STAGE: 

st1 <- subset(tactile, Developmental.Stage == "1", na.rm=T) #0
st1.NH <- subset(st1, st1$Response==0) #0
st1.H <- subset(st1, st1$Response==1) #0
st2 <- subset(tactile, Developmental.Stage == "2", na.rm=T) #4
st2.NH <- subset(st2, st2$Response==0) #4
st2.H <- subset(st2, st2$Response==1) #0
st3 <- subset(tactile, Developmental.Stage == "3", na.rm=T) #29
st3.NH <- subset(st3, st3$Response==0) #24
st3.H <- subset(st3, st3$Response==1) #5
st4 <- subset(tactile, Developmental.Stage == "4", na.rm=T) #19
st4.NH <- subset(st4, st4$Response==0) #10
st4.H <- subset(st4, st4$Response==1) #9
st5 <- subset(tactile, Developmental.Stage == "5", na.rm=T) #24
st5.NH <- subset(st5, st5$Response==0) #6
st5.H <- subset(st5, st5$Response==1) #18
st6 <- subset(tactile, Developmental.Stage == "6", na.rm=T) #33
st6.NH <- subset(st6, st6$Response==0) #7
st6.H <- subset(st6, st6$Response==1) #26
st7 <- subset(tactile, Developmental.Stage == "7", na.rm=T) #3
st7.NH <- subset(st7, st7$Response==0) #0
st7.H <- subset(st7, st7$Response==1) #3

hat <- subset(tactile, Response == "1") #61
nohat <- subset(tactile, Response == "0") #51

sum(st3$Response==0)
sum(st3$Response==1)
#propH (of st3) = 5/(5+24) = 0.2083333

sum(st4$Response==1)
sum(st4$Response==2)
#propH (of st4) = 9/(9+10) = 0.4736842

sum(st5$Response==1)
sum(st5$Response==2)
#propH (of st5) = 18/(18+6) = 0.75

sum(st6$Response==1)
sum(st6$Response==2)
#propH (of st6) = 26/(26+7) = 0.7878788

sum(st7$Response==1)
sum(st7$Response==2)
#propH (of st7) = 3/(3+0) = 1

stage<-c(2, 3, 4, 5, 6, 7)
N<-c(4, 29, 19, 24, 33, 3)
propH<-c(0, 0.2083333, 0.4736842, 0.75, 0.7878788, 1)
df<- data.frame (stage, N, propH)

#plot of propH x stage
ggplot(df, aes(x=as.factor(stage), y=propH)) + 
  geom_boxplot(data=df, size=1) +
  ylab("Proportion hatched\n")+
  theme_bw(20) +
  xlab("\n Developmental stage")


#plot of VOR x stage (H vs. NH **NOT** separated)
tactile$Response<-as.factor(tactile$Response)
ggplot(tactile, aes(x=as.factor(Developmental.Stage), y=Average.Amplitude)) + 
  geom_boxplot(data=tactile, size=1) +
  ylab("VOR amplitude (°)\n")+
  theme_bw(20) +
  xlab("\n Developmental stage")+
  theme(legend.position="none")

##### TUKEY TEST #####
# Tukey test VOR
tactile$Developmental.Stage<-as.factor(tactile$Developmental.Stage)
lm1<-lm(Average.Amplitude~Developmental.Stage, data=tactile)
lm2<-glht(lm1, linfct=mcp(Developmental.Stage ="Tukey"))
summary(lm2)
cld(lm2) 


# in stage 4, 5, 6
stgs456 <-subset (tactile, tactile$Developmental.Stage==4|tactile$Developmental.Stage==5|tactile$Developmental.Stage==6)

# variation in stage 5 vs variation in 4 and 6. 

stgs456$Average.Amplitude<- as.numeric(stgs456$Average.Amplitude)
stgs456$Developmental.Stage<- as.numeric(stgs456$Developmental.Stage)

VOR_asymptotic_test <- 
  with(stgs456, 
       asymptotic_test(Average.Amplitude, 
                       Developmental.Stage))

VOR_mlrt_test <- 
  with(stgs456, 
       mslr_test(nr = 1e4, 
                 Average.Amplitude, 
                 Developmental.Stage))


# it looks like across stages 5-6 the hatched ones have higher VOR but not stage 3-4

# Tukey test at stage 3, H vs. NH. 
st3$Response<-as.factor(st3$Response)
lm1<-lm(Average.Amplitude~Response, data=st3)
lm2<-glht(lm1, linfct=mcp(Response ="Tukey"))
summary(lm2)
cld(lm2) 

# Tukey test at stage 4, H vs. NH. 
st4$Response<-as.factor(st4$Response)
lm1<-lm(Average.Amplitude~Response, data=st4)
lm2<-glht(lm1, linfct=mcp(Response ="Tukey"))
summary(lm2)
cld(lm2) 

# Tukey test at stage 5, H vs. NH. 
st5$Response<-as.factor(st5$Response)
lm1<-lm(Average.Amplitude~Response, data=st5)
lm2<-glht(lm1, linfct=mcp(Response ="Tukey"))
summary(lm2)
cld(lm2) 

# Tukey test at stage 6, H vs. NH. 
st6$Response<-as.factor(st6$Response)
lm1<-lm(Average.Amplitude~Response, data=st6)
lm2<-glht(lm1, linfct=mcp(Response ="Tukey"))
summary(lm2)
cld(lm2) 



hat <- subset(tactile, Response == "1") #61
nohat <- subset(tactile, Response == "0") #51

# test for stage & VOR effects & interactions?
# an analysis of hatching looking at stage & VOR as predictors 

#stage as predictor of hatching
hist(tactile$Developmental.Stage)
boxplot(hat$Developmental.Stage, nohat$Developmental.Stage, xlab="Response", ylab="Developmental Stage")
axis(1, at=1:2, labels=c("Hatched", "Not Hatched"))
t.test(as.numeric(hat$Developmental.Stage), as.numeric(nohat$Developmental.Stage), mu = 0, alt="two.sided", paired = F, conf.int=T, conf.level=0.99)
## wilcox.test(as.numeric(hat$Developmental.Stage), as.numeric(nohat$Developmental.Stage), mu = 0, alt="two.sided", paired = F, conf.int=T, conf.level=0.99)

#VOR as predictor of hatching
hist(tactile$Average.Amplitude)
boxplot(hat$Average.Amplitude, nohat$Average.Amplitude, xlab="Response", ylab="VOR Amplitude")
axis(1, at=1:2, labels=c("Hatched", "Not Hatched"))
t.test(as.numeric(hat$Developmental.Stage), as.numeric(nohat$Developmental.Stage), mu = 0, alt="two.sided", paired = F, conf.int=T, conf.level=0.99)
wtest<-wilcox.test(as.numeric(hat$Average.Amplitude), as.numeric(nohat$Average.Amplitude), mu = 0, alt="two.sided", paired = F, conf.int=T, conf.level=0.99)
qnorm(wtest$p.value) #z-value to report

## We should look at the rlp between latency & VOR. Have you done that?]

#plot of latency x stage
ggplot(tactile, aes(x=as.factor(Developmental.Stage), y=Latency.to.Hatch.in.Minutes)) + 
  geom_boxplot(data=tactile, size=1) +
  ylab("Latency to hatch (mins)\n")+
  theme_bw(20) +
  xlab("\n Developmental stage")

#plot of latency x VOR 

# redefine st1 etc. 
st1 <- subset(tactile, Developmental.Stage == "1", na.rm=T) #0
st2 <- subset(tactile, Developmental.Stage == "2", na.rm=T) #4
st3 <- subset(tactile, Developmental.Stage == "3", na.rm=T) #29
st4 <- subset(tactile, Developmental.Stage == "4", na.rm=T) #19
st5 <- subset(tactile, Developmental.Stage == "5", na.rm=T) #24
st6 <- subset(tactile, Developmental.Stage == "6", na.rm=T) #33
st7 <- subset(tactile, Developmental.Stage == "7", na.rm=T) #3

colour = c("black", "red", "orange", "yellow", "green", "blue", "purple")

str(tactile$Average.Amplitude)
str(tactile$Latency.to.Hatch.in.Minutes)
range(tactile$Latency.to.Hatch.in.Minutes, na.rm=T)
range(tactile$Average.Amplitude, na.rm=T)

ggplot(tactile, aes(x=Average.Amplitude, y=Latency.to.Hatch.in.Minutes, color=Developmental.Stage)) + 
  geom_point(data=st1, colour = colour[1], size=3) +
  geom_point(data=st2, colour = colour[2], size=3) +
  geom_point(data=st3, colour = colour[3], size=3) +
  geom_point(data=st4, colour = colour[4], size=3) +
  geom_point(data=st5, colour = colour[5], size=3) +
  geom_point(data=st6, colour = colour[6], size=3) +
  geom_point(data=st7, colour = colour[7], size=3) +
  scale_x_continuous(breaks= c(0.0,10, 20.0,30.0,40.0,50.0), labels = c("0", "10","20","30", "40","50"))+
  ylab("Latency to hatch (mins)\n")+
  theme_bw(20) +
  xlab("\n VOR amplitude (°)")+
  theme(legend.position="none")

cor.test(tactile$Average.Amplitude, tactile$Latency.to.Hatch.in.Minutes, method="pearson")
cor.test(tactile$Developmental.Stage, tactile$Latency.to.Hatch.in.Minutes, method="pearson")

### start here







glm1 <- glm(Response ~ Developmental.Stage, data=tactile)
glm2 <- glm(Response ~ Average.Amplitude, data=tactile)
glm3 <- glm(Average.Amplitude ~ Clutch, data=tactile)
glm4 <- glm(Average.Amplitude ~ Developmental.Stage + Age.Block, data=tactile)
glm5 <- glm(Average.Amplitude ~ Developmental.Stage * Age.Block, data=tactile)

glms<-list(glm1, glm2, glm3, glm4, glm5)
aictab(glms)

library(car)
Anova(glm1)
Anova(glm2)
Anova(glm3)
Anova(glm4)
Anova(glm5)


colour = c("red", "orange", "yellow", "green", "blue", "purple", "brown")
tactile$Response<-as.factor(tactile$Response)
str(tactile$Response)
ggplot(tactile, aes(x=Developmental.Stage, y=Average.Amplitude, colour=Developmental.Stage, shape=Response)) + 
  geom_point(data=st1, colour = colour[1], size=5) +
  geom_point(data=st2, colour = colour[2], size=5) +
  geom_point(data=st3, colour = colour[3], size=5) +
  geom_point(data=st4, colour = colour[4], size=5) +
  geom_point(data=st5, colour = colour[5], size=5) +
  geom_point(data=st6, colour = colour[6], size=5) +
  geom_point(data=st7, colour = colour[7], size=5) +
  #geom_point(data=hat, size=5, pch=1) +
  #geom_point(data=nohat, size=5, pch=19) + 
  theme_bw(20)+
  #theme (panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())+
  #geom_errorbar(data = SibshipsSonia.df, aes(ymin=AverageAmp-SE, ymax=AverageAmp+SE), width = 0.05, colour = Clutch)+  
  ylab("VOR amplitude (°)\n")+ 
  xlab("\nHours since first sibling hatched") 


# BIN into YES and NO VOR: 
noVOR <- subset(tactile, Average.Amplitude < 5)
yesVOR <- subset(tactile, Average.Amplitude >= 5)

sum(noVOR$Response=="1") #hatch
sum(noVOR$Response=="2") #nh
#propH (of noVOR) = 4/(4+17) = 0.1904762
# 1-0.1904762 =0.8095238

sum(yesVOR$Response=="1")
sum(yesVOR$Response=="2")
#propH (of yesVOR) = 57/(57+32) = 0.6404494

tactile$Response<-as.factor(tactile$Response)
#plot of VOR x stage (H vs. NH separated)
ggplot(tactile, aes(x=as.factor(Developmental.Stage), y=Average.Amplitude, color=Response)) + 
  geom_boxplot(data=tactile, size=1) +
  ylab("VOR amplitude (°)\n")+
  theme_bw(20) +
  xlab("\n Developmental stage")+
  theme(legend.position="none")

ggplot(tactile, aes(x=Developmental.Stage, y=Average.Amplitude, colour=Developmental.Stage, shape=Response)) + 
  geom_point(data=st1, colour = colour[1], size=5) +
  geom_point(data=st2, colour = colour[2], size=5) +
  geom_point(data=st3, colour = colour[3], size=5) +
  geom_point(data=st4, colour = colour[4], size=5) +
  geom_point(data=st5, colour = colour[5], size=5) +
  geom_point(data=st6, colour = colour[6], size=5) +
  geom_point(data=st7, colour = colour[7], size=5) +
  #geom_point(data=hat, size=5, pch=1) +
  #geom_point(data=nohat, size=5, pch=19) + 
  theme_bw(20)+
  #theme (panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())+
  #geom_errorbar(data = SibshipsSonia.df, aes(ymin=AverageAmp-SE, ymax=AverageAmp+SE), width = 0.05, colour = Clutch)+  
  ylab("VOR amplitude (°)\n")+ 
  xlab("\nDevelopmental Stage") 


library(ggplot2)

# ggplot(tactile, aes(x=TtoH, y=Average.Amp)) + geom_point(shape=1) + xlim(c(0,0.1))
#   #as VOR increases, they take a shorter time to hatch. 
# 
# # hatched vs. nonhatched #as 
# ggplot(tactile, aes(x=AgeBlock, y=Average.Amp, color=Response)) + geom_point(shape=1) + scale_colour_hue(l=50)
# #as VOR increases, they tend to hatch (instead of not hatch)
# 
# # trendline through all points
# ggplot(tactile, aes(x=AgeBlock, y=Average.Amp)) + geom_point(shape=1) + scale_colour_hue(l=50) + geom_smooth() 
# 
# # by clutch, corrected for time hatched
# ggplot(tactile, aes(x=HsinceH, y=Average.Amp, color=Clutch)) + geom_point(shape=1) + geom_smooth() +
#   ylab("Average Amplitude") +
#   xlab("Hours Since First Hatch")
# 

#############  MAIN FIGURE !!! ####################
# hatched bs. nonhatched, corrected for time hatched 
ggplot(tactile, aes(x=Hours.Since.First.Hatch, y=Average.Amplitude, color=Response)) + 
  geom_point(shape=16, size=4) +
  scale_shape_discrete(solid=T) +
  ylab("VOR amplitude (°)\n") +
  xlab("\nHours since first sibling hatched") + 
  #          theme(axis.title.x = element_text(vjust=-0.8)) +
  #          theme(axis.title.y = element_text(vjust=0.8)) +
  theme_bw(20) +
  theme(panel.grid.major = element_line(colour = "#808080"))+
  scale_colour_manual(values= c("forestgreen", "goldenrod"))

mean(tactile$Age.Block)
min(tactile$Age.Block)
max(tactile$Age.Block)

st1 <- subset(tactile, Developmental.Stage == "1", na.rm=T)
st2 <- subset(tactile, Developmental.Stage == "2", na.rm=T)
st3 <- subset(tactile, Developmental.Stage == "3", na.rm=T)
st4 <- subset(tactile, Developmental.Stage == "4", na.rm=T)
st5 <- subset(tactile, Developmental.Stage == "5", na.rm=T)
st6 <- subset(tactile, Developmental.Stage == "6", na.rm=T)
st7 <- subset(tactile, Developmental.Stage == "7", na.rm=T)

hat <- subset(tactile, Response == "Hatched")
nohat <- subset(tactile, Response == "Not Hatched")

colour = c("red", "orange", "yellow", "green", "blue", "purple", "brown")
tactile$Response<-as.factor(tactile$Response)
str(tactile$Response)
ggplot(tactile, aes(x=Hours.Since.First.Hatch, y=Average.Amplitude, colour=Developmental.Stage, shape=Response)) + 
  geom_point(data=st1, colour = colour[1], size=5) +
  geom_point(data=st2, colour = colour[2], size=5) +
  geom_point(data=st3, colour = colour[3], size=5) +
  geom_point(data=st4, colour = colour[4], size=5) +
  geom_point(data=st5, colour = colour[5], size=5) +
  geom_point(data=st6, colour = colour[6], size=5) +
  geom_point(data=st7, colour = colour[7], size=5) +
  #geom_point(data=hat, size=5, pch=1) +
  #geom_point(data=nohat, size=5, pch=19) + 
  theme_bw(20)+
  #theme (panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())+
  #geom_errorbar(data = SibshipsSonia.df, aes(ymin=AverageAmp-SE, ymax=AverageAmp+SE), width = 0.05, colour = Clutch)+  
  ylab("VOR amplitude (°)\n")+ 
  xlab("\nHours since first sibling hatched") 

############ ANOVA Stats ###############
hist(tactile$Average.Amplitude)
hist()

m_1 <- aov(Average.Amplitude ~ 1, data=tactile)
m_2 <-aov(Average.Amplitude ~ Response, data = tactile)
m_3 <- aov(Average.Amplitude ~ Clutch, data = tactile)
m_4 <- aov(Average.Amplitude ~ Clutch*Response, data = tactile)
m_5 <- aov(Average.Amplitude ~ Response*Hours.Since.First.Hatch, data=tactile)
#m_6 <- aov(Average.Amplitude ~ Response*DiffRandL, data=tactile)

AIC(m_1) 
AIC(m_2) #3rd
AIC(m_3) 
AIC(m_4)
AIC(m_5)#lowest
AIC(m_6)#2nd

summary(m_1) 
summary(m_2) #significant 2.86e-08
summary(m_3) #significant 0.0075
summary(m_4) #significant
summary(m_5) #significant

anova(m_2) #p=2.948e-08 ***
anova(m_4)
anova(m_5) #p=1.585e-10 ***
anova(m_6) #3.398e-09 ***


# hatched bs. nonhatched, corrected for time hatched 
ggplot(tactile, aes(x=Hours.Since.First.Hatch, y=Average.Amplitude, color=Response)) + 
  geom_point(shape=16, size=3) +
  scale_shape_discrete(solid=T) +
  ylab("VOR amplitude (°)\n")+  
  xlab("\nHours since first sibling hatched") + 
  theme(axis.title.x = element_text(vjust=-0.8)) +
  theme(axis.title.y = element_text(vjust=0.8)) +
  theme_bw(30) +
  theme(panel.grid.major = element_line(colour = "#808080"))+
  scale_colour_manual(values= c("forestgreen", "goldenrod"))

ggplot(tactile, aes(x=Developmental.Stage, y=Average.Amplitude, colour=Response)) + 
  geom_point(size=4, position = "jitter") +
  theme_bw(20)+
  #theme (panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())+
  #geom_errorbar(data = SibshipsSonia.df, aes(ymin=AverageAmp-SE, ymax=AverageAmp+SE), width = 0.05, colour = Clutch)+  
  ylab("VOR amplitude (°)\n")+  
  xlab("\nDevelopmental stage")

ggplot(tactile, aes(x=Developmental.Stage, y=Average.Amplitude, colour=Response)) + 
  geom_point(size=3) +
  theme_bw(20)+
  #theme (panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())+
  #geom_errorbar(data = SibshipsSonia.df, aes(ymin=AverageAmp-SE, ymax=AverageAmp+SE), width = 0.05, colour = Clutch)+  
  ylab("VOR amplitude (°)\n")+  
  xlab("\nDevelopmental stage")

# omit-rows-containing-specific-column-of-na
tactile<-tactile[!is.na(tactile$Developmental.Stage),]

ggplot(tactile, aes(x=as.factor(Developmental.Stage), y=Average.Amplitude, colour = Response)) + 
  geom_boxplot(data=tactile, size=1, na.rm=T) +
  geom_point(data=tactile, size=2, position="jitter")+
  ylab("VOR amplitude (°)\n")+
  theme_bw(20) +
  xlab("\n Developmental stage")

tactile$Response<-as.numeric(tactile$Response)
cor.test(tactile$Average.Amplitude, tactile$Developmental.Stage, method="pearson")
cor.test(tactile$Average.Amplitude, tactile$Response, method="pearson")
cor.test(tactile$Average.Amplitude, tactile$Hours.Since.First.Hatch, method="pearson")



# ############ ANOVA Stats ###############
# 
# m_1 <- aov(Average.Amp ~ 1, data=tactile)
# m_2 <-aov(Average.Amp ~ Response, data = tactile)
# m_3 <- aov(Average.Amp ~ Clutch, data = tactile)
# m_4 <- aov(Average.Amp ~ Clutch*Response, data = tactile)
# m_5 <- aov(Average.Amp ~ Response*HsinceH, data=tactile)
# m_6 <- aov(Average.Amp ~ Response*DiffRandL, data=tactile)
# 
# AIC(m_1) 
# AIC(m_2) #3rd
# AIC(m_3) 
# AIC(m_4)
# AIC(m_5)#
# AIC(m_6)#lowest
# 
# summary(m_1) 
# summary(m_2) #significant
# summary(m_3) 
# summary(m_4)
# summary(m_5)
# 
# anova(m_2) #***
# anova(m_4)
# anova(m_5) # ***
# anova(m_6) #3.393e-07 ***


# # #DF[!is.na(DF$y),]
# # #tactile$Average.Amp[!is.na(tactile$AverageR2),]
# # #tactileomitR2<- subset(tactile, !is.na(AverageR2))
# # # hatched vs. nonhatched, corrected for time hatched
# # ggplot(tactile, aes(x=HsinceH, y=Average.Amp, color=Hatched)) + geom_point(shape=1)
# # model_3 <- aov(Average.Amp ~ Hatched*HsinceH, data=tactile)
# # anova(model_3)
# # 
# # ggplot(tactile, aes(x=Average.Amp, y=AverageR2, color=Hatched)) + geom_point(shape=1)
# # 
# # 
# # #subset out low R2 values to look at them individually (portraits)
# # lowR2 <- subset(tactile, Average.Amp == 0)
# 
# # grid.locator()
# 
# # identifyPch <- function(x, y = NULL, n = length(x), pch = 19)
# # {
# #   xy <- xy.coords(x, y); x <- xy$x; y <- xy$y
# #   sel <- rep(FALSE, length(x)); res <- integer(0)
# #   while(sum(sel) < n) {
# #     ans <- identify(x[!sel], y[!sel], n = 1, plot = FALSE, ...)
# #     if(!length(ans)) break
# #     ans <- which(!sel)[ans]
# #     points(x[ans], y[ans], pch = pch)
# #     sel[ans] <- TRUE
# #     res <- c(res, ans)
# #   }
# #   res
# # }
# # identifyPch(tactile$Hatched=="H" && x<-5)
# 
# # by clutch, corrected for time hatched
# #ggplot(tactile, aes(x=HsinceH, y=Average.Amp, color=Hatched)) + geom_point(shape=1) + geom_smooth() 
# 
# # subset by clutch (11 total clutches)
# c254 <- subset(tactile, Clutch == "254")
# c255 <- subset(tactile, Clutch == "255")
# c256 <- subset(tactile, Clutch == "256")
# c257 <- subset(tactile, Clutch == "257")
# c258 <- subset(tactile, Clutch == "258")
# c259 <- subset(tactile, Clutch == "259")
# c262 <- subset(tactile, Clutch == "262")
# c263 <- subset(tactile, Clutch == "263")
# c264 <- subset(tactile, Clutch == "264")
# c265 <- subset(tactile, Clutch == "265")
# c266 <- subset(tactile, Clutch == "266")
# 
# # individual plots by clutch
# ggplot(c254, aes(x=HsinceH, y=Average.Amp)) + geom_point(shape=1) + scale_colour_hue(l=50) + geom_smooth(method=lm) 
# ggplot(c255, aes(x=HsinceH, y=Average.Amp)) + geom_point(shape=1) + scale_colour_hue(l=50) + geom_smooth() 
# ggplot(c256, aes(x=HsinceH, y=Average.Amp)) + geom_point(shape=1) + scale_colour_hue(l=50) + geom_smooth(method=lm) 
# ggplot(c257, aes(x=HsinceH, y=Average.Amp)) + geom_point(shape=1) + scale_colour_hue(l=50) + geom_smooth() 
# ggplot(c258, aes(x=HsinceH, y=Average.Amp)) + geom_point(shape=1) + scale_colour_hue(l=50) + geom_smooth() 
# ggplot(c259, aes(x=HsinceH, y=Average.Amp)) + geom_point(shape=1) + scale_colour_hue(l=50) + geom_smooth() 
# ggplot(c262, aes(x=HsinceH, y=Average.Amp)) + geom_point(shape=1) + scale_colour_hue(l=50) + geom_smooth() 
# ggplot(c263, aes(x=HsinceH, y=Average.Amp)) + geom_point(shape=1) + scale_colour_hue(l=50) + geom_smooth() 
# ggplot(c264, aes(x=HsinceH, y=Average.Amp)) + geom_point(shape=1) + scale_colour_hue(l=50) + geom_smooth() 
# ggplot(c265, aes(x=HsinceH, y=Average.Amp)) + geom_point(shape=1) + scale_colour_hue(l=50) + geom_smooth(method=lm) 
# ggplot(c266, aes(x=HsinceH, y=Average.Amp)) + geom_point(shape=1) + scale_colour_hue(l=50) + geom_smooth() 
# 
# # all together now
# dat <- rbind(c254,c255, c256,c257, c258, c259, c262, c263, c264, c265, c266)
# #dat <- rbind(df1,df2,df3)
# dat$grp <- rep(factor(1:11),times=c(nrow(c254),nrow(c255),nrow(c256),nrow(c257),nrow(c258),nrow(c259),nrow(c262),nrow(c263),nrow(c264),nrow(c265),nrow(c266)))
# cols <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "red", "green", "gray")
# ggplot(data = dat, aes(x=HsinceH, y=Average.Amp, colour = grp)) + 
#   geom_point(shape=1) + 
#   scale_colour_hue(l=50) + 
#   geom_smooth(method=lm,   # Add linear regression lines
#               se=FALSE,    # Don't add shaded confidence region
#               fullrange=TRUE) +# Extend regression lines
#   ylab("Average Amplitude") +
#   xlab("Hours Since First Hatch")

############ Standard Error Stats ###############

# errorstats <- summarySE(tactile, measurevar="Average.Amp", groupvars=c("AgeBlock"))
# errorstats
# 
# # Standard error of the mean
# ggplot(errorstats, aes(x=AgeBlock, y=Average.Amp)) + 
#   geom_errorbar(aes(ymin=Average.Amp-se, ymax=Average.Amp+se), width=.1) +
#   geom_line() +
#   geom_point()
# 
# # Use 95% confidence interval instead of SEM
# ggplot(errorstats, aes(x=AgeBlock, y=Average.Amp)) + 
#   geom_errorbar(aes(ymin=Average.Amp-ci, ymax=Average.Amp+ci), width=.1) +
#   geom_line() +
#   geom_point()
#

ShakeRoll.df$proportionNontest<-ShakeRoll.df$Nontest/(ShakeRoll.df$Test+ShakeRoll.df$Nontest)
propNontest<-na.omit(ShakeRoll.df$proportionNontest)
se <- function(x) sqrt(var(x)/length(x))
mean(propNontest)
se(propNontest)

# 
# ####### DiffRandL instead of HsinceH ########
# 
# coef(lm(DiffRandL~HsinceH, data=tactile))
# ggplot(tactile, aes(x=HsinceH, y=DiffRandL, color=Hatched)) + geom_point(shape=1) + geom_abline(intercept=-10.1924927)
# 
# model_1 <- aov(DiffRandL ~ HsinceH, data=tactile)
# anova(model_1)
# 
# # by clutch, corrected for pigeon-eyed developmental marker
# ggplot(tactile, aes(x=DiffRandL, y=Average.Amp, color=Clutch)) + geom_point(shape=1) + geom_smooth() 
# 
# # hatched bs. nonhatched, corrected for age/pigeon-eyed-ness
# ggplot(tactile, aes(x=DiffRandL, y=Average.Amp, color=Hatched)) + geom_point(shape=1)
# 
# model_2 <- aov(DiffRandL ~ Average.Amp, data=tactile)
# anova(model)
# 
# #########################################################
# ################### ONTOGENY FIGURE #####################
# #########################################################
# 
# hyp <-subset(VOR.df, Stimulus == "H", na.rm=T, select=c(Stimulus, Clutch, AgeBlock, NumH, Individual, Hatched, HatchTime, HatchAge, HsinceH, corrHsinceH, AverageR2, Average.Amp, DiffRandL))
# mech <- subset(VOR.df, Stimulus == "T", na.rm=T, select=c(Stimulus, Clutch, AgeBlock, NumH, Individual, Hatched, HatchTime, HatchAge, HsinceH, corrHsinceH, AverageR2, Average.Amp, DiffRandL))
# 
# hyp$corrHsinceH[is.na(hyp$NumH)] <- NA
# mech$corrHsinceH[is.na(mech$NumH)] <- NA
# 
# mech$NumH[mech$NumH == 0] <- 3
# mech$NumH[mech$NumH == 1] <- 4
# mech$NumH[mech$NumH == 2] <- 5
# 
# c254h <- subset(hyp, Clutch == "254")
# c255h <- subset(hyp, Clutch == "255")
# c256h <- subset(hyp, Clutch == "256")
# c257h <- subset(hyp, Clutch == "257")
# c258h <- subset(hyp, Clutch == "258")
# c259h <- subset(hyp, Clutch == "259")
# c262h <- subset(hyp, Clutch == "262")
# c263h <- subset(hyp, Clutch == "263")
# c264h <- subset(hyp, Clutch == "264")
# c265h <- subset(hyp, Clutch == "265")
# c266h <- subset(hyp, Clutch == "266")
# 
# c254m <- subset(mech, Clutch == "254")
# c255m <- subset(mech, Clutch == "255")
# c256m <- subset(mech, Clutch == "256")
# c257m <- subset(mech, Clutch == "257")
# c258m <- subset(mech, Clutch == "258")
# c259m <- subset(mech, Clutch == "259")
# c262m <- subset(mech, Clutch == "262")
# c263m <- subset(mech, Clutch == "263")
# c264m <- subset(mech, Clutch == "264")
# c265m <- subset(mech, Clutch == "265")
# c266m <- subset(mech, Clutch == "266")
# 
# clutch_colors <- rainbow(11)
# 
# plot(c254h$corrHsinceH, c254h$NumH, pch=16, ylim=c(0,5), xlim=c(0,25), axes=FALSE, xlab="", ylab="", col=clutch_colors[1])
# #type="b" WHY WON"T THIS WORK???
# axis(2, ylim=c(0,5), lab=F, col="black", las=1)
# text(labels=c(0, 1, 2, 0, 1, 2))
# mtext("Hypoxia Hatch               Mechanical Hatch",side=2,line=2.5)
# box()
# 
# par(new=TRUE)## Allow a second plot on the same graph
# 
# ## Plot the second plot and put axis scale on right
# lines(c254m$corrHsinceH, c254m$NumH, pch=15,  xlab="", ylab="", ylim=c(0,5), xlim=c(0,25), type="b", ann=T, col=clutch_colors[1])
# ## a little farther out (line=4) to make room for labels
# #mtext("Mechanical Hatch", side=4, line=1.5) 
# #axis(4, ylim=c(0,5), las=1)
# 
# ### Draw the time axis
# #axis(1,pretty(range(corrHsinceH),1))
# mtext("Time (Hours since first hatch)",side=1,col="black",line=2.5)  
# axis(1, xlim=c(0,25), las=1)
# ## Add Legend
# legend("topleft",legend=c("Hypoxia Hatch","Mechanical Hatch"),pch=c(16,15))
# 
# par(new=TRUE)
# plot(c255h$corrHsinceH, c255h$NumH, type="l", pch=16, ylim=c(0,5), xlim=c(0,25), axes=FALSE, xlab="", ylab="", type="b", col=clutch_colors[2])
# lines(c255m$corrHsinceH, c255m$NumH, pch=15,  xlab="", ylab="", ylim=c(0,5), xlim=c(0,25), type="b", col=clutch_colors[2])
# 
# par(new=TRUE)
# plot(c256h$corrHsinceH, c256h$NumH, pch=16, ylim=c(0,5), xlim=c(0,25), axes=FALSE, xlab="", ylab="", type="b", col=clutch_colors[3])
# lines(c256m$corrHsinceH, c256m$NumH, pch=15,  xlab="", ylab="", ylim=c(0,5), xlim=c(0,25), type="b", col=clutch_colors[3])
# 
# par(new=TRUE)
# plot(c257h$corrHsinceH, c257h$NumH, pch=16, ylim=c(0,5), xlim=c(0,25), axes=FALSE, xlab="", ylab="", type="b", col=clutch_colors[4])
# lines(c257m$corrHsinceH, c257m$NumH, pch=15,  xlab="", ylab="", ylim=c(0,5), xlim=c(0,25), type="b", col=clutch_colors[4])
# 
# par(new=TRUE)
# plot(c258h$corrHsinceH, c258h$NumH, pch=16, ylim=c(0,5), xlim=c(0,25), axes=FALSE, xlab="", ylab="", type="b", col=clutch_colors[5])
# lines(c258m$corrHsinceH, c258m$NumH, pch=15,  xlab="", ylab="", ylim=c(0,5), xlim=c(0,25), type="b", col=clutch_colors[5])
# 
# par(new=TRUE)
# plot(c259h$corrHsinceH, c259h$NumH, pch=16, ylim=c(0,5), xlim=c(0,25), axes=FALSE, xlab="", ylab="", type="b", col=clutch_colors[6])
# lines(c259m$corrHsinceH, c259m$NumH, pch=15,  xlab="", ylab="", ylim=c(0,5), xlim=c(0,25), type="b", col=clutch_colors[6])
# 
# par(new=TRUE)
# plot(c262h$corrHsinceH, c262h$NumH, pch=16, ylim=c(0,5), xlim=c(0,25), axes=FALSE, xlab="", ylab="", type="b", col=clutch_colors[7])
# lines(c262m$corrHsinceH, c262m$NumH, pch=15,  xlab="", ylab="", ylim=c(0,5), xlim=c(0,25), type="b", col=clutch_colors[7])
# 
# par(new=TRUE)
# plot(c263h$corrHsinceH, c263h$NumH, pch=16, ylim=c(0,5), xlim=c(0,25), axes=FALSE, xlab="", ylab="", type="b", col=clutch_colors[8])
# lines(c263m$corrHsinceH, c263m$NumH, pch=15,  xlab="", ylab="", ylim=c(0,5), xlim=c(0,25), type="b", col=clutch_colors[8])
# 
# par(new=TRUE)
# plot(c264h$corrHsinceH, c264h$NumH, pch=16, ylim=c(0,5), xlim=c(0,25), axes=FALSE, xlab="", ylab="", type="b", col=clutch_colors[9])
# lines(c264m$corrHsinceH, c264m$NumH, pch=15,  xlab="", ylab="", ylim=c(0,5), xlim=c(0,25), type="b", col=clutch_colors[9])
# 
# par(new=TRUE)
# plot(c265h$corrHsinceH, c265h$NumH, pch=16, ylim=c(0,5), xlim=c(0,25), axes=FALSE, xlab="", ylab="", type="b", col=clutch_colors[10])
# lines(c265m$corrHsinceH, c265m$NumH, pch=15,  xlab="", ylab="", ylim=c(0,5), xlim=c(0,25), type="b", col=clutch_colors[10])
# 
# par(new=TRUE)
# plot(c266h$corrHsinceH, c266h$NumH, pch=16, ylim=c(0,5), xlim=c(0,25), axes=FALSE, xlab="", ylab="", type="b", col=clutch_colors[11])
# lines(c266m$corrHsinceH, c266m$NumH, pch=15,  xlab="", ylab="", ylim=c(0,5), xlim=c(0,25), type="b", col=clutch_colors[11])
# 
# # # add lines 
# # for (i in 1:11) { 
# #   Clutch <- subset(Orange, Tree==i) 
# #   lines(tree$age, tree$circumference, type="b", lwd=1.5,
# #         lty=linetype[i], col=colors[i], pch=plotchar[i]) 
# # } 
# 
# c254 <- rbind(c254h, c254m)
# c255 <- rbind(c255h, c255m)
# c256 <- rbind(c256h, c256m)
# c257 <- rbind(c257h, c257m)
# c258 <- rbind(c258h, c258m)
# c259 <- rbind(c259h, c259m)
# c262 <- rbind(c262h, c262m)
# c263 <- rbind(c263h, c263m)
# c264 <- rbind(c264h, c264m)
# c265 <- rbind(c265h, c265m)
# c266 <- rbind(c266h, c266m)
# 
# # ggplot(c255, aes(x=corrHsinceH, y=NumH, shape=Stimulus))+ geom_point(shape=1) + geom_line() +
# #   ylab("Hypoxia Hatch                                                 Mechanical Hatch") +
# #   xlab("Hours Since First Hatch")
# # par(new=TRUE)
# # ggplot(c256, aes(x=corrHsinceH, y=NumH, shape=Stimulus))+ geom_point(shape=2) + geom_line() +
# #   ylab("Hypoxia Hatch                                                 Mechanical Hatch") +
# #   xlab("Hours Since First Hatch") 






# ggplot(ShakeRoll.df, aes(x=AverageAmp, y=PropH, color=HorNH)) + 
#   geom_point(shape=16, size=4, position="jitter") +
#   ylab("Proportion of clutch hatched\n") +  
#   theme_bw(20) +
#   xlab("VOR amplitude (°)\n")
# 
# # ShakeRoll.df$HorNH<-as.factor(ShakeRoll.df$HorNH)
# # ShakeRoll.df$AverageR2<-as.numeric(ShakeRoll.df$AverageR2)
# # ShakeRoll.df$AverageAmp<-as.numeric(ShakeRoll.df$AverageAmp)
# # ShakeRoll.df$SUM.of.trait.values<-as.numeric(ShakeRoll.df$SUM.of.trait.values)
# 
# 
# ggplot(ShakeRoll.df, aes(x=SUM.of.trait.values, y=PropH, color = HorNH)) + 
#   geom_point(shape=16, size=3, position = "jitter") +
#   ylab("Proportion of clutch hatched\n") +  
#   theme_bw(20) +
#   xlab("\n Developmental stage")
# 
# 
# col = rainbow(14)
# ggplot(ShakeRoll.df, aes(x=AverageAmp, y=PropH, color=SUM.of.trait.values)) + 
#   geom_point(shape=16, size=4, position="jitter") +
#   ylab("Proportion of clutch hatched\n") +  
#   theme_bw(20) +
#   xlab("VOR amplitude (°)\n")
# 
# yesVOR <- subset(ShakeRoll.df, AverageAmp > 5)
# yesVORhat <- subset(yesVOR, HorNH=="Hatched")
# yesVORnohat <- subset(yesVOR, HorNH=="Not Hatched")
# noVOR <- subset(ShakeRoll.df, AverageAmp < 5)
# noVORhat <- subset(noVOR, HorNH=="Hatched")
# noVORnohat <- subset(noVOR, HorNH=="Not Hatched")









################################################################################
################################################################################
################################################################################
################################################################################
################################################################################


#geom_errorbar(data=error, aes(ymin=Average-se, ymax=Average+se), width=.1)
#scale_fill_continuous(guide=FALSE)
# ggplot(ShakeRoll.df, aes(x=AverageAmp, y=PropH, color = Hcontext)) + 
#   geom_point(shape=16, size=4) +
#   ylab("Proportion Hatched\n") +  
#   theme_bw(20) +
#   xlab("\nIndividual VOR Amplitude (degrees)")
# 
HatchingContext.df<-read.csv(file="HatchingContext.csv")
# ggplot(HatchingContext.df, aes(x=AverageAmp, y=proportionH, color = Hcontext)) + 
#   geom_point(shape=16, size=4) +
#   ylab("Proportion Hatched\n") +  
#   theme_bw(20) +
#   xlab("\nIndividual VOR Amplitude (degrees)")
# 
# ggplot(HatchingContext.df, aes(x=Time, y=AverageAmp, color = SUM.of.trait.values)) + 
#   geom_point(shape=16, size=4) +
#   ylab("Proportion Hatched\n") +  
#   theme_bw(20) +
#   xlab("\nIndividual VOR Amplitude (degrees)")
# 
# m_1 <- aov(AverageAmp ~ 1, data=HatchingContext.df)
# m_2 <-aov(AverageAmp ~ Hcontext, data = HatchingContext.df)
# AIC(m_1) 
# AIC(m_2) 
# summary(m_1) 
# summary(m_2)
# anova(m_1) 
# anova(m_2)

######START #In clutches with partial vibration cued hatching, hatched embryos did have greater VOR than unhatched (hypoxia hatched) siblings. 

# # partial = clutches where not all H or not all NH.
# # take out clutches where sumH (hatched in response to playback or hatched inresponse to tactile) is 0
# HatchingContext.df$sumH[HatchingContext.df$sumH==0] <- NA
# # same take out hypoxiaH == 0
# HatchingContext.df$hypoxiaH[HatchingContext.df$hypoxiaH==0] <- NA
# #left with only clutches with partial vibration cued hatching

HatchingContext.df$propH[HatchingContext.df$propH==0] <- NA
HatchingContext.df$propH[HatchingContext.df$propH==1] <- NA
HatchingContext.df$AverageAmp[HatchingContext.df$propH==0] <- NA
HatchingContext.df$AverageAmp[HatchingContext.df$propH==1] <- NA
HatchingContext.df$contxt[HatchingContext.df$propH==0] <- NA
HatchingContext.df$contxt[HatchingContext.df$propH==1] <- NA

hist(HatchingContext.df$AverageAmp) #normal--> parametric test

#results1 = lm(AverageAmp ~ contxt, data=HatchingContext.df)
#anova(results1)

#incorporate clutch as a random effect
library(lme4)
results0 = lmer(AverageAmp ~ 1 + (1|Clutch), data=HatchingContext.df)
#summary(results0)
#qqnorm(resid(results0)); qqline (resid(results0))
#plot(results0$fitted, results0$res, xlab="Fitted", ylab="Residuals")
results1 = lmer(AverageAmp ~ contxt + (1|Clutch), data=HatchingContext.df)
#results2 = lmer(AverageAmp ~ propH + contxt + (1|Clutch), data=HatchingContext.df)
#results3 = lmer(AverageAmp ~ propH * contxt + (1|Clutch), data=HatchingContext.df)

anova(results0,results1)

#anova(results1,results2) #propH has a significant effect on Average Amp
#anova(results0,results2) #hatching context has a nearly significant effect on Average Amp
#anova(results2,results3) #significant interaction effect between propH and hatching context


#HatchingContext.df$AverageAmp[HatchingContext.df$AverageAmp<5] <- NA
#HatchingContext.df$AverageAmp <- replace(HatchingContext.df$AverageAmp, is.na(HatchingContext.df$AverageAmp), 0)
# ggplot(HatchingContext.df, aes(x=AverageAmp, y=proportionH, color = Hcontext), na.rm=T) + 
#   geom_point(shape=16, size=4) +
#   ylab("Proportion Hatched\n") +  
#   theme_bw(20) +
#   xlab("\nIndividual VOR Amplitude (degrees)")



#subset H vs. NH individuals
Hatched <- subset(HatchingContext.df, Hcontext == 1 , na.rm=T)
NotHatched <- subset(HatchingContext.df, Hcontext ==2 | Hcontext==3, na.rm=T)

######################################################
boxplot(Hatched$AverageAmp, NotHatched$AverageAmp, xlab="Hatching Context", ylab="VOR Amplitude")
axis(1, at=1:2, labels=c("Vibration-cued", "Tactile/Hyp"))
#####################################################

t.test(Hatched$AverageAmp, NotHatched$AverageAmp)


#Maximum likelihood analysis from Moczek & Nijhout 2003
#use new package: 
#"grofit: fitting biological growth curves with R"



# # ###########take out low R2 values ########
ggplot(tactile, aes(x=Average.Amp, y=AverageR2, color=Response)) +
  geom_point(shape=16, size=3) +
  scale_shape_discrete(solid=T) +
  ylab("Average R2") +  
  xlab("VOR amplitude (°)\n")+
  theme(axis.title.x = element_text(vjust=-0.8)) +
  theme(axis.title.y = element_text(vjust=0.8)) +
  theme_bw(30) +
  theme(panel.grid.major = element_line(colour = "#808080"))+
  scale_colour_manual(values= c("forestgreen", "goldenrod"))
################according to this figure - should make all VOR =0 all individuals with average R2<0.6
# ## put these as VOR = 0 instead. 
## not a relationship --> no response. 
#

subsettactile <- tactile
subsettactile$AverageR2[subsettactile$AverageR2<0.6] <- NA
subsettactile$Average.Amp <- replace(subsettactile$Average.Amp, is.na(subsettactile$AverageR2), 0)

ggplot(subsettactile, aes(x=HsinceH, y=Average.Amp, color=Response)) + 
  geom_point(shape=16, size=3) +
  geom_jitter(position = "jitter")+        
  scale_shape_discrete(solid=T) +
  ylab("VOR amplitude (°)\n")+ 
  xlab("\nHours since first sibling hatched") + 
  theme(axis.title.x = element_text(vjust=-0.8)) +
  theme(axis.title.y = element_text(vjust=0.8)) +
  theme_bw(30) +
  theme(panel.grid.major = element_line(colour = "#808080"))+
  scale_colour_manual(values= c("forestgreen", "goldenrod"))


str(ShakeRoll.df)
ShakeRoll.df$HorNH<-as.factor(ShakeRoll.df$HorNH)
ggplot(ShakeRoll.df, aes(x=AverageAmp, y=AverageR2, color=HorNH)) +
  geom_point(shape=16, size=3) +
  scale_shape_discrete(solid=T) +
  ylab("Average R2") +  
  xlab("VOR amplitude (°)\n")+ 
  theme(axis.title.x = element_text(vjust=-0.8)) +
  theme(axis.title.y = element_text(vjust=0.8)) +
  theme_bw(30) +
  theme(panel.grid.major = element_line(colour = "#808080"))+
  scale_colour_manual(values= c("forestgreen", "goldenrod"))

subsetshakeroll <- ShakeRoll.df
subsetshakeroll$AverageR2[subsetshakeroll$AverageR2<0.7] <- NA
subsetshakeroll$Average.Amp <- replace(subsetshakeroll$Average.Amp, is.na(subsetshakeroll$AverageR2), 0)

ggplot(subsetshakeroll, aes(x=AverageAmp, y=AverageR2, color=HorNH)) +
  geom_point(shape=16, size=3) +
  scale_shape_discrete(solid=T) +
  ylab("Average R2") +  
  xlab("VOR amplitude (°)\n")+
  theme(axis.title.x = element_text(vjust=-0.8)) +
  theme(axis.title.y = element_text(vjust=0.8)) +
  theme_bw(30) +
  theme(panel.grid.major = element_line(colour = "#808080"))+
  scale_colour_manual(values= c("forestgreen", "goldenrod"))

ggplot(subsetshakeroll, aes(x=AverageAmp, y=PropH, color=HorNH)) +
  geom_point(shape=16, size=3) +
  scale_shape_discrete(solid=T) +
  ylab("Proportion of clutch hatched") +  
  xlab("VOR amplitude (°)\n")+
  theme(axis.title.x = element_text(vjust=-0.8)) +
  theme(axis.title.y = element_text(vjust=0.8)) +
  theme_bw(30) +
  theme(panel.grid.major = element_line(colour = "#808080"))+
  scale_colour_manual(values= c("forestgreen", "goldenrod"))



