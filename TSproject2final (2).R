##################################################################################
#                        Time Series Econometrics: Lab #3                        #
#                                 October 2, 2019                                #
#                                 Julieta Yung                                   #
##################################################################################

# preliminaries
install.packages("xts")
install.packages("pdfetch")
install.packages("stargazer")
install.packages("forecast")
install.packages("vars")
install.packages("corrplot")
install.packages("tseries")
install.packages("zoo")
install.packages("fUnitRoots")
install.packages("hydroTSM")
install.packages("tsbox")

library(xts)       # time series
library(pdfetch)   # loads data from the internet
library(stargazer) # nicely formatted tables
library(forecast)  # has Acf and Pacf functions
library(vars)      # run VARs and diagnostics
library(corrplot)  # correlation matrix
library(tseries)   # for adf.test()
library(zoo)       # idk
library(fUnitRoots)#for unit roots
library(hydroTSM)  #handles converting dates
library(tsbox)     #converts zoo into ts
##################################################################################
# STEP 1: Data Formatting                                                        #
##################################################################################

#note that dates are formatted in YYYY-MM-DD
#designated frequency is monthly

# load data from FRED (xts)
Swiss.XPV.FRED <- pdfetch_FRED("XTEXVA01CHM664S")  # Value of Exports, start=1960-01-01, end=2019-06-01
Swiss.CPI.FRED <- pdfetch_FRED("CHECPIALLMINMEI")  # CPI for Switzerland, Seasonally Adjusted, start=1960-01-01, end=2019-08-01
Swiss.1M.UNCLEAN.Daily.FRED <- pdfetch_FRED("CHF1MTD156N")       # Swiss 1 month start=1989-01-03, end=2019-10-01
Swiss.US.ER.FRED <- pdfetch_FRED("EXSZUS")        # US/Franc Exchange rate, start=1971-01-01, end=2019-09-01

Swiss.1M.FRED = daily2monthly(Swiss.1M.UNCLEAN.Daily.FRED, mean ,na.rm=TRUE)

Swiss.1M.FRED = xts(Swiss.1M.FRED)
####################

#Swiss.XPV <- ts(Swiss.XPV.FRED, frequency=12, class = "ts")
#Swiss.CPI <- ts(Swiss.CPI.FRED, frequency=12, class = "ts")
#Swiss.1M <- ts(Swiss.1M.FRED, frequency=12, class = "ts")
#Swiss.US.ERb <- ts(Swiss.US.ER.FRED, frequency=12, class = "ts")
#Swiss.US.ER <- lag(Swiss.US.ERb, k=215)

Swiss.XPV <- ts(Swiss.XPV.FRED, start=c(1960,1), end = c(2019,8), frequency=12, class = "ts")
Swiss.CPI <- ts(Swiss.CPI.FRED, start=c(1960,1), end = c(2019,9), frequency=12, class = "ts")
Swiss.1M <- ts(Swiss.1M.FRED, start=c(1989,1), end = c(2019,10), frequency=12, class = "ts")
Swiss.US.ERb <- ts(Swiss.US.ER.FRED, start=c(1971,1), end = c(2019,9), frequency=12, class = "ts")

plot(Swiss.1M)
plot(Swiss.XPV)
plot(Swiss.CPI)
plot(Swiss.US.ERb)
#################### transform variables

Swiss.Inflation <- diff(log(Swiss.CPI),lag=12)*100 # annualized inflation 
#Swiss.Inflation <- lag(Swiss.Inflation, k=347)
adfTest(Swiss.Inflation, lags = 3, type=c("c","ct","nc"))

#adfTest(Swiss.PercentD.XPV, lags = 3, type=c("c","ct","nc"))
Swiss.PercentD.XPV <- diff(log(Swiss.XPV), lag=12)*100 #annualized %change in export value
#Swiss.PercentD.XPV <- lag(Swiss.PercentD.XPV, k=346)
adfTest(Swiss.PercentD.XPV, lags = 3, type=c("c","ct","nc"))

#adfTest(Swiss.US.ER, lags = 3, type=c("c","ct","nc"))
Swiss.PercentD.ER <- diff(log(Swiss.US.ERb), lag=12)*100 #annualized %change in export value
#Swiss.PercentD.ER <- lag(Swiss.PercentD.ER, k=215)
adfTest(Swiss.PercentD.ER, lags = 3, type=c("c","ct","nc"))

#adfTest(Swiss.1M, lags = 3, type=c("c","ct","nc"))
Swiss.1M.Diff <- diff(Swiss.1M, lag=12)
adfTest(Swiss.1M.Diff, lags = 3, type=c("c","ct","nc"))

###################### bind all data together and give variables a name


#VAR OPTION 1: normal exchange rate
#VAR.Swiss.data.all <- cbind(XPV.Delta=Swiss.PercentD.XPV, Inflation=Swiss.Inflation, One.M.Diff=Swiss.1M.Diff, Swiss.Us.XRate=Swiss.US.ER)

#VAR OPTION 2: %change ER
#VAR.Swiss.data.all <- cbind(export=Swiss.PercentD.XPV, inflation=Swiss.Inflation, policyrate=Swiss.1M.Diff, exchangerate=Swiss.PercentD.ER)

#VAR OPTION 3: Inflation as first variable, 1M diff as second, ER as third, XPV Percent change as last
VAR.Swiss.data.all <- cbind(inflation=Swiss.Inflation, export=Swiss.PercentD.XPV, policyrate=Swiss.1M.Diff, exchangerate=Swiss.PercentD.ER)

plot(VAR.Swiss.data.all)

VAR.Swiss.Window.data <- na.omit(VAR.Swiss.data.all)

plot(VAR.Swiss.Window.data)

##################################################################################
# Preliminary data analysis                                                       #
##################################################################################

# Plot the auto and cross-correlation
# first column: GDP and the lags of GDP; lags of inf, lags of pol
Acf( VAR.Swiss.Window.data, lag=24, las=1, cex=1.5)
Pacf(VAR.Swiss.Window.data, lag=24, las=1, cex=1.5)

# Correlation matrix (not very pretty)
heatmap(x = cor(VAR.Swiss.Window.data), symm = TRUE)

# but we can format it (thanks Bennett 2018!)
col        <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
cex.before <- par("cex")
par(cex = 1.5) #Set text size for matrix

corrplot(cor(VAR.Swiss.Window.data), 
         method="color", col=col(200),  
         type="upper", order="FPC", 
         addCoef.col = "gray15", # Add coefficient of correlation
         tl.cex = 1.5/par("cex"),
         cl.cex = 1.5/par("cex"),
         tl.col="black", tl.srt=25, #Text label color and rotation
         diag=TRUE # show correlation coefficient on the principal diagonal
)

##################################################################################
# STEP 2: Decide how many lags to include in the VAR                             #
##################################################################################

# select VAR lags
VARselect(VAR.Swiss.Window.data, type="both") #We seemingly receive a rec for very high number of lags

##################################################################################
# STEP 3: Estimate the reduced form VAR to get coefficients and residuals        #
##################################################################################

# VAR model
VAR.Swiss <- VAR(VAR.Swiss.Window.data, p=5)

##################################################################################
# STEP 4: Evaluate the model(s) and run diagnostics                              #
##################################################################################

summary(VAR.Swiss)
stargazer(VAR.Swiss$varresult, type="text")

max(roots(VAR.Swiss)) # must be less than 1, we get greater than 1 for current set of variables

# Test for serial correlation in the residuals using the Portmanteau test
# Null: No autocorrelation of residuals (want to fail to reject, p-value large)
serial.test(VAR.Swiss, lags.pt=(300))

# Check the serial correlation of residuals using Ljung-Box test 
# (Null: independence, want to fail to reject, large p-value)
Box.test(VAR.Swiss$varresult$export$residuals)
Box.test(VAR.Swiss$varresult$inflation$residuals)
Box.test(VAR.Swiss$varresult$policyrate$residuals)
Box.test(VAR.Swiss$varresult$exchangerate$residuals)

# Check root-stability and stability of the empirical fluctuation processes
plot(stability(VAR.Swiss))

# test whether residuals are normally distributed
# Null: the data are normally distributed (you want to fail to reject, p-value large)
# multivariate.only = FALSE: Does joint testing (with Cholesky), otherwise, individual
normality.test(VAR.Swiss) 

# Plot the histogram and a normal curve
x    <- VAR.Swiss$varresult$export$residuals
h    <- hist(x, breaks=30, col="firebrick3", xlab="%",
             main="Histogram with Normal Curve")
xfit <- seq(min(x),max(x),length=40)
yfit <- dnorm(xfit,mean=mean(x),sd=sd(x))
yfit <- yfit*diff(h$mids[1:2])*length(x)
lines(xfit, yfit, col="blue", lwd=2) 

x    <- VAR.Swiss$varresult$exchangerate$residuals
h    <- hist(x, breaks=30, col="firebrick3", xlab="%",
             main="Histogram with Normal Curve")
xfit <- seq(min(x),max(x),length=40)
yfit <- dnorm(xfit,mean=mean(x),sd=sd(x))
yfit <- yfit*diff(h$mids[1:2])*length(x)
lines(xfit, yfit, col="blue", lwd=2)

x    <- VAR.Swiss$varresult$policyrate$residuals
h    <- hist(x, breaks=30, col="firebrick3", xlab="%",
             main="Histogram with Normal Curve")
xfit <- seq(min(x),max(x),length=40)
yfit <- dnorm(xfit,mean=mean(x),sd=sd(x))
yfit <- yfit*diff(h$mids[1:2])*length(x)
lines(xfit, yfit, col="blue", lwd=2)

x    <- VAR.Swiss$varresult$inflation$residuals
h    <- hist(x, breaks=30, col="firebrick3", xlab="%",
             main="Histogram with Normal Curve")
xfit <- seq(min(x),max(x),length=40)
yfit <- dnorm(xfit,mean=mean(x),sd=sd(x))
yfit <- yfit*diff(h$mids[1:2])*length(x)
lines(xfit, yfit, col="blue", lwd=2)

# Test for heteroskedasticity in the residuals
# Null: squared residuals are a sequence of white noise (you want to fail to reject, p-value large)
arch.test(VAR.Swiss)


##################################################################################
# STEP 5: Estimate the IRF for each structural shock                             #
##################################################################################

# Check out all IRFs
VAR.Swiss.IRF <-irf(VAR.Swiss)
plot(VAR.Swiss.IRF) #Note for the response of xpv due to Swiss us exrate


# Pick only a specific impulse and a specific response variable
var.3.IRF.y.i <-irf(VAR.Swiss, impulse="policyrate", response="export")
plot(var.3.IRF.y.i)

# Adjust preferences (n.ahead= periodS after the shock, ci= confidence interval, runs=how many times to run to get C.I.
var.3.IRF.y.j <-irf(VAR.Swiss, impulse="exchangerate", response="export", n.ahead=48, ci= 0.90, runs=10000, cumulative=FALSE)
plot(var.3.IRF.y.j)

var.3.IRF.y.j <-irf(VAR.Swiss, impulse="policyrate", response="exchangerate", n.ahead=48, ci= 0.90, runs=1000, cumulative=FALSE)
plot(var.3.IRF.y.j)

var.3.IRF.y.j <-irf(VAR.Swiss, impulse="exchangerate", response="export", n.ahead=48, ci= 0.90, runs=1000, cumulative=FALSE)
plot(var.3.IRF.y.j)

var.3.IRF.y.k <-irf(VAR.Swiss, impulse="exchangerate", response="policyrate", n.ahead=48, ci= 0.90, runs=1000, cumulative=FALSE)
plot(var.3.IRF.y.k)

# more information for the command: https://www.rdocumentation.org/packages/vars/versions/1.5-3/topics/irf

##################################################################################
# beautiful IRF (my code to make the IRFs look better for my papers)             #
##################################################################################
n <- 20 # how many periods you want to plot?
z        <- var.3.IRF.y.j       # name of the one response you want to plot
xx       <- 0:n                 # horizon
Exp      <- xx*0                # zero line
Con.High <- z$Upper$exchangerate  # upper CB
Con.Low  <- z$Lower$exchangerate  # lower CB
One.Run  <- z$irf$exchangerate    # IRF

# create plot (fix the small plot box problem!)
par(mar=c(4,6,1,1))                 # set the size of the plot
lim.high=max(Con.High)  # set the highest limit for the plot
lim.low =min(Con.Low)   # set the lowest limit for the plot
plot(Exp~xx,type="l",
     xlab= "month",
     ylab= "",
     lty="dashed",
     lwd=3,
     ylim=range(One.Run,Exp,lim.high,lim.low),
     cex.lab=2,las=1,cex.axis=2)

# use polygon to shade the area between confidence intervals and plot the IRF
polygon(c(xx,rev(xx)),c(Con.Low,rev(Con.High)),col="slategray2",border=NA)
lines(xx,Exp,lty="dashed",lwd=3)
lines(xx,Exp,lty="solid",col="darkblue",lwd=1.5)
lines(xx,One.Run, lwd=3)

##################################################################################
# Exercise: Change the order of the variables? How do the IRFs change?           #
##################################################################################

# Example: Policy Rate first, Inflation, GDP growth
VAR.Swiss2 <- cbind(data.var[,4], data.var[,3], data.var[,1], data.var[,2])
plot(VAR.Swiss2)
var.3.irf.y <- VAR(data.var.i.pi.y, p=5)

var.3.i.pi.y.IRF <- irf(var.3.i.pi.y)
plot(var.3.i.pi.y.IRF)

# Include more lags in the VAR, how do the IRF change?

##################################################################################
# STEP 6: Estimate the FEVD for each structural shock                            #
##################################################################################
VAR.Swiss.FEVD <- fevd(VAR.Swiss, n.ahead=30)
plot(VAR.Swiss.FEVD)

# How much does each variable (in the legend) contribute to explaining the variance 
# in the forecast error of each variable (in the subplot)?
# Example: For GDP growth, GDP growth accounts for most of the variance... (>80%)

# more information for the command: https://www.rdocumentation.org/packages/vars/versions/1.5-3/topics/fevd

##################################################################################
# STEP 7: Estimate the HVD for each structural shock                             #
##################################################################################

# Run all of this first, these are some auxiliary functions that are not included 
# in the packages. Based on Cesa-Bianchi's Matlab toolbox
##################################################################################
VARhd <- function(Estimation){
  ## make X and Y
  nlag    <- Estimation$p   # number of lags
  DATA    <- Estimation$y   # data
  QQ      <- VARmakexy(DATA,nlag,1)
  ## Retrieve and initialize variables 
  invA    <- t(chol(as.matrix(summary(Estimation)$covres)))   # inverse of the A matrix
  Fcomp   <- companionmatrix(Estimation)                      # Companion matrix
  #det     <- c_case                                           # constant and/or trends
  F1      <- t(QQ$Ft)                                         # make comparable to notes
  eps     <- ginv(invA) %*% t(residuals(Estimation))          # structural errors 
  nvar    <- Estimation$K                                     # number of endogenous variables
  nvarXeq <- nvar * nlag                                      # number of lagged endogenous per equation
  nvar_ex <- 0                                                # number of exogenous (excluding constant and trend)
  Y       <- QQ$Y                                             # left-hand side
  #X       <- QQ$X[,(1+det):(nvarXeq+det)]                    # right-hand side (no exogenous)
  nobs    <- nrow(Y)                                          # number of observations
  
  ## Compute historical decompositions
  # Contribution of each shock
  invA_big <- matrix(0,nvarXeq,nvar)
  invA_big[1:nvar,] <- invA
  Icomp <- cbind(diag(nvar), matrix(0,nvar,(nlag-1)*nvar))
  HDshock_big <- array(0, dim=c(nlag*nvar,nobs+1,nvar))
  HDshock <- array(0, dim=c(nvar,(nobs+1),nvar))
  
  for (j in 1:nvar){  # for each variable
    eps_big <- matrix(0,nvar,(nobs+1)) # matrix of shocks conformable with companion
    eps_big[j,2:ncol(eps_big)] <- eps[j,]
    for (i in 2:(nobs+1)){
      HDshock_big[,i,j] <- invA_big %*% eps_big[,i] + Fcomp %*% HDshock_big[,(i-1),j]
      HDshock[,i,j] <-  Icomp %*% HDshock_big[,i,j]
    } 
  } 
  HD.shock <- array(0, dim=c((nobs+nlag),nvar,nvar))   # [nobs x shock x var]
  for (i in 1:nvar){
    for (j in 1:nvar){
      HD.shock[,j,i] <- c(rep(NA,nlag), HDshock[i,(2:dim(HDshock)[2]),j])
    }
  }
  return(HD.shock)
}
VARmakexy <- function(DATA,lags,c_case){
  nobs <- nrow(DATA)
  #Y matrix 
  Y <- DATA[(lags+1):nrow(DATA),]
  Y <- DATA[-c(1:lags),]
  #X-matrix 
  if (c_case==0){
    X <- NA
    for (jj in 0:(lags-1)){
      X <- rbind(DATA[(jj+1):(nobs-lags+jj),])
    } 
  } else if(c_case==1){ #constant
    X <- NA
    for (jj in 0:(lags-1)){
      X <- rbind(DATA[(jj+1):(nobs-lags+jj),])
    }
    X <- cbind(matrix(1,(nobs-lags),1), X) 
  } else if(c_case==2){ # time trend and constant
    X <- NA
    for (jj in 0:(lags-1)){
      X <- rbind(DATA[(jj+1):(nobs-lags+jj),])
    }
    trend <- c(1:nrow(X))
    X <-cbind(matrix(1,(nobs-lags),1), t(trend))
  }
  A <- (t(X) %*% as.matrix(X)) 
  B <- (as.matrix(t(X)) %*% as.matrix(Y))
  Ft <- ginv(A) %*% B
  retu <- list(X=X,Y=Y, Ft=Ft)
  return(retu)
}
companionmatrix <- function (x) 
{
  if (!(class(x) == "varest")) {
    stop("\nPlease provide an object of class 'varest', generated by 'VAR()'.\n")
  }
  K <- x$K
  p <- x$p
  A <- unlist(Acoef(x))
  companion <- matrix(0, nrow = K * p, ncol = K * p)
  companion[1:K, 1:(K * p)] <- A
  if (p > 1) {
    j <- 0
    for (i in (K + 1):(K * p)) {
      j <- j + 1
      companion[i, j] <- 1
    }
  }
  return(companion)
}
##################################################################################


# Now run the command we just created
var.3.HVD <- VARhd(Estimation=VAR.Swiss) # this gives you the HVD

# plot results (but R Studio does not plot results well with barplot)
barplot(t(var.3.HVD[,,4]),
        main="HVD of the Policy Rate", col=rainbow(4),
        legend=colnames(VAR.Swiss.Window.data)) 

# option 1: you can look at the HVD one shock at a time
barplot(t(var.3.HVD[,4,1]),
        main="HVD of the Export", 
        legend="Exchange rate Shock") 
barplot(t(var.3.HVD[,4,2]), 
        main="HVD of the Inflation", 
        legend="Exchange rate shock") 
barplot(t(var.3.HVD[,4,3]), 
        main="HVD of the Policy Rate", 
        legend="Exchange Rate shock") 


barplot(t(var.3.HVD[,1,3]),
        main="HVD of the policy Rate", 
        legend="export") 
barplot(t(var.3.HVD[,2,3]), 
        main="HVD of the policy Rate", 
        legend="inflation shock") 
barplot(t(var.3.HVD[,3,3]), 
        main="HVD of the policy Rate", 
        legend="policy rate") 
barplot(t(var.3.HVD[,4,3]),
        main="HVD of the policy Rate", 
        legend="exchange rate") 

# option 2: you can export the data in Excel and plot the bar chart
write.csv(var.3.HVD[,,3], "var.3.HVD.csv")

# option 3: you can use ggplot (but you will have to convert the data into a data frame)
# https://blog.rstudio.com/2016/11/14/ggplot2-2-2-0/

##################################################################################
# STEP 8: OTHER TYPES OF RESTRICTIONS?                                           #
##################################################################################
# Zero long run restrictions: BQ command
var.3.BQ <- BQ(VAR.Swiss)
summary(var.3.BQ)
plot(irf(var.3.BQ))
plot(fevd(var.3.BQ))

# more information for the command: https://www.rdocumentation.org/packages/vars/versions/1.5-3/topics/BQ
# IMPOSE YOUR RESTRICTIONS
amat <- diag(4)
diag(amat) <- NA
amat[2, 1] <- NA
amat[3, 1] <- NA
amat[4, 1] <- NA

## Estimation method scoring
var.3.SVAR <- SVAR(x= VAR.Swiss, estmethod = "scoring", Amat = amat, Bmat = NULL,
                   max.iter = 100, maxls = 1000, conv.crit = 1.0e-8) 
plot(irf(var.3.SVAR))
plot(fevd(var.3.SVAR))

# more information for the command: https://www.rdocumentation.org/packages/vars/versions/1.5-3/topics/SVAR

#################################################################################
# Figures for Slides 
# choose colors: http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
# control the size of all pictures
resize.win <- function(Width=8, Height=6){windows(record=TRUE, width=Width, height=Height)}

resize.win(8,3) 
par(mar=c(3.1,3.1,1.7,1.7), cex=1.5) # control the margins and the font size
plot(VAR.Swiss[,1], lwd=3, xlab="", ylab="", las=1, col="dodgerblue3")
abline(h=mean(data.var[,1]), col="black", lwd=3, lty=3, ylab.las=2)
mtext("%", side=2, line=2.1, las=1, cex=1.5)

resize.win(8,3) 
par(mar=c(3.1,3.1,1.7,1.7), cex=1.5) # control the margins and the font size
plot(data.var[,2], lwd=3, xlab="", ylab="", las=1, col="dodgerblue3", ylim=c(0,20))
abline(h=mean(data.var[,2]), col="black", lwd=3, lty=3, ylab.las=2)
mtext("%", side=2, line=2.1, las=1, cex=1.5)

resize.win(8,3) 
par(mar=c(3.1,3.1,1.7,1.7), cex=1.5) # control the margins and the font size
plot(data.var[,3], lwd=3, xlab="", ylab="", las=1, col="dodgerblue3", ylim=c(0,20))
abline(h=mean(data.var[,3]), col="black", lwd=3, lty=3, ylab.las=2)
mtext("%", side=2, line=2.1, las=1, cex=1.5)


###############################################################################
# Simulation of a VAR
set.seed(1000) # so the simulation result can be duplicated
n   <- 200 # sample size = 200
z   <- as.matrix(cbind(rep(0, n),rep(0, n)))
w   <- as.matrix(cbind(rnorm(n), rnorm(n)))
phi <- as.matrix(cbind(c(0.3, 0.5), c(0, 0.6)))
for (i in 2:n) {
  z[i,] <- phi %*% z[i-1,] + w[i,]
}

adf.test(z[,1], k = 1)
adf.test(z[,1], k = 4)
adf.test(z[,2], k = 1)
adf.test(z[,2], k = 4)
# In this case, unit root is rejected (with p-value less than 0.05), so VAR can be applied. 
# If unit roots cannot be rejected then
# 1. Run VAR using differenced data if series are not cointegrated
# 2. Run Error Correction Model is series are cointegrated

# In theory we should include sufficient lagged values so that the error
# term is serially uncorrelated. In practice, we can choose p by minimizing AIC
VARselect(z)
#In this case, the VAR(1) has the smallest AIC -0.06176601, and so is chosen.

# 1. We start with a VAR(1) with both intercept term and trend term. Then we see both terms are insignificant.
var.1c <- VAR(z, p=1, type = "both")
summary(var.1c)

# 2. Next we run VAR(1) without intercept term or trend.
var.1c <- VAR(z, p=1, type = "none")
summary(var.1c)

# We can get the same result by using these commands to estimate OLS standard regression
y     <- z[,1]
ylag1 <- c(NA, y[1:n-1])
x     <- z[,2]
xlag1 <- c(NA, x[1:n-1])
eqy   <- lm(y~ylag1+xlag1-1) # without intercept
summary(eqy)

# To get omega hat we need to keep the residuals and then compute two variances and one covariance 
# (with adjustment of degree of freedom)
resy <- var.1c$varresult$y1$residuals
resx <- var.1c$varresult$y2$residuals

sigmau2  <- var(resy)*(n-1)/(n-2)
sigmav2  <- var(resx)*(n-1)/(n-2)
sigmauv  <- cov(resy,resx)*(n-1)/(n-2)
omegahat <- as.matrix(cbind(c(sigmau2, sigmauv), c(sigmauv, sigmav2)))

Aprime <- chol(omegahat)
A      <- t(Aprime)
#The matrix A will be used to construct impulse response to orthogonal (or structural form) errors

# The responses of y and x to the one-unit impulse of the orthogonal (structural-form) the error for y can be obtained as
phihat     <- matrix(0,2,2)
phihat[1,] <- var.1c$varresult$y1$coefficients
phihat[2,] <- var.1c$varresult$y2$coefficients

wtilde <- as.matrix(c(1,0))
A1c    <- A%*%wtilde # or A1c = A[,1]
if1    <- A1c
if2    <- phihat%*%if1
if3    <- phihat%*%if2
if4    <- phihat%*%if3

# You can get the same result using the R command
var1c.11 <- irf(var.1c, impulse = "y1", response="y1", boot=TRUE)
var1c.11
plot(var1c.11)
var1c.21 <- irf(var.1c, impulse = "y1", response="y2", boot=TRUE)
var1c.21
plot(var1c.21)

# The R command to test for Granger causality is
causality(var.1c, cause = c("y2"))
causality(var.1c, cause = c("y1"))