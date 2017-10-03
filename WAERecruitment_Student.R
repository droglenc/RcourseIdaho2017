########################################################################
##
##  Initial Preparations
##
########################################################################
# Clear console and global environment
rm(list = ls())
cat("\014")     # or ctrl-L in RStudio
# Load packages
library(FSA)
library(FSAdata)
library(dplyr)
library(magrittr)
library(minpack.lm)
library(AICcmodavg)



########################################################################
##
##  Initial Data Wrangling (and Quick Summary)
##
########################################################################
# Load data (from the FSAdata package)
data(WalleyeEL)
str(WalleyeEL)
WalleyeEL %<>% mutate(logage0=log(age0))

plot(age0~age5,data=WalleyeEL)


########################################################################
##
##  Fit Stock-Recruit Models
##
########################################################################
# Create R functions for each stock-recruitment function
SRi <- srFuns("independence")
SRbh <- srFuns("BevertonHolt")
SRr <- srFuns("Ricker")
SRs <- srFuns("Shepherd")
SRsl <- srFuns("SailaLorda")

# Get initial values
SVi <- srStarts(age0~age5,data=WalleyeEL,type="independence",plot=TRUE)
SVbh <- srStarts(age0~age5,data=WalleyeEL,type="BevertonHolt",plot=TRUE)
SVbh <- srStarts(age0~age5,data=WalleyeEL,type="BevertonHolt",
                 fixed=list(a=30,b=0.001),plot=TRUE)
SVr <- srStarts(age0~age5,data=WalleyeEL,type="Ricker",plot=TRUE)
SVs <- srStarts(age0~age5,data=WalleyeEL,type="Shepherd",plot=TRUE)
SVs <- srStarts(age0~age5,data=WalleyeEL,type="Shepherd",
                 fixed=list(a=30,b=0.001,c=1),plot=TRUE)
SVsl <- srStarts(age0~age5,data=WalleyeEL,type="SailaLorda",plot=TRUE)

# Fit using Levenberg-Marquardt algorithm (less sensitive to starts)
SRi.fit <- nlsLM(logage0~log(SRi(age5,a)),data=WalleyeEL,start=SVi)
SRbh.fit <- nlsLM(logage0~log(SRbh(age5,a,b)),data=WalleyeEL,start=SVbh)
SRr.fit <- nlsLM(logage0~log(SRr(age5,a,b)),data=WalleyeEL,start=SVr)
SRs.fit <- nlsLM(logage0~log(SRs(age5,a,b,c)),data=WalleyeEL,start=SVs)
SRsl.fit <- nlsLM(logage0~log(SRsl(age5,a,b,c)),data=WalleyeEL,start=SVsl)

# Compare models
ms <- list(SRi.fit,SRbh.fit,SRr.fit,SRs.fit,SRsl.fit)
mnames <- c("Independence","Beverton-Holt","Ricker","Shepherd","Saila-Lorda")
aictab(ms,modnames=mnames)

# Examine model fits ... ugh
clrs <- c("black","black","red","blue","orange")
plot(age0~age5,data=WalleyeEL,pch=19,
     xlim=c(0,3500),ylim=c(0,40000),
     xlab="Number of Age-5 and Older",ylab="Number of Age-0")
curve(SRi(x,coef(SRi.fit)),from=0,to=4000,add=TRUE,lty=2,col=clrs[1])
curve(SRbh(x,coef(SRbh.fit)),from=0,to=4000,add=TRUE,lwd=2,col=clrs[2])
curve(SRr(x,coef(SRr.fit)),from=0,to=4000,add=TRUE,lwd=2,col=clrs[3])
curve(SRs(x,coef(SRs.fit)),from=0,to=4000,add=TRUE,lwd=2,col=clrs[4])
curve(SRsl(x,coef(SRsl.fit)),from=0,to=4000,add=TRUE,lwd=2,col=clrs[5])
legend("topright",mnames,col=clrs,lty=c(2,1,1,1,1),lwd=c(1,2,2,2,2),
       bty="n",cex=0.9)


# --- ASSIGNMENT -- Fit models without the year with most age-5 fish----




