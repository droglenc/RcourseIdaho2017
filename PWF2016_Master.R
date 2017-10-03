########################################################################
##
##  Initial Preparation
##
########################################################################
# Clear console and global environment
rm(list = ls())
cat("\014")     # or ctrl-L in RStudio
# Load packages
library(FSA)
library(FSAdata)
library(nlstools)
library(AICcmodavg)
library(dplyr)
library(magrittr)



########################################################################
##
##  Initial Data Wrangling (and Quick Summaries)
##
########################################################################
# Load data
setwd("C:/aaaWork/Web/GitHub/RCourseIdaho2017")
pwf <- read.csv("PWF2016.csv")
str(pwf)
headtail(pwf)
levels(pwf$region)

# === DEMO ONLY -- 4 WAYS TO ADD VARIABLE TO DATA.FRAME ================
tmp <- pwf
tmp$logtl <- log(tmp$tl)
tmp <- mutate(tmp,logtl=log(tl))
tmp <- tmp %>% mutate(logtl=log(tl))
tmp %<>% mutate(logtl=log(tl))
rm(tmp)  # cleanup memory
# --- DEMO ONLY END ====================================================

# Modifiy data.frame
#   Re-order levels of region variable
#   Add country and length bins variables
#   Add logs of tl and wt
#   Remove towID and date variables (for clean-up)
#   Sort individuals by country, then region, then age, then tl
pwf %<>% mutate(region=factor(region,levels=c("Ontario-East","Ontario-West",
                                              "Michigan","Minnesota","Wisconsin")),
                country=mapvalues(region,from=levels(region),
                                  to=c("Canada","Canada","USA","USA","USA")),
                lcat=lencat(tl,w=10),
                logtl=log(tl),logwt=log(wt)) %>%
  select(-towID,-date) %>%
  arrange(country,region,age,tl)

headtail(pwf,n=6)
levels(pwf$country)

## Very quick summaries
Summarize(~tl,data=pwf,digits=1)
Summarize(age~region,data=pwf,digits=1)
xtabs(~region,data=pwf)
xtabs(~region+age,data=pwf)



########################################################################
##
##  Create Age-Length Key
##
########################################################################
# Get aged USA fish
USA.aged <- filterD(pwf,country=="USA",!is.na(age))
headtail(USA.aged)
# Compute two-way frequency table of length bins (rows) and ages
USA.raw <- xtabs(~lcat+age,data=USA.aged)
addmargins(USA.raw)
# Compute row proportions table (this is the ALK)
USA.alk <- prop.table(USA.raw,margin=1)
round(USA.alk*100,1)
alkPlot(USA.alk)
alkPlot(USA.alk,type="lines")
alkPlot(USA.alk,type="bubble")


# --- ASSIGNMENT -- Make a similar ALK for Canadian fish ---------------
CAN.aged <- filterD(pwf,country=="Canada",!is.na(age))
head(CAN.aged)
CAN.raw <- xtabs(~lcat+age,data=CAN.aged)
addmargins(CAN.raw)
CAN.alk <- prop.table(CAN.raw,margin=1)
round(CAN.alk*100,1)
alkPlot(CAN.alk)
# --- ASSIGNMENT END ---------------------------------------------------



########################################################################
##
##  Apply Age-Length Key (assign ages to unaged fish)
##
########################################################################
# Get unaged USA fish
USA.unaged <- filterD(pwf,country=="USA",is.na(age))
headtail(USA.unaged)
# Examine length-frequency of unaged fish
xtabs(~lcat,data=USA.unaged)
# Assign ages to unaged fish
USA.unaged <- alkIndivAge(USA.alk,age~tl,data=USA.unaged)
headtail(USA.unaged)
# Combine previously aged (from strux) and newly aged (from ALK) fish
USA <- rbind(USA.aged,USA.unaged)
headtail(USA,n=5)
# clean-up temporary, intermediate, or special use objects
rm(USA.aged,USA.unaged,USA.raw)


# -- ASSIGNMENT -- Construct data.frame of ALL Canadian fish with ages
CAN.unaged <- filterD(pwf,country=="Canada",is.na(age))
head(CAN.unaged)
xtabs(~lcat,data=CAN.unaged)
CAN.unaged <- alkIndivAge(CAN.alk,age~tl,data=CAN.unaged)
CAN <- rbind(CAN.aged,CAN.unaged)
# clean-up temporary, intermediate, or special use objects
rm(CAN.aged,CAN.unaged,CAN.raw)
# --- ASSIGNMENT END ---------------------------------------------------

# -- ASSIGNMENT -- Construct data.frame of ALL fish with ages
pwf <- rbind(CAN,USA)
# --- ASSIGNMENT END ---------------------------------------------------



########################################################################
##
##  Estimate Mortality Rate
##
########################################################################
# Compute age frequency (and include the log frequency)
USA.af <- group_by(USA,age) %>%
  summarise(freq=n()) %>%
  mutate(logfreq=log(freq)) %>%
  as.data.frame()
USA.af
str(USA.af)

# plot catch-curve
plot(logfreq~age,data=USA.af)

# Catch-curve analysis (using "first principles")
USA.af.gte2 <- filterD(USA.af,age>=2)
USA.cc1 <- lm(logfreq~age,data=USA.af.gte2)
coef(USA.cc1)
confint(USA.cc1)
USA.af.gte2 %<>% mutate(wts=predict(USA.cc1))
USA.af.gte2
USA.cc2 <- lm(logfreq~age,data=USA.af.gte2,weights=wts)
cbind(Est=coef(USA.cc2),confint(USA.cc2))

# Catch-curve analysis (using convenience function)
USA.cc1 <- catchCurve(freq~age,data=USA.af,ages2use=2:7)
cbind(Est=coef(USA.cc1),confint(USA.cc1))
plot(USA.cc1)

USA.cc2 <- catchCurve(freq~age,data=USA.af,ages2use=2:7,weighted=TRUE)
cbind(Est=coef(USA.cc2),confint(USA.cc2))

# Chapman-Robson analysis (using convenience function)
USA.cr <- chapmanRobson(freq~age,data=USA.af,ages2use=2:7)
cbind(Est=coef(USA.cr),confint(USA.cr))
plot(USA.cr)


# -- ASSIGNMENT -- Compute Z/A for ALL Canadian fish (use FSA functions)
CAN.af <- group_by(CAN,age) %>%
  summarise(freq=n()) %>%
  mutate(logfreq=log(freq)) %>%
  as.data.frame()
CAN.af

plot(logfreq~age,data=CAN.af)

CAN.cc <- catchCurve(freq~age,data=CAN.af,ages2use=2:9,weighted=TRUE)
cbind(Est=coef(CAN.cc),confint(CAN.cc))
plot(CAN.cc)

CAN.cr <- chapmanRobson(freq~age,data=CAN.af,ages2use=2:9)
cbind(Est=coef(CAN.cr),confint(CAN.cr))
plot(CAN.cr)
# --- ASSIGNMENT END ---------------------------------------------------



########################################################################
##
##  Compare Mortality Rates
##
########################################################################
# Compute age frequency (separated by locations)
ALL.af <- group_by(pwf,country,age) %>%
  summarise(freq=n()) %>%
  mutate(logfreq=log(freq)) %>%
  as.data.frame()
ALL.af

# Plot catch-curve
clrs <- c("black","blue")
plot(logfreq~age,data=ALL.af,pch=19,col=clrs[as.numeric(country)])

# Fit weighted IVR model to compare slopes
ALL.af.gte2 <- filterD(ALL.af,age>=2)
ALL.cc1 <- lm(logfreq~age*country,data=ALL.af.gte2)
ALL.af.gte2 %<>% mutate(wts=predict(ALL.cc1))
ALL.af.gte2
ALL.cc2 <- lm(logfreq~age*country,data=ALL.af.gte2,weights=wts)
anova(ALL.cc2)
cbind(Est=coef(ALL.cc2),confint(ALL.cc2))

# Fancy plot
plot(logfreq~age,data=ALL.af,col=clrs[as.numeric(country)],
     xlab="Age (years)",ylab="log(Frequency of Fish)")
points(logfreq~age,data=filterD(ALL.af,age>=2),
       pch=19,col=clrs[as.numeric(country)])
tmp <- c(2,9)
lines(tmp,predict(ALL.cc2,data.frame(age=tmp,country="Canada")),
      col=clrs[1],lwd=2)
tmp <- c(2,7)
lines(tmp,predict(ALL.cc2,data.frame(age=tmp,country="USA")),
      col=clrs[2],lwd=2)
legend("topright",levels(pwf$country),col=clrs,pch=19,lwd=2,bty="n")

# clean-up temporary, intermediate, or special use objects
rm(tmp,ALL.af.gte2,USA.af.gte2,USA.cc1,ALL.cc1)



########################################################################
##
##  Fit Growth Model
##
########################################################################
# Create axis label objects to save typing
xlbl <- "Age (yrs)"
ylbl <- "Total Length (mm)"
# Create colors that are transparent for help with overplotting
clrs2 <- col2rgbt(clrs,1/20)

# Examine lenght-at-age plot
plot(tl~age,data=USA,pch=19,col=clrs2[2],xlab=xlbl,ylab=ylbl)

# Create function with "typical" von Bertalanffy growth function
vb <- vbFuns("Typical",msg=TRUE)
vb

# === DEMO ONLY -- HOW TO USE vb() =====================================
vb(4,Linf=150,K=0.3,t0=-1)
vb(4,c(150,0.3,-1))
vb(1:5,c(150,0.3,-1))
# === DEMO ONLY END ====================================================

# Generate starting values ... automatic
USA.vbs <- vbStarts(tl~age,data=USA,type="Typical",plot=TRUE)
# Generate starting values ... manual iteration
USA.vbs <- vbStarts(tl~age,data=USA,type="Typical",plot=TRUE,
                    fixed=list(Linf=160,K=0.3,t0=0))
USA.vbs

# Fit the VBGF to data
USA.vbf <- nls(tl~vb(age,Linf,K,t0),data=USA,start=USA.vbs)
# Examine residual plot
residPlot(USA.vbf)

# Extract results
summary(USA.vbf,correlation=TRUE)
( USA.vbc <- coef(USA.vbf) )
cbind(Est=USA.vbc,confint(USA.vbf))

# Bootstrap estimates of uncertainty
USA.vbb <- nlsBoot(USA.vbf,niter=999)
str(USA.vbb)
headtail(USA.vbb$coefboot)
cbind(EST=USA.vbc,confint(USA.vbb,plot=TRUE,rows=1,cols=3))

# Uncertainty on predictions
ageX <- 4
predict(USA.vbf,data.frame(age=ageX))
USA.vbbp <- apply(USA.vbb$coefboot,MARGIN=1,FUN=vb,t=ageX)
c(pred=predict(USA.vbf,data.frame(age=ageX)),
  quantile(USA.vbbp,c(0.025,0.975)))

plot(tl~age,data=USA,xlab=xlbl,ylab=ylbl,pch=19,col=clrs2[2])
curve(vb(x,USA.vbc),from=1,to=8,n=500,lwd=2,col=clrs[2],add=TRUE)


# - ASSIGNMENT -- Fit "typical" VBGF to ALL Canadian fish --------------
( CAN.vbs <- vbStarts(tl~age,data=CAN,type="Typical",plot=TRUE) )

CAN.vbf <- nls(tl~vb(age,Linf,K,t0),data=CAN,start=CAN.vbs)
residPlot(CAN.vbf)

summary(CAN.vbf,correlation=TRUE)
( CAN.vbc <- coef(CAN.vbf) )
CAN.vbb <- nlsBoot(CAN.vbf,niter=999)
cbind(EST=CAN.vbc,confint(CAN.vbb))

CAN.vbbp <- apply(CAN.vbb$coefboot,MARGIN=1,FUN=vb,t=ageX)
c(pred=predict(CAN.vbf,data.frame(age=ageX)),
  quantile(CAN.vbbp,c(0.025,0.975)))

plot(tl~age,data=CAN,xlab=xlbl,ylab=ylbl,pch=19,col=clrs2[1],ylim=c(20,130))
curve(vb(x,CAN.vbc),from=1,to=9,n=500,lwd=2,col=clrs[1],add=TRUE)
# --- ASSIGNMENT END ---------------------------------------------------



########################################################################
##
##  Compare Growth Model Parameters
##
########################################################################
plot(tl~age,data=pwf,pch=19,col=clrs2[country],xlab=xlbl,ylab=ylbl)

( svOm <- vbStarts(tl~age,data=pwf) )
( svLKt <- Map(rep,svOm,c(2,2,2)) )

vbLKt <- tl~Linf[country]*(1-exp(-K[country]*(age-t0[country])))
fitLKt <- nls(vbLKt,data=pwf,start=svLKt)
residPlot(fitLKt,col=col2rgbt("black",1/3))

vbOm <- tl~Linf*(1-exp(-K*(age-t0)))
fitOm <- nls(vbOm,data=pwf,start=svOm)

extraSS(fitOm,com=fitLKt,sim.name="{Omega}",com.name="{Linf,K,t0}")
lrt(fitOm,com=fitLKt,sim.name="{Omega}",com.name="{Linf,K,t0}")

vbLK <- tl~Linf[country]*(1-exp(-K[country]*(age-t0)))
( svLK <- Map(rep,svOm,c(2,2,1)) )
fitLK <- nls(vbLK,data=pwf,start=svLK)
vbLt <- tl~Linf[country]*(1-exp(-K*(age-t0[country])))
svLt <- Map(rep,svOm,c(2,1,2))
fitLt <- nls(vbLt,data=pwf,start=svLt)
vbKt <- tl~Linf*(1-exp(-K[country]*(age-t0[country])))
svKt <- Map(rep,svOm,c(1,2,2))
fitKt <- nls(vbKt,data=pwf,start=svKt)
extraSS(fitLK,fitLt,fitKt,com=fitLKt,com.name="{Linf,K,t0}",
        sim.names=c("{Linf,K}","{Linf,t0}","{K,t0}"))

vbL <- tl~Linf[country]*(1-exp(-K*(age-t0)))
( svL <- Map(rep,svOm,c(2,1,1)) )
fitL <- nls(vbL,data=pwf,start=svL)
vbt <- tl~Linf*(1-exp(-K*(age-t0[country])))
svt <- Map(rep,svOm,c(1,1,2))
fitt <- nls(vbt,data=pwf,start=svt)
extraSS(fitL,fitt,com=fitLt,com.name="{Linf,t0}",sim.names=c("{Linf}","{t0}"))

summary(fitLt,correlation=TRUE)
round(cbind(Est=coef(fitLt),confint(fitLt)),3)

vbK <- tl~Linf*(1-exp(-K[country]*(age-t0)))
svK <- Map(rep,svOm,c(1,2,1))
fitK <- nls(vbK,data=pwf,start=svK)

ms <- list(fitOm,fitL,fitK,fitt,fitLK,fitLt,fitKt,fitLKt)
mnames <- c("{Omega}","{Linf}","{K}","{t0}","{Linf,K}","{Linf,t0}","{K,t0}","{Linf,K,t0}")
aictab(ms,mnames)

( cfLt <- coef(fitLt) )
jit <- 0.05
plot(tl~I(age-jit),data=filterD(pwf,country=="Canada"),
     pch=19,col=clrs2[1],xlab=xlbl,ylab=ylbl,ylim=c(0,160),xlim=c(0,10))
points(tl~I(age+jit),data=filterD(pwf,country=="USA"),
       pch=19,col=clrs2[2])
curve(vb(x,cfLt[c("Linf1","K","t01")]),from=0,to=10,add=TRUE,col=clrs[1],lwd=2)
curve(vb(x,cfLt[c("Linf2","K","t02")]),from=0,to=10,add=TRUE,col=clrs[2],lwd=2)
legend("topleft",levels(pwf$country),col=clrs,pch=19,lwd=2,bty="n")



########################################################################
##
##  Compute Weight-Length Relationship
##
########################################################################
clrs2 <- col2rgbt(clrs,1/3)
plot(wt~tl,data=USA,pch=19,col=clrs2[1])
plot(logwt~logtl,data=USA,pch=19,col=clrs2[1])

USA.lw <- lm(logwt~logtl,data=USA)
cbind(Est=coef(USA.lw),confint(USA.lw))

lens <- seq(30,150,length.out=100)
USA.pwt <- exp(predict(USA.lw,data.frame(logtl=log(lens)),
                       interval="prediction"))
head(USA.pwt)
plot(wt~tl,data=USA,pch=19,col=clrs2[2],
     xlab="Total Length (mm)",ylab="Weight (g)")
lines(USA.pwt[,"fit"]~lens,lwd=2,col=clrs[2])
lines(USA.pwt[,"lwr"]~lens,lwd=2,lty=2,col=clrs[2])
lines(USA.pwt[,"upr"]~lens,lwd=2,lty=2,col=clrs[2])


# - ASSIGNMENT -- Fit weight-length relationship to Canadian fish ------
CAN.lw <- lm(logwt~logtl,data=CAN)
cbind(Est=coef(CAN.lw),confint(CAN.lw))

lens <- seq(30,150,length.out=100)
CAN.pwt <- exp(predict(CAN.lw,data.frame(logtl=log(lens)),
                       interval="prediction"))
head(CAN.pwt)
plot(wt~tl,data=CAN,pch=19,col=clrs2[1],
     xlab="Total Length (mm)",ylab="Weight (g)")
lines(CAN.pwt[,"fit"]~lens,lwd=2)
lines(CAN.pwt[,"lwr"]~lens,lwd=2,lty=2)
lines(CAN.pwt[,"upr"]~lens,lwd=2,lty=2)
# --- ASSIGNMENT END ---------------------------------------------------



########################################################################
##
##  Compare Weight-Length Model Parameters
##
########################################################################
plot(logwt~logtl,data=pwf,pch=19,col=clrs2[country])

ALL.lw <- lm(logwt~logtl*country,data=pwf)
anova(ALL.lw)
cbind(Est=coef(ALL.lw),confint(ALL.lw))

plot(wt~tl,data=pwf,pch=19,col=clrs2[country])
CAN.pwt <- exp(predict(CAN.lw,data.frame(logtl=log(lens),country="Canada"),
                       interval="prediction"))
lines(CAN.pwt[,"fit"]~lens,lwd=2,col=clrs[1])
lines(CAN.pwt[,"lwr"]~lens,lwd=2,lty=2,col=clrs[1])
lines(CAN.pwt[,"upr"]~lens,lwd=2,lty=2,col=clrs[1])
USA.pwt <- exp(predict(USA.lw,data.frame(logtl=log(lens),country="USA"),
                       interval="prediction"))
lines(USA.pwt[,"fit"]~lens,lwd=2,col=clrs[2])
lines(USA.pwt[,"lwr"]~lens,lwd=2,lty=2,col=clrs[2])
lines(USA.pwt[,"upr"]~lens,lwd=2,lty=2,col=clrs[2])

lwCompPreds(ALL.lw,lens=c(40,70,100,130))



## Use makeStudentScript() to make student script from the master
source("zzzDontShare/Workshop_Helpers.R")
makeStudentScript("PWF2016")
