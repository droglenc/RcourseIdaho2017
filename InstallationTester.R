## Run this script to test the installation of packages needed for the
## workshop. If all packates are installed correctly this entire script
## should run without error (you may get some warnings and messages, but
## there should be no errors). Contact me (derek@derekogle.com) if you
## have any questions.

library(FSA)
library(FSAdata)

binCI(7,20)
poiCI(12)
chooseColors("rich",5)
data(BlackDrum2001)
dunnTest(tl~sex,data=BlackDrum2001)
tmp <- filterD(BlackDrum2001,sex %in% c("female","male"))
lm1 <- lm(weight~tl,data=BlackDrum2001)
fitPlot(lm1)
residPlot(lm1)
lm2 <- lm(tl~sex,data=BlackDrum2001)
fitPlot(lm2)
residPlot(lm2)
nls1 <- nls(tl~Linf*(1-exp(-K*(otoage-t0))),data=tmp,
            start=list(Linf=1200,K=0.15,t0=-1.5))
fitPlot(nls1)
residPlot(nls1)
nls2 <- nls(tl~Linf[sex]*(1-exp(-K[sex]*(otoage-t0[sex]))),data=tmp,
            start=list(Linf=c(1200,1200),K=c(0.15,0.15),t0=c(-1.5,-1.5)))
lrt(nls1,com=nls2)

library(nlstools)
boot1 <- nlsBoot(nls1,niter=10)
confint(boot1,plot=TRUE)

library(AICcmodavg)
aictab(list(nls1,nls2),modnames=c("omega","ultimate full"))

library(minpack.lm)
nls3 <- nlsLM(tl~Linf*(1-exp(-K*(otoage-t0))),data=tmp,
              start=list(Linf=1200,K=0.15,t0=-1.5))

library(car)
library(dplyr)
library(dunn.test)
library(epitools)
library(gplots)
library(lmtest)
library(plotrix)
library(plyr)
library(sciplot)