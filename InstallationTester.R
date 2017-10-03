## Run this script to test the installation of packages needed for the
## workshop. If all packages are installed correctly, then this entire
## script should run without error (you may get some warnings and
## messages, but there should be no errors). Contact me
## (derek@derekogle.com) if you have any questions.

## Checks FSA install and all packages that it depends on
source(system.file("helpers","InstallTester.R",package="FSA"),echo=TRUE)

## Check extra packages required for the workshop
data(SpotVA1)

## Test magrittr functions
library(magrittr)
SpotVA1 %<>% mutate(sex=factor(sample(c("F","M"),nrow(SpotVA1),replace=TRUE)))

## Test nlstools functions
nls1 <- nls(tl~Linf*(1-exp(-K*(age-t0))),data=SpotVA1,
            start=list(Linf=16,K=0.22,t0=-2))
nls2 <- nls(tl~Linf[sex]*(1-exp(-K[sex]*(age-t0[sex]))),data=SpotVA1,
            start=list(Linf=c(16,16),K=c(0.22,0.22),t0=c(-2,-2)))
library(nlstools)
boot1 <- nlsBoot(nls1,niter=10)
confint(boot1,plot=TRUE)

## Test AICcmodavg functions
library(AICcmodavg)
aictab(list(nls1,nls2),modnames=c("omega","ultimate full"))

## Test minpack.lm functions
library(minpack.lm)
nls3 <- nlsLM(tl~Linf*(1-exp(-K*(age-t0))),data=SpotVA1,
              start=list(Linf=16,K=0.22,t0=-2))
