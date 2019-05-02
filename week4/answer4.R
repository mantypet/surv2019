##########################################################################
# Exercise 4: Proportional hazards models and model check
#   
##########################################################################

# Download the required R packages
library(survival)
library(KMsurv)
library(eha)


###############################################################
# Q 3, 4, 5
###############################################################
library(timereg)
data(melanoma, package="boot")
melanoma$log2thick <- log2(melanoma$thickness)

# ulcer as the only covariate and plot Nelson-Aalen plots for patients with and without ulcer
fit.su=coxph(Surv(time,status==1)~strata(ulcer),data=melanoma)
surv.su=survfit(fit.su)
plot(surv.su,fun="cumhaz", mark.time=F, xlab="Years since operation",ylab="Cumulative hazard",lty=1:2)
legend("topleft",c("No ulcer","Ulcer"),lty=1:2)
# Then fit a Cox model with ulcer as the only covariate and plot the model based estmates of the cumulative hazards in the same plot:
fit.u=coxph(Surv(time,status==1)~ulcer,data=melanoma)
surv.u=survfit(fit.u,newdata=data.frame(ulcer=c(0,1)))
lines(surv.u,fun="cumhaz", mark.time=F,conf.int=F, lty=1:2,col="red")


###############################################################
# Q 6
###############################################################
# Consider the model with ulcer and log-thickness
fit.ut=coxph(Surv(time,status==1)~ulcer+log2(thickness),data=melanoma)
summary(fit.ut)
# cumulative hazards for the four covariate combinations
# 1) ulcer=0, thickness=1
# 2) ulcer=0, thickness=4
# 3) ulcer=1, thickness=4
# 3) ulcer=1, thickn=10
new.covariates=data.frame(ulcer=c(0,0,1,1),thickness=c(1,4,4,10))
surv.ut=survfit(fit.ut,newdata= new.covariates)
plot(surv.ut,fun="cumhaz", mark.time=F, xlab="Years since operation",ylab="Cumulative hazard",lty=1:4, main="Strata-specific cumulative hazards (ulcer, thickness)")
legend("topleft",c("(0,1)","(0,4)","(1,4)","(1,10)"), lty=1:4)
# To plot the survival functions for the same combinations
plot(surv.ut,mark.time=F, xlab="Years since peration",lty=1:4)
legend("bottomleft",c("(0,1)","(0,4)","(1,4)","(1,10)"), lty=1:4)

###############################################################
# Q 7 Stratified analysis
###############################################################
# Fit a model with strata defined by ulceration and log-thickness as covariate
fit.strat=coxph(Surv(time,status==1)~log2(thickness)+strata(ulcer), data=melanoma)
summary(fit.strat)

# Plot the cumuative baseline hazards for the two ulcer strata:
baseline.covar=data.frame(thickness=3)
surv.strat=survfit(fit.strat,newdata=baseline.covar)
plot(surv.strat,fun="cumhaz", mark.time=F,xlab="Years since operation",ylab="Cumulative hazard",lty=1:2)
legend("topleft",c("No ulcer","Ulcer"),lty=1:2)


###############################################################
# Q 8 Graphical check of proportionality
###############################################################

## these functions are additional to get output for all time points
cumhaz.cox=function(cox.surv)
{
no.strata=length(cox.surv$strata)
strata.ind=NULL
for (i in 1:no.strata){
strata.ind=c(strata.ind,rep(i,cox.surv$strata[i]))}
cumhaz.cox=vector('list',length=no.strata)
for (g in 1:no.strata) {
cumhaz.cox[[g]]=cox.surv$cumhaz[strata.ind==g]}
return(cumhaz.cox)
}
time.cox=function(cox.surv)
{
no.strata=length(cox.surv$strata)
strata.ind=NULL
for (i in 1:no.strata){
strata.ind=c(strata.ind,rep(i,cox.surv$strata[i]))}
time.cox=vector('list',length=no.strata)
for (g in 1:no.strata) {
time.cox[[g]]=cox.surv$time[strata.ind==g]}
return(time.cox)
}
##

# Then we check proportional hazards for ulceration:
par(mfrow=c(1,1))
fit.logtstu=coxph(Surv(time,status==1)~strata(ulcer)+log2thick, data=melanoma)
surv.logtstu=survfit(fit.logtstu,newdata=data.frame(log2thick=0))
plot(surv.logtstu[1]$time, log(-log(surv.logtstu[1]$surv)), type="l",ylim=c(-6,-1),xlab="Years since operation",ylab="log cumulative hazard",lty=1)
lines(surv.logtstu[2]$time, log(-log(surv.logtstu[2]$surv)),lty=2)
legend("topleft",c("No ulcer","Ulcer"),lty=1:2)

#cumhaz.logtstu=cumhaz.cox(surv.logtstu)
#time.logtstu=time.cox(surv.logtstu)
#plot(time.logtstu[[1]],log(cumhaz.logtstu[[1]]),type='s',ylim=c(-6,0),xlab="Years since operation",ylab="log cumulative hazard", main="Ulcer")
#lines(time.logtstu[[2]],log(cumhaz.logtstu[[2]]),type='s',lty=2)
# the cruves are fairly parallel

# proportionality for log2thick
# create groups from thickness
# look at
summary(melanoma$log2thick)
melanoma$grthick = 0
melanoma$grthick[melanoma$log2thick > 0.37 & melanoma$log2thick<=1.37] = 1
melanoma$grthick[melanoma$log2thick>1.37] = 2


fit.sttu=coxph(Surv(time,status==1)~factor(ulcer)+strata(grthick), data=melanoma)
surv.sttu=survfit(fit.sttu,newdata=data.frame(ulcer=0))
plot(surv.sttu[1]$time, log(-log(surv.sttu[1]$surv)), type="l",ylim=c(-4,0),xlab="Years since operation",ylab="log cumulative hazard",lty=1)
lines(surv.sttu[2]$time, log(-log(surv.sttu[2]$surv)),lty=2)
lines(surv.sttu[3]$time, log(-log(surv.sttu[3]$surv)),lty=3)
legend("topleft",c("[0,1]","(1,3]","(3,-]"),lty=1:3)

#cumhaz.sttu=cumhaz.cox(surv.sttu)
#time.sttu=time.cox(surv.sttu)
#plot(time.sttu[[1]],log(cumhaz.sttu[[1]]),type='s',xlim=c(0,10),ylim=c(-4,0.5),
#xlab="Years since operation", ylab="log cumulative hazard", main="Tumor thickness")
#lines(time.sttu[[2]],log(cumhaz.sttu[[2]]),type='s',lty=2)
#lines(time.sttu[[3]],log(cumhaz.sttu[[3]]),type='s',lty=3)

## There is a tendency that the curves are closer together for large values
## of t than for smaller ones, indicating a non-proportional effect of thickness

###############################################################
# Q 9 Formal test of proportionality
###############################################################
# This is done by adding time-dependent covariate x*log(t) for each covariate x,
# and testing whether the time-dependent covariates are significant:
fit.logtu=coxph(Surv(time,status==1)~factor(ulcer)+log2thick, data=melanoma)
cox.zph(fit.logtu,transform='log')

###### Not to be discussed
# We also make plots that give nonparametric estimates of the
# (possible) time dependent effect of the covariates:
par(mfrow=c(1,2))
plot(cox.zph(fit.logtu))
##############

