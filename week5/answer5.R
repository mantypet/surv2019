#########################################################################
# Exercise 5: Competing risks as a multi-state model
#########################################################################

library(survival)
library(eha)

###########################################
# (1) Simulatation of (competing risk) data
###########################################

# Parameters in a competing risk model with two absorbing states
N      = 1000; # Number of individuals
alpha1 = 0.2;  # constant hazard for transition 0 -> 1
alpha2 = 0.4;  # constant hazard for transition 0 -> 2
theta  = 1.5   # the effect of treatment (relative rate)

########################################################
# Simulate the event times and the event in the controls
########################################################

t1_control = rexp(N,alpha1); # "latent" time for event 0 -> 1
t2_control = rexp(N,alpha2); # "latent" time for event 0 -> 2
T          = 3;              # censoring time (end of follow-up)

# Select the smallest and determine the appropriate indicator value
t_control  = pmin(t1_control,t2_control,T);
D1         = 1*(t_control==t1_control); # event 0->1
D2         = 2*(t_control==t2_control); # event 0->2
D3         = 3*(t_control==T);          # censoring
D_control  = pmax(D1,D2,D3)             # the event indicator

# Repeat the above for the treated group
t1_treat = rexp(N,theta*alpha1); # "latent" time for event 0 -> 1
t2_treat = rexp(N,theta*alpha2); # "latent" time for event 0 -> 2
T        = 3;                    # censoring time (end of follow-up)

# Select the smallest and determine the appropriate indicator value
t_treat    = pmin(t1_treat,t2_treat,T);
D1         = 1*(t_treat==t1_treat); # event 0->1
D2         = 2*(t_treat==t2_treat); # event 0->2
D3         = 3*(t_treat==T);        # censoring
D_treat    = pmax(D1,D2,D3)         # the event indicator

###############################################
# Form the data matrix (one row per individual)
data_control   = cbind(rep(0,N),t_control,D_control,rep(0,N))
data_treat     = cbind(rep(0,N),t_treat,D_treat,rep(1,N));

data           = rbind(data_control,data_treat)
data           = as.data.frame(data)
colnames(data) = c("entry","exit","status","treatment")

##################################################################
# Analysis option (a): different baselines for the two transitions,
#                      different effects of treatment on both transitions
##################################################################

# Stratify the data according to the outcome being '1' or '2'

dataA1           = cbind(data$exit,data$status==1,data$treatment)
dataA1           = as.data.frame(dataA1)
colnames(dataA1) = c("exit","status","treatment")

dataA2           = cbind(data$exit,data$status==2,data$treatment)
dataA2           = as.data.frame(dataA2)
colnames(dataA2) = c("exit","status","treatment")

cA1 = weibreg(Surv(exit,status)~treatment,data=dataA1,shape=1);
cA2 = weibreg(Surv(exit,status)~treatment,data=dataA2,shape=1);

######################################################################
# Analysis option(b): different baselines for the two transitions,
#                     a common effect of treatment on both transitions
######################################################################
dataB1          = cbind(data$exit,data$status==1,data$treatment,1);
dataB2          = cbind(data$exit,data$status==2,data$treatment,2);
dataB           = rbind(dataB1,dataB2);
dataB           = as.data.frame(dataB);
colnames(dataB) = c("exit","status","treatment","stratum")

cB             = weibreg(Surv(exit,status)~treatment+as.factor(stratum),data=dataB,shape=1)
cB_alternative = coxph(Surv(exit,status)~treatment+strata(stratum),data=dataB);


##################################################################
# Analysis (c): proportional baseline hazards,
#               different effects of treatment on both transitions
##################################################################
dataC1          = cbind(data$exit,data$status==1,data$treatment,0,1);
dataC2          = cbind(data$exit,data$status==2,0,data$treatment,2);
dataC           = rbind(dataC1,dataC2)
dataC           = as.data.frame(dataC)
colnames(dataC) = c("exit","status","treatment1","treatment2","stratum")

cC            = weibreg(Surv(exit,status)~treatment1+treatment2+as.factor(stratum),data=dataC,shape=1);
cC_alternativ =  coxph(Surv(exit,status)~treatment1+treatment2+as.factor(stratum),data=dataC);







