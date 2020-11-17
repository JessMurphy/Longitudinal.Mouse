########## LOAD LIBRARIES ##########

library(lme4)
library(lmerTest)
library(MASS) #mvnorm
library(binom) #binom.confint


########## READ DATA ##########

setwd("C:/Repositories/Longitudinal.Mouse/")

# load in data and necessary models
data_long = read.csv("./EasyLME/data/data_long.csv")
load("./models/mouse_slope_int.R")
load("./models/mouse_intercept.R")
load("./models/no_RE.R")


########## POWER ANALYSIS ##########

# define study design
n = 20  #number of subjects per group (2 groups)
n_time = 12 #number of time points
n_mice = 1 #number of mice per donor
nsim0 = 10000 #number of simulations
alpha = 0.05 #significance level

# extract time points
Time = unique(data_long$Time) #(0,1,3,4,7,11,14,18,21,25,28,32)

# define simulation parameters for random mouse intercept and slope model (lmod4)
mu_slope = summary(lmod4)$coefficients[,1] #(96.731, -0.189, 0.99, -0.425)
var_slope = as.data.frame(VarCorr(lmod4))
mouseInt = var_slope[1,4] #4.667
mouseSlope = var_slope[2,4] #0.1808
mouseCov = var_slope[3,4] #0.078
resid_slope = var_slope[4,4] #8.36

# calculate variance matrix of random slope & intercept model (lmod4)
G_slope = matrix(data=mouseCov, nrow=2, ncol=2)
G_slope[1,1] = mouseInt
G_slope[2,2] = mouseSlope
Z_slope = cbind(1, Time)
R_slope = diag(resid_slope, n_time)
V_slope = Z_slope%*%G_slope%*%t(Z_slope)+R_slope

# define simulation parameters for random mouse intercept model (lmod7)
mu_int = summary(lmod7)$coefficients[,1] #(95.948, 0.244, 1.024, -0.4)
var_int = as.data.frame(VarCorr(lmod7))
mouseINDonor_RE = var_int[1,4] #35.0577
resid_int = var_int[2,4] #27.65

# calculate variance matrix of random intercept model (lmod7)
R_int = diag(resid_int, n_time)
V_int = matrix(mouseINDonor_RE, ncol = n_time, nrow = n_time) 
V_int = V_int + R_int

# define simulation parameters for no random effects model (null)
mu_noRE = summary(null)$coefficients[,1] #(96.65, -0.814, 0.994, -0.325)
R_noRE = diag(sigma(null), n_time)

# combine the variances for each model into one list
V = list(V_slope, V_int, R_noRE) 

# build similuation design matrix
Subject = factor(sort(rep(1:(2*n), each=n_time*n_mice)))
Mice = rep(rep(1:n_mice, each=n_time), n*2)
Mice2 = factor(paste(Subject, Mice, sep="."))
group = rep(c("Healthy", "Undernourished"), each=n*n_time*n_mice)  
theData = data.frame(Time, Mice2, Subject)
names(theData) = c("Time", "Mice", "Subject")
theData$Group = factor(group)

Healthy = data.frame(Intercept=1, Group=rep(0, each=n_time*n_mice), Time=Time)
Healthy$Int = Healthy$Group * Healthy$Time
Undernourished = data.frame(Intercept=1, Group=rep(1, each=n_time*n_mice), Time=Time)
Undernourished$Int = Undernourished$Group * Undernourished$Time

# define mean structure
mu_slope0 = mu_slope
mu_slope0[4] = 0

# define lists to store outputs from simulation
power = list()
coefficients = list()
errors = list()
stats = list()

# set seed for reproducible results
set.seed(0)

# loop through the different models to simulate data from
# model 1: random mouse slope & intercept (lmod4)
# model 2: random mouse intercept (lmod7)
# model 3: no random effects (null)

for(j in 1:3){
  
  # use the appropriate sigma matrix for each model
  V_loop = V[[j]]
  
  # define vectors to store outputs from model loop
  pval_out_slope = pval_out_int = pval_out_noRE = c()
  slope_effect = int_effect = noRE_effect = c()
  slope_error = int_error = noRE_error = c()
  slope_stat = int_stat = noRE_stat = c()
  
  set.seed(0)
  
  # perform simulations
  for(i in 1:nsim0){
    
    # simulate new resposne values for each group based on the specified effect size and variance matrix
    y_healthy = as.vector(t(mvrnorm(n*n_mice, mu=t(as.vector(mu_slope0))%*%t(as.matrix(Healthy[1:n_time,])), Sigma=V_loop)))
    y_undernourished = as.vector(t(mvrnorm(n*n_mice, mu=t(as.vector(mu_slope0))%*%t(as.matrix(Undernourished[1:n_time,])), Sigma=V_loop)))
    
    # combine responses with the design matrix
    theData$growth = c(y_healthy, y_undernourished)
    
    # fit models 1-3 to the simulated data
    tmp_slope = lmer(growth ~ Group*Time + (Time|Group:Subject),  data = theData)
    tmp_int = lmer(growth ~ Group*Time + (1|Group:Subject),  data = theData)
    tmp_noRE = lm(growth ~ Group*Time,  data = theData)
    
    # extract the p-values of the interaction term for each model
    pval_slope = summary(tmp_slope)$coef[4,5]
    pval_out_slope = c(pval_out_slope, pval_slope)
    pval_int = summary(tmp_int)$coef[4,5]
    pval_out_int = c(pval_out_int, pval_int)
    pval_noRE = summary(tmp_noRE)$coef[4,4]
    pval_out_noRE = c(pval_out_noRE, pval_noRE)
    
    # extract the coefficients of the interaction term for each model
    interaction_slope = summary(tmp_slope)$coef[4,1]
    slope_effect = c(slope_effect, interaction_slope)
    interaction_int = summary(tmp_int)$coef[4,1]
    int_effect = c(int_effect, interaction_int)
    interaction_noRE = summary(tmp_noRE)$coef[4,1]
    noRE_effect = c(noRE_effect, interaction_noRE)
    
    # extract the standard error of the interaction term for each model
    se_slope = summary(tmp_slope)$coef[4,2]
    slope_error = c(slope_error, se_slope)
    se_int = summary(tmp_int)$coef[4,2]
    int_error = c(int_error, se_int)
    se_noRE = summary(tmp_noRE)$coef[4,2]
    noRE_error = c(noRE_error, se_noRE)
    
    # extract the t-statistic of the interaction term for each model
    t_slope = summary(tmp_slope)$coef[4,4]
    slope_stat = c(slope_stat, t_slope)
    t_int = summary(tmp_int)$coef[4,4]
    int_stat = c(int_stat, t_int)
    t_noRE = summary(tmp_noRE)$coef[4,3]
    noRE_stat = c(noRE_stat, t_noRE)
    
  } # end sim loop
  
  # calculate the power for the designated effect size per model
  power_slope = length(pval_out_slope[pval_out_slope <= alpha])/length(pval_out_slope)
  power_int = length(pval_out_int[pval_out_int <= alpha])/length(pval_out_int)
  power_noRE = length(pval_out_noRE[pval_out_noRE <= alpha])/length(pval_out_noRE)
  
  # combine all the power calculations simulated under the same dataset 
  power[[j]] = cbind(power_slope, power_int, power_noRE) 
  colnames(power[[j]]) = c("Slope & Intercept", "Intercept", "None")
  
  # combine the interaction coefficients simulated under the same dataset 
  coefficients[[j]] = cbind(slope_effect, int_effect, noRE_effect) 
  colnames(coefficients[[j]]) = c("Slope & Intercept", "Intercept", "None")
  
  # combine the standard errors simulated under the same dataset 
  errors[[j]] = cbind(slope_error, int_error, noRE_error) 
  colnames(errors[[j]]) = c("Slope & Intercept", "Intercept", "None")
  
  # combine the t-statistics simulated under the same dataset 
  stats[[j]] = cbind(slope_stat, int_stat, noRE_stat) 
  colnames(stats[[j]]) = c("Slope & Intercept", "Intercept", "None")
  
  # print the model number
  print(j)
  
} # end model loop

# warning (2 and 3): boundary (singular) fit: see ?isSingular
# warning message: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#                  Model failed to converge with max|grad| = 0.00313721 (tol = 0.002, component 1)

names(power) = c("Slope & Intercept", "Intercept", "None")
names(coefficients) = c("Slope & Intercept", "Intercept", "None")
names(errors) = c("Slope & Intercept", "Intercept", "None")
names(stats) = c("Slope & Intercept", "Intercept", "None")

# confidence interval function
confint <- function(x){
  # calculate confidence interval
  interval = binom.confint(x*nsim0, nsim0, method="exact")
  # round to two decimal places
  paste0("[", round(interval[1,5], 3), ", ", round(interval[1,6], 3), "]")
}

# confidence intervals for each simulation scenario
# c(1,2) means apply the function over rows and columns
apply(power[[1]], c(1,2), confint) 
apply(power[[2]], c(1,2), confint)
apply(power[[3]], c(1,2), confint)
