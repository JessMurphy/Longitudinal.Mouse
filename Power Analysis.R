#################### LOAD LIBRARIES ####################

library(lme4)
library(lmerTest)
library(MASS) #mvnorm
library(binom) #binom.confint
library(tidyverse)
library(ggpubr)


#################### READ DATA ####################

#setwd("C:/Repositories/Longitudinal.Mouse/")

# load in data and necessary models
data_long = read.csv("C:/Repositories/EasyLME/data/data_long.csv")
load("./models/mouse_slope_int.R")
load("./models/mouse_intercept.R")
load("./models/no_RE.R")


#################### POWER ANALYSIS ####################

# define study design
n = 4  #number of subjects per group (2 groups)
n_time = 12 #number of time points
n_mice = 5 #number of mice per donor
nsim = 1000 #number of simulations
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
mouseINDonor_RE = var_slope[1,4] #var_int[1,4] #35.0577
resid_int = var_slope[4,4] #var_int[2,4] #27.65

# calculate variance matrix of random intercept model (lmod7)
R_int = diag(resid_int, n_time)
V_int = matrix(mouseINDonor_RE, ncol = n_time, nrow = n_time) 
V_int = V_int + R_int

# define simulation parameters for no random effects model (null)
mu_noRE = summary(null)$coefficients[,1] #(96.65, -0.814, 0.994, -0.325)
R_noRE = diag(sigma(null), n_time) #7.85

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

# define the interaction effect sizes to loop through
# use the mean effects from lmod4 (mouse random slope & intercept) 
effect_size = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.5)*mu_slope[4]

# define lists to store outputs from simulation
power = list()

# set seed for reproducible results
set.seed(0)

# loop through the three different models to simulate data from
# model 1: random mouse slope & intercept (lmod4)
# model 2: random mouse intercept (lmod7)
# model 3: no random effects (null)

for(j in 1:3){
  
  # use the appropriate sigma matrix for each model
  V_loop = V[[j]]
  
  # define vectors to store outputs from model loop
  power_out_slope = power_out_int = power_out_noRE = c()
  
  # loop through different effect sizes for the interaction term
  for(int in effect_size){
    
    # set the interaction term equal to the designated effect size
    mu_slope[4] = int
    
    # define vectors to store outputs from effect size loop
    pval_out_slope = pval_out_int = pval_out_noRE = c()
    
    # perform simulations
    for(i in 1:nsim){
      
      # simulate new response values for each group based on the specified effect size and variance matrix
      y_healthy = as.vector(t(mvrnorm(n*n_mice, mu=t(as.vector(mu_slope))%*%t(as.matrix(Healthy[1:n_time,])), Sigma=V_loop)))
      y_undernourished = as.vector(t(mvrnorm(n*n_mice, mu=t(as.vector(mu_slope))%*%t(as.matrix(Undernourished[1:n_time,])), Sigma=V_loop)))
      
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
      
    } # end sim loop
    
    # calculate the power for the designated effect size per model
    power_slope = length(pval_out_slope[pval_out_slope <= alpha])/length(pval_out_slope)
    power_int = length(pval_out_int[pval_out_int <= alpha])/length(pval_out_int)
    power_noRE = length(pval_out_noRE[pval_out_noRE <= alpha])/length(pval_out_noRE)
    
    # combine the power for each effect size fit to the same model
    power_out_slope = rbind(power_out_slope, power_slope)
    power_out_int = rbind(power_out_int, power_int)
    power_out_noRE = rbind(power_out_noRE, power_noRE)
    
    # print the effect_size 
    print(int)
    
  } # end mu loop
  
  # combine all the power calculations simulated under the same dataset 
  power[[j]] = cbind(power_out_slope, power_out_int, power_out_noRE) 
  colnames(power[[j]]) = c("Mouse Slope/Int", "Mouse Intercept", "No Random Effects") 
  rownames(power[[j]]) = c("0%", "5%", "10%", "15%", "20%", "25%", "50%")
  
  # print the model number
  print(j)
  
} # end model loop

# warning (all models): boundary (singular) fit: see ?isSingular

models = c("Mouse Intercept/Slope", "Mouse Intercept", "No Random Effects")
names(power) = models

# confidence interval function
confint <- function(x){
  # calculate confidence interval
  interval = binom.confint(x*nsim, nsim, method="exact")
  # round to two decimal places
  paste0("[", round(interval[1,5], 2), ", ", round(interval[1,6], 2), "]")
}

# confidence intervals for each simulation scenario
# c(1,2) means apply the function over rows and columns
apply(power[[1]], c(1,2), confint)
apply(power[[2]], c(1,2), confint)
apply(power[[3]], c(1,2), confint)


#################### POWER CURVES ####################

# add the effect sizes as a variable
power = lapply(power, as.data.frame)
power2 = lapply(power, tibble::rownames_to_column, var="Effect.Size")

# convert to long format
power.long = lapply(power2, pivot_longer, cols=all_of(models), names_to="Model", values_to="Power")

# combine all simulation scenarios together & add simulation scenario variable
data.all = rbind(power.long[[1]], power.long[[2]], power.long[[3]])
data.all$Scenario = rep(models, each=length(models)*length(effect_size))

# ensure correct ordering of the factor variables
data.all2$Scenario = factor(data.all2$Scenario, levels=models)
data.all2$Model = factor(data.all2$Model, levels=models)

# remove the percent sign from the effect sizes
data.all$Effect.Size = as.numeric(sapply(strsplit(data.all$Effect.Size, "%", fixed=T), head, 1))

# define confidence interval function
confint.data <- function(data, nsim){
  
  # calculate confidence interval
  interval = binom.confint(data$Power*nsim, nsim, method="exact")
  
  # add lower and upper bounds to data frame
  data$Lower = round(interval$lower, 3)
  data$Upper = round(interval$upper, 3)
  
  return(data)
}

# call the confidence interval function
data.all2 = confint.data(data.all, 1000)

# save power results
#write.table(data.all2, "power_results_long.txt", quote=F, sep="\t")

# read in power results
#data.all2 = read.table("power_results_long.txt", header=T, sep="\t")

# remove the inflated power to plot
plot.data = data.all2 %>% dplyr::filter(!(Scenario=="Mouse Intercept/Slope" & Model!="Mouse Intercept/Slope" & Effect.Size>0))

# plot of power faceted by simulation scenarios
ggplot() +
  geom_point(data=plot.data, aes(x=Effect.Size, y=Power, col=Model, shape=Model), size=1.7) + 
  geom_line(data=plot.data %>% filter(Effect.Size!=0), aes(x=Effect.Size, y=Power, col=Model), size=0.5) + 
  geom_ribbon(data=plot.data %>% filter(Effect.Size!=0), aes(x=Effect.Size, y=Power, ymin=Lower, ymax=Upper, fill=Model), alpha=0.25) +
  geom_abline(data=plot.data, aes(intercept=0.05, slope=0), lty=2, col="darkgrey") +
  labs(x="Effect Size (%)") + facet_wrap(~Scenario) +
  theme_pubr(border=T)
