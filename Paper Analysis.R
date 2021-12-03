#################### LOAD LIBRARIES ####################

library(lme4)
library(lmerTest)
library(nlme)
library(tidyverse)
library(ggpubr)
library(huxtable) #quick_docx
library(pbkrtest) #PBmodcomp
library(kableExtra) #knitr::kable
library(broom.mixed) #tidy


#################### READ DATA ####################

#setwd("C:/Repositories/Longitudinal.Mouse/")
data = read.csv("undernourished_study_science.csv", as.is=T)

# helper script with functions for plotting trend lines
source("C:/Repositories/EasyLME/Functions.R")


#################### FORMAT DATA ####################

# change missing or "dead" entries to NA
data = na_if(data, "nd")
data = na_if(data, "Dead")

# remove mice with NA values
#data2 = drop_na(data, X0:X32)

# Current setup of data has all time points for a mouse in one row. But we need
# the data to be in long format, where each time point is a separate observation (row).
data_long = gather(data, time, perc_weight, X0:At.sacrifice)

# remove time point of "At.sacrifice"
data_long = data_long[!data_long$time=="At.sacrifice",]

# remove the "X" in front of the time points and treat them as numeric
data_long$time = as.numeric(sapply(data_long$time, function(x) gsub("X", "", x)))

# change donor age to a binary variable (1 for 6 months, 0 for 18 months)
data_long$Donor.Age = ifelse(as.numeric(data_long$Donor.Age..months.) <= 11 , 1, 0)

# consider only observations with an efficiency above 50 percent (stated in paper)
data_long = data_long[data_long$Transplantation.Efficiency=="> 50%",]

# create a new mouse ID number to indicate the correct nested nature of mice within donors
data_long$Mouse = paste(data_long$Donor, data_long$Mouse.ID.Number, sep=".")

# treat percent weight change as numeric
data_long$perc_weight = as.numeric(data_long$perc_weight)

# select only necessary variables
data_long = data_long %>% select(Donor, Donor.Status, Time=time, Perc.Weight=perc_weight, Donor.Age, Mouse)

# relevel donor variable by decreasing Perc.Weight
max_donors = data_long %>% group_by(Donor) %>% summarize(max=max(Perc.Weight, na.rm=T)) %>% arrange(desc(max))
data_long$Donor = factor(data_long$Donor, levels=max_donors$Donor)

# remove % weight na values (necessary for PBmodcomp) & arrange by donor
data_long2 = data_long %>% drop_na(Perc.Weight) %>% arrange(Donor)

# treat character variables as factors
data_long2$Mouse = factor(data_long2$Mouse, levels=unique(data_long2$Mouse))
data_long2$Donor.Status = as.factor(data_long2$Donor.Status)

# remove mice missing more than one time point
#data_long3 = data_long2 %>% filter(Donor != "185A.2", Mouse != "4092C.9", Mouse != "3114C.11")

# remove the first few time points (get rid of the dip)
data_long_nodip = data_long2 %>% filter(Time>4)

# save processed data
#save(data_long, file="./EasyLME/data/data_long.R")

# create a data summary table (summary_table function in Functions.R script)
summary_table(data_long2, "Perc.Weight", "Time", "Donor.Status", "Donor", "Mouse") %>% 
  knitr::kable(format="html") %>% kable_styling(full_width=F) %>%
  add_header_above(c(" "=3, "Time"=2, "Perc.Weight"=2))


#################### EXPLORE DATA ####################

# histogram of response variable (% weight)
ggplot(data_long2, aes(x=Perc.Weight)) +
  geom_histogram() + theme_bw(base_size=13) + 
  labs(x="% Weight", title="Histogram of % Weight")

# scatterplot of overall data trends
ggplot(data_long2, aes(x=Time, y=Perc.Weight, color=Donor)) +
  geom_point() + theme_pubr(legend="right", border=T) + 
  labs(x="Time", y="Perc.Weight")

# trendlines per donor
ggplot(data_long2, aes(x=Time, y=Perc.Weight, colour=Donor.Status)) +
  geom_line(aes(group=Mouse), lwd=1) + facet_wrap(~Donor) + 
  theme_pubr(border=T)

# boxplots based on mouse
ggplot(data_long2, aes(x=Mouse, y=Perc.Weight, color=Donor)) +
  geom_boxplot() + theme_bw(base_size=13) + 
  labs(x="Mouse", y="% Weight", title="Growth per Mouse") +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())

# boxplots based on donor
ggplot(data_long2, aes(x=Donor, y=Perc.Weight, color=Donor.Status)) +
  geom_boxplot() + theme_bw(base_size=13) + 
  labs(x="Donor", y="% Weight", title="Growth per Donor")

# boxplots based on donor status
ggplot(data_long2, aes(x=Donor.Status, y=Perc.Weight)) +
  geom_boxplot() + theme_bw(base_size=13) + 
  labs(x="Donor Status", y="% Weight", title="Growth per Status")

# average trendlines per donor
avg_donors = data_long2 %>% group_by(Donor, Time) %>% 
  summarize(Perc.Weight=mean(Perc.Weight), Donor.Status=first(Donor.Status))

ggplot() +
  geom_line(avg_donors, mapping=aes(x=Time, y=Perc.Weight, color=Donor, linetype=Donor.Status), lwd=0.75) +
  theme_pubr(legend="right", border=T) 


##################### FIT MODELS ####################

# fit complete mixed effects model (donor slope/intercept + mouse slope/intercept)
complete = lmerTest::lmer(Perc.Weight ~ Donor.Status*Time + (Time|Donor) + 
                        (Time|Donor:Mouse), data_long2)

# check convergence warnings across optimizers (default optimizer - nloptwrap)
summary(allFit(complete))

# the following optimizers produced convergence warnings:
#     bobyqa: boundary (singular) fit: see ?isSingular
#     Nelder_Mead: unable to evaluate scaled gradient
#                  Model failed to converge: degenerate  Hessian with 1 negative eigenvalues
#     optimx.L-BFGS-B: boundary (singular) fit: see ?isSingular
# - nlminbwrap, nmkbw, nloptwrap.NLOPT_LN_NELDERMEAD, and nloptwrap.NLOPT_LN_BOBYQA 
#   did not produce convergence warnings

# try making random the effect terms uncorrelated
complete2 = lmerTest::lmer(Perc.Weight ~ Donor.Status*Time + (Time||Donor) + 
                            (Time||Donor:Mouse), data_long2)

# check the convergence warnings again
summary(allFit(complete2))

# only nloptwrap.NLOPT_LN_BOBYQA (the default) produced a warning:
#     Model failed to converge with max|grad| = 0.00708574 (tol = 0.002, component 1)

# fit nested mixed effects models
lmod2 = update(complete, .~. - (Time|Donor) + (1|Donor)) #mouse slope/intercept + donor intercept
lmod3 = update(complete, .~. - (Time|Donor:Mouse) + (1|Donor:Mouse)) #donor slope/intercept + mouse intercept
lmod4 = update(lmod2, .~. - (1|Donor)) #mouse slope/intercept
lmod5 = update(lmod3, .~. - (Time|Donor) + (1|Donor)) #mouse intercept + donor intercept
lmod6 = update(lmod5, .~. - (1|Donor:Mouse)) #donor intercept
lmod7 = update(lmod5, .~. - (1|Donor)) #mouse intercept
lmod8 = update(lmod3, .~. - (1|Donor:Mouse)) #donor slope/intercept

# lmod3 produced the following warning: 
#     In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#     Model failed to converge with max|grad| = 0.00295291 (tol = 0.002, component 1)
# - nmkbw and optimx.L-BFGS-B were the only optimizers that did not produce warnings for lmod3

# lmod8 produced the following warning: 
#     optimx.L-BFGS-B: boundary (singular) fit: see ?isSingular

# fit null model (no random effects)
null = lm(Perc.Weight ~ Donor.Status*Time, data_long2)

# fit repeated measures anova
rm = aov(formula = Perc.Weight ~ Donor.Status*Time + Error(Mouse), data_long2)

# fit autocorrelation model
ac = gls(Perc.Weight ~ Donor.Status*Time, data=data_long2, 
         correlation=corCAR1(form=~1|Mouse)) # logLik = -1175.597

# fit autocorrelation model with a random slope for mouse
re.ac = lme(Perc.Weight~Donor.Status*Time, data=data_long2, random=~0+Time|Mouse, 
             correlation=corCAR1(form=~1|Mouse)) # logLik = -1133.994

# add age as a covariate
re.ac2 = update(re.ac, .~. + Donor.Age) # logLik = -1130.725

# save models for power analysis
#save(lmod4, file="mouse_slope_int.R")
#save(lmod7, file="mouse_intercept.R")
#save(null, file="no_RE.R")

# refit the models with the first few time points removed (no convergence warnings)
complete.nodip = update(complete, data=data_long_nodip)
lmod2.nodip = update(lmod2, data=data_long_nodip)
lmod4.nodip = update(lmod4, data=data_long_nodip)
lmod7.nodip = update(lmod7, data=data_long_nodip) # interaction term much less significant, sign for group effect changed
null.nodip = update(null, data=data_long_nodip) # interaction term much less significant


#################### RESIDUAL PLOTS ####################

# define model names
names = c("Donor Intercept/Slope + Mouse Intercept/Slope",
          "Donor Intercept + Mouse Intercept/Slope",
          "Mouse Intercept/Slope", "Mouse Intercept", "No Random Effects")

# create list of models
models = list(complete, lmod2, lmod4, lmod7, null)

# define variables for data frame
Time = model.frame(complete)[["Time"]]
Model = rep(names, each=length(Time))
RE = model.frame(complete)[["Mouse"]]

# extract the residuals for each model
Residuals = lapply(models, FUN=residuals)
Residuals = unlist(Residuals)

# create a data frame with information for the combined residual plots
resid_data = data.frame(Time, RE, Residuals, Model)
resid_data$Model = factor(resid_data$Model, levels=rev(names))
                        
# plot time vs residuals faceted by model
ggplot(resid_data, aes(x=Time, y=Residuals)) +
  geom_line(aes(group=RE)) + facet_wrap(~Model, ncol=2) + 
  theme_pubr(border=T)
# No Random Effects - variance increases over time suggesting need for random effects
# Mouse Intercept - variance still not constant across time - add random slope


#################### COEFFICIENT PLOT ####################

# create data frame of the model coefficient estimates, standard errors, and confidence intervals
# (coef_data function in Functions.R script)
terms = c("Group", "Time", "Interaction")
model.data = coef_data(models, names) %>% arrange(desc(Model))

# reorder the factor variables
model.data$Term = reorder(model.data$Term, model.data$Estimate)
model.data$Model = reorder(model.data$Model, desc(model.data$Model))

# dot and whisker plot of the coefficient estimates and 95% confidence intervals
ggplot(model.data) + geom_point(aes(x=Estimate, y=Term, color=Model), position=position_dodge(width = 1/2)) +
  geom_linerange(aes(y=Term, xmin=ci.low, xmax=ci.high, color=Model), position=position_dodge(width = 1/2), show.legend=F, lwd=1) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
  scale_color_manual(name="Model", values=gg_color_hue(5), breaks=names) + 
  labs(x="Estimate", y=element_blank()) +
  theme_pubr(border=T, legend="right")


#################### TEST MODELS ####################

# likelihood ratio tests (lrt) of the random effects
# (use refit=FALSE for restricted likelihood ratio test)
# extract chi-squared test statistic (column 6) & p-value (column 8)
lrt1v2 = anova(complete, lmod2)[2,c(6,8)] #donor slope
lrt1v3 = anova(complete, lmod3)[2,c(6,8)] #mouse slope
lrt2v4 = anova(lmod2, lmod4)[2,c(6,8)] #donor intercept
lrt2v5 = anova(lmod2, lmod5)[2,c(6,8)] #mouse slope
lrt3v5 = anova(lmod3, lmod5)[2,c(6,8)] #donor slope
lrt3v8 = anova(lmod3, lmod8)[2,c(6,8)] #mouse intercept
lrt4v7 = anova(lmod4, lmod7)[2,c(6,8)] #mouse slope
lrt5v6 = anova(lmod5, lmod6)[2,c(6,8)] #mouse intercept
lrt5v7 = anova(lmod5, lmod7)[2,c(6,8)] #donor intercept
lrt8v6 = anova(lmod8, lmod6)[2,c(6,8)] #donor slope
lrt6vnull = anova(lmod6, null)[2,c(6,8)] #donor intercept
lrt7vnull = anova(lmod7, null)[2,c(6,8)] #mouse intercept

# parametric bootstrap tests of the random effects (these take a little while to run)
#pb1v2 = pbkrtest::PBmodcomp(complete, lmod2, seed=1)
#pb2v4 = pbkrtest::PBmodcomp(lmod2, lmod4, seed=1)
#pb4v7 = pbkrtest::PBmodcomp(lmod4, lmod7, seed=1)
#pb7v8 = pbkrtest::PBmodcomp(lmod7, null, seed=1) #doesn't work with lm model


#################### RESULTS ####################

# extract the lrt p-values
pvalues = as.numeric(c(lrt1v2[,2], lrt2v4[,2], lrt4v7[,2], lrt7vnull[,2], ""))
pvalues = formatC(pvalues, format="e", digits=2)
pvalues[length(pvalues)] = ""

# make a model comparison table for the five models in paper 
# (results_table function in Functions.R script)
results = results_table(models, pvalues, names)

# bold the estimates with a p-value < 0.05
results2 = results[[1]] %>% knitr::kable(format="html", format.args=list(big.mark = ','), escape=F) %>%
  column_spec(2, bold=ifelse(results[[2]][,2]>0.05, FALSE, TRUE)) %>% 
  column_spec(3, bold=ifelse(results[[2]][,3]>0.05, FALSE, TRUE)) %>%
  column_spec(4, bold=ifelse(results[[2]][,4]>0.05, FALSE, TRUE))
#quick_docx(results2)

# make a test results table for the five models in paper
tests = rbind(lrt1v2, lrt2v4, lrt4v7, lrt7vnull)
models = c("1 vs 2", "2 vs 3", "3 vs 4", "4 vs 5")
tests = cbind(models, tests)
colnames(tests) = c("Models", "Test Statistic", "P-value")
tests = theme_plain(hux(tests, add_colnames = TRUE))
#quick_docx(tests)


#################### FITTED LINES ####################

# donor_lines and mouse_lines functions from Functions.R script

# donor intercept/slope + mouse intercept/slope model (complete)
donor_lines(complete, data_long2, "Perc.Weight", "Time", "Donor.Status", "Donor")
mouse_lines(complete, data_long2, "Perc.Weight", "Time", "Donor.Status", "Donor", "Mouse")

# donor intercept + mouse intercept/slope model (lmod2)
donor_lines(lmod2, data_long2, "Perc.Weight", "Time", "Donor.Status", "Donor")
mouse_lines(lmod2, data_long2, "Perc.Weight", "Time", "Donor.Status", "Donor", "Mouse")

# mouse intercept/slope model (lmod4)
donor_lines(lmod4, data_long2, "Perc.Weight", "Time", "Donor.Status", "Donor")
mouse_lines(lmod4, data_long2, "Perc.Weight", "Time", "Donor.Status", "Donor", "Mouse")

# mouse intercept model (lmod7)
donor_lines(lmod7, data_long2, "Perc.Weight", "Time", "Donor.Status", "Donor")
mouse_lines(lmod7, data_long2, "Perc.Weight", "Time", "Donor.Status", "Donor", "Mouse")

# no random effects model (null)
donor_lines(null, data_long2, "Perc.Weight", "Time", "Donor.Status", "Donor")
