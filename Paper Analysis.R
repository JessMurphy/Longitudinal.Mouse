########## LOAD LIBRARIES ##########

library(lme4)
library(lmerTest)
library(ggplot2)
library(gridExtra) #grid.arrange
library(tidyr) #gather
library(huxtable) #theme_plain, hux
library(jtools) #export_summs
library(pbkrtest) #PBmodcomp
library(dplyr) #na_if


########## READ DATA ##########

#setwd("C:/Repositories/Longitudinal.Mouse/")
data = read.csv("undernourished_study_science.csv", as.is=T)

# helper script with functions for plotting trend lines
source("./EasyLME/Functions.R")


########## FORMAT DATA ##########

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

# treat character variables as factors
data_long$Donor.Status = as.factor(data_long$Donor.Status)
data_long$Donor = as.factor(data_long$Donor)
data_long$Mouse = as.factor(data_long$Mouse)

# select only necessary variables
data_long = data_long %>% select(Donor, Donor.Status, Time=time, Perc.Weight=perc_weight, Donor.Age, Mouse)

# remove % weight na values (necessary for PBmodcomp) &
# reorder Donor & Mouse variables by decreasing Perc.Weight
data_long2 = data_long %>% drop_na(Perc.Weight) %>%
  mutate(across(Donor, ~reorder(factor(.), Perc.Weight, FUN=max))) %>%
  mutate(across(Mouse, ~factor(., levels=unique(Mouse[order(Donor)]))))

# remove mice missing more than one time point
#data_long3 = data_long2 %>% filter(Donor != "185A.2", Mouse != "4092C.9", Mouse != "3114C.11")

# save processed data
#save(data_long, file="./EasyLME/data/data_long.R")


########## EXPLORE DATA ##########

# histogram of response variable (% weight)
ggplot(data_long2, aes(x=Perc.Weight)) +
  geom_histogram() + theme_bw(base_size=13) + 
  labs(x="% Weight", title="Histogram of % Weight")

# scatterplot of overall data trends
ggplot(data_long2, aes(x=Time, y=Perc.Weight, color=Donor)) +
  geom_point() + theme_bw(base_size=13) + 
  labs(x="Days", y="% Weight", title="Overall Growth Trends")

# trendlines per donor
ggplot(data_long2, aes(x=Time, y=Perc.Weight, colour=Donor.Status)) +
  geom_point() + geom_line(aes(group=Mouse)) + facet_wrap(~Donor) + 
  theme_bw(base_size=13)

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


########### FIT MODELS ##########

# fit complete mixed effects model (donor slope/intercept + mouse slope/intercept)
complete = lmerTest::lmer(Perc.Weight ~ Donor.Status*Time + (Time|Donor) + 
                        (Time|Donor:Mouse), data_long2)

# check convergence warnings across optimizers (default optimizer - nloptwrap)
summary(allFit(complete))

# the following optimizers produced convergence warnings:
# bobyqa: boundary (singular) fit: see ?isSingular
# Nelder_Mead: unable to evaluate scaled gradient
#              Model failed to converge: degenerate  Hessian with 1 negative eigenvalues
# optimx.L-BFGS-B: boundary (singular) fit: see ?isSingular

# nlminbwrap, nmkbw, nloptwrap.NLOPT_LN_NELDERMEAD, and nloptwrap.NLOPT_LN_BOBYQA 
# did not produce convergence warnings

# fit nested mixed effects models
lmod2 = update(complete, .~. - (Time|Donor) + (1|Donor)) #mouse slope/intercept + donor intercept
lmod3 = update(complete, .~. - (Time|Donor:Mouse) + (1|Donor:Mouse)) #donor slope/intercept + mouse intercept
lmod4 = update(lmod2, .~. - (1|Donor)) #mouse slope/intercept
lmod5 = update(lmod3, .~. - (Time|Donor) + (1|Donor)) #mouse intercept + donor intercept
lmod6 = update(lmod5, .~. - (1|Donor:Mouse)) #donor intercept
lmod7 = update(lmod5, .~. - (1|Donor)) #mouse intercept
lmod8 = update(lmod3, .~. - (1|Donor:Mouse)) #donor slope/intercept

# lmod3 produced the following warning: boundary (singular) fit: see ?isSingular
# - nlminbwrap and nmkbw were the only optimizers that did not produce warnings for lmod3

# fit null model (no random effects)
null = lm(Perc.Weight ~ Donor.Status*Time, data_long2)

# save models for power analysis
#save(lmod4, file="mouse_slope_int.R")
#save(lmod7, file="mouse_intercept.R")
#save(null, file="no_RE.R")


########## DIAGNOSTIC PLOTS ##########

# time vs residuals (null model) to determine need for random effects
ggplot(data_long2, aes(x=Time, y=residuals(null), group=Mouse)) +
  geom_line() + facet_wrap(~Donor) +
  labs(x="Days", y="Residuals", title="No Random Effects") + theme_bw()
# variance increases over time suggesting need for random effects

# time vs residuals (mouse intercept model) to determine need for random slope
ggplot(data_long2, aes(x=Time, y=residuals(lmod7), group=Mouse)) +
  geom_line() + facet_wrap(~Donor) +
  scale_y_continuous(limits=c(-20,20)) + theme_bw() +
  labs(x="Days", y="Residuals", title="Mouse Intercept") 
# variance still not constant across time - add random slope

# time vs residuals (mouse intercept & slope model)
ggplot(data_long2, aes(x=Time, y=residuals(lmod4), group=Mouse)) +
  geom_line() + facet_wrap(~Donor) +
  scale_y_continuous(limits=c(-20,20)) + theme_bw() +
  labs(x="Days", y="Residuals", title="Mouse Intercept/Slope")

# time vs residuals (donor intercept + mouse intercept & slope model)
ggplot(data_long2, aes(x=Time, y=residuals(lmod2), group=Mouse)) +
  geom_line() + facet_wrap(~Donor) +
  scale_y_continuous(limits=c(-20,20)) + theme_bw() +
  labs(x="Days", y="Residuals", title="Donor Intercept + Mouse Int/Slope")

# time vs residuals (donor intercept & slope + mouse intercept & slope model)
ggplot(data_long2, aes(x=Time, y=residuals(complete), group=Mouse)) +
  geom_line() + facet_wrap(~Donor) +
  scale_y_continuous(limits=c(-20,20)) + theme_bw() +
  labs(x="Days", y="Residuals", title="Donor Int/Slope + Mouse Int/Slope")


########## TEST MODELS ##########

# likelihood ratio tests (lrt)
# (use refit=FALSE for restricted liklihood ratio test)
# extract chi-squared test statistic (column 6) & p-value (column 8)
lrt1v2 = anova(complete, lmod2)[2,c(6,8)]
lrt1v3 = anova(complete, lmod3)[2,c(6,8)]
lrt2v4 = anova(lmod2, lmod4)[2,c(6,8)]
lrt2v5 = anova(lmod2, lmod5)[2,c(6,8)]
lrt3v5 = anova(lmod3, lmod5)[2,c(6,8)]
lrt4v7 = anova(lmod4, lmod7)[2,c(6,8)]
lrt5v6 = anova(lmod5, lmod6)[2,c(6,8)]
lrt5v7 = anova(lmod5, lmod7)[2,c(6,8)]
lrt7v8 = anova(lmod7, null)[2,c(6,8)]

# parametric bootstrap tests (these take a little while to run)
#pb1v2 = pbkrtest::PBmodcomp(complete, lmod2, seed=1)
#pb2v4 = pbkrtest::PBmodcomp(lmod2, lmod4, seed=1)
#pb4v7 = pbkrtest::PBmodcomp(lmod4, lmod7, seed=1)
#pb7v8 = pbkrtest::PBmodcomp(lmod7, null, seed=1) #doesn't work with lm model


########## RESULTS ##########

# make a model comparison table for all models
results = export_summs(complete, lmod2, lmod3, lmod4, lmod5, lmod6, lmod7, null,
                       statistics = c("logLik"), error_pos = 'same', align = "center",
                       bold_signif = 0.05, stars = NULL, number_format = 2,
                       model.names = c("Donor Intercept/Slope + Mouse Intercept/Slope",
                                       "Donor Intercept + Mouse Intercept/Slope",
                                       "Donor Intercept/Slope + Mouse Intercept",
                                       "Mouse Intercept/Slope", "Donor Intercept + Mouse Intercept",
                                       "Donor Intercept", "Mouse Intercept",
                                       "No Random Effects"),
                       coefs = c("Group" = "Donor.StatusUndernourished",
                                 "Time" = "time",
                                 "Interaction" = "Donor.StatusUndernourished:time"))
#quick_docx(theme_plain(t(results)))

# make a model comparison table for just the five models in paper
results2 = export_summs(complete, lmod2, lmod4, lmod7, null,
             statistics = c("logLik"), error_pos = 'same', align = "center",
             bold_signif = 0.05, stars = NULL, number_format = 2,
             model.names = c("Donor Intercept/Slope + Mouse Intercept/Slope",
                             "Donor Intercept + Mouse Intercept/Slope",
                             "Mouse Intercept/Slope", "Mouse Intercept",
                             "No Random Effects"),
             coefs = c("Group" = "Donor.StatusUndernourished",
                       "Time" = "time",
                       "Interaction" = "Donor.StatusUndernourished:time"))
#quick_docx(theme_plain(t(results2)))

# make a test results table for the five models in paper
tests = rbind(lrt1v2, lrt2v4, lrt4v7, lrt7v8)
models = c("1 vs 2", "2 vs 3", "3 vs 4", "4 vs 5")
tests = cbind(models, tests)
colnames(tests) = c("Models", "Test Statistic", "P-value")
tests = theme_plain(hux(tests, add_colnames = TRUE))
#quick_docx(tests)


########## FITTED LINES ##########

# donor_lines and mouse_lines from Functions.R script

# donor intercept/slope + mouse intercept/slope model (complete)
donor_lines(complete, data_long2, "perc_weight", "time", "Donor.Status", "Donor")
mouse_lines(complete, data_long2, "perc_weight", "time", "Donor.Status", "Donor", "Mouse")

# donor intercept + mouse intercept/slope model (lmod2)
donor_lines(lmod2, data_long2, "perc_weight", "time", "Donor.Status", "Donor")
mouse_lines(lmod2, data_long2, "perc_weight", "time", "Donor.Status", "Donor", "Mouse")

# mouse intercept/slope model (lmod4)
donor_lines(lmod4, data_long2, "perc_weight", "time", "Donor.Status", "Donor")
mouse_lines(lmod4, data_long2, "perc_weight", "time", "Donor.Status", "Donor", "Mouse")

# mouse intercept model (lmod7)
donor_lines(lmod7, data_long2, "perc_weight", "time", "Donor.Status", "Donor")
mouse_lines(lmod7, data_long2, "perc_weight", "time", "Donor.Status", "Donor", "Mouse")

# no random effects model (null)
donor_lines(null, data_long2, "perc_weight", "time", "Donor.Status", "Donor")
