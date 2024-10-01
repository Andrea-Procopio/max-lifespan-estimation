#### R SCRIPT FOR STATS PROJECT
library(ggplot2)
library(maxLik)
library(splines)
library(survival)
library(survminer)
library(readxl)
library(ggthemes)
library(longevity)

#1. IDL DATASET

idl_path = "IDL_dataset_fil_path.csv"
idl = read.csv(idl_path, sep = ";")

# Kaplan Meier estimation of the Hazard function
h <- c()
for (i in 105:122) {
  d_t <- sum(idl$AGEYEARS == i)
  r_t <- sum(idl$AGEYEARS >= i)
  h <- c(h, d_t / r_t)
}

# plot the hazard function
hazard <- data.frame(Age = 105:122, Hazard = h)
ggplot(hazard, aes(x = Age, y = Hazard)) +
  geom_point(color = "blue", size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(x = "Age",
       y = "Hazard rate") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

### Fit an exponential and Generalized Pareto (with interval truncation and interval censoring)

# look at IDL data between calendar times c1, c2
c1 = 1995; c2 = 2010
idl_s = idl[(idl$DEATH_YEAR < c2) & (idl$DEATH_YEAR >= c1),]

# define the variables
t = idl_s$AGEYEARS - 105
x = idl_s$BIRTH_YEAR + 105
C1 = pmax(0, c1 - x)
C2 = c2 - x

# likelihood functions
llk_exp_trun = function(par) {
  b = par[1]
  llk_value = log((pexp(t + 1, rate=b) - pexp(t, rate=b)) / (pexp(C2, rate=b) - pexp(C1, rate=b)))
  sum(llk_value)
}
llk_gp_trun = function(par) {
  g = par[1]
  s = par[2]
  llk = log((pgpd(t + 1, shape=g, scale=s) - pgpd(t, shape=g, scale=s)) / (pgpd(C2, shape=g, scale=s) - pgpd(C1, shape=g, scale=s)))
  sum(llk)
}

# ml estimation
ml_exp = maxLik(llk_exp_trun, start=c(1)) #EXP
ml_gp = maxLik(llk_gp_trun, start=c(1, 1)) #GP


# plot
barplot(prop.table(table(idl_s$AGEYEARS)), ylim = 0:1)
lines((pexp(1:17, ml_exp$estimate) - pexp(0:16, ml_exp$estimate)), type="l", col="blue")
lines((pgpd(1:17, shape=ml_gp$estimate[1], scale=ml_gp$estimate[2]) - pgpd(0:16, shape=ml_gp$estimate[1], scale=ml_gp$estimate[2])), type="l", col="red")


# this analysis allows us to formulate the null hypothesis that the exponential and GP fit equally well
# in order to test such hypothesis, we'll use the likelihood ratio test

# likelihood functions (no log)
lik_exp_trun = function(par) {
  b = par[1]
  lik = (pexp(t + 1, rate=b) - pexp(t, rate=b)) / (pexp(C2, rate=b) - pexp(C1, rate=b))
  sum(lik)
}
lik_gp_trun = function(par) {
  g = par[1]
  s = par[2]
  lik = (pgpd(t + 1, shape=g, scale=s) - pgpd(t, shape=g, scale=s)) / (pgpd(C2, shape=g, scale=s) - pgpd(C1, shape=g, scale=s))
  sum(lik)
}

# likelihood ratio
ratio = lik_gp_trun(ml_gp$estimate) / lik_exp_trun(c(ml_exp$estimate))

# the lambda which follows the chi squared distribution
lambda = -2 * log(ratio)

# the p-value: probability that such ratio is higher than the one observed (so that GP is even better relatively to exp)
p_value <- pchisq(lambda, 1, lower.tail = FALSE)

# since the p-value is higher (much higher) than 0.05 --> I cannot reject the NULL HYPOTHESIS (the difference in fit between the two models is not statistically significant) 
# thus I conclude that the models fit equally well




##########################2. LTL DATASET
path_data = "dataset_fil_path.xlsx"
data = read_excel(path_data)
path_followup = "dataset_fil_path_followup.xlsx"
followup = read_excel(path_followup)
#2. LTL DATASET
# Here, we report the code for the analysis of the dataset about LTL.
# There are two datasets: 
# "data", which includes over 500 observations, with a single measure of LTL and the survival time (in months)
# "followup", which includes 36 observations, with two measures of LTL for each subject and the survival time after the second LTL measure (in months)
              # Such measures were carried out at a mean distance of 7 years

# AVERAGE LTL for dead people VS alive people (full dataset)

# average ltl for alive people
ltl_alive = data$TL[data$`event (0=alive, 1=dead)`==0]
ltl_alive = na.omit(ltl_alive)
avg_ltl_alive = mean(ltl_alive)

# average ltl for dead people
ltl_dead = data$TL[data$`event (0=alive, 1=dead)`==1]
ltl_dead = na.omit(ltl_dead)
avg_ltl_dead = mean(ltl_dead)


# Z-TEST TO ADDRESS STATISTICAL SIGNIFICANCE

# compute sd and sample sizes
sd_ltl_dead = sd(ltl_dead)
sd_ltl_alive = sd(ltl_alive)
n_dead = length(ltl_dead)
n_alive = length(ltl_alive)

# calculate the z-statistic
z = ((avg_ltl_alive-avg_ltl_dead) / sqrt((sd_ltl_alive^2 / n_alive) + (sd_ltl_dead^2 / n_dead)))

# compute p-value
p_value = 2 * (1 - pnorm(abs(z))) # since the p-value is much smaller than 0.05, the difference is statistically significant


# STANDARD DEVIATION OF LTL BY QUARTILE OF AGE

# sort data by age
data_sorted = data[order(data$Age), ]
# split by quartile of Age
data_sorted$Quartile = cut(data_sorted$Age, breaks = quantile(data_sorted$Age, probs = seq(0, 1, 0.25)), include.lowest = TRUE, labels = FALSE)
#mean LTL for each quartile
mean_ltl_quartile = tapply(data_sorted$TL, data_sorted$Quartile, mean)
# sd of LTL by quartile of age
std_dev_ltl_quartile = tapply(data_sorted$TL, data_sorted$Quartile, sd)


# MEAN ATTRITION FOR PEOPLE WHO DIED VS PEOPLE WHO SURVIVED (follow-up dataset)

# mean attrition of alive people = -0.06
att_fu_alive = followup$LTL2[followup$`event (0=alive, 1=dead)`==0] - followup$LTL1[followup$`event (0=alive, 1=dead)`==0]
mean(att_fu_alive)

# mean attrition of dead people = -0.23
att_fu_dead = followup$LTL2[followup$`event (0=alive, 1=dead)`==1] - followup$LTL1[followup$`event (0=alive, 1=dead)`==1]
mean(att_fu_dead)


#### PLOTS
# define variables
ltl = data$TL
months = data$months
death = data$`event (0=alive, 1=dead)`
sex = data$`Sex (1=male, 2=female)`
age = data$Age
age_death = data$`age at death`


## DENSITY PLOT OF T/S BY AGE QUANTILES 
quantile_breaks = quantile(age, probs = c(0, 0.25, 0.5, 0.75, 1),na.rm = TRUE)
`Age Group` = cut(age, breaks = quantile_breaks,
                   include.lowest = TRUE,
                   right = FALSE,
                   labels = c("Q1 [65-70]", "Q2 [71-78]",
                              "Q3 [79-89]", "Q4 [90-106]"))

ggplot(data, aes(x = ltl, fill = `Age Group`)) + 
  geom_density(alpha = 0.5) + 
  labs(x = "LTL (T/S)", y = "Density") +
  theme_linedraw()

#SAME PLOT IN LOG SCALE FOR LTL
ggplot(data, aes(x = log(ltl), fill = `Age Group`)) + 
  geom_density(alpha = 0.5) + 
  labs(x = "log(T/S)", y = "Density") +
  theme_linedraw()

## SCATTER PLOT OF T/S AND AGE

##LOESS REGRESSION
ggplot(data, aes(x = Age, y = ltl)) + 
  geom_point() +
  geom_smooth(method = "loess", span = 1) + #more stable, less overfitting
  labs(x = "Age", y = "LTL (T/S)") + theme_linedraw()


## KAPLAN MEIER LTL AND AGE
surv_object = Surv(time = months, event = death)
km_fit = survfit(surv_object ~ 1)

# Plotting the overall survival curve
ggsurvplot(km_fit,
           data = data,
           conf.int = TRUE,
           conf.int.alpha = 0.5,
           palette = "blue",
           xlab = "Time (months)",
           ylab = "Survival probability",
           title = "Overrall survival curve",
           legend = "none")


## KM ESTIMATE BY LTL QUARTILES
summary(ltl)
breaks_ltl =  c(0.156, 0.4928, 0.7269, 1.1875, 5.4958)

data$LTL = cut(ltl, breaks = breaks_ltl,include.lowest = TRUE)

surv_object = Surv(time = months, event = death)
km_fit_ltl = survfit(surv_object ~ LTL, data = data)


ggsurvplot(km_fit_ltl,
           data = data,
           conf.int = TRUE,
           conf.int.alpha = 0.1,
           pval = TRUE,
           xlab = "Time (months)",
           ylab = "Survival probability",
           legend.title = "",
           ggtheme = theme_linedraw())
