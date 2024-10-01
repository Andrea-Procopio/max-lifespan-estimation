## CODE DATAGP - GETTING THE DATASET AND THE VARIABLES

datagp = read_excel("data_path_here.xlsx")
ltl = datagp$TL
months = datagp$months
death = datagp$`event (0=alive, 1=dead)`
sex = datagp$`Sex (1=male, 2=female)`
age = datagp$Age
age_death = datagp$`age at death`

library(ggplot2)
library(splines)
library(survival)
library(survminer)

## DENSITY PLOT OF T/S BY AGE QUANTILES 
quantile_breaks <- quantile(age, probs = c(0, 0.25, 0.5, 0.75, 1),na.rm = TRUE)
`Age Group` <- cut(age, breaks = quantile_breaks,
                   include.lowest = TRUE,
                   right = FALSE,
                   labels = c("Q1 [65-70]", "Q2 [71-78]",
                              "Q3 [79-89]", "Q4 [90-106]"))

ggplot(datagp, aes(x = ltl, fill = `Age Group`)) + 
  geom_density(alpha = 0.5) + 
  labs(x = "T/S (LTL)", y = "Density") +
  theme_minimal()

#SAME PLOT IN LOG SCALE FOR LTL
ggplot(datagp, aes(x = log(ltl), fill = `Age Group`)) + 
  geom_density(alpha = 0.5) + 
  labs(x = "log(T/S)", y = "Density") +
  theme_minimal()


## DENSITY PLOT LTL BY SEX
ggplot(data = datagp, aes(x = log(ltl), fill = factor(sex))) +
  geom_density(alpha = 0.5) +
  labs(title = "Density Plot of log(ltl) by Sex", x = "Attrition", y = "Density") +
  scale_fill_manual(values = c("blue", "pink"), labels = c("Male", "Female"))


## SCATTER PLOT OF T/S AND AGE

##LOESS REGRESSION
ggplot(datagp, aes(x = Age, y = ltl)) + 
  geom_point() +
  geom_smooth(method = "loess", span = 1) + #more stable, less overfitting
  labs(x = "Age", y = "T/S (LTL)") + theme_minimal()

##SPLINE REGRESSION
model <- lm(ltl ~ ns(Age, df = 3), data = datagp)
age_seq <- seq(min(datagp$Age), max(datagp$Age), length.out = 300)
pred_data <- data.frame(Age = age_seq)
pred_data$ltl_pred <- predict(model, newdata = pred_data)

ggplot(datagp, aes(x = Age, y = ltl)) +
  geom_point() +
  geom_line(data = pred_data, aes(x = Age, y = ltl_pred), col = "blue") +
  labs(x = "Age", y = "T/S (LTL)") + theme_minimal()


## KAPLAN MEIER LTL AND AGE
library(wesanderson)
names(wes_palettes)
palette1 <- wes_palette("GrandBudapest2", 4, type = "discrete")
palette2 <- wes_palette("FantasticFox1", 4, type = "discrete")

surv_object <- Surv(time = months, event = death)
km_fit <- survfit(surv_object ~ 1)

# Plotting the overall survival curve
ggsurvplot(km_fit,
           data = datagp,
           conf.int = TRUE,
           conf.int.alpha = 0.5,
           palette = "blue",
           xlab = "Time (months)",
           ylab = "Survival probability",
           title = "Overrall survival curve",
           legend = "none")

##STRATIFIED BY AGE QUARTILES
summary(age)
breaks_age <- c(65, 71, 79, 90, 107)
labels_age <- c("Q1 [65-70]", "Q2 [71-78]", "Q3 [79-89]", "Q4 [90-106]")

datagp$AgeGroup <- cut(age,
                   breaks = breaks_age,
                   include.lowest = TRUE,
                   right = FALSE,
                   labels = labels_age)

km_fit_age <- survfit(surv_object ~ AgeGroup, data = datagp)

ggsurvplot(km_fit_age,
           data = datagp,
           conf.int = TRUE,
           conf.int.alpha = 0.1,
           pval = TRUE,
           xlab = "Time (months)",
           ylab = "Survival probability",
           legend.title = "")


## STRATIFIED BY LTL QUARTILES
summary(ltl)
breaks_ltl =  c(0.156, 0.4928, 0.7269, 1.1875, 5.4958)

datagp$LTL_Quartiles <- cut(ltl, 
                            breaks = breaks_ltl,
                            include.lowest = TRUE)

surv_object <- Surv(time = months, event = death)
km_fit_ltl <- survfit(surv_object ~ LTL_Quartiles, data = datagp)

library(ggthemes)
ggsurvplot(km_fit_ltl,
           data = datagp,
           conf.int = TRUE,
           conf.int.alpha = 0.1,
           pval = TRUE,
           xlab = "Time (months)",
           ylab = "Survival probability",
           legend.title = "",
           ggtheme = theme_linedraw())


## INTERACTIVE PLOT ON KAPLAN MEIER STRATIFIED FOR QUARTILES OF LTL
library(plotly)

ggsurvplotly <- ggsurvplot(km_fit_ltl, data = datagp,
                           conf.int = TRUE, pval = TRUE, risk.table = TRUE)

ggplotly(ggsurvplotly$plot)


