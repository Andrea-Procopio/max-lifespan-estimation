## CODE DGP36 - GETTING THE DATASET AND THE VARIABLES

dgp36 = read_excel("attrition_dataset_path_here.xlsx")
ltl1 = dgp36$LTL1
ltl2 = dgp36$LTL2
months = dgp36$months
death = dgp36$`event (0=alive, 1=dead)`
sex = dgp36$`Sex (1=male, 2=female)`
age = dgp36$Age
attrition = ltl1-ltl2

library(ggplot2)
library(splines)
library(survival)
library(survminer)

ggplot(data = data.frame(attrition), aes(x = attrition)) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(title = "Density Plot of Attrition", x = "Attrition (LTL1 - LTL2)", y = "Density") +
  theme_linedraw()

ggplot(data = data.frame(attrition), aes(x = attrition)) +
  geom_histogram(bins = 30, fill = "blue", color = "black") +
  labs(title = "Histogram of Attrition", x = "Attrition", y = "Count")

ggplot(data = dgp36, aes(x = attrition, fill = factor(sex))) +
  geom_density(alpha = 0.5) +
  labs(title = "Density Plot of Attrition by Sex", x = "Attrition", y = "Density") +
  scale_fill_manual(values = c("blue", "pink"), labels = c("Male", "Female"))

#same density plot with adjusted bandwith
ggplot(data = data.frame(attrition), aes(x = attrition, fill = factor(sex))) +
  geom_density(alpha = 0.5, adjust = 0.3) + 
  labs(title = "Density Plot of Attrition with Adjusted Bandwidth", x = "Attrition", y = "Density")

#MALE ATTRITION, very small sample of 10 obs. 20% dead at followup 
mean(attrition[sex==1]) #0.53557
median(attrition[sex==1]) #0.5299

#FEMALE ATTRITION, small sample of 26 obs. ≈42% dead at followup
mean(attrition[sex==2]) #-0.03769615
median(attrition[sex==2]) #0.0735

table(sex) ## small sample, not conclusive
prop.table(table(death, sex), 2)


## KAPLAN MEIER FOR ATTRITION
summary(attrition)
attrition_quartiles= c(-1.534, -0.07032, 0.0918, 0.40877, 1.4889)
breaks_attirtion =  c(-Inf, 0.0918, Inf)

dgp36$Attrition_Quartiles <- cut(attrition, 
                            breaks = breaks_attirtion,
                            include.lowest = TRUE)

surv_object <- Surv(time = months/12, event = death)
km_fit_ltl <- survfit(surv_object ~ Attrition_Quartiles, data = dgp36)


library(ggthemes)
ggsurvplot(km_fit_ltl,
           data = dgp36,
           conf.int = TRUE,
           conf.int.alpha = 0.1,
           pval = TRUE,
           xlab = "Time (months)",
           ylab = "Survival probability",
           legend.title = "",
           ggtheme = theme_linedraw())

#coherent with the hypotesis that lower attrition enhance survival probability,
#but p = 0.66 -> not statistically significant
