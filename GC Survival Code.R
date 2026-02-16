library(readr)
library(tidyverse)
library(survival)
library(asaur)

# Loading Data

PBC276 <- read_csv("~/UTH Survival Analysis/PH1831_Project/PH1831 Project Data/PBC276.csv")
names(PBC276)

attach(PBC276)

null.model.pbc <- coxph(Surv(futime, status) ~ 1 , data=PBC276)
null.model.pbc

### Step 2: Find the residuals for the Null Model

rr <- residuals(null.model.pbc)

# Numeric Variables

# need to look at bili copper albumin protime age

# Variable = age
plot(rr ~ age, xlab="Age at Time of Diagnosis", ylab="Residual")
lines(lowess(rr ~ age), col="blue")
# try log transformation and spline - no idea what might work

plot(rr ~ bili, xlab="Bili", ylab="Residual")
lines(lowess(rr ~ bili), col="blue")
# looks like log transformation would work

plot(rr ~ copper, xlab="Copper", ylab="Residual")
lines(lowess(rr ~ copper), col="blue")
# looks like log transformation would work

plot(rr ~ albumin, xlab="Albumin", ylab="Residual")
lines(lowess(rr ~ albumin), col="blue")
# no idea what would work here

plot(rr ~ protime, xlab="Protime", ylab="Residual")
lines(lowess(rr ~ protime), col="blue")
# looks like log transformation would work

# log hazard transformation as discussed in class
plot(rr ~ age, lty=2, log="x", xlab="Age", ylab="Martingale Residuals")
lines(lowess(rr ~ age), col="red")
title("Martingale Residuals versus log age")
# no change keep as age

plot(rr ~ bili, lty=2, log="x", xlab="Bili", ylab="Martingale Residuals")
lines(lowess(rr ~ bili), col="red")
title("Martingale Residuals versus log bili")
# there is change - think about log of bili

plot(rr ~ copper, lty=2, log="x", xlab="Copper", ylab="Martingale Residuals")
lines(lowess(rr ~ copper), col="red")
title("Martingale Residuals versus log bili")
# there is a change think about log of copper

plot(rr ~ albumin, lty=2, log="x", xlab="Albunim", ylab="Martingale Residuals")
lines(lowess(rr ~ albumin), col="red")
title("Martingale Residuals versus log bili")
# no change - keep as albumin

plot(rr ~ protime, lty=2, log="x", xlab="Protime", ylab="Martingale Residuals")
lines(lowess(rr ~ protime), col="red")
title("Martingale Residuals versus log bili")
# no change - keep as albumin


# Individual Spline transformations
(fit.spline1 <- coxph(Surv(futime,status) ~ pspline(age) ) )
termplot(fit.spline1, term=1, se=TRUE, col.term=1, col.se="blue")
# keep as linear

(fit.spline2 <- coxph(Surv(futime,status) ~ pspline(bili) ) )
termplot(fit.spline2, term=1, se=TRUE, col.term=1, col.se="blue")
# either linear or nonlinear works - try log of bili

(fit.spline3 <- coxph(Surv(futime,status) ~ pspline(copper) ) )
termplot(fit.spline3, term=1, se=TRUE, col.term=1, col.se="blue")
# either way works but linear might be better and easier to work with
# looks like a cubed transformation or log of copper

(fit.spline4 <- coxph(Surv(futime,status) ~ pspline(albumin) ) )
termplot(fit.spline4, term=1, se=TRUE, col.term=1, col.se="blue")
# linear

(fit.spline5 <- coxph(Surv(futime,status) ~ pspline(protime) ) )
termplot(fit.spline5, term=1, se=TRUE, col.term=1, col.se="blue")
# either way works but linear might be better and easier to work with
# hard to tell what type of transformation here


# Checking transformations in R 
finalmodel <- coxph(Surv(futime, status) ~ bili + copper + albumin + protime 
                        + age + as.factor(stage) + as.factor(drug))

summary(finalmodel)
step(finalmodel, direction = "both", criterion = "AIC")

# After transformation models test variabe selection
finalmodel2 <- coxph(Surv(futime, status) ~ bili + copper + albumin + protime + 
                       age + as.factor(stage) + as.factor(drug) + log(bili) + 
                       log(copper))
step(finalmodel2, direction = "both", criterion = "AIC")

finalmodel3 <- coxph(Surv(futime, status) ~ albumin + protime + 
                       age + as.factor(stage) + as.factor(drug) + log(bili) + 
                       log(copper))
step(finalmodel3, direction = "both", criterion = "AIC")

# only variables left are albumin, protime, age, log(bili) and log(copper)

# Final Model from SAS Output 

finalmodel.pbc1 <- coxph(Surv(futime, status) ~ albumin + protime + age + 
                          log(bili) + log(copper) + as.factor(drug))
summary(finalmodel.pbc1)

finalmodel.pbc2 <- coxph(Surv(futime, status) ~ albumin + protime + age + 
                           log(bili) + log(copper))
summary(finalmodel.pbc2)

## Check proportional hazards assumption

(cox.prop.assump1 <- cox.zph(finalmodel.pbc1) )
(cox.prop.assump2 <- cox.zph(finalmodel.pbc2) )

# protime has a p-value less than 0.05 in both models

## Checking the DfBeta Residuals for Final Model 

# DFbetas for Cox Model 

# for Variable 1 - albunim

rr.dfbeta <- residuals(finalmodel.pbc1, type="dfbeta")
subject.index <- 1:nrow(rr.dfbeta)
plot(rr.dfbeta[,1] ~ subject.index, axes=F)
axis(1, at=seq(from=1, to=length(subject.index), by=5))
axis(2)
title("DFbeta for Final Model for Albunim")

# for Variable 2 - protime

rr.dfbeta <- residuals(finalmodel.pbc1, type="dfbeta")
subject.index <- 1:nrow(rr.dfbeta)
plot(rr.dfbeta[,2] ~ subject.index, axes=F)
axis(1, at=seq(from=1, to=length(subject.index), by=5))
axis(2)
title("DFbeta for Final Model for Protime")

# for Variable 3 - age

rr.dfbeta <- residuals(finalmodel.pbc1, type="dfbeta")
subject.index <- 1:nrow(rr.dfbeta)
plot(rr.dfbeta[,3] ~ subject.index, axes=F)
axis(1, at=seq(from=1, to=length(subject.index), by=5))
axis(2)
title("DFbeta for Final Model for Age")

# for Variable 4 - log(bili)

rr.dfbeta <- residuals(finalmodel.pbc1, type="dfbeta")
subject.index <- 1:nrow(rr.dfbeta)
plot(rr.dfbeta[,5] ~ subject.index, axes=F)
axis(1, at=seq(from=1, to=length(subject.index), by=5))
axis(2)
title("DFbeta for Final Model for Log of Bili")

# for Variable 5 - log(copper)

rr.dfbeta <- residuals(finalmodel.pbc1, type="dfbeta")
subject.index <- 1:nrow(rr.dfbeta)
plot(rr.dfbeta[,5] ~ subject.index, axes=F)
axis(1, at=seq(from=1, to=length(subject.index), by=5))
axis(2)
title("DFbeta for Final Model for Log of Copper")

# for Variable 6 - drug

rr.dfbeta <- residuals(finalmodel.pbc1, type="dfbeta")
subject.index <- 1:nrow(rr.dfbeta)
plot(rr.dfbeta[,6] ~ subject.index, axes=F)
axis(1, at=seq(from=1, to=length(subject.index), by=5))
axis(2)
title("DFbeta for Final Model for Treatment")

## Deviance Residual Plots

resid.deviance <- residuals(finalmodel.pbc1, type = "deviance")
plot(resid.deviance ~ age, xlab = "Age at Diagnosis", ylab = "Deviance Residuals", 
     main = "Deviance Residuals versus Age at Diagnosis")

plot(resid.deviance ~ protime, xlab = "Protime", ylab = "Deviance Residuals", 
     main = "Deviance Residuals versus Protime")

plot(resid.deviance ~ albumin, xlab = "Albumin", ylab = "Deviance Residuals", 
     main = "Deviance Residuals versus Albunim")

plot(resid.deviance ~ log(copper), xlab = "Log(Copper)", ylab = "Deviance Residuals", 
     main = "Deviance Residuals versus Log(Copper")

plot(resid.deviance ~ log(bili), xlab = "Log(Bili)", ylab = "Deviance Residuals", 
     main = "Deviance Residuals versus log(Bili")

## Case-Deletion Residual Plot for Final Model and identifying outliers
resid.dfbeta <- residuals(finalmodel.pbc1, type = "dfbeta")
n.obs <- length(futime)
index.obs <- 1:n.obs
plot(resid.dfbeta[,1] ~ index.obs, type = "h", xlab="Observation Number", ylab="Change in Albumin Coefficent", 
     main = "Case Deletion Plot for the Variable Albumin")
abline(h=0)

resid.dfbeta <- residuals(finalmodel.pbc1, type = "dfbeta")
n.obs <- length(futime)
index.obs <- 1:n.obs
plot(resid.dfbeta[,2] ~ index.obs, type = "h", xlab="Observation Number", ylab="Change in Protime Coefficent", 
     main = "Case Deletion Plot for the Variable Protime")
abline(h=0)

resid.dfbeta <- residuals(finalmodel.pbc1, type = "dfbeta")
n.obs <- length(futime)
index.obs <- 1:n.obs
plot(resid.dfbeta[,3] ~ index.obs, type = "h", xlab="Observation Number", ylab="Change in Age Coefficent", 
     main = "Case Deletion Plot for the Variable Age")
abline(h=0)

resid.dfbeta <- residuals(finalmodel.pbc1, type = "dfbeta")
n.obs <- length(futime)
index.obs <- 1:n.obs
plot(resid.dfbeta[,4] ~ index.obs, type = "h", xlab="Observation Number", ylab="Change in Log(Bili) Coefficent", 
     main = "Case Deletion Plot for the Variable Log(Bili)")
abline(h=0)

resid.dfbeta <- residuals(finalmodel.pbc1, type = "dfbeta")
n.obs <- length(futime)
index.obs <- 1:n.obs
plot(resid.dfbeta[,5] ~ index.obs, type = "h", xlab="Observation Number", ylab="Change in Log(Copper) Coefficent", 
     main = "Case Deletion Plot for the Variable Log(Copper)")
abline(h=0)

resid.dfbeta <- residuals(finalmodel.pbc1, type = "dfbeta")
n.obs <- length(futime)
index.obs <- 1:n.obs
plot(resid.dfbeta[,6] ~ index.obs, type = "h", xlab="Observation Number", ylab="Change in Drug Coefficent", 
     main = "Case Deletion Plot for the Variable Drug")
abline(h=0)

## Figures for Project Report

PBC276v2 <- PBC276 %>% 
  mutate(drug = as.numeric(drug)) %>%
  mutate(sex = as.numeric(sex))

View(PBC276v2)
attach(PBC276v2)

final.model.pbc <- coxph(Surv(futime, status) ~ drug + sex 
                         + age + log(bili) + albumin + copper + sgot + protime 
                         + as.factor(stage))
summary(final.model.pbc)

PBC276v2 <- PBC276
PBC276v2$drug[PBC276v2$drug == "0"] <- "Placebo"
PBC276v2$drug[PBC276v2$drug == "1"] <- "Treatment"



library(ggplot2)
library(survminer)

resid.dfbeta <- residuals(final.model.pbc, type = "dfbeta")
n.obs <- length(futime)
index.obs <- 1:n.obs

# Potential Figures

# DfBeta plot

# add this plot to google 
y = c("Drug", "Sex", "Age", "Log(Bili)", "Albumin", "Copper","SGOT","Prothrombin Time")
par(mfrow=c(2,4))
for (i in c(1:8)){
  rr.dfbeta <- residuals(final.model.pbc, type="dfbeta")
  subject.index <- 1:nrow(rr.dfbeta)
  plot(rr.dfbeta[,i] ~ subject.index, axes=F, ylab=y[i])
  axis(1, at=seq(from=1, to=length(subject.index), by=5))
  axis(2)
}



rr.dfbeta <- residuals(final.model.pbc, type="dfbeta")
subject.index <- 1:nrow(rr.dfbeta)
plot(rr.dfbeta[,1] ~ subject.index, axes=F)
axis(1, at=seq(from=1, to=length(subject.index), by=5))
axis(2)
title("DFbeta for Final Model for Albumin")

# Residual Deviance Plot
ggplot(PBC276, aes(y=resid.deviance, x=log(bili))) + 
  geom_point(aes(color = as.factor(drug))) +
  facet_wrap(~as.factor(drug)) + 
  labs(title = "Case Deletion Plot for the Variable Albunim", 
       subtitle = "Need subtitle?",
       y = "Deviance Residuals",
       x = "Log(Bili)",
       color = "Treatment",
       main = "Deviance Residuals versus log(Bili") + 
  theme_bw() + 
  theme(legend.position = "bottom")

# Case Deletion Plot
ggplot(PBC276, aes(y=resid.dfbeta[,1], x=index.obs)) + 
  geom_point(aes(color = as.factor(drug))) +
  facet_wrap(~as.factor(drug)) + 
  labs(title = "Case Deletion Plot for the Variable Albunim", 
       subtitle = "Need subtitle?",
       y = "Change in Drug Coefficent",
       x = "Observation Number",
       color = "Treatment") + 
  theme_bw() + 
  theme(legend.position = "bottom")


# plots final model 
par(mfrow=c(1,1))
plot(survfit(final.model.pbc), 
     ylab = "Probability of Survival", 
     xlab = "Time", col = c("black", "black", "black"),
     main = "Cox Proportional Hazards Model",
     conf.int = FALSE) 

par(mfrow=c(2,4))
plot(cox.zph(final.model.pbc)) # add to Google doc
ggcoxzph(cox.zph(final.model.pbc))

# Fit (complete) survival curves
#++++++++++++++++++++++++++++++++++++
require("survival")
fit2 <- survfit( Surv(futime, status) ~ as.factor(drug) + as.factor(sex),
                 data = PBC276v2 )
# Visualize
#++++++++++++++++++++++++++++++++++++
# Visualize: add p-value, change y limits
# change color using brewer palette
ggsurvplot(fit2, pval = TRUE, 
           break.time.by = 400,
           risk.table = FALSE,
           risk.table.height = 0.5#Useful when you have multiple groups
)


# estimate survival based on here- not working like hoped 
# have to add more variables to work
est <- survfit(final.model.pbc, 
               newdata = data.frame(albumin=albumin,
                                    protime=protime,
                                    age=age,
                                    bili=2000,
                                    copper=1000,
                                    drug=="Placebo"))
lines(est$time, est$surv, col = 'blue', type = 's')


coef <- c(0.02607,-0.10620,0.02547,0.73970,-0.84222,0.00192,0.00254,0.31222,1.30263,1.59215,1.82491)
sd <- c(0.20146,0.30254,0.01073,0.12581,0.25978,0.00118,0.00190,0.10669,1.07317,1.03864,1.04185)

coef-1.96*sd
coef+1.968*sd

# Identify outliers in variable albumin
identify(rr.dfbeta[,1] ~ subject.index)

identify(resid.dfbeta[,1] ~ index.obs)
which(resid.dfbeta[,1] > 0.05) # outliers identified here
which(resid.dfbeta[,1] < -0.05)

# Identify outliers in variable protime
identify(rr.dfbeta[,2] ~ subject.index)

identify(resid.dfbeta[,2] ~ index.obs)
which(resid.dfbeta[,2] > 0.03) # outliers identified here
which(resid.dfbeta[,2] < -0.03)

# Identify outliers in variable age
identify(rr.dfbeta[,3] ~ subject.index)

identify(resid.dfbeta[,3] ~ index.obs)
which(resid.dfbeta[,3] > 0.004) # outliers identified here
which(resid.dfbeta[,3] < -0.004)

# Identify outliers in variable log(bili)
identify(rr.dfbeta[,4] ~ subject.index)

identify(resid.dfbeta[,4] ~ index.obs)
which(resid.dfbeta[,4] > 0.04) # outliers identified here
which(resid.dfbeta[,4] < -0.04)

# Identify outliers in variable log(copper)
identify(rr.dfbeta[,5] ~ subject.index)

identify(resid.dfbeta[,5] ~ index.obs)
which(resid.dfbeta[,5] > 0.05) # outliers identified here
which(resid.dfbeta[,5] < -0.05)

# Identify outliers in variable drug
identify(rr.dfbeta[,6] ~ subject.index)

identify(resid.dfbeta[,6] ~ index.obs)
which(resid.dfbeta[,6] > 0.04) # outliers identified here
which(resid.dfbeta[,6] < -0.04)
