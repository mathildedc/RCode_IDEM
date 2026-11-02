# required packages
library(MEMSS)
library(survival)
library(MASS)
library(fitdistrplus)
library(sandwich)
library(wCorr)
library(tableone)
library(survey)
library(dplyr)
library(table1)

data <- Phenobarb

data$t1[1] <- 0
data$t2[1] <- 1

for (j in 2:length(data$Subject)){
  if (data$Subject[j]==data$Subject[j-1]){
    data$t1[j]<-data$t2[j-1] 
    data$t2[j]<-data$time[j]
  }
  else {
    data$t1[j] <- 0
    data$t2[j] <- 1
  }
}


#subject 24 : 63+63 at baseline (instead of t=0 and t=0.7)
data[340,6] <- 63+63
data[342,8] <- 1

#subject 31 : observation Y at t=1 (t1=1 and t2=1)
data[421,7] <- data[422,7]

#subject 40 : observation Y at t=1 (t1=1 and t2=1)
data[527,7] <- data[528,7]

#subject 46 : 11+11 at baseline (instead of t=0 and t=0.5)
data[581,6] <- 11+11
data[583,8] <- 1

#subject 57 : 20+20 at baseline (instead of t=0 and t=0.5)
data[702,6] <- 20+20
data[704,8] <- 1

#subject 58 : 14+14 at baseline (instead of t=0 and t=0.5)
data[714,6] <- 14+14
data[716,8] <- 1

newData <- data[-c(341,422,528,582,703,715),]




newData$cumuldose[1] <- newData$dose[1]
for (j in 2:length(newData$Subject)){
  if (newData$Subject[j]==newData$Subject[j-1]){
    if (is.na(newData$dose[j]) == TRUE){
      newData$cumuldose[j] <- newData$cumuldose[j-1]
    }
    else newData$cumuldose[j] <- newData$cumuldose[j-1] + newData$dose[j] 
  }
  else newData$cumuldose[j] <- newData$dose[j]
}




newData$visit <- ifelse(is.na(newData$conc) == TRUE,0,1)


newData$doset1[1] <- newData$dose[1]
for (j in 2:length(newData$Subject)){
  if (newData$Subject[j] == newData$Subject[j-1]){
    newData$doset1[j] <- newData$doset1[j-1]}
  else newData$doset1[j] <- newData$dose[j]
}


# goodness of fit of normal, gamma and log-normal
fit_norm  <- fitdist(newData$doset1, "norm")
fit_gamma <- fitdist(newData$doset1, "gamma")
fit_lnorm <- fitdist(newData$doset1, "lnorm")
AIC_fit_norm  <- fitdist(newData$doset1, "norm")$aic
AIC_fit_gamma <- fitdist(newData$doset1, "gamma")$aic
AIC_fit_lnorm <- fitdist(newData$doset1, "lnorm")$aic


plot.legend <- c("normal (AIC=6503)", "Gamma (AIC=5911)", "log-normal (AIC=5736)")
denscomp(list(fit_norm, fit_gamma, fit_lnorm), legendtext = plot.legend,
         xlab = "Initial dose", main="")


coxfit <- coxph(Surv(t1,t2, visit) ~ cumuldose+Wt+as.numeric(Apgar), data=newData)

gamma  <- coxfit$coef 

newData$rho_i<- exp(gamma[1]*newData$cumuldose + gamma[2]*newData$Wt + gamma[3]*as.numeric(newData$Apgar))    



dataIDT<- newData

dataIDT$pred_dose <-predict(glm(log(doset1) ~ Wt + as.numeric(Apgar), data=dataIDT),type='response')

dataIDT$residuals <- log(dataIDT$doset1) - dataIDT$pred_dose
sigma2<- var(c(dataIDT$residuals))

dataIDT$pred_doset1_stab <- predict(glm(log(doset1) ~ 1, data=dataIDT), type='response')
dataIDT$residuals_stab <- log(dataIDT$doset1) - dataIDT$pred_doset1_stab
sigma2_stab <- var(dataIDT$residuals_stab)


dataIDT$density <- 1/(dataIDT$doset1*sqrt(2*pi))*1/sqrt(sigma2)*exp(-dataIDT$residuals^2/(2*sigma2))
dataIDT$density_stab <- 1/(dataIDT$doset1*sqrt(2*pi))*1/sqrt(sigma2_stab)*exp(-dataIDT$residuals_stab^2/(2*sigma2_stab))

data_final<- cbind(newData, dataIDT[,c('density','density_stab')])

data_final$idt <- data_final$density_stab/(data_final$density)



OLS <- lm(conc ~ 1 + cumuldose, data=data_final)$coef[2] 

IDT <- lm(conc ~ 1 + cumuldose, weights=idt, data=data_final)$coef[2]

IIV <- lm(conc ~ 1 + cumuldose, weights=1/rho_i, data=data_final)$coef[2]

IDEM <- lm(conc~ 1 + cumuldose, weight=idt*1/rho_i, data=data_final)$coef[2]


estimates <- matrix(c(OLS,IDT,IIV,IDEM), ncol=4)
colnames(estimates) <- c('OLS','IDT','IIV','IDEM')


# Results in tables

#table 1
data_final$dose_median <- ifelse(is.na(data_final$dose) == TRUE,3, ifelse(data_final$dose > 24 ,0,1))
table1(~Wt + as.numeric(Apgar) + ApgarInd + conc | time*dose_median
       , data = data_final)


#Balance before/after IDT weighting
dat_tab<- data_final
dat_tab$Apgar<- as.numeric(data_final$Apgar)
myVars<- c('Wt', 'Apgar')
dat_tab <- dat_tab |> mutate(dose_strat = ntile(doset1, 4))
tab <- CreateTableOne(vars = myVars, data = dat_tab , strata='dose_strat')
tab
nhanes.IDT.lr <- svydesign(ids=~0, data= dat_tab, weights= ~idt)
tabIDT <- svyCreateTableOne(vars = myVars, strata = "dose_strat", data = nhanes.IDT.lr, test = FALSE)
tabIDT


#Log Prescribed dose linear model
summary(glm(log(doset1) ~ Wt + as.numeric(Apgar), data=dataIDT))
(c(0.5937-1.96*0.0173,0.5937+1.96*0.0173 ))
(c(-0.0336-1.96*0.0048,-0.0336+1.96*0.0048))



#Proportional rate model for IIV
summary(coxfit)


# sandwich variance
sd_OLS <- sqrt(sandwich(lm(conc ~ 1 + cumuldose, data=data_final))[2,2])
sd_IDT <- sqrt(sandwich(lm(conc ~ 1 + cumuldose, weights=idt, data=data_final))[2,2])
sd_IIV <- sqrt(sandwich(lm(conc ~ 1 + cumuldose, weights=1/rho_i, data=data_final))[2,2])
sd_IDEM <- sqrt(sandwich(lm(conc~ 1 + cumuldose, weight=idt*1/rho_i, data=data_final))[2,2])

CI_OLS <- c(estimates[1]-1.96*sd_OLS, estimates[1]+1.96*sd_OLS)
CI_IDT <- c(estimates[2]-1.96*sd_IDT, estimates[2]+1.96*sd_IDT)
CI_IIV <- c(estimates[3]-1.96*sd_IIV, estimates[3]+1.96*sd_IIV)
CI_IDEM <- c(estimates[4]-1.96*sd_IDEM, estimates[4]+1.96*sd_IDEM)

