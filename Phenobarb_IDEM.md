Phenobarb Dataset Application
================

Required R packages

``` r
library(MEMSS)
library(survival)
```

Phenobarb dataset : We are interested in estimating the causal effect of
dose on conc

``` r
data <- Phenobarb
```

Creation of times t1 and t2 to define the beginning and end of time
point

``` r
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
```

Manual modifications for specific subjects : additional doses given
within one hour of the first dose (time $\leq$ 1) were combined with the
initial dose (at time = 0). For subjects with a concentration response
available within the first hour, we considered that this was observed at
the time of the initial dose.

``` r
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
```

Creation of the cumulative dose variable (cumuldose)

``` r
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
```

Creation of visits indicator variable

``` r
newData$visit <- ifelse(is.na(newData$conc) == TRUE,0,1)
```

Creation of initial prescribed dose (doset1)

``` r
newData$doset1[1] <- newData$dose[1]
for (j in 2:length(newData$Subject)){
  if (newData$Subject[j] == newData$Subject[j-1]){
      newData$doset1[j] <- newData$doset1[j-1]}
  else newData$doset1[j] <- newData$dose[j]
}
```

IIV weight : we assume a proportional rate model

``` r
coxfit <- coxph(Surv(t1,t2, visit) ~ cumuldose+Wt+as.numeric(Apgar), data=newData)

gamma  <- coxfit$coef 

newData$rho_i<- exp(gamma[1]*newData$cumuldose + gamma[2]*newData$Wt + gamma[3]*as.numeric(newData$Apgar))    
```

IDT weight : we assume a normal distribution

``` r
dataIDT<- newData

dataIDT$pred_dose <-predict(glm(doset1 ~ Wt + as.numeric(Apgar), data=dataIDT),type='response')
dataIDT$residuals <- dataIDT$doset1 - dataIDT$pred_dose
sigma2<- var(c(dataIDT$residuals))

dataIDT$pred_doset1_stab <- predict(glm(doset1 ~ 1, data=dataIDT), type='response')
dataIDT$residuals_stab <- dataIDT$doset1 - dataIDT$pred_doset1_stab
sigma2_stab <- var(dataIDT$residuals_stab)


dataIDT$density <- 1/sqrt(2*pi)*1/sqrt(sigma2)*exp(-dataIDT$residuals^2/(2*sigma2))
dataIDT$density_stab <- 1/sqrt(2*pi)*1/sqrt(sigma2_stab)*exp(-dataIDT$residuals_stab^2/(2*sigma2_stab))

data_final<- cbind(newData, dataIDT[,c('Subject','density','density_stab')])

data_final$idt <- data_final$density_stab/(data_final$density)
```

ATE estimation with OLS, IDT, IIV and IDEM

``` r
OLS <- lm(conc ~ 1 + cumuldose, data=data_final)$coef[2] 

IDT <- lm(conc ~ 1 + cumuldose, weights=idt, data=data_final)$coef[2]

IIV <- lm(conc ~ 1 + cumuldose, weights=1/rho_i, data=data_final)$coef[2]

IDEM <- lm(conc~ 1 + cumuldose, weight=idt*1/rho_i, data=data_final)$coef[2]


estimates <- matrix(c(OLS,IDT,IIV,IDEM), ncol=4)
colnames(estimates) <- c('OLS','IDT','IIV','IDEM')
```
