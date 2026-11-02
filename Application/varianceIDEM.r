#required packages
library(MEMSS)
library(survival)
library(MASS)
library(fitdistrplus)
library(sandwich)
library(wCorr)



## data analysis ###############################################################
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

data_final<- cbind(newData, dataIDT[,c('Subject','density','density_stab')])

data_final$idt <- data_final$density_stab/(data_final$density)



OLS <- lm(conc ~ 1 + cumuldose, data=data_final)$coef[2] 

IDT <- lm(conc ~ 1 + cumuldose, weights=idt, data=data_final)$coef[2]

IIV <- lm(conc ~ 1 + cumuldose, weights=1/rho_i, data=data_final)$coef[2]

IDEM <- lm(conc~ 1 + cumuldose, weight=idt*1/rho_i, data=data_final)$coef[2]


estimates <- matrix(c(OLS,IDT,IIV,IDEM), ncol=4)
colnames(estimates) <- c('OLS','IDT','IIV','IDEM')


## Two-step variance estimation #################################################

betas<-matrix(lm(conc~ 1 + cumuldose, weight=idt*1/rho_i , data=data_final)$coef,ncol=1)


allrows<-as.data.frame(cbind(c(rep(1,length(data_final$Subject))),data_final))


allvisits<-allrows[allrows$visit==1,]

colnames(allvisits)[1]<-'intercept'
colnames(allrows)[1]<-'intercept'


Xvisits<-allvisits[,c(1,11)]


Kvisits<-allvisits[,c(1,3:4)]
Kvisits$Apgar <- as.numeric(Kvisits$Apgar)

alpha<- matrix(glm(log(data_final$doset1)~data_final$Wt+as.numeric(data_final$Apgar))$coef,ncol=1)
gamma<-matrix(glm(data_final$visit~data_final$cumuldose+data_final$Wt + as.numeric(data_final$Apgar)
                    ,family='binomial')$coef)
expgammavisits<- exp(gamma[1]+ gamma[2]*allvisits$cumuldose+gamma[3]*allvisits$Wt+
                       gamma[4]*as.numeric(allvisits$Apgar))
oneoverexpgammavisits<- 1/expgammavisits
expomegavisits<- matrix(exp( t(alpha)%*%t(Kvisits)),ncol=1)
betasvisits<- matrix(t(t(betas)%*%t(Xvisits)),ncol=1)
yminusbetasvisits<- allvisits$conc-betasvisits



nbpersons<-length(data_final[!duplicated(data_final$Subject),]$Subject) 
nit<- length(data_final[is.na(data_final$conc)==0,]$conc) 
nbrows<-length(data_final[,1]) 


########################## 
## 1. compute the g(o): ##
##########################

tt2<-glm(conc~ 1 + cumuldose, weight=idt*1/rho_i , data=data_final)

estfun.glm2<-function(glm.obj){
  if (is.matrix(glm.obj$x))
    xmat<-glm.obj$x
  else {
    mf<-model.frame(glm.obj)
    xmat<-model.matrix(terms(glm.obj),mf)
  }
  residuals(glm.obj,"working")*glm.obj$weights*xmat   }

aa5<-estfun.glm2(tt2)
go<- t(aa5)

vec<-matrix(NA,nrow=2,ncol=nbrows)
indic<-ifelse(is.na(data_final$conc),0,1)
ll<-1

for(ik in 1:nbrows){   
  if(indic[ik]==1){vec[,ik]<-go[,ll]; ll<-ll+1}
  if(indic[ik]==0){vec[,ik]<-c(rep(0,2)) }
}


################### 
## 2. compute mo ##
###################

resfoismat<-function(glm.obj){
  if (is.matrix(glm.obj$x))
    xmat<-glm.obj$x
  else {
    mf<-model.frame(glm.obj)
    xmat<-model.matrix(terms(glm.obj),mf)
  }
  residuals(glm.obj,"working")*glm.obj$weights*xmat 
}

dg <- data_final[!duplicated(data_final$Subject),]

m1<- glm(log(dg$doset1)~dg$Wt + as.numeric(dg$Apgar)) 
resm1<-resfoismat(m1)

inda<-ifelse(is.na(data_final$conc),0,1)
resm1m<-matrix(NA, nrow=nbrows,ncol=3)

resm1m[1,]<- resm1[1,]

jindex<-2

for(i in 2:nbrows){
  if(data_final$Subject[i]==data_final$Subject[i-1] ){  resm1m[i,] <-c(0,0,0) }
  if(data_final$Subject[i]!=data_final$Subject[i-1] ){  resm1m[i,] <- resm1[jindex,] ; jindex<-jindex+1} 
}

m2<- glm(data_final$visit~data_final$cumuldose+data_final$Wt + as.numeric(data_final$Apgar),
         family='binomial')
resm2<-resfoismat(m2)
mo<- rbind(t(resm1m[,1]), t(resm1m[,2]), t(resm1m[,3]), t(resm2[,1]), t(resm2[,2]), t(resm2[,3]), t(resm2[,4]))


###################
## 3. compute M  ##
###################


crossprod<- cbind( c(rep(1,length(dg[,1]) ) ), dg$Wt,as.numeric(dg$Apgar))

Mrow1<- c(1/nbrows*sum(-( t(alpha) %*% t(crossprod))),
          1/nbrows*sum(-( t(alpha) %*% t(crossprod))*dg$Wt) ,
          1/nbrows*sum(-( t(alpha) %*% t(crossprod))*as.numeric(dg$Apgar))
          ,0,0,0,0 )

Mrow2<- c(1/nbrows*sum(-dg$Wt*(t(alpha) %*% t(crossprod))),
          1/nbrows*sum(-( t(alpha) %*% t(crossprod))*dg$Wt^2) ,
          1/nbrows*sum(-( t(alpha) %*% t(crossprod))*dg$Wt*as.numeric(dg$Apgar))
          , 0,0,0 ,0 )

Mrow3<- c(1/nbrows*sum(as.numeric(dg$Apgar)*( t(alpha) %*% t(crossprod))),
          1/nbrows*sum(-( t(alpha) %*% t(crossprod))*dg$Wt*as.numeric(dg$Apgar)) ,
          1/nbrows*sum(-( t(alpha) %*% t(crossprod))*as.numeric(dg$Apgar)^2)
          , 0,0,0,0 )


gammav2<-matrix(glm(data_final$visit~data_final$cumuldose+data_final$Wt+as.numeric(data_final$Apgar),
                    family='binomial')$coef)

double<-cbind(c(rep(1,length(data_final[,1]))) ,data_final$cumuldose,data_final$Wt,as.numeric(data_final$Apgar))

Mrow4<- c(0,0,0, 1/nbrows*sum( -exp( t(gammav2)%*%t(double) ) /(1+exp( t(gammav2)%*%t(double)) )^2  ), 
          1/nbrows*sum(- data_final$cumuldose*exp( t(gammav2)%*%t(double) ) /(1+exp( t(gammav2)%*%t(double)) )^2  ),
          1/nbrows*sum(- data_final$Wt*exp( t(gammav2)%*%t(double) ) /(1+exp( t(gammav2)%*%t(double)) )^2  ),
          1/nbrows*sum(- as.numeric(data_final$Apgar)*exp( t(gammav2)%*%t(double) ) /(1+exp( t(gammav2)%*%t(double)) )^2  )) 

Mrow5<- c(0,0,0, 1/nbrows*sum( -data_final$cumuldose*exp( t(gammav2)%*%t(double) ) /(1+exp( t(gammav2)%*%t(double)) )^2  ), 
          1/nbrows*sum(-data_final$cumuldose*data_final$cumuldose*exp( t(gammav2)%*%t(double) ) /(1+exp( t(gammav2)%*%t(double)) )^2  ),
          1/nbrows*sum(-data_final$cumuldose*data_final$Wt*exp( t(gammav2)%*%t(double) ) /(1+exp( t(gammav2)%*%t(double)) )^2  ),
          1/nbrows*sum(-data_final$cumuldose*as.numeric(data_final$Apgar)*exp( t(gammav2)%*%t(double) ) /(1+exp( t(gammav2)%*%t(double)) )^2  )) 

Mrow6<- c(0,0,0, 1/nbrows*sum( -data_final$Wt*exp( t(gammav2)%*%t(double) ) /(1+exp( t(gammav2)%*%t(double)) )^2  ), 
          1/nbrows*sum(-data_final$Wt*data_final$cumuldose*exp( t(gammav2)%*%t(double) ) /(1+exp( t(gammav2)%*%t(double)) )^2  ),
          1/nbrows*sum(-data_final$Wt*data_final$Wt*exp( t(gammav2)%*%t(double) ) /(1+exp( t(gammav2)%*%t(double)) )^2  ),
          1/nbrows*sum(-data_final$Wt*as.numeric(data_final$Apgar)*exp( t(gammav2)%*%t(double) ) /(1+exp( t(gammav2)%*%t(double)) )^2)) 

Mrow7<- c(0,0,0, 1/nbrows*sum( -as.numeric(data_final$Apgar)*exp( t(gammav2)%*%t(double) ) /(1+exp( t(gammav2)%*%t(double)) )^2  ), 
          1/nbrows*sum(-as.numeric(data_final$Apgar)*data_final$cumuldose*exp( t(gammav2)%*%t(double) ) /(1+exp( t(gammav2)%*%t(double)) )^2  ),
          1/nbrows*sum(-as.numeric(data_final$Apgar)*data_final$Wt*exp( t(gammav2)%*%t(double) ) /(1+exp( t(gammav2)%*%t(double)) )^2  ),
          1/nbrows*sum(-as.numeric(data_final$Apgar)*as.numeric(data_final$Apgar)*exp( t(gammav2)%*%t(double) ) /(1+exp( t(gammav2)%*%t(double)) )^2)) 

M<- rbind(Mrow1,Mrow2,Mrow3,Mrow4,Mrow5,Mrow6,Mrow7)


################
## 4. Gbeta  ###
################


g1<- c(1/nbrows*sum(allvisits$idt*1/allvisits$rho_i*-allvisits[,1]*allvisits[,1]) , 1/nbrows*sum(allvisits$idt*1/allvisits$rho_i*-allvisits[,11]*allvisits[,1])) 
g2<- c(1/nbrows*sum(allvisits$idt*1/allvisits$rho_i*-allvisits[,1]*allvisits[,11]) , 1/nbrows*sum(allvisits$idt*1/allvisits$rho_i*-allvisits[,11]*allvisits[,11])) 

gbeta<- rbind(g1,g2)


################
### 5. G psi ###
################

ag<-allvisits[,c(1,11)]
allk<-allvisits[,c(1,3:4)]
allk$Apgar <- as.numeric(allk$Apgar)

gp1<- c(1/nbrows*sum(allvisits$idt/sigma2*(log(allvisits$doset1)-t(t(alpha)%*%t(allk)))*1/allvisits$rho_i*(allvisits$conc-t(t(betas)%*%t(ag)))*ag[,1]),
        1/nbrows*sum(allvisits$idt/sigma2*(log(allvisits$doset1)-t(t(alpha)%*%t(allk)))*1/allvisits$rho_i*(allvisits$conc-t(t(betas)%*%t(ag)))*ag[,1]* allk[,2]), 
        1/nbrows*sum(allvisits$idt/sigma2*(log(allvisits$doset1)-t(t(alpha)%*%t(allk)))*1/allvisits$rho_i*(allvisits$conc-t(t(betas)%*%t(ag)))*ag[,1] * allk[,3]), 
        1/nbrows*sum(allvisits$idt*(allvisits$conc-t(t(betas)%*%t(ag)))*ag[,1]*  -1*1/allvisits$rho_i ), 
        1/nbrows*sum(allvisits$idt*(allvisits$conc-t(t(betas)%*%t(ag)))*ag[,1]* -ag[,2] *1/allvisits$rho_i ), 
        1/nbrows*sum(allvisits$idt*(allvisits$conc-t(t(betas)%*%t(ag)))*ag[,1]* -allvisits$Wt * 1/allvisits$rho_i ),
        1/nbrows*sum(allvisits$idt*(allvisits$conc-t(t(betas)%*%t(ag)))*ag[,1]* -as.numeric(allvisits$Apgar) * 1/allvisits$rho_i ))

gp2<- c(1/nbrows*sum(allvisits$idt/sigma2*(log(allvisits$doset1)-t(t(alpha)%*%t(allk)))*1/allvisits$rho_i*(allvisits$conc-t(t(betas)%*%t(ag)))*ag[,2]),
        1/nbrows*sum(allvisits$idt/sigma2*(log(allvisits$doset1)-t(t(alpha)%*%t(allk)))*1/allvisits$rho_i*(allvisits$conc-t(t(betas)%*%t(ag)))*ag[,2]* allk[,2]), 
        1/nbrows*sum(allvisits$idt/sigma2*(log(allvisits$doset1)-t(t(alpha)%*%t(allk)))*1/allvisits$rho_i*(allvisits$conc-t(t(betas)%*%t(ag)))*ag[,2] * allk[,3]), 
        1/nbrows*sum(allvisits$idt*(allvisits$conc-t(t(betas)%*%t(ag)))*ag[,2]*  -1*1/allvisits$rho_i ), 
        1/nbrows*sum(allvisits$idt*(allvisits$conc-t(t(betas)%*%t(ag)))*ag[,2]* -ag[,2] *1/allvisits$rho_i ), 
        1/nbrows*sum(allvisits$idt*(allvisits$conc-t(t(betas)%*%t(ag)))*ag[,2]* -allvisits$Wt * 1/allvisits$rho_i ),
        1/nbrows*sum(allvisits$idt*(allvisits$conc-t(t(betas)%*%t(ag)))*ag[,2]* -as.numeric(allvisits$Apgar) * 1/allvisits$rho_i ))



Gpsi<-rbind(gp1,gp2) 


##############
## variance ##
##############

psipart<-Gpsi %*% ( -solve(M)%*%mo ) 
var<- 1/nbrows^2 * solve(gbeta) %*%( ( vec +    psipart ) %*% t( vec + psipart)) %*%  solve(gbeta)
sandwich <- var[2,2]
sqrt(sandwich)

(CD_IDEM <- c(estimates[4]-1.96*sqrt(sandwich), estimates[4]+1.96*sqrt(sandwich)))



 