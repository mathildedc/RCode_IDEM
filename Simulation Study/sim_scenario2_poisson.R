#required packages
library(survival)
library(rootSolve)
library(geex)
library(ggplot2)
library(reshape2)
library(dplyr)


set.seed(123)

n_sample <-  1000 # either 500 or 1000
n_simul <- 1000
tau <- 200


#parameters for prescribed dose
alpha0 <- 5.9; alphaK1 <- 0.5; alphaK2 <- -0.5; alphaK3 <- -0.04


#parameters for confounders
K1_mean<-5;  K1_sd<-0.5;  K2_prob<-0.55; K3_mean<-100; K3_sd<-5

#parameters for mediator
Z1_mean<-4; Z1_sd<-sqrt(2); Z2_mean<-2; Z2_sd<-sqrt(1)

#parameters for response model
beta0<-2; betaK1<-1; betaK2<-0.25; betaK3<--2; betaQ<-0; betaA<-1.4;betaZ<-3;
epsilon_sd<-0.1; phi_mean<-0; phi_sd<-0.2


#parameters for visits model 
gamma0<--3.5; gammaA<-0.01;  gammaZ<--1.2; gammaK1<-0.5 ; gammaK2<-0.5; gammaK3<-0.03 ; gammaQ<--3 


ATE_Mat <- data.frame(matrix(nr=n_simul, nc=1*10))
names(ATE_Mat) <- c("OLS","adjusted OLS", "IDT","IIV","IDEM", 
                    "IDEM with 0-97.5 trunc. IIV", "IDEM with 0-99 trunc. IIV ", 
                    "IDEM with 20 str. FSW-IIV ", "IDEM with 500 str. FSW-IIV", "IDEM with NB IDT")




S <- 1 
for(S in 1  : n_simul){
  
  ## Dataset simulation
  
  # Initial matrix for the simulation of variables
  n_patients <- n_sample
  n_total <- n_patients*tau
  mat <- data.frame(matrix(NA, nrow=n_total, ncol= 11))
  
  names(mat)<- c('ID','time','D','Y','K1','K2','K3','Q','visit','Z', 'ratei')
  
  mat[,'ID'] <- rep(1:n_patients, each=tau)  
  mat[,'time'] <- rep(1:tau/100,  n_patients) 
  
  # Confounders
  mat[,'K1'] <- rep(rnorm(mean=K1_mean, sd=K1_sd, n=n_patients),each=tau)
  mat[,'K2'] <- rep(rbinom(p=K2_prob, size=1, n=n_patients),each=tau)
  mat[,'K3'] <- rep(rnorm(mean=K3_mean, sd=K3_sd, n=n_patients),each=tau)
  
  # Prescribed dose
  mat_oneID <- mat[!duplicated(mat$ID),] 
  lambda_dose <- exp(alpha0 + alphaK1*mat_oneID$K1+ alphaK2*mat_oneID$K2+ alphaK3*mat_oneID$K3)
  mat_oneID[,'doset1'] <- rpois(n=length(mat_oneID$ID), lambda = lambda_dose) 
  mat_final <- merge(mat, mat_oneID[, c('ID', 'doset1')], by = 'ID', all.x = TRUE)
  
  # Cumulative exposure
  mat_final$cumuldose <- ifelse(mat_final$doset1>100*mat_final$time,100*mat_final$time, floor(mat_final$doset1))
  
  # Mediator
  Z1 <- rnorm(mean=Z1_mean, sd=Z1_sd, n=n_total) 
  Z2 <- rnorm(mean=Z2_mean, sd=Z2_sd, n=n_total) 
  mat_final[,'Z'] <- ifelse(mat_final$cumuldose>90,Z2,Z1)                      
  centerZ <- ifelse(mat_final$cumuldose>90, Z2_mean, Z1_mean)
  
  mat_final[,'Q'] <- rnorm(mean=0.5, sd=0.1, n=n_total) 
  
  # Response variable
  phi <- rnorm(mean=phi_mean, sd=phi_sd,n=n_total)
  epsilon <- rnorm(mean=phi, sd=epsilon_sd, n=n_total)
  
  mat_final[,'Y'] <- beta0 + betaA*mat_final$cumuldose + betaK1*mat_final$K1 + betaK2*mat_final$K2 +
    betaK3*mat_final$K3 + betaZ*(mat_final$Z - centerZ) + betaQ*mat_final$Q+epsilon
  
  
  # Visits rate
  ratei <- exp(gamma0 + gammaA*mat_final[,'cumuldose'] + gammaZ*mat_final[,'Z'] + gammaK1*mat_final[,'K1'] +
                 gammaK2*mat_final[,'K2'] + gammaK3*mat_final[,'K3'] + gammaQ*mat_final[,'Q'] ) 
  
  ratei[ratei>1] <- 0.9999999999999 
  mat_final[,'ratei'] <- ratei
  
  
  # Visits indicator
  mat_final[,'visit'] <- rbinom(prob= ratei, size=1, n=n_total)  
  
  
  
  ## ATE Estimation
  
  Data <- mat_final
  Data[Data$visit==0,]$Y <- NA 
  
  
  # IIV weight: we assume a proportional rate model
  
  Data$t1 <- Data$time-0.01
  Data$t2 <- Data$time
  coxfit <- coxph(Surv(t1,t2,visit)~cumuldose+Z+K1+K2+K3+Q,data=Data)
  gamma <- coxfit$coef 
  
  
  # Baseline hazard
  baseh <- basehaz(coxfit, centered=FALSE)
  colnames(baseh)[2] <- 'time2'
  haz <- baseh[,1][-1] - baseh[,1][-(length(baseh[,1]))]
  baseh$col3 <- round(baseh[,2],2) 
  baseh$haz <- c(baseh$hazard[1],haz)
  Data$col3<- round(Data$time,2)
  dat1 <- merge( baseh, Data, all.y=TRUE, by='col3') 
  Data <- dat1[order(dat1$ID, dat1$time),]
  
  
  Data$w_iiv <- Data$haz*exp(gamma[1]* Data$cumuldose + gamma[2]*Data$Z + gamma[3]*Data$K1 + gamma[4]*Data$K2+ 
                               gamma[5]*Data$K3 + gamma[6]*Data$Q) 
  
  
  # Truncation 0-97.5 and 0-99:
  
  q_97.5 <- quantile(Data$w_iiv, probs = c(0, 0.975), na.rm = TRUE)
  Data$w_iiv_trunc_97.5 <- pmin(pmax(Data$w_iiv, q_97.5[1]), q_97.5[2])
  
  q_99 <- quantile(Data$w_iiv, probs = c(0, 0.99), na.rm = TRUE)
  Data$w_iiv_trunc_99 <- pmin(pmax(Data$w_iiv, q_99[1]), q_99[2]) 
  
  # FSW:
  
  # 20 strata
  Data <- Data %>%
    dplyr::group_by(time) %>%
    dplyr::mutate(q_20 = ntile(w_iiv, 20),
                  total_visits = sum(visit),
                  total_n = n()) %>%
    dplyr::group_by(time, q_20) %>%
    dplyr::mutate(q20_visits  = sum(visit),
                  q20_n = n(),
                  w_iiv_fsw_20 = (((q20_visits/q20_n)/(total_visits/total_n))^(-1))) %>%
    dplyr::ungroup() %>%
    as.data.frame()
  
  
  # 500 strata
  Data <- Data %>%
    dplyr::group_by(time) %>%
    dplyr::mutate(q_500 = ntile(w_iiv, 499)) %>%
    dplyr::group_by(time, q_500) %>%
    dplyr::mutate(q500_visits  = sum(visit),
                  q500_n = n(),
                  w_iiv_fsw_500 = (((q500_visits/q500_n)/(total_visits/total_n))^(-1))) %>%
    dplyr::ungroup() %>%
    as.data.frame() 
  
  
  
  # IDT weight: we assume a poisson distribution for prescribed dose
  
  DataIDT <- Data[!duplicated(Data$ID),]
  
  DataIDT$pred_doset1 <- predict( glm(doset1 ~ K1 + K2 + K3,data=DataIDT, 
                                      family=poisson), type='response')
  
  DataIDT$pred_doset1_stab <-predict( glm(doset1 ~ 1, data=DataIDT, 
                                          family=poisson), type='response')
  
  
  DataIDT$density <- dpois(DataIDT$doset1, lambda=DataIDT$pred_doset1)
  DataIDT$density_stab <- dpois(DataIDT$doset1, lambda=DataIDT$pred_doset1_stab)
  
  Datafinal <- merge(Data, DataIDT[,c('ID','density','density_stab')], by = "ID")
  
  
  Datafinal$w_idt <- Datafinal$density_stab/(Datafinal$density)
  
  # IDT with negative binomial: 
  
  glm_nb <- glm.nb(doset1 ~ K1 + K2 + K3, data = DataIDT)
  glm_nb_stab <- glm.nb(doset1 ~ 1, data = DataIDT)
  
  mu_nb <- predict(glm_nb, type = "response")
  mu_nb_stab <- predict(glm_nb_stab, type = "response")
  
  DataIDT$density_nb <- dnbinom(DataIDT$doset1, size = glm_nb$theta, mu = mu_nb)
  DataIDT$density_nb_stab <- dnbinom(DataIDT$doset1, size = glm_nb_stab$theta, mu = mu_nb_stab)
  
  DataIDT$w_idt_nb <- DataIDT$density_nb_stab/ DataIDT$density_nb
  
  Datafinal <- merge(Datafinal, DataIDT[,c('ID','w_idt_nb')], by = "ID")
  
  # OLS, IDT, IIV and IDEM estimators
  
  OLS <-lm(Y~ 1+ cumuldose, data=Datafinal)$coef[2] 
  OLS_adj <- lm(Y ~ 1 + cumuldose + K1 + K2 + K3, data = Datafinal)$coef[2]
  IDT <-lm(Y~ 1+ cumuldose, weights=w_idt, data=Datafinal)$coef[2]
  IIV <- lm(Y~ 1 + cumuldose, weights=1/w_iiv, data=Datafinal)$coef[2]
  IDEM <- lm(Y~ 1 + cumuldose, weight=w_idt*1/w_iiv, data=Datafinal)$coef[2]
  IDEM_trunc_97.5 <- lm(Y~ 1 + cumuldose, weights=w_idt/w_iiv_trunc_97.5, data=Datafinal)$coef[2]
  IDEM_trunc_99 <- lm(Y~ 1 + cumuldose, weights=w_idt/w_iiv_trunc_99, data=Datafinal)$coef[2]
  IDEM_fsw_20 <- lm(Y~ 1 + cumuldose, weights=w_idt/w_iiv_fsw_20, data=Datafinal)$coef[2]
  IDEM_fsw_500 <- lm(Y~ 1 + cumuldose, weights=w_idt/w_iiv_fsw_500, data=Datafinal)$coef[2]
  IDEM_nb <- lm(Y~ 1 + cumuldose, weights=w_idt_nb/w_iiv, data=Datafinal)$coef[2]
  
  ATE_Mat[S,] <- c(OLS,OLS_adj,IDT,IIV,IDEM, IDEM_trunc_97.5, IDEM_trunc_99,
                   IDEM_fsw_20, IDEM_fsw_500, IDEM_nb) 

  
  
}

# Results plot
ggplot(data = reshape2::melt(ATE_Mat, id.vars = NULL), aes(x=variable, y=value)) +
  geom_boxplot(aes(fill=variable)) + labs(x = "estimators", y = "ATE estimates") + 
  geom_hline(yintercept =1.4, linetype = "dashed", color = "blue") +
  scale_x_discrete(labels = scales::label_wrap(8)) +
  scale_fill_manual(values=c("#CCCCCC","azure4" , "#FF9999", "#66CC99","mediumpurple1","yellow2",
                             "darkorange2", "#99CCFF", "blue3", "red3"))+ ylim(c(1.1, 1.8))





