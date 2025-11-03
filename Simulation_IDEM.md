Simulation Study for IDEM estimator
================

### Simulation Setup

Required R packages

``` r
library(survival)
library(rootSolve)
library(geex)
```

Sample size ($n$), number of simulations and maximum follow-up time
($\tau$)

``` r
n_sample <-  100
n_simul<- 100
tau <- 200
```

Prescribed dose $D_i$ model parameters
($\alpha_{0}, \alpha_{K_1}, \alpha_{K_2}, \alpha_{K_3}$)

``` r
alpha0= 120; alphaK1= 1.5; alphaK2= -2 ; alphaK3= -0.8 
```

Confounders ($K_1, K_2, K_3$) parameters

``` r
K1_mean=5;  K1_sd=1;  K2_prob=0.55; K3_mean=100; K3_sd=5
```

Mediator ($Z_1$ or $Z_2$) parameters

``` r
Z1_mean=4; Z1_sd=2; Z2_mean=2; Z2_sd=1
```

Response $Y_i(t)$ model parameters
($\beta_{0}, \beta_{K_1}, \beta_{K_2}$,$\beta_{K_3},\beta_{A(t)},\beta_{Q}$,$\beta_Z$,$\epsilon$,$\phi$)

``` r
beta0 <- 2; betaK1<-1; betaK2 <-0.25; betaK3<- -2; betaQ<- 0.; betaA<- 1.4; betaZ<- 3;
epsilon_sd <- 0.1; phi_mean =0; phi_sd =0.2
```

Visits indicator $dN_i(t)$ model parameters
($\gamma_0, \gamma_{A(t)}, \gamma_{Z}, \gamma_{K_1}$,$\gamma_{K_2},\gamma_{K_3},\gamma_{Q}$)

``` r
gamma0 <- -3.5; gammaA<-  0.05;  gammaZ<- 0.5; gammaK1<- 0.1 ; gammaK2<- -0.5;
gammaK3<- -0.05 ; gammaQ<- 2 
```
