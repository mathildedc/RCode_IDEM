# Estimating the Causal Effect of a Cumulative Exposure on an Outcome in Studies Prone to Confounding and Irregular Visits
**Authors**: Mathilde Dicaire-Cartier and Janie Coulombe

## Overview 
R code for the simulation study and real-data application of the preprint *Estimating the Causal Effect of a Cumulative Exposure on an Outcome in Studies Prone to Confounding and Irregular Visits*. [(OSF preprint)](https://osf.io/43ndz/overview?view_only=0e1440c147a24e95a081ae7c6dda8295)

**Simulation Study/** contains the four R scripts for the different treatment and visit scenarios covered.

**Application/** contains the main application script, as well as the two-step variance estimation scripts for IDT, IIV and IDEM estimators. 

## Repository structure
````
RCode_IDEM/
├── Simulation Study/
│   ├── sim_scenario1_normal.R
│   ├── sim_scenario1_poisson.R
│   ├── sim_scenario2_normal.R
│   └── sim_scenario2_poisson.R
├── Application/
│   ├── application.R
│   ├── varianceIDT.R
│   ├── varianceIIV.R
│   └── varianceIDEM.R
└── README.md
````
