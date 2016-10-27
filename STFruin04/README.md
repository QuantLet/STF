
[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="880" alt="Visit QuantNet">](http://quantlet.de/index.php?p=info)

## [<img src="https://github.com/QuantLet/Styleguide-and-Validation-procedure/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **STFruin04** [<img src="https://github.com/QuantLet/Styleguide-and-Validation-procedure/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/d3/ia)

```yaml

Name of QuantLet : STFruin04

Published in : Statistical Tools for Finance and Insurance

Description : 'Produces the exact ruin probability in infinite time for mixture of 2 exponentials
claims. Needs "ruinmix2exps.m" function.'

Keywords : 'heavy-tailed, pdf, probability, poisson process, approximation, simulation,
exponential, empirical'

See also : 'STFruin02, STFruin06, STFruin07, STFruin08, STFruin09, STFruin10, STFruin12, STFruin13,
STFruin14, STFruin17, ruinmix2exps'

Author : Zografia Anastasiadou

Submitted : Fri, May 04 2012 by Dedy Dwi Prastyo

Example : psi = ruinmix2exps(u1, theta1, dparameters1)

```


### R Code:
```r
rm(list = ls(all = TRUE))
# setwd('C:/...')

# produces the exact ruin probability in infinite time for insurance collective risk model with mixture of 2 exponentials
# distribution claims
ruinmix2exps <- function(u, theta, dparameters) {
    # u: initial capital for risk process theta: security loading in insurance collective risk model dparameters: list, composed
    # of 2 vectors containing the parameters of loss distribution, exponential parameters (first vector) and weights (second
    # vector)
    p1 <- dparameters[[1]]  # exponential parameters
    p2 <- dparameters[[2]]  # weights
    p <- p2[1]/p1[1]/(p2[1]/p1[1] + (1 - p2[1])/p1[2])
    psii <- p1[1] * (1 - p) + p1[2] * p
    r1 <- (psii + theta * sum(p1) - sqrt((psii + theta * sum(p1))^2 - 4 * prod(p1) * theta * (1 + theta)))/(2 * (1 + theta))
    r2 <- (psii + theta * sum(p1) + sqrt((psii + theta * sum(p1))^2 - 4 * prod(p1) * theta * (1 + theta)))/(2 * (1 + theta))
    y <- 1/((1 + theta) * (r2 - r1)) * ((psii - r1) * exp(-r1 * u) + (r2 - psii) * exp(-r2 * u))  # ruin probability using the Laplace transform inversion
    return(y)
}

u1 <- c(0, 10^9, 5 * 10^9, 10^10, 2 * 10^10, 5 * 10^10)  # initial capital of insurance company (in USD)
theta1 <- 0.3  # relative safety loading
dparameters1 <- list(c(3.59e-10, 7.5088e-09), c(0.0584, 0.9416))  # exponential parameters (first vector) and weights (second vector)

# ruin probability for mixture of 2 exponentials claims with exponetial parameters \fbeta1=3.5900e-10, beta2=7.5088e-09,
# alpha=0.0584 and theta=0.3 (u in USD)
psi <- ruinmix2exps(u1, theta1, dparameters1) 

```

### MATLAB Code:
```matlab
clear all;
close all,
clc;
format long;

% produces the exact ruin probability in infinite time for insurance
% collective risk model with mixture of 2 exponentials distribution claims

u1=[0;10^9;5*10^9;10^10;2*10^10;5*10^10]; % initial capital of insurance company (in USD)
theta1=0.3;                               % relative safety loading
dparameters1=[3.5900e-10 , 0.0584 ; 7.5088e-09 , 0.9416]; % exponential parameters (first column) and weights (second column)

% ruin probability for mixture of 2 exponentials claims with exponetial parameters 
%beta1=3.5900e-10, beta2=7.5088e-09, alpha=0.0584 and theta=0.3 (u in USD)

psi=ruinmix2exps(u1,theta1,dparameters1)

```
