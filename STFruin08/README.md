
[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="880" alt="Visit QuantNet">](http://quantlet.de/index.php?p=info)

## [<img src="https://github.com/QuantLet/Styleguide-and-Validation-procedure/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **STFruin08** [<img src="https://github.com/QuantLet/Styleguide-and-Validation-procedure/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/d3/ia)

```yaml

Name of QuantLet : STFruin08

Published in : Statistical Tools for Finance and Insurance

Description : 'Produces the ruin probability in infinite time for mixture of 2 exponentials
distribution claims given by Beekman-Bowers approximation. Needs the "moments.m" function.'

Keywords : 'heavy-tailed, pdf, probability, poisson process, approximation, simulation,
exponential, empirical'

See also : 'STFruin02, STFruin03, STFruin04, STFruin05, STFruin06, STFruin07, STFruin09, STFruin09,
STFruin10, STFruin12, STFruin13, STFruin14, STFruin17, moments'

Author : Zografia Anastasiadou

Submitted : Fri, May 04 2012 by Dedy Dwi Prastyo

```


### R Code:
```r
rm(list = ls(all = TRUE))
# setwd('C:/...')

# returns the k-th moment (up to fourth) of the mixture of 2 exponentials distribution claims
moments <- function(k, dparameters) {
    # k: order of moment to calculate dparameters: list, composed of 2 vectors containing the parameters of loss distribution,
    # weights (first vector) and exponential parameters (second vector)
    p1 <- dparameters[[1]]  # weights
    p2 <- dparameters[[2]]  # exponential parameters
    if (k == 1) {
        mk <- sum(p1/p2)
    } else {
        if (k == 2) {
            mk <- 2 * sum(p1/p2^2)
        } else {
            if (k == 3) {
                mk <- 6 * sum(p1/p2^3)
            } else {
                if (k == 4) {
                  mk <- 24 * sum(p1/p2^4)
                }
            }
        }
    }
    return(mk)  # k-th raw moment of the mixture of 2 exponentials distribution claims
}

u <- c(0, 10^9, 5 * 10^9, 10^10, 2 * 10^10, 5 * 10^10)  # initial capital of insurance company (in USD)
theta <- 0.3  # relative safety loading
dparameters1 <- list(c(0.0584, 0.9416), c(3.59e-10, 7.5088e-09))  # weights (first vector) and exponential parameters (second vector)

m <- moments(1, dparameters1)  # 1st raw moment
m2 <- moments(2, dparameters1)  # 2nd raw moment
m3 <- moments(3, dparameters1)  # 3nd raw moment


delta1 <- (2 * m * theta)/(m2 + ((4 * m * m3/3/m2) - m2) * theta)  # scale parameter of gamma distribution function
delta2 <- (1 + theta)/(1 + ((4 * m * m3/3/m2^2) - 1) * theta)  # shape parameter of gamma distribution function

# the Beekman-Bowers approximation for mixture of 2 exponentials claims with \fbeta1=3.5900e-10, beta2=7.5088e-09,
# alpha=0.0584 and theta=0.3 (u in USD)
psi <- (1 - pgamma(u, shape = delta2, scale = 1/delta1))/(1 + theta) 

```

### MATLAB Code:
```matlab
clear all;
close all,
clc;
format long;

u=[0;10^9;5*10^9;10^10;2*10^10;5*10^10]; % initial capital of insurance company (in USD)
theta=0.3;                               % relative safety loading
dparameters1=[0.0584 , 3.5900e-10 ; 0.9416 , 7.5088e-09]; % weights (first column) and exponential parameters (second column)

m=moments(1,dparameters1);  % 1st raw moment
m2=moments(2,dparameters1); % 2nd raw moment
m3=moments(3,dparameters1); % 3nd raw moment

delta1=(2*m*theta)/(m2+((4*m*m3/3/m2)-m2)*theta) % scale parameter of gamma distribution function
delta2=(1+theta)/(1+((4*m*m3/3/m2^2)-1)*theta)   % shape parameter of gamma distribution function

% the Beekman-Bowers approximation for mixture of 2 exponentials claims with beta1=3.5900e-10, beta2=7.5088e-09, alpha=0.0584 and theta=0.3 (u in USD) 
psi=(1-gamcdf(u,delta2,1/delta1))/(1+theta)
```
