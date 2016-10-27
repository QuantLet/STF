
[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="888" alt="Visit QuantNet">](http://quantlet.de/)

## [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **STFruin02** [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/)

```yaml

Name of QuantLet : STFruin02

Published in : Statistical Tools for Finance and Insurance

Description : Produces the exact ruin probability in infinite time for exponential claims.

Keywords : heavy-tailed, pdf, probability, poisson process, exponential, empirical

See also : 'STFruin04, STFruin06, STFruin07, STFruin08, STFruin09, STFruin10, STFruin12, STFruin13,
STFruin14, STFruin17'

Author : Zografia Anastasiadou

Submitted : Thu, July 28 2011 by Dedy Dwi Prastyo

```


### R Code:
```r
rm(list = ls(all = TRUE))
# setwd('C:/...')

# produces the exact ruin probability in infinite time for insurance collective risk model with exponential claims
ruinexp <- function(u, theta, beta) {
    # u: initial capital for risk process theta: security loading in insurance collective risk model beta: parameter for
    # exponential loss distribution
    y <- (1/(1 + theta)) * exp(-(theta * beta * u/(1 + theta)))  # ruin probability
    return(y)
}

u1 <- c(0, 10^9, 2 * 10^9, 3 * 10^9, 4 * 10^9, 5 * 10^9)  # initial capital of insurance company (in USD)
theta1 <- 0.3  # relative safety loading                      
beta1 <- 6.378e-09  # parameter for exponential loss distribution

# ruin probability in infinite time for exponential claims with \fbeta=6.3780e-09 and theta=0.3 (u in USD)
psi <- ruinexp(u1, theta1, beta1)
psi 

```
