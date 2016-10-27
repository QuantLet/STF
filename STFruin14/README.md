
[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="888" alt="Visit QuantNet">](http://quantlet.de/)

## [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **STFruin14** [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/)

```yaml

Name of QuantLet : STFruin14

Published in : Statistical Tools for Finance and Insurance

Description : 'Produces the ruin probability in infinite time for mixture of 2 exponentials
distribution claims given by heavy-light traffic approximation. Needs the "moments.m" function.'

Keywords : 'heavy-tailed, pdf, probability, poisson process, approximation, simulation,
exponential, empirical'

See also : 'STFruin02, STFruin03, STFruin04, STFruin05, STFruin06, STFruin07, STFruin08, STFruin09,
STFruin10, STFruin11, STFruin12, STFruin13, STFruin17, moments'

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
    p1 <- dparameters[[1]]
    p2 <- dparameters[[2]]
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

psi1 <- exp((-2 * theta * m * u)/m2)  # the heavy traffic approximation

u1 <- (1 - 1/(1 + theta)) * u  # capital for the light traffic approximation

pa <- matrix(dparameters1[[1]])  # 2 x 1 matrix of weights
pb <- matrix(dparameters1[[2]])  # 2 x 1 matrix of exponential parameters

paa <- matrix(pa, nrow = 2, ncol = length(u1))
pbb <- matrix(pb, nrow = 2, ncol = length(u1))

d <- rowSums(t(paa/pbb * exp(-pb %*% u1)))
c <- matrix(m, nrow = length(d)) - d

psi2 <- (m - c)/(1 + theta)/m  # the light traffic approximation

# the heavy-light traffic approximation for mixture of 2 exponentials claims with \fbeta1=3.5900e-10, beta2=7.5088e-09,
# alpha=0.0584 and theta=0.3 (u in USD)
psi <- (1/(1 + theta)^2) * psi1 + (1 - 1/(1 + theta)) * psi2 

```

### MATLAB Code:
```matlab
close all
clear all
clc

format long;
u = [0,10.^9,5.*10.^9,10.^10,2.*10.^10,5.*10.^10]'; % initial capital of insurance company (in USD)
theta = 0.3;                                % relative safety loading
dparameters1 = [[0.0584,0.9416]',[3.5900e-10,7.5088e-09]']; % weights (first vector) and exponential parameters (second vector)
m = moments(1,dparameters1);    % 1st raw moment
m2 = moments(2,dparameters1);   % 2nd raw moment
psi1 = exp((-2*theta*m*u)/m2);  % the heavy traffic approximation
u1 = (1-1/(1+theta))*u;         % capital for the light traffic approximation
pa = dparameters1(:,1); % 2 x 1 matrix of weights
pb = dparameters1(:,2); %2 x 1 matrix of exponential parameters
paa = repmat(pa,1,6);
pbb = repmat(pb,1,6);

n = paa./pbb.*exp(-pb*u1');
d = sum(n);
c = meshgrid(m,d)-d';
psi2 = (m-c)/(1+theta)/m       % the light traffic approximation
%the heavy-light traffic approximation for mixture of 2 exponentials claims
%with beta1=3.5900e-10, beta2=7.5088e-09, alpha=0.0584 and theta=0.3 (u in USD)
psi = (1/(1+theta)^2)*psi1 + (1-1/(1+theta))*psi2

```
