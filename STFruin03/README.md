
[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="880" alt="Visit QuantNet">](http://quantlet.de/index.php?p=info)

## [<img src="https://github.com/QuantLet/Styleguide-and-Validation-procedure/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **STFruin03** [<img src="https://github.com/QuantLet/Styleguide-and-Validation-procedure/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/d3/ia)

```yaml

Name of QuantLet : STFruin03

Published in : Statistical Tools for Finance and Insurance

Description : Produces the exact ruin probability in infinite time for gamma claims.

Keywords : monte-carlo, heavy-tailed, pdf, probability, poisson process, exponential, empirical

See also : STFruin07, STFruin08, STFruin09, STFruin10, STFruin12, STFruin13, STFruin14, STFruin17

Author : Zografia Anastasiadou

Submitted : Fri, September 16 2011 by Dedy Dwi Prastyo

```


### R Code:
```r
rm(list = ls(all = TRUE))
# setwd('C:/...')

# install.packages('MASS')
library(MASS)

# returns the k-th moment (up to fourth) of the mixture of 2 gamma distribution claims
moments <- function(k, dparameters) {
    # k: order of moment to calculate dparameters: list of scalars, parameters of gamma distribution
    p1 <- dparameters[[1]]
    p2 <- dparameters[[2]]
    if (k == 1) {
        mk <- p1/p2
    } else {
        if (k == 2) {
            mk <- (p1^2 + p1)/(p2^2)
        } else {
            if (k == 3) {
                mk <- (p1^3 + 3 * p1^2 + 2 * p1)/(p2^3)
            } else {
                if (k == 4) {
                  mk <- p1 * (p1 + 1) * (p1 + 2) * (p1 + 3)/(p2^4)
                }
            }
        }
    }
    return(mk)  # k-th raw moment of gamma distribution claims 
}

# returns the moment generating function or its k-th derivative (up to third) for gamma distribution claims
mgfs <- function(x, k, dparameters) {
    # x: scalar, n x 1 vector or m x n matrix, argument of the moment generating function k: scalar, integer, 0 =< k <= 3, order
    # of the derivative dparameters: list of scalars, parameters of gamma distribution
    p1 <- dparameters[[1]]
    p2 <- dparameters[[2]]
    if (k == 0) {
        y <- p2^p1/((p2 - x)^p1)
    } else {
        if (k == 1) {
            y <- (p1/p2) * (p2/(p2 - x))^(p1 + 1)
        } else {
            if (k == 2) {
                y <- p1 * (p1 + 1) * p2^p1/(p2 - x)^(p1 + 2)
            } else {
                if (k == 3) {
                  y <- p1 * (p1 + 1) * (p1 + 2) * p2^p1/(p2 - x)^(p1 + 3)
                }
            }
        }
    }
    return(y)
}

u1 <- c(0, 10^9, 2 * 10^9, 3 * 10^9, 4 * 10^9, 5 * 10^9)  # initial capital of insurance company (in USD)
theta1 <- 0.3  # relative safety loading
dparameters1 <- list(0.9185, 6.1662e-09)  # gamma parameters

m <- moments(1, dparameters1)  # 1st raw moment
m2 <- moments(2, dparameters1)  # 2nd raw moment
m3 <- moments(3, dparameters1)  # 3nd raw moment

# returns the adjustment coefficient R for gamma distribution claims
adjR <- function(theta, dparameters) {
    # theta: security loading in insurance collective risk model dparameters: list of scalars, parameters of gamma distribution
    p1 <- dparameters[[1]]
    p2 <- dparameters[[2]]
    R0 <- 0.99999999 * p2
    R0 <- c(12 * theta * m/(3 * m2 + sqrt(9 * m2^2 + 24 * m * m3 * theta)), R0)
    R0 <- min(R0)
    r <- R0
    err <- 1
    while (err > 1e-09) {
        D1 <- 1 + (1 + theta) * m * r - mgfs(r, 0, dparameters)
        D2 <- (1 + theta) * m - mgfs(r, 1, dparameters)
        err <- r
        r <- r - D1/D2
        err <- abs(err - r)/r
        R <- r
    }
    return(R)  # adjustment coefficient R 
}

R <- adjR(theta1, dparameters1)  # adjustment coefficient R for gamma distribution claims
mgfprim <- mgfs(R, 1, dparameters1)  # moment generating function for gamma distribution claims 

C <- (theta1 * m)/(mgfprim - m * (1 + theta1))
Cram <- C * exp(-R * u1)  # the Cramer-Lundberg approximation for gamma claims 

p1 <- dparameters1[[1]]
p2 <- dparameters1[[2]]
u1 <- u1 * p2/p1
b <- 1/p1
n <- length(u1)

# the function to be integrated
exactgamint <- function(x) {
    j <- 1
    while (j < n + 1) {
        uj <- u1[j]
        L <- x^(1/b) * exp(-(x + 1) * uj/b)
        M <- (x^(1/b) * (1 + (1 + theta1) * (x + 1)/b) - cos(pi/b))^2 + sin(pi/b)^2
        j <- j + 1
    }
    y <- L/M
    return(y)
}

# integrates exactgamint function using the Simpson's method
d <- area(exactgamint, 0, 0.001)
d <- d + area(exactgamint, 0.001, 1)
d <- rbind(1, d)

err <- 1e-05
int <- matrix(1, n)
j <- 1
while (j < n + 1) {
    i <- 2
    while (abs((d[i - 1] - d[i])/d[i]) > err) {
        v <- area(exactgamint, i - 1, i)
        d <- rbind(d, (v + d[i]))
        i <- i + 1
    }
    endd <- length(d)
    int[j] <- d[endd]
    j <- j + 1
}

# the ruin probability for gamma claims with alpha=0.9185, beta=6.1662e-09 and theta=0.3 (u in USD)
psi <- Cram + as.vector((theta1 * sin(pi/b)/pi/b) * int)
psi 

```
