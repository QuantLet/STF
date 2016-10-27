rm(list = ls(all = TRUE))
# setwd('C:/...')

# install.packages('MASS')
library(MASS)

# returns the k-th moment (up to fourth) of the mixture of 2 exponentials distribution claims
moments <- function(k, dparameters) {
    # k: order of moment to calculate dparameters: list, composed of 2 vectors containing the parameters of the loss
    # distribution, weights (first vector) and exponential parameters (second vector)
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

u1 <- c(0, 10^9, 5 * 10^9, 10^10, 2 * 10^10, 5 * 10^10)  # initial capital of insurance company (in USD)
theta1 <- 0.3  # relative safety loading
dparameters1 <- list(c(0.0584, 0.9416), c(3.59e-10, 7.5088e-09))  # weights (first vector) and exponential parameters (second vector)

# moments for mixture of 2 exponentials distribution claims
m <- moments(1, dparameters1)
m2 <- moments(2, dparameters1)
m3 <- moments(3, dparameters1)
m4 <- moments(4, dparameters1)

if (m2 * m4 < 3/2 * m3^2) {
    thetanew <- (theta1 * m * (2 * m3^2 - m2 * m4)/m2^2/m3)
    mnew <- ((3 * m3^2 - 2 * m2 * m4)/m2/m3)
    m2new <- ((m2 * m4 - 2 * m3^2) * (2 * m2 * m4 - 3 * m3^2)/m3^2/m2^2)
} else {
    thetanew <- 1/2 * theta1/m2^2 * m * (m3 + m2 * m)
    mnew <- m
    m2new <- 1/2/m2 * m * (m3 + m2 * m)
}

p1 <- mnew^2/(m2new - mnew^2)
p2 <- mnew/(m2new - mnew^2)
dparametersnew <- list(p1, p2)  # gamma parameters

# returns the k-th moment (up to fourth) of the mixture of 2 gamma distribution claims
momentsgam <- function(k, dparameters) {
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
mgfsg <- function(x, k, dparameters) {
    # x: scalar, argument of the moment generating function k: scalar, integer, 0 =< k <= 3, order of the derivative
    # dparameters: list of scalars, parameters of gamma distribution
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

# moments for gamma distribution claims
mg <- momentsgam(1, dparametersnew)  # 1st raw moment
mg2 <- momentsgam(2, dparametersnew)  # 2nd raw moment
mg3 <- momentsgam(3, dparametersnew)  # 3nd raw moment

# returns the adjustment coefficient R for gamma distribution claims
adjR <- function(theta, dparameters) {
    # theta: security loading in insurance collective risk model dparameters: list of scalars, parameters of gamma distribution
    p1 <- dparameters[[1]]
    p2 <- dparameters[[2]]
    R0 <- 0.99999999 * p2
    R0 <- c(12 * theta * mg/(3 * mg2 + sqrt(9 * mg2^2 + 24 * mg * mg3 * theta)), R0)
    R0 <- min(R0)
    r <- R0
    err <- 1
    while (err > 1e-09) {
        D1 <- 1 + (1 + theta) * mg * r - mgfsg(r, 0, dparameters)
        D2 <- (1 + theta) * mg - mgfsg(r, 1, dparameters)
        err <- r
        r <- r - D1/D2
        err <- abs(err - r)/r
        R <- r
    }
    return(R)  # adjustment coefficient R 
}

R <- adjR(thetanew, dparametersnew)  # adjustment coefficient R for gamma distribution claims
mgfprim <- mgfsg(R, 1, dparametersnew)  # moment generating function for gamma distribution claims 

C <- (thetanew * mg)/(mgfprim - mg * (1 + thetanew))
Cram <- C * exp(-R * u1)  # the Cramer-Lundberg approximation for gamma claims 

u1 <- u1 * p2/p1
b <- 1/p1
n <- length(u1)

# the function to be integrated
exactgamint <- function(x) {
    j <- 1
    while (j < n + 1) {
        uj <- u1[j]
        L <- x^(1/b) * exp(-(x + 1) * uj/b)
        M <- (x^(1/b) * (1 + (1 + thetanew) * (x + 1)/b) - cos(pi/b))^2 + sin(pi/b)^2
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

# the 4-moment gamma De Vylder approximation for mixture of 2 exponentials claims with \fbeta1=3.5900e-10,
# beta2=7.5088e-09, alpha=0.0584 and theta=0.3 (u in USD)
psi <- Cram + as.vector((thetanew * sin(pi/b)/pi/b) * int)
psi 
