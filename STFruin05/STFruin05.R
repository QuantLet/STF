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

# returns the moment generating function or its k-th derivative (up to third) for mixture of 2 exponentials distribution
# claims
mgfs <- function(x, k, dparameters) {
    # x: scalar, n x 1 vector or m x n matrix, argument of the moment generating function k: scalar, integer, 0 =< k <= 3, order
    # of the derivative dparameters: list, composed of 2 vectors containing the parameters of the loss distribution, weights
    # (first vector) and exponential parameters (second vector)
    p1 <- dparameters[[1]]  # weights
    p2 <- dparameters[[2]]  # exponential parameters
    if (k == 0) {
        y <- sum((p1 * p2)/(p2 - t(x)))
    } else {
        if (k == 1) {
            y <- sum((p1 * p2)/(p2 - t(x))^2)
        } else {
            if (k == 2) {
                y <- 2 * sum((p1 * p2)/(p2 - t(x))^3)
            } else {
                if (k == 3) {
                  y <- 6 * sum((p1 * p2)/(p2 - t(x))^4)
                }
            }
        }
    }
    return(y)
}

u <- c(0, 10^9, 5 * 10^9, 10^10, 2 * 10^10, 5 * 10^10)  # initial capital of insurance company (in USD)
theta1 <- 0.3  # relative safety loading
dparameters1 <- list(c(0.0584, 0.9416), c(3.59e-10, 7.5088e-09))  # weights (first vector) and exponential parameters (second vector)

m <- moments(1, dparameters1)  # 1st raw moment
m2 <- moments(2, dparameters1)  # 2nd raw moment
m3 <- moments(3, dparameters1)  # 3nd raw moment

# returns the adjustment coefficient R for mixture of 2 exponentials distribution claims
adjR <- function(theta, dparameters) {
    # theta: security loading in insurance collective risk model dparameters: list, composed of 2 vectors containing the
    # parameters of the loss distribution, weights (first vector) and exponential parameters (second vector)
    p1 <- dparameters[[1]]  # weights
    p2 <- dparameters[[2]]  # exponential parameters
    R0 <- min(p2)
    R0 <- c(12 * theta * m/(3 * m2 + sqrt(9 * m2^2 + 24 * m * m3 * theta)), R0)
    R0 <- min(R0)
    r <- R0
    err = 1
    while (err > 1e-09) {
        D1 <- 1 + (1 + theta) * m * r - mgfs(r, 0, dparameters1)
        D2 <- (1 + theta) * m - mgfs(r, 1, dparameters1)
        err <- r
        r <- r - D1/D2
        err <- abs(err - r)/r
        R <- r
    }
    return(R)  # adjustment coefficient R 
}

R <- adjR(theta1, dparameters1)  # adjustment coefficient R for mixture of 2 exponentials distribution claims
mgfprim <- mgfs(R, 1, dparameters1)  # moment generating function for mixture of 2 exponentials distribution claims 

C <- (theta1 * m)/(mgfprim - m * (1 + theta1))

# the Cramer-Lundberg approximation for mixture of 2 exponentials claims with \fbeta1=3.5900e-10, beta2=7.5088e-09,
# alpha=0.0584 and theta=0.3 (u in USD)
psi <- C * exp(-R * u)
psi 
