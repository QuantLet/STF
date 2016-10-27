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
