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

# the exponential approximation for mixture of 2 exponentials claims with \fbeta1=3.5900e-10, beta2=7.5088e-09,
# alpha=0.0584 and theta=0.3 (u in USD)
psi <- exp(-1 - (2 * m * theta * u - m2)/sqrt(m2^2 + 4 * theta * m * m3/3)) 
