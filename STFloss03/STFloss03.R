# clear all variables and close windows
rm(list = ls(all = TRUE))
graphics.off()

mixexppdf = function(x, alpha, beta1, beta2) {
    # MIXEXPPDF Mixed exponential probability density function (pdf).  Y = MIEXEXPPDF(X,ALPHA,BETA1,BETA2) returns the pdf of
    # the mixed exponential distribution with mixing probability A and distributions parameters BETA1, BETA2, evaluated at the
    # values in X.  For CONTROL=0 the error message is displayed, if the parmeters are negative or a>1.  The default values for
    # A, BETA1, BETA2 are 0.5, 1, 2.
    if (missing(alpha) == T) 
        alpha = 0.5
    if (missing(beta1) == T) 
        beta1 = 1
    if (missing(beta2) == T) 
        beta2 = 2
    if (beta1 <= 0) {
        stop("Non-positive beta1! Please insert a positive beta1!")
    }
    if (beta2 <= 0) {
        stop("Non-positive beta2! Please insert a positive beta2!")
    }
    if (alpha <= 0) {
        stop("Alpha lesser or equal 0! Please insert alpha between 0 and 1!")
    }
    if (alpha >= 1) {
        stop("Alpha greater or equal 1! Please insert alpha between 0 and 1!")
    }
    
    y = matrix(0, dim(x)[1], dim(x)[2])
    pos = x > 0
    y[pos] = alpha * beta1 * exp(-beta1 * x[pos]) + (1 - alpha) * beta2 * exp(-beta2 * x[pos])
    return(cbind(y))
}

step = 10

x = cbind((1:(8 * step))/step)

y1 = dexp(x, 1/3)
y2 = dexp(x, 1)
y3 = mixexppdf(x, 0.5, 0.3, 1)

plot(x, y1, col = "black", type = "l", lwd = 2, xlab = "x", ylab = "PDF(x)", ylim = c(-0.01, 0.8))
title("Mixture of two exponential densities")
lines(x, y2, col = "red3", lty = 3, lwd = 2)
lines(x, y3, col = "blue3", lty = 2, lwd = 2)

dev.new()

plot(x, y1, log = "y", col = "black", type = "l", lwd = 2, xlab = "x", ylab = "PDF(x)", frame = T, ylim = c(0.01, 1))
par(new = T)
plot(x, y2, log = "y", axes = F, frame = F, col = "red3", type = "l", lty = 3, lwd = 2, ylab = "", xlab = "", ylim = c(0.01, 
    1))
par(new = T)
plot(x, y3, log = "y", axes = F, frame = F, col = "blue3", type = "l", lty = 2, lwd = 2, ylab = "", xlab = "", ylim = c(0.01, 
    1))
title("Mixture of two exponential semi-log densities") 
