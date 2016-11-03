# clear variables and close windows
rm(list = ls(all = TRUE))
graphics.off()

Burrpdf = function(x, alpha, lambda, tau) {
    # BURRPDF Burr probability density function (pdf).  Y = BURRPDF(X,ALPHA,LAMBDA,TAU) returns the pdf of the Burr distribution
    # with parameters ALPHA,LAMBDA,TAU, evaluated at the values in X.  For CONTROL=0 the error message is displayed, if the
    # parmeters are negative.  The default values for the parameters ALPHA, LAMBDA, TAU, CONTROL are 1, 1, 2, 0, respectively.
    if (missing(tau) == TRUE) {
        tau = 2
    }
    if (missing(lambda) == TRUE) {
        lambda = 1
    }
    if (missing(alpha) == TRUE) {
        alpha = 1
    }
    if (missing(x) == TRUE) {
        stop("stats:normpdf:TooFewInputs! Input argument X is undefined.")
    }
    if (lambda <= 0) {
        stop("Non-positive lambda!")
    }
    if (alpha <= 0) {
        stop("Non-positive alpha!")
    }
    
    x = cbind(x)
    y = matrix(0, dim(x)[1], dim(x)[2])
    pos = x > 0
    y[pos] = tau * alpha * lambda^alpha * x[pos]^(tau - 1) * (lambda + x[pos]^tau)^(-alpha - 1)
}

step = 20

# Burr densities
x = (1:(144 * step))/step
y1 = Burrpdf(x, 0.5, 2, 1.5)
y2 = Burrpdf(x, 0.5, 0.5, 5)
y3 = Burrpdf(x, 2, 1, 0.5)


# Burr linear plot
plot(x, y1, col = "black", type = "l", lwd = 2, xlab = "x", ylab = "PDF(x)", ylim = c(0, 1.2), xlim = c(0, 8))
title("Burr densities")
lines(x, y2, col = "red3", lty = 3, lwd = 2)
lines(x, y3, col = "blue3", lty = 2, lwd = 2)

# Burr double-logarithmic densities
dev.new()
plot(x, y1, log = "xy", col = "black", type = "l", lwd = 2, xlab = "x", ylab = "PDF(x)", ylim = c(1e-04, 1), xlim = c((0.1), 
    (100)))
par(new = T)
plot(x, y2, type = "l", log = "xy", axes = F, frame = F, ylab = "", xlab = "", col = "red3", lty = 3, lwd = 2, ylim = c(1e-04, 
    1), xlim = c((0.1), (100)))
par(new = T)
plot(x, y3, type = "l", log = "xy", axes = F, frame = F, ylab = "", xlab = "", col = "blue3", lty = 2, , lwd = 2, ylim = c(1e-04, 
    1), xlim = c((0.1), (100)))
title("Burr densities")

# Weibull densities
y1 = dweibull(x, shape = 0.5, scale = (1^(-1/0.5)))
y2 = dweibull(x, shape = 2, scale = (1^(-1/2)))
y3 = dweibull(x, shape = 6, scale = (0.01^(-1/6)))

# Weibull linear plot
dev.new()
plot(x, y1, col = "black", type = "l", lwd = 2, xlab = "x", ylab = "PDF(x)", ylim = c(0, 1.2), xlim = c(0, 5))
title("Weibull densities")
lines(x, y2, col = "red3", lty = 3, lwd = 2)
lines(x, y3, col = "blue3", lty = 2, lwd = 2)

# Weibull semi-logarithmic plot
dev.new()
plot(x, y1, log = "y", col = "black", type = "l", lwd = 2, xlab = "x", ylab = "PDF(x)", ylim = c(0.001, 1), xlim = c(0, 5))
par(new = T)
plot(x, y2, type = "l", log = "y", axes = F, frame = F, ylab = "", xlab = "", col = "red3", lty = 3, lwd = 2, ylim = c(0.001, 
    1), xlim = c(0, 5))
par(new = T)
plot(x, y3, type = "l", log = "y", axes = F, frame = F, ylab = "", xlab = "", col = "blue3", lty = 2, , lwd = 2, ylim = c(0.001, 
    1), xlim = c(0, 5))
title("Weibull densities") 
