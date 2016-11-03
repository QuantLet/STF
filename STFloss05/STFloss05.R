
# clear variables and close windows
rm(list = ls(all = TRUE))
graphics.off()


Paretopdf = function(x, alpha, lambda) {
    # PARETOPDF Pareto probability density function (pdf).  Y = PARETOPDF(X,ALPHA,LAMBDA) returns the pdf of the Pareto
    # distribution with parameters ALPHA, LAMBDA, evaluated at the values in X.  For CONTROL=0 the error message is displayed,
    # if the parmeters are negative.  The default values for ALPHA, LAMBDA are 1, 1 respectively.
    
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
    y[pos] = alpha * lambda^alpha/(lambda + x[pos])^(alpha + 1)
    return(y)
}

# Pareto densities
step = 10

x = (1:(144 * step))/step
y1 = Paretopdf(x, 0.5, 2)
y2 = Paretopdf(x, 2, 0.5)
y3 = Paretopdf(x, 2, 1)

# linear plot
plot(x, y1, col = "black", type = "l", lwd = 2, xlab = "x", ylab = "PDF(x)", ylim = c(0, 1), xlim = c(0, 5))
title("Pareto densities")
lines(x, y2, col = "red3", lty = 3, lwd = 2)
lines(x, y3, col = "blue3", lty = 2, lwd = 2)

dev.new()
plot(x, y1, log = "xy", col = "black", type = "l", lwd = 2, xlab = "x", ylab = "PDF(x)", ylim = c(1e-04, 1), xlim = c((0.1), 
    (100)))
par(new = T)
plot(x, y2, type = "l", log = "xy", axes = F, frame = F, ylab = "", xlab = "", col = "red3", lty = 3, lwd = 2, ylim = c(1e-04, 
    1), xlim = c((0.1), (100)))
par(new = T)
plot(x, y3, type = "l", log = "xy", axes = F, frame = F, ylab = "", xlab = "", col = "blue3", lty = 2, , lwd = 2, ylim = c(1e-04, 
    1), xlim = c((0.1), (100)))
title("Pareto densities") 
