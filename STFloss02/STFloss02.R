
# clear variables and close windows
rm(list = ls(all = TRUE))
graphics.off()

# Pseudo random numbers
set.seed(25)

r = rlnorm(10, 0.5, 0.5)

edf = function(x) {
    n = length(x)
    x = sort(x)
    a = (1:n - 1)/n
    aa = cbind(x[2:n], a[2:n])
    bb = cbind(x[1:n - 1], a[2:n])
    cc = rbind(aa, bb)
    
    edf = apply(aa, 2, sort)
    end = dim(edf)[1]
    edf = rbind(c(edf[1, 1], 0), edf, c(edf[end, 1], 1))
    return(edf)
}

y = edf(r)

x = y[, 1]
y = y[, 2]

xw = (10:500)/100
w = plnorm(xw, 0.5, 0.5)

plot(x, y, col = "black", type = "s", lty = 1, lwd = 1.5, xlab = "x", ylab = "CDF(x)")
title("Empirical distribution function")

dev.new()
plot(x, y, type = "l", col = "black", lwd = 1.5, xlab = "x", xlim = c(0, 5), ylab = "CDF(x)")
lines(xw, w, col = "red3", lty = 3, lwd = 1.5)
title("Empirical and lognormal distributions")
 
