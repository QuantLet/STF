# clear variables and close windows
rm(list = ls(all = TRUE))
graphics.off()

# gamma densities
step = 10

x = (1:(8 * step))/step
y1 = dgamma(x, shape = 1, scale = 2)
y2 = dgamma(x, shape = 2, scale = 1)
y3 = dgamma(x, shape = 3, scale = 0.5)

# linear plot
plot(x, y1, col = "black", type = "l", lwd = 2, xlab = "x", ylab = "PDF(x)", ylim = c(-0.01, 0.6))
title("Gamma densities")
lines(x, y2, col = "red3", lty = 3, lwd = 2)
lines(x, y3, col = "blue3", lty = 2, lwd = 2)

# semi-logarithmic plot
dev.new()
plot(x, y1, log = "y", col = "black", type = "l", lwd = 2, xlab = "x", ylab = "PDF(x)", ylim = c(0.001, 1))
par(new = T)
plot(x, y2, type = "l", log = "y", axes = F, frame = F, ylab = "", xlab = "", col = "red3", lty = 3, lwd = 2, ylim = c(0.001, 
    1))
par(new = T)
plot(x, y3, log = "y", type = "l", axes = F, frame = F, ylab = "", xlab = "", col = "blue3", lty = 2, lwd = 2, ylim = c(0.001, 
    1))
title("Gamma densities")

# log-normal densities
x = (1:(25 * step))/step

y1 = dlnorm(x, 2, 1)
y2 = dlnorm(x, 2, 0.1)
y3 = dlnorm(x, 0.5, 2)

# linear plot
dev.new()
plot(x, y1, col = "black", type = "l", lwd = 2, xlab = "x", ylab = "PDF(x)", ylim = c(-0.01, 0.6))
title("Log-normal densities")
lines(x, y2, col = "red3", lty = 3, lwd = 2)
lines(x, y3, col = "blue3", lty = 2, lwd = 2)

# semi-logarithmic plot

dev.new()
plot(x, y1, log = "y", col = "black", type = "l", lwd = 2, xlab = "x", ylab = "PDF(x)", ylim = c(0.001, 1))
par(new = T)
plot(x, y2, type = "l", log = "y", axes = F, frame = F, ylab = "", xlab = "", col = "red3", lty = 3, lwd = 2, ylim = c(0.001, 
    1))
par(new = T)
plot(x, y3, log = "y", type = "l", axes = F, frame = F, ylab = "", xlab = "", col = "blue3", lty = 2, lwd = 2, ylim = c(0.001, 
    1))
title("Log-normal densities") 
