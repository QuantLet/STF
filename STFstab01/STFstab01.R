graphics.off()
rm(list = ls(all = TRUE))
# setwd('C:/...')

# install.packages('fBasics') install.packages('stabledist')
library(stabledist)
library(fBasics)

x1 <- c(-30:30)/5
x2 <- c(-50:50)/5

alpha1 <- c(2, 1.8, 1.5, 1)

# symmetric stable pdfs
w1 <- dstable(x1, alpha1[1], beta = 0)
w2 <- dstable(x2, alpha1[2], beta = 0)
w3 <- dstable(x2, alpha1[3], beta = 0)
w4 <- dstable(x2, alpha1[4], beta = 0)

plot(x2, log(w2), type = "l", main = "Dependence on alpha", xlab = "x", ylab = "log(PDF(x))", cex.axis = 2, cex.lab = 1.4, cex.main = 2, 
    lwd = 3, col = "red", lty = 3)
lines(x1, log(w1), lwd = 3)
lines(x2, log(w3), col = "blue", lwd = 3, lty = 2)
lines(x2, log(w4), col = "green", lwd = 3, lty = 5)

########################################### 

x3 <- c(5:60)/10
x4 <- c(1:200)/20 + 60/10

x <- c(x3, x4)

alpha2 <- c(2, 1.95, 1.8, 1.5)

# symmetric stable cdfs
w5 <- pstable(x3, alpha2[1], beta = 0)
w6 <- pstable(x, alpha2[2], beta = 0)
w7 <- pstable(x, alpha2[3], beta = 0)
w8 <- pstable(x, alpha2[4], beta = 0)

s1 <- cbind(x3, 1 - w5)
s2 <- cbind(x, 1 - w6)
s3 <- cbind(x, 1 - w7)
s4 <- cbind(x, 1 - w8)
dev.new()
plot(log(s2), type = "l", ylim = c(-11.5, -0.5), main = "Tails of stable laws", xlab = "log(x)", ylab = "log(1-CDF(x))", cex.axis = 2, 
    cex.lab = 1.4, cex.main = 2, lwd = 3, col = "red", lty = 3)
lines(log(s1), lwd = 3)
lines(log(s3), col = "blue", lwd = 3, lty = 2)
lines(log(s4), col = "green", lwd = 3, lty = 5) 
