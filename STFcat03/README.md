
[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="880" alt="Visit QuantNet">](http://quantlet.de/index.php?p=info)

## [<img src="https://github.com/QuantLet/Styleguide-and-Validation-procedure/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **STFcat03** [<img src="https://github.com/QuantLet/Styleguide-and-Validation-procedure/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/d3/ia)

```yaml

Name of QuantLet : STFcat03

Published in : Statistical Tools for Finance and Insurance

Description : 'Produces the plot of the aggregate quarterly number of losses for the PCS data
together with different mean value functions.'

Keywords : calibration, plot, mean, currency, empirical, mean, CAT bond

See also : STFcat01, STFcat01, STFcat02

Author : Zografia Anastasiadou

Submitted : Fri, May 04 2012 by Dedy Dwi Prastyo

Datafile : ncl.dat

Example : Plot of the aggregate number of losses / Mean value function over time.

```

![Picture1](plotR.png)

![Picture2](plotm.png)


### MATLAB Code:
```matlab
rm(list = ls(all = TRUE))
# setwd('C:/...')

m <- read.table("ncl.dat")

lambdaHPP <- 34.2
lambdaNHPP1 <- c(35.32, 2.32 * 2 * pi, -0.2)
lambdaNHPP2 <- c(35.22, 0.224, -0.16)

A <- 4
B <- 6
Delta <- 0.01

n <- table(m[, 5])
n <- as.data.frame(n)$Freq
t <- (1:length(n))/4
tn <- cbind(t, n)
ncum <- cumsum(n)
tncum <- cbind(t, ncum)

t <- seq(A, B, Delta)

HP <- lambdaHPP * t
NHP1 <- lambdaNHPP1[1] * t - lambdaNHPP1[2]/(2 * pi) * cos(2 * pi * (t + lambdaNHPP1[3])) + lambdaNHPP1[2]/(2 * pi) * cos(2 * 
    pi * lambdaNHPP1[3])
NHP2 <- lambdaNHPP2[1] * t + lambdaNHPP2[2] * (0.5 * t - 1/(8 * pi) * sin(4 * pi * (t + lambdaNHPP2[3])) + 1/(8 * pi) * sin(4 * 
    pi * lambdaNHPP2[3]))

plot(tncum[16:24, ], type = "l", col = "blue", lwd = 3, xlim = c(4, 6), ylim = c(139, 219), xlab = "Years", ylab = "Aggregate number of losses / Mean value function", 
    cex.lab = 1.4, cex.axis = 1.4)
lines(t, HP, col = "red", lwd = 3, lty = 3)
lines(t, NHP1, col = "black", lwd = 3, lty = 2)
lines(t, NHP2, col = "green", lwd = 3, lty = 4) 

```

### R Code:
```r
NA
```
