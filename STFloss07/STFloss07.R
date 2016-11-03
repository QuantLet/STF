
# clear variables and close windows
rm(list = ls(all = TRUE))
graphics.off()


# install packages for gamma and beta functions
install.packages("zipfR")
library(zipfR)

install.packages("fAsianOptions")
library(fAsianOptions)

################################## subroutine betacore(x,a,b) ###

betacore = function(x, a, b) {
    # BETACORE Core algorithm for the incomplete beta function.  BETACORE(x,a,b) is used twice by BETAINC(X,A,B).
    
    x = as.matrix(x)
    
    y = x
    qab = a + b
    qap = a + 1
    qam = a - 1
    am = matrix(1, nrow(x), ncol(x))
    bm = am
    y = am
    bz = 1 - qab * x/qap
    d = matrix(0, nrow(x), ncol(x))
    app = d
    ap = d
    bpp = d
    bp = d
    yold = d
    m = 1
    eps = 2^(-52)
    while (any(abs(y - yold) > 10 * eps * abs(y))) {
        tem = 2 * m
        d = m * (b - m) * x/((qam + tem) * (a + tem))
        ap = y + d * am
        bp = bz + d * bm
        d = -(a + m) * (qab + m) * x/((a + tem) * (qap + tem))
        app = ap + d * y
        bpp = bp + d * bz
        yold = y
        am = ap/bpp
        bm = bp/bpp
        y = app/bpp
        if (m == 1) {
            bz = 1  # only need to do this first time through
        }
        m = m + 1
    }
    return(y)
}

################################# subroutine betainc(x,a,b) ###

betainc = function(x, a, b) {
    # BETAINC Incomplete beta function.  The elements of X must be in the interval [0,1].  The arguments X, A and B must all be
    # the same size except that scalar arguments function as constant matrices of the common size of the other arguments.
    
    x = as.vector(x)
    a = c(a)
    b = c(b)
    
    m = max(nrow(x), max(length(a), length(b)))
    n = max(ncol(x), max(1, 1))
    
    y = matrix(0, m, n)
    bt = exp(lgamma(a + b) - lgamma(a) - lgamma(b) + a * c(log(x)) + b * c(log(1 - x)))
    
    k = which(x < (c(a) + 1)/(c(a) + c(b) + 2))
    if (length(k) > 0) {
        if (length(x) == 1) {
            xk = x
        } else {
            xk = x[k]
        }
        if (length(a) == 1) {
            ak = a
        } else {
            ak = a[k]
        }
        if (length(b) == 1) {
            bk = b
        } else {
            bk = b[k]
        }
        y[k] = bt[k] * betacore(xk, ak, bk)/ak
    }
    
    k = which(x >= (a + 1)/(a + b + 2))
    if (length(k) > 0) {
        if (length(x) == 1) {
            xk = x
        } else {
            xk = x[k]
        }
        if (length(a) == 1) {
            ak = a
        } else {
            ak = a[k]
        }
        if (length(b) == 1) {
            bk = b
        } else {
            bk = b[k]
        }
        y[k] = 1 - bt[k] * betacore(1 - xk, bk, ak)/bk
    }
    return(y)
}


################################################## subroutine mef(param1,param2,param3,x,ind) ###

mef = function(param1, param2, param3, x, ind) {
    # MEF Mean excess function.  y = mef (param1, param2, param3, x, ind) returns the values of the mean excess function for the
    # distribution specified by IND: 1-lognormal, 2-gamma, 3-weibull, 4-pareto, 5-Burr, 6-mixture of exponential distributions
    # with parameters PARAM1, PARAM2, PARAM3.
    
    if (ind == 1) {
        # Lognormal
        mu = param1
        sigma = param2
        FP = exp(mu + sigma^2/2) * (1 - pnorm((log(x) - mu - sigma^2)/sigma))
        FI = pnorm((log(x) - mu)/sigma)
        y = FP/(1 - FI) - x
        return(y)
    } else if (ind == 2) {
        # Gamma
        alpha = param1
        beta = param2
        FP = alpha/beta * (1 - pgamma(x, shape = alpha + 1, scale = 1/beta))
        FI = pgamma(x, shape = alpha, scale = 1/beta)
        y = FP/(1 - FI) - x
        return(y)
    } else if (ind == 3) {
        # Weibull, c.d.f = 1-exp(-t^alpha * beta)
        alpha = param1
        beta = param2
        # y = beta^(-1/alpha)* gamma(1+1/alpha)*(1-Igamma(1+1/alpha,x^alpha*beta,lower=F))*exp(beta*x^alpha)-x
        y = beta^(-1/alpha) * gamma(1 + 1/alpha) * Igamma(1 + 1/alpha, x^alpha * beta, lower = F) * exp(beta * x^alpha) - x
        return(y)
    } else if (ind == 4) {
        # Pareto
        alpha = param1
        lambda = param2
        y = (x + lambda)/(alpha - 1)
        return(y)
    } else if (ind == 5) {
        # Burr
        alpha = param1
        lambda = param2
        tau = param3
        FP1 = lambda^(1/tau) * gamma(alpha - 1/tau) * gamma(1 + 1/tau)/gamma(alpha) * (lambda/(lambda + x^tau))^(-alpha)
        # ibeta = function(x,a,b){ pbeta(x,a,b)*beta(a,b) }
        FP2 = 1 - betainc(x^tau/(lambda + x^tau), 1 + 1/tau, alpha - 1/tau)
        y = FP1 * FP2 - x
        return(y)
    } else if (ind == 6) {
        # mixture of exponentials
        beta1 = param1
        beta2 = param2
        ce = param3
        FP1 = ce/beta1 * exp(-beta1 * x) + (1 - ce)/beta2 * exp(-beta2 * x)
        FP2 = ce * exp(-beta1 * x) + (1 - ce) * exp(-beta2 * x)
        y = FP1/FP2
        return(y)
    }
}

######################## main calculation ###

xaxis = (1:2000)/100

par(mfrow = c(1, 2))

# log
p1 = mef(0.001, 1, -1, xaxis, 1)
# gamma 1
p2 = mef(0.5, 1/1.6351, -1, xaxis, 2)
# gamma 2
p3 = mef(1.5, 1/1.6351, -1, xaxis, 2)
# mix
p4 = mef(1/4, 1/0.8, 0.1, xaxis, 6)

plot(xaxis, p1, col = "black", type = "l", lwd = 2, ylim = c(0, max(p1, p2, p3, p4)), ylab = "e(x)", xlab = "x")
lines(xaxis, p2, col = "blue3", lty = 2, lwd = 2)
lines(xaxis, p3, col = "green4", lty = 4, lwd = 2)
lines(xaxis, p4, col = "red3", lty = 3, lwd = 2)

# MEF for Pareto, Burr, 2 x Weibull

# Weibull
p1 = mef(1.1, 0.5, -1, xaxis, 3)
# Weibull
p2 = mef(0.9, 0.5, -1, xaxis, 3)
# Pareto
p3 = mef(5, 3.03, -1, xaxis, 4)
# Burr
p4 = mef(2.39, 3.03, 3, xaxis, 5)


plot(xaxis, p1, col = "black", type = "l", lwd = 2, ylab = "e(x)", xlab = "x", xlim = c(0, 20), ylim = c(0, 6))
lines(xaxis, p2, col = "blue3", lty = 2, lwd = 2)
lines(xaxis, p3, col = "green4", lty = 4, lwd = 2)
lines(xaxis, p4, col = "red3", lty = 3, lwd = 2)
 
