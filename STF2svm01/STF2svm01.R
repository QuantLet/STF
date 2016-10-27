graphics.off()
rm(list = ls(all = TRUE))

# load required packages install.packages('kernlab') install.packages('quadprog') install.packages('zoo')
# install.packages('tseries')
library(kernlab)
library(quadprog)
library(zoo)
library(tseries)

# setwd('C:/Users/...') #set your working direktory

x1 = 28
x2 = 7
y = 3
C = c(1, 1, 1, 200)
sig = c(100, 2, 0.5, 2)

G = read.matrix("Training100by100noNA.txt", header = TRUE, sep = "")

stf2svm01 = function(G, x1, x2, y, C, sig) {
    # g = data.frame(x2=G[,x2], x1=G[,x1],y=G[,y])
    
    # Calculate the score values for each company
    scores <- ksvm(G[, c(x2, x1)], G[, y], type = "C-svc", kernel = "rbfdot", kpar = list(sigma = 1/sig), C = C)
    
    # Make a 2-Dim plot
    
    plot(scores, data = G[, c(x2, x1, y)])
}


for (i in 1:4) {
    dev.new()
    stf2svm01(G = G, x1 = x1, x2 = x2, y = y, C = C[i], sig = sig[i])
} 
