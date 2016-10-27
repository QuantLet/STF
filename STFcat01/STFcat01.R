rm(list = ls(all = TRUE))
# setwd('C:/...')

d <- read.table("ncl.dat")

plot(1990 + d[, 2], d[, 3]/1e+09, type = "l", col = "blue", xlab = "Years", ylab = "Adjusted PCS catastrophe claims (USD billion)", 
    lwd = 2, cex.lab = 1.4, cex.axis = 1.4) 
