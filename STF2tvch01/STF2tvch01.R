# Clear variables and close windows
rm(list = ls(all = TRUE))
graphics.off()

# Please change working directory setwd('C:/...')

data <- read.delim2("SP1997-2005s.txt")

time <- (1:length(data[, 1]))
dat0 <- data[, 1] - c(mean(data[, 1]))
dat0 <- dat0/sd(dat0)

timet <- (time - 1078)/250 + 2001
plot(timet[time >= 1075], dat0[time >= 1075], xaxp = c(2001, 2005, 4), xlab = "Time", ylab = "Log-returns", type = "l") 
