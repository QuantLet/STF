
# clear variables and close windows
rm(list = ls(all = TRUE))
graphics.off()

#################################### Subroutine getQnumber(x) #####

getQnumber = function(x) {
    # GETQNUMBER Auxiliary function for STF2loss9 y=getQnumber(x) returns the quarterly number of losses.
    
    maxx = max(x)
    l = length(x)
    y = matrix(0, maxx, 1)
    i = 1
    q = 1
    
    while (i <= l) {
        if (x[i] == q) {
            y[q] = y[q] + 1
            i = i + 1
        } else {
            q = q + 1
        }
    }
    return(y)
}

############################ main calculation #####


m = read.table("dfl.txt")
v = m[which(m[, 5] != 0), 5:6]

n = getQnumber(v[, 2])

# Qnumbers
t = (1:length(n))/4

plot(t, n, type = "l", col = "blue3", lwd = 2, xlab = "Time (years)", ylab = "Number of events", xlim = c(0, 23))

# Mean value function

ncum = cumsum(n)
dev.new()
plot(t, ncum, type = "l", lty = 2, col = "blue3", lwd = 2, xlab = "Time (years)", ylab = "Aggregate number of losses / Mean value function", 
    xlim = c(0, 23))
lines(t, (100.2394 * t), col = "black", lwd = 2)
lines(t, 17.9937 * t + 3.5759 * t^2, col = "red3", lty = 3, lwd = 2)

 
