# clear variables and close windows
rm(list = ls(all = TRUE))
graphics.off()

samplemef = function(data, xaxis) {
    # SAMPLEMEF Sample mean excess function.  RES = SAMPLEMEF (DATA, XAXIS) returns the value of the sample mean excess function
    # for the vector DATA in the points from vector XAXIS.
    
    
    sorteddata = sort(data)
    dataLength = length(sorteddata)
    resLength = length(xaxis)
    res = xaxis
    
    
    i = 0
    valuesOnTheLeft = 0
    smef = mean(sorteddata[1:dataLength])
    
    while (i < resLength) {
        i = i + 1
        while (sorteddata[valuesOnTheLeft + 1] < xaxis[i]) {
            valuesOnTheLeft = valuesOnTheLeft + 1
        }
        smef = mean(sorteddata[(valuesOnTheLeft + 1):dataLength])
        res[i] = smef - xaxis[i]
    }
    return(res)
}


# load the data
data = read.table("dfl.txt")
data = data[, 5]
xaxis = sort(data)

# calculate empirical mean excess function on dfl data
sampleMef = samplemef(data, xaxis)
sampleMef = sampleMef/1e+06

xaxis = xaxis/1e+06

plot(xaxis, sampleMef, col = "black", xlab = "x (DKK million)", ylab = "e_n(x) (DKK million)", xlim = c(0, 18), ylim = c(0, 
    max(sampleMef))) 
