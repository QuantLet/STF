# clear variables and close windows
rm(list = ls(all = TRUE))
graphics.off()

########################################## Subroutine manh(x) ############


manh = function(x) {
    # Manhattan distence between time series normalised by the time series length.  x - time series
    x = as.matrix(x)
    h = nrow(x)
    k = ncol(x)
    result = matrix(0, k, k)
    for (j in 1:k) {
        for (i in 1:k) {
            result[i, j] = abs(mean(x[, i] - x[, j]))
        }
    }
    retval = result
    return(retval)
}


########################################### Subroutine mst(x) ############


mst = function(x) {
    # Algorithm generates minimum spanning tree The rsult is presentes as a set of links between nodes
    n = nrow(x)
    m = ncol(x)
    true = upper.tri(x)
    x = true * x
    net = matrix(0, n - 1, 3)
    onnet = rep(as.integer(0), n)
    klaster = rep(as.integer(0), n)
    klast = 0L
    licz = 0L
    # check if the matrics is symmetric and positive
    maxx = max(apply(x, 2, max))
    smax = 10 * abs(maxx)
    x[x == 0] = smax
    while (licz < n - 1) {
        
        minx = min(apply(x, 2, min))
        d = which(x <= minx, arr.ind = T)
        i = d[, 1]
        j = d[, 2]
        if (length(i) > 1) {
            ii = i[1]
            jj = j[1]
            i = 0
            j = 0
            i = ii
            j = jj
        }
        
        if (onnet[i] == 0 & onnet[j] == 0) {
            licz = licz + 1L
            net[licz, 1] = i
            net[licz, 2] = j
            klast = klast + 1L
            klaster[i] = klast
            klaster[j] = klast
            net[licz, 3] = min(x[i, j], x[j, i])
            onnet[i] = 1
            onnet[j] = 1
            x[i, j] = smax
            x[j, i] = smax
            
        } else if (onnet[i] == 0 & onnet[j] == 1) {
            licz = licz + 1
            net[licz, 1] = i
            net[licz, 2] = j
            net[licz, 3] = min(x[i, j], x[j, i])
            onnet[i] = 1
            klaster[i] = klaster[j]
            x[i, j] = smax
            x[j, i] = smax
        } else if (onnet[i] == 1 & onnet[j] == 0) {
            licz = licz + 1L
            net[licz, 1] = i
            net[licz, 2] = j
            net[licz, 3] = min(x[i, j], x[j, i])
            onnet[j] = 1
            klaster[j] = klaster[i]
            x[i, j] = smax
            x[j, i] = smax
        } else if (onnet[i] == 1 & onnet[j] == 1 & klaster[i] == klaster[j]) {
            x[i, j] = smax
            x[j, i] = smax
        } else if (onnet[i] == 1 & onnet[j] == 1 & klaster[i] != klaster[j]) {
            licz = licz + 1L
            net[licz, 1] = i
            net[licz, 2] = j
            net[licz, 3] = min(x[i, j], x[j, i])
            klaster[klaster == klaster[i]] = klaster[j]
        }
    }
    retval = net
    return(retval)
}

############################################# Subroutine theil(x) ############


theil = function(x, n) {
    # Converts the given time series into time series of Theil index n - size of moving window x - analysed time series
    x = as.matrix(x)
    n_pocz = nrow(x)
    k = ncol(x)
    
    rozm = n_pocz - n + 1
    
    retval = matrix(0, rozm, k)
    for (i in 1:rozm) {
        sr = apply(x[i:(i + n - 1), ], 2, mean)
        temp = x[i:(i + n - 1), ]/(matrix(1, n, 1) %*% sr)
        temp = temp * log(temp)
        retval[i, ] = apply(temp, 2, mean)
    }
    return(retval)
}



########################################### Main calculation #############




dataSP = read.table("close.csv", header = T)
dataSP = as.matrix(dataSP)
dataSP = abs(diff(log(dataSP)))  # abs log return

dl_szer = nrow(dataSP)
podmioty = ncol(dataSP)

dataSP[dataSP == 0] = 1e-07

window1 = 50
window2 = 50
theil_data = theil(dataSP, window1)

wynikSP = matrix(0, dl_szer - window1 - 1 - window2, 4)
czes = matrix(0, podmioty, podmioty)

for (t in 1:(dl_szer - window1 - 1 - window2)) {
    window_data = theil_data[t:(t + window2), ]
    wind_dist = manh(window_data)
    wind_mst = mst(wind_dist)
    for (i in 1:(podmioty - 1)) {
        if (wind_mst[i, 1] < wind_mst[i, 2]) {
            czes[wind_mst[i, 1], wind_mst[i, 2]] = czes[wind_mst[i, 1], wind_mst[i, 2]] + 1
        } else {
            czes[wind_mst[i, 2], wind_mst[i, 1]] = czes[wind_mst[i, 2], wind_mst[i, 1]] + 1
        }
    }
    wind_mst = numeric()
    wind_dist = numeric()
}

companies = c("ABB", "AAPL", "BA", "KO", "EMR", "GE", "HPQ", "HIT", "IBM", "INTC", "JNJ", "LMT", "MSFT", "NOC", "NVS", "CL", 
    "PEP", "PG", "TSEM", "WEC")  # company names

freq = round(1e+05 * czes/sum(sum(czes)))/1000

freqTEMP = freq
rownames(freqTEMP) = companies
colnames(freqTEMP) = companies
freqTEMP

 
