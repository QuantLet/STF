
[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="888" alt="Visit QuantNet">](http://quantlet.de/)

## [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **STFdmm10** [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/)

```yaml

Name of QuantLet : STFdmm10

Published in : Statistical Tools for Finance and Insurance

Description : 'Presents the analysis of the network evolution. The procedure generate networks in
moving time window and calculate the frequency of connections between elements of the chosen S&P
500 companies (close.csv). The minimum spanning tree is constructed. Requires ultra.m, mst.m to run
the program.'

Keywords : financial, distance, tree, portfolio, asset, historical moving window

See also : STFdmm11, mst, ultra

Author : Janusz Miskiewicz, Awdesch Melzer

Submitted : Wed, November 21 2012 by Dedy Dwi Prastyo

Datafile : close.csv

```


### R Code:
```r
# --------------------------------------------------------------------- Book: STF2
# --------------------------------------------------------------------- Quantlet: STFdmm10
# --------------------------------------------------------------------- Description: STFdmm10 presents the analysis of the
# network evolution.  The procedure generate networks in moving time window and calculate the frequency of connections
# between elements of the chosen S&P 500 companies (close.csv).  The minimum spanning tree is constructed.
# --------------------------------------------------------------------- Usage: -
# --------------------------------------------------------------------- See also: STFdmm01, STFdmm02, STFdmm03, STFdmm04,
# STFdmm05, STFdmm06, STFdmm07, STFdmm08, STFdmm09 ---------------------------------------------------------------------
# Inputs: none --------------------------------------------------------------------- Output: The minimum spanning tree is
# constructed.  --------------------------------------------------------------------- Example: -
# --------------------------------------------------------------------- Author: Matlab: Janusz Miskiewicz R: Awdesch Melzer
# 20121121 ---------------------------------------------------------------------

rm(list = ls(all = TRUE))
graphics.off()


########################################### Subroutine ultra(x) ############


ultra = function(x) {
    # Ultrametric distance between time series.  x - time series matrix
    h = nrow(x)
    k = ncol(x)
    retval = sqrt(abs(0.5 * (matrix(1, k, k) - cor(x))))
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


########################################### Main calculation #############



data = read.table("close.csv", header = T)  # load data
data = as.matrix(data)
data = diff(log(data))  # log return
dl_szer = nrow(data)
podmioty = ncol(data)
czes = matrix(0, podmioty, podmioty)
window = 100  #time window size



# moving time window
wynik = matrix(0, dl_szer - window - 1, 5)

for (t in 1:(dl_szer - window - 1)) {
    window_data = data[t:(t + window - 1), ]
    wind_dist = ultra(window_data)
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
    window_data = numeric()
}



companies = c("ABB", "AAPL", "BA", "KO", "EMR", "GE", "HPQ", "HIT", "IBM", "INTC", "JNJ", "LMT", "MSFT", "NOC", "NVS", "CL", 
    "PEP", "PG", "TSEM", "WEC")  # company names

freq = round(1e+05 * czes/sum(sum(czes)))/1000

freqTEMP = freq
rownames(freqTEMP) = companies
colnames(freqTEMP) = companies
freqTEMP
 

```

### MATLAB Code:
```matlab
clear all
close all
clc

data    = load('close.csv'); % load data
data    = diff(log(data)); % log return
[dl_szer,podmioty] = size(data);
czes    = zeros(podmioty);
window  = 100; %time window size

% moving time window
  wynik = zeros(dl_szer - window-1,5);
  
for t=1:(dl_szer - window-1)
    window_data = data(t:(t+window-1),:);
    wind_dist   = ultra(window_data);
    wind_mst    = mst(wind_dist);
  for i=1:podmioty-1
        if (wind_mst(i,1)<wind_mst(i,2))
            czes(wind_mst(i,1),wind_mst(i,2))=czes(wind_mst(i,1),wind_mst(i,2))+1;
        else
            czes(wind_mst(i,2),wind_mst(i,1))=czes(wind_mst(i,2),wind_mst(i,1))+1;
        end;
  end;
    wind_mst    = [];
    wind_dist   = [];
    window_data = [];
end;


companies       = ['   ','ABB  ','AAPL ','BA   ','KO   ','EMR   ','GE   ','HPQ  ','HIT  ','IBM  ','INTC ','JNJ  ','LMT  ','MSFT ','NOC  ','NVS  ','CL    ','PEP  ','PG    ','TSEM ','WEC  '];
companies1      = char('ABB  ','AAPL ','BA   ','KO   ','EMR   ','GE   ','HPQ  ','HIT  ','IBM  ','INTC ','JNJ  ','LMT  ','MSFT ','NOC  ','NVS  ','CL    ','PEP  ','PG    ','TSEM ','WEC  ');
freq            = round(100000*czes/sum(sum(czes)))/1000;

freqTEMP        = [companies1(1:20,:),num2str(freq(:,:))];

disp('      ABB        AAPL        BA           KO          EMR         GE      HPQ         HIT         IBM         INTC        JNJ         LMT         MSFT        NOC         NVS             CL          PEP         PG      TSEM        WEC')
disp(freqTEMP)


```
