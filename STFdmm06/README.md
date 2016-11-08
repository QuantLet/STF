
[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="888" alt="Visit QuantNet">](http://quantlet.de/)

## [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **STFdmm06** [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/)

```yaml

Name of QuantLet : STFdmm06

Published in : Statistical Tools for Finance and Insurance

Description : 'Presents the analysis of the network evolution. The procedure generate networks in
moving time window and calculate the frequency of connections between elements of the chosen WIG 20
companies (gwp.csv). The bidirectional minimum length path is constructed. Requires ultra.m, bmlp.m
to run the program.'

Keywords : financial, distance, tree, portfolio, asset, historical moving window

See also : STFdmm11, bmlp, ultra

Author : Janusz Miskiewicz, Awdesch Melzer

Submitted : Wed, November 21 2012 by Dedy Dwi Prastyo

Datafile : gwp.csv

```


### R Code:
```r
# clear variables and close windows
rm(list = ls(all = TRUE))
graphics.off()



########################################### Subroutine bmlp(x) ############

bmlp = function(x) {
    # bidirextional minimum length path algoritm x - the distance matrix The result is presentes as a set of links between nodes
    n = nrow(x)
    m = ncol(x)
    net = matrix(0, n - 1, 3)
    onnet = matrix(0, n, 1)
    end1 = 0
    end2 = 0
    licz = 0
    
    # the distance matrics should be symmetric and positive
    maxx = 10 * max(x)
    smax = maxx * diag(nrow = nrow(x), ncol = ncol(x))
    x = x + smax
    
    # the first pair
    minx = min(x)
    ij = which(x == minx, arr.ind = T)
    i = ij[, 1]
    j = ij[, 2]
    
    if (length(i) == 1) {
        end1 = i
        end2 = j
        onnet[end1] = 1
        onnet[end2] = 1
        net[1, 1] = end1
        net[1, 2] = end2
        net[1, 3] = minx
        licz = 1
        x[end1, end2] = maxx
        x[end2, end1] = maxx
    } else {
        end1 = i[1]
        end2 = j[1]
        onnet[end1] = 1
        onnet[end2] = 1
        net[1, 1] = end1
        net[1, 2] = end2
        net[1, 3] = minx
        licz = 1
        x[end1, end2] = maxx
        x[end2, end1] = maxx
    }
    
    while (licz < n - 1) {
        minx1 = min(x[end1, ])
        minx2 = min(x[end2, ])
        if (minx1 < minx2) {
            y = end1
            minx = minx1
        } else {
            y = end2
            minx = minx2
        }
        i = which(x[y, ] == minx, arr.ind = T)
        if (length(i) > 1) {
            tmp = 1
            while ((onnet[i[tmp]] == 1) && (tmp < length(i))) {
                tmp = tmp + 1
            }
            if (onnet[i[tmp]] == 0) {
                ii = i(tmp)
                i = NULL
                i = ii
                tmp = 0
            } else {
                ii = i[1]
                i = NULL
                i = ii
            }
        }
        if (onnet[i] == 0) {
            licz = licz + 1
            net[licz, 1] = y
            net[licz, 2] = i
            net[licz, 3] = x[i, y]
            onnet[i] = 1
            x[i, y] = maxx
            x[y, i] = maxx
            y = i
        }
        if ((onnet[i] == 1) && (onnet[y] == 1)) {
            x[i, y] = maxx
            x[y, i] = maxx
        }
    }
    retval = net
    return(retval)
}

########################################### Subroutine ultra(x) ############


ultra = function(x) {
    # Ultrametric distance between time series.  x - time series matrix
    h = nrow(x)
    k = ncol(x)
    retval = sqrt(abs(0.5 * (matrix(1, k, k) - cor(x))))
    return(retval)
}


########################################### Main calculation #############




data = read.table("gwp.csv", header = T)  # load data
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
    wind_bmlp = bmlp(wind_dist)
    
    for (i in 1:(podmioty - 1)) {
        if (wind_bmlp[i, 1] < wind_bmlp[i, 2]) {
            czes[wind_bmlp[i, 1], wind_bmlp[i, 2]] = czes[wind_bmlp[i, 1], wind_bmlp[i, 2]] + 1
        } else {
            czes[wind_bmlp[i, 2], wind_bmlp[i, 1]] = czes[wind_bmlp[i, 2], wind_bmlp[i, 1]] + 1
        }
    }
    wind_bmlp = numeric()
    wind_dist = numeric()
    window_data = numeric()
}

companies = c("WIG20 ", "ACP ", "BIO ", "BRE ", "BZW  ", "CEZ ", "CPS ", "GTN ", "GTC ", "KGH ", "LTS ", "PBG ", "PEO ", "PKN ", 
    "PKO ", "PXM ", "TPS ", "TVN ")
freq = round(1e+05 * czes/sum(sum(czes)))/1000

freqTEMP = freq
rownames(freqTEMP) = companies
colnames(freqTEMP) = companies
freqTEMP
 

```

### MATLAB Code:
```matlab
% clear variables and close windows
clear all
close all
clc


data    = load('gwp.csv'); % load data
data    = diff(log(data)); % log return
[dl_szer,podmioty] = size(data);
czes    = zeros(podmioty);
window  = 100; %time window size

% moving time window
  wynik = zeros(dl_szer - window-1,5);
for t=1:(dl_szer - window-1)
    window_data = data(t:(t+window-1),:);
    wind_dist   = ultra(window_data);
    wind_bmlp   = bmlp(wind_dist);
  for i=1:podmioty-1
        if (wind_bmlp(i,1)<wind_bmlp(i,2))
            czes(wind_bmlp(i,1),wind_bmlp(i,2))=czes(wind_bmlp(i,1),wind_bmlp(i,2))+1;
        else
            czes(wind_bmlp(i,2),wind_bmlp(i,1))=czes(wind_bmlp(i,2),wind_bmlp(i,1))+1;
        end;
  end;
    wind_bmlp   = [];
    wind_dist   = [];
    window_data = [];
end;

companies       = ['   ','WIG20 ','ACP ','BIO ','BRE ','BZW  ','CEZ ','CPS ','GTN ','GTC ','KGH ','LTS ','PBG ','PEO ','PKN ','PKO ','PXM ','TPS ','TVN '];
companies1      = char('WIG20 ','ACP ','BIO ','BRE ','BZW  ','CEZ ','CPS ','GTN ','GTC ','KGH ','LTS ','PBG ','PEO ','PKN ','PKO ','PXM ','TPS ','TVN ');
freq            = round(100000*czes/sum(sum(czes)))/1000;

freqTEMP        = [companies1(1:18,:),num2str(freq(:,:))];

disp('   WIG20      ACP         BIO         BRE         BZW         CEZ         CPS         GTN         GTC         KGH         LTS         PBG         PEO         PKN         PKO         PXM         TPS         TVN ')
disp(freqTEMP)


```
