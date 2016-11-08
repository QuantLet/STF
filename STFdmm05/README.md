
[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="888" alt="Visit QuantNet">](http://quantlet.de/)

## [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **STFdmm05** [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/)

```yaml

Name of QuantLet : STFdmm05

Published in : Statistical Tools for Finance and Insurance

Description : 'Presents the analysis of the network evolution. The procedure generate networks in
moving time window and calculate the frequency of connections between elements of the chosen WIG 20
companies (gwp.csv). The unidirectional minimum length path is constructed. Requires umlp.m,
ultra.m and {Statistical Toolbox} from Matlab to run the quantlet.'

Keywords : financial, distance, tree, portfolio, asset, historical moving window

See also : STFdmm11, ultra, umlp

Author : Janusz Miskiewicz, Awdesch Melzer

Submitted : Fri, November 09 2012 by Dedy Dwi Prastyo

Datafile : gwp.csv

```


### R Code:
```r
# clear history
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

########################################### Subroutine umlp(x,y) ##########

umlp = function(x, y) {
    # unidirextional minimum length path algorithm x - the distance matrix y - root of the chain (the number of column) The
    # rsult is presentes as a set of links between nodes
    
    n = nrow(x)
    m = ncol(x)
    net = matrix(0, n - 1, 3)
    onnet = matrix(0, n, 1)
    licz = 0
    onnet[y] = 1
    maxx = 10 * max(x)
    smax = maxx * diag(nrow = nrow(x), ncol = ncol(x))
    x = x + smax
    
    while (licz < n - 1) {
        minx = min(x[y, ])
        it = which(x[y, ] == minx, arr.ind = T)
        
        if (length(it) > 1) {
            tmp = 1
            
            while ((onnet[it[tmp]] == 1) && (tmp < length(it))) {
                tmp = tmp + 1
            }
            
            if (onnet[it[tmp]] == 0) {
                ii = it(tmp)
                it = NULL
                it = ii
                tmp = 0
            } else {
                ii = it[1]
                it = NULL
                it = ii
            }
        }
        if (onnet[it] == 0) {
            licz = licz + 1
            net[licz, 1] = y
            net[licz, 2] = it
            net[licz, 3] = x[it, y]
            onnet[it] = 1
            x[it, y] = maxx
            x[y, it] = maxx
            y = it
        }
        
        if ((onnet[it] == 1) && (onnet[y] == 1)) {
            x[it, y] = maxx
            x[y, it] = maxx
        }
    }
    retval = net
    return(retval)
}

########################################### Main calculation #############



data = read.table("gwp.csv", header = T)  # load data
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
    wind_umlp = umlp(wind_dist, 5)
    for (i in 1L:(podmioty - 1)) {
        if (wind_umlp[i, 1] < wind_umlp[i, 2]) {
            czes[wind_umlp[i, 1], wind_umlp[i, 2]] = czes[wind_umlp[i, 1], wind_umlp[i, 2]] + 1
        } else {
            czes[wind_umlp[i, 2], wind_umlp[i, 1]] = czes[wind_umlp[i, 2], wind_umlp[i, 1]] + 1
        }
    }
    wind_umlp = numeric()
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
% clear history
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
    wind_umlp   = umlp(wind_dist,5);
  for i=1:podmioty-1
        if (wind_umlp(i,1)<wind_umlp(i,2))
            czes(wind_umlp(i,1),wind_umlp(i,2))=czes(wind_umlp(i,1),wind_umlp(i,2))+1;
        else
            czes(wind_umlp(i,2),wind_umlp(i,1))=czes(wind_umlp(i,2),wind_umlp(i,1))+1;
        end;
  end;
    wind_umlp   = [];
    wind_dist   = [];
    window_data = [];
end;

companies   = ['   ','WIG20 ','ACP ','BIO ','BRE ','BZW  ','CEZ ','CPS ','GTN ','GTC ','KGH ','LTS ','PBG ','PEO ','PKN ','PKO ','PXM ','TPS ','TVN '];
companies1  = char('WIG20 ','ACP ','BIO ','BRE ','BZW  ','CEZ ','CPS ','GTN ','GTC ','KGH ','LTS ','PBG ','PEO ','PKN ','PKO ','PXM ','TPS ','TVN ');
freq        = round(100000*czes/sum(sum(czes)))/1000;

freqTEMP    = [companies1(1:18,:),num2str(freq(:,:))];

disp('   WIG20      ACP         BIO         BRE         BZW         CEZ         CPS         GTN         GTC         KGH         LTS         PBG         PEO         PKN         PKO         PXM         TPS         TVN ')
disp(freqTEMP)


```
