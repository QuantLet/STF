clear all
close all 
clc


dataSP  = load('close.csv');
dataWIG = load('gwp.csv');

dataSP  = diff(log(dataSP)); % log return
dataWIG = diff(log(dataWIG));

[dl_szerSP,podmiotySP]   = size(dataSP);
[dl_szerWIG,podmiotyWIG] = size(dataWIG);

window  = 100; % size of the time window

% analysis of S&P 500
wynikSP = zeros(dl_szerSP - window-1,4);
for t=1:(dl_szerSP - window-1)
    window_data  = dataSP(t:(t+window-1),:);
    wind_dist    = ultra(window_data);
    wind_umlp    = umlp(wind_dist,1);
    wind_bmlp    = bmlp(wind_dist);
    wind_mst     = mst(wind_dist);
    wynikSP(t,1) = 233-t;
    wynikSP(t,2) = mean(wind_umlp(:,3));
    wynikSP(t,3) = mean(wind_bmlp(:,3));
    wynikSP(t,4) = mean(wind_mst(:,3));
    wind_umlp    = [];
    wind_bmlp    = [];
    wind_dist    = [];
    wind_mst     = [];
    window_data  = [];
end;

  
 % analysis of WIG 20
 wynikWIG = zeros(dl_szerWIG - window-1,4);
 for t=1:(dl_szerWIG - window-1)
    window_data   = dataWIG(t:(t+window-1),:);
    wind_dist     = ultra(window_data);
    wind_umlp     = umlp(wind_dist,1);
    wind_bmlp     = bmlp(wind_dist);
    wind_mst      = mst(wind_dist);
    wynikWIG(t,1) = 233-t;
    wynikWIG(t,2) = mean(wind_umlp(:,3));
    wynikWIG(t,3) = mean(wind_bmlp(:,3));
    wynikWIG(t,4) = mean(wind_mst(:,3));
    wind_umlp     = [];
    wind_bmlp     = [];
    wind_dist     = [];
    wind_mst      = [];
    window_data   = [];
end;

subplot(2,1,1), plot(wynikWIG(:,1), wynikWIG(:,4),'*',wynikWIG(:,1), wynikWIG(:,2),'o',wynikWIG(:,1), wynikWIG(:,3),'+')
legend('MST','UMLP','BMLP',2)
xlabel('time');
ylabel('mean distance');
title('UD WIG 20, time window =100d');
subplot(2,1,2), plot(wynikSP(:,1), wynikSP(:,4),'*',wynikSP(:,1), wynikSP(:,2),'o',wynikSP(:,1), wynikSP(:,3),'+')
legend('MST','UMLP','BMLP')
xlabel('time');
ylabel('mean distance');
title('UD S&P 500, time window =100d');