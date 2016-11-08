clear all
close all
clc

dataSP     = load('close.csv');

dataSP     = abs(diff(log(dataSP))); % abs log return
[dl_szer,podmioty] = size(dataSP);

dataSP(dataSP==0)  = 0.0000001;
window1    = 50;
window2    = 50;
theil_data = theil(dataSP,window1);

wynikSP    = zeros(dl_szer - window1-1-window2,4);
czes       = zeros(podmioty);

  for t=1:(dl_szer - window1-1-window2)
    window_data  = theil_data(t:(t+window2),:);
    wind_dist    = manh(window_data);
    wind_mst     = mst(wind_dist);
    wind_umlp    = umlp(wind_dist,5);
    wind_bmlp    = bmlp(wind_dist);
    wynikSP(t,1) = dl_szer-t;
    wynikSP(t,2) = mean(wind_mst(:,3));
    wynikSP(t,3) = mean(wind_umlp(:,3));
    wynikSP(t,4) = mean(wind_bmlp(:,3));
    wind_umlp    = [];
    wind_bmlp    = [];
    wind_mst     = [];
    wind_dist    = [];
  end;
  
dataWIG    = load('gwp.csv');

dataWIG    = abs(diff(log(dataWIG))); % abs log return
[dl_szer,podmioty]  = size(dataWIG);

dataWIG(dataWIG==0) = 0.0000001;
window1    = 50;
window2    = 50;
theil_data = theil(dataWIG,window1);

wynikWIG   = zeros(dl_szer - window1-1-window2,4);
czes       = zeros(podmioty);

  for t=1:(dl_szer - window1-1-window2)
    window_data   = theil_data(t:(t+window2),:);
    wind_dist     = manh(window_data);
    wind_mst      = mst(wind_dist);
    wind_umlp     = umlp(wind_dist,5);
    wind_bmlp     = bmlp(wind_dist);
    wynikWIG(t,1) = dl_szer-t;
    wynikWIG(t,2) = mean(wind_mst(:,3));
    wynikWIG(t,3) = mean(wind_umlp(:,3));
    wynikWIG(t,4) = mean(wind_bmlp(:,3));
    wind_umlp     = [];
    wind_bmlp     = [];
    wind_mst      = [];
    wind_dist     = [];
  end;  
  
subplot(2,1,1), plot(wynikWIG(:,1), wynikWIG(:,2),'*',wynikWIG(:,1), wynikWIG(:,3),'o',wynikWIG(:,1), wynikWIG(:,4),'+')
legend('MST','UMLP','BMLP')
xlabel('time');
ylabel('mean distance');
title('ED WIG 20, T_1 =50d, T_2=50d');
subplot(2,1,2), plot(wynikSP(:,1), wynikSP(:,2),'*',wynikSP(:,1), wynikSP(:,3),'o',wynikSP(:,1), wynikSP(:,4),'+')
legend('MST','UMLP','BMLP')
xlabel('time');
ylabel('mean distance');
title('ED S&P 500, T_1 =50d, T_2=50d');