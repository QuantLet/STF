% clear variables and close windows
clear all
close all
clc

data = load('close.csv'); % load data
data = diff(log(data));   % log return
[dl_szer,podmioty]=size(data);
okno = 10;

% time window loop
for window=10:(dl_szer-okno-1) 
  wynik_mst=[];
  wynik_umlp=[];
  wynik_bmlp=[];
  for t=1:(dl_szer - window-1)
    window_data=data(t:(t+window),:);
    wind_dist=ultra(window_data);
    wind_mst=mst(wind_dist);
    wynik_mst=[wynik_mst;wind_mst(:,3)];
    wind_umlp=umlp(wind_dist,1);
    wynik_umlp=[wynik_umlp;wind_umlp(:,3)];
    wind_bmlp=bmlp(wind_dist);
    wynik_bmlp=[wynik_bmlp;wind_bmlp(:,3)];
    wind_bmlp=[];
    wind_umlp=[];
    wind_mst=[];
    wind_dist=[];
  end;
  wynik_wind(window-9,1)=window;
  wynik_wind_mst(window-9,1)=mean(wynik_mst);
  wynik_wind_mst(window-9,2)=std(wynik_mst);
  wynik_wind_mst(window-9,3)=skewness(wynik_mst);
  wynik_wind_mst(window-9,4)=kurtosis(wynik_mst);
  wynik_wind_bmlp(window-9,1)=mean(wynik_bmlp);
  wynik_wind_bmlp(window-9,2)=std(wynik_bmlp);
  wynik_wind_bmlp(window-9,3)=skewness(wynik_bmlp);
  wynik_wind_bmlp(window-9,4)=kurtosis(wynik_bmlp);
  wynik_wind_umlp(window-9,1)=mean(wynik_umlp);
  wynik_wind_umlp(window-9,2)=std(wynik_umlp);
  wynik_wind_umlp(window-9,3)=skewness(wynik_umlp);
  wynik_wind_umlp(window-9,4)=kurtosis(wynik_umlp);
end;

subplot(2,2,1), plot(wynik_wind(:,1), wynik_wind_mst(:,1),wynik_wind(:,1), wynik_wind_umlp(:,1),wynik_wind(:,1), wynik_wind_bmlp(:,1),'LineWidth',2)
legend('MST','UMLP','BMLP')
xlabel('window');
ylabel('mean');
title('S&P 500');
subplot(2,2,2),plot(wynik_wind(:,1), wynik_wind_mst(:,2),wynik_wind(:,1), wynik_wind_umlp(:,2),wynik_wind(:,1), wynik_wind_bmlp(:,2),'LineWidth',2) 
legend('MST','UMLP','BMLP')
xlabel('window');
ylabel('std');
title('S&P 500');
subplot(2,2,3),plot(wynik_wind(:,1), wynik_wind_mst(:,3),wynik_wind(:,1), wynik_wind_umlp(:,3),wynik_wind(:,1), wynik_wind_bmlp(:,3),'LineWidth',2)
legend('MST','UMLP','BMLP')
xlabel('window');
ylabel('skewness');
title('S&P 500');
subplot(2,2,4),plot(wynik_wind(:,1), wynik_wind_mst(:,4),wynik_wind(:,1), wynik_wind_umlp(:,4),wynik_wind(:,1), wynik_wind_bmlp(:,4),'LineWidth',2)
legend('MST','UMLP','BMLP')
xlabel('window');
ylabel('kurtosis');
title('S&P 500');
