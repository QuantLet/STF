clear all
close all
clc

data  = load('close.csv');

data  = abs(diff(log(data))); % abs log return

[dl_szer,podmioty] = size(data);
data(data==0)      = 0.0000001;
disp('The procedure may last several hours since 9025 networks are constructed.')

% time window loop
for window1=5:100
    theil_data  = theil(data,window1);
    for window2 = 5:100
    
% moving time window
        wynik   = [];
        for t=1:(dl_szer - window1-1-window2)
            window_data = theil_data(t:(t+window2),:);
            wind_dist   = manh(window_data);
            wind_mst    = mst(wind_dist);
            wynik       = [wynik;wind_mst(:,3)];
            wind_mst    = [];
            wind_dist   = [];
        end;
        wynik_wind_mean(window1-4,window2-4) = mean(wynik);
        wynik_wind_std(window1-4,window2-4)  = std(wynik);
        wynik_wind_skew(window1-4,window2-4) = skewness(wynik);
        wynik_wind_kurt(window1-4,window2-4) = kurtosis(wynik);
    end;
end;

subplot(2,2,1), mesh(wynik_wind_mean,'DisplayName','mean distance, MST');figure(gcf)
xlabel('T_1');
ylabel('T_2');
zlabel('mean');
title('mean distance, MST');
subplot(2,2,2),mesh(wynik_wind_std,'DisplayName','Std, MST');figure(gcf) 
xlabel('T_1');
ylabel('T_2');
zlabel('std');
title('Std, MST');
subplot(2,2,3),mesh(wynik_wind_skew,'DisplayName','Skewness, MST');figure(gcf)
xlabel('T_1');
ylabel('T_2');
zlabel('skewness');
title('Skewness, MST');
subplot(2,2,4),mesh(wynik_wind_kurt,'DisplayName','Kurtosis, MST');figure(gcf)
xlabel('T_1');
ylabel('T_2');
zlabel('kurtosis');
title('Kurtosis, MST');