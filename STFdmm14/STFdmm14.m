clear all
close all
clc



dataSP      = load('close.csv');

dataSP      = abs(diff(log(dataSP))); % abs log return
[dl_szer,podmioty] = size(dataSP);

dataSP(dataSP==0)  = 0.0000001;
window1     = 50;
window2     = 50;
theil_data  = theil(dataSP,window1);

wynikSP     = zeros(dl_szer - window1-1-window2,4);
czes        = zeros(podmioty);
for t=1:(dl_szer - window1-1-window2)
    window_data     = theil_data(t:(t+window2),:);
    wind_dist       = manh(window_data);
    wind_mst        = mst(wind_dist);
    for i=1:podmioty-1
        if (wind_mst(i,1)<wind_mst(i,2))
            czes(wind_mst(i,1),wind_mst(i,2)) = czes(wind_mst(i,1),wind_mst(i,2))+1;
        else
            czes(wind_mst(i,2),wind_mst(i,1)) = czes(wind_mst(i,2),wind_mst(i,1))+1;
        end;
    end;
    wind_mst        = [];
    wind_dist       = [];
end;
freq        = round(100000*czes/sum(sum(czes)))/1000;
companies   = ['   ','ABB  ','AAPL ','BA   ','KO   ','EMR   ','GE   ','HPQ  ','HIT  ','IBM  ','INTC ','JNJ  ','LMT  ','MSFT ','NOC  ','NVS  ','CL    ','PEP  ','PG    ','TSEM ','WEC  '];
companies1  = char('ABB  ','AAPL ','BA   ','KO   ','EMR   ','GE   ','HPQ  ','HIT  ','IBM  ','INTC ','JNJ  ','LMT  ','MSFT ','NOC  ','NVS  ','CL    ','PEP  ','PG    ','TSEM ','WEC  ');

freqTEMP    = [companies1(1:20,:),num2str(freq(:,:))];

disp('      ABB        AAPL        BA           KO          EMR         GE      HPQ         HIT         IBM         INTC        JNJ         LMT         MSFT        NOC         NVS             CL          PEP         PG      TSEM        WEC')
disp(freqTEMP)
  
