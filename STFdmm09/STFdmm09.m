clear all
close all
clc


data    = load('close.csv');
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


companies       = ['   ','ABB  ','AAPL ','BA   ','KO   ','EMR   ','GE   ','HPQ  ','HIT  ','IBM  ','INTC ','JNJ  ','LMT  ','MSFT ','NOC  ','NVS  ','CL    ','PEP  ','PG    ','TSEM ','WEC  '];
companies1      = char('ABB  ','AAPL ','BA   ','KO   ','EMR   ','GE   ','HPQ  ','HIT  ','IBM  ','INTC ','JNJ  ','LMT  ','MSFT ','NOC  ','NVS  ','CL    ','PEP  ','PG    ','TSEM ','WEC  ');
freq            = round(100000*czes/sum(sum(czes)))/1000;

freqTEMP        = [companies1(1:20,:),num2str(freq(:,:))];

disp('      ABB        AAPL        BA           KO          EMR         GE      HPQ         HIT         IBM         INTC        JNJ         LMT         MSFT        NOC         NVS             CL          PEP         PG      TSEM        WEC')
disp(freqTEMP)

