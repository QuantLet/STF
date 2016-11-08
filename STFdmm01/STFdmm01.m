% clear variables and close windows
clear all
close all
clc


data      = load('close.csv'); % load data
data      = diff(log(data));   % log returns
data      = corr(data);        % correlations of log returns
odl       = sqrt(0.5*(1-data));
mst_odl   = mst(odl);          % minimum span tree

format short;
companies = char('ABB  ','AAPL ','BA   ','KO   ','EMR   ','GE   ','HPQ  ','HIT  ','IBM  ','INTC ','JNJ  ','LMT  ','MSFT ','NOC  ','NVS  ','CL    ','PEP  ','PG    ','TSEM ','WEC  ');
odl       = triu(odl,1);
MST_UD    = [companies(mst_odl(:,1),:),companies(mst_odl(:,2),:),num2str(mst_odl(:,3),4)];
MST_UD