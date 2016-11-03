% clear variables and close windows
clear all
close all
clc


% Pseudo random numbers
RandStream.setGlobalStream(RandStream('mt19937ar','seed',1));

nsim = 100 ;
T=23;
u = 4e8;
theta = 0.5;
a = 17.9937;
b = 3.5759*2;
q = [0.001,0.01,0.05,0.25,0.5,0.75,0.95,0.99,0.999]; 

data = load('dfl.dat');
data = data(find(data(:,5)~=0),[2,5]);
dantime = data(:,1);
danval = data(:,2);

% lognormal
mulogn = 12.5247;
sigmalogn = 1.5384;

[t1,y1]=simNHPPRP(u,theta,1,[a,b],'lognormal',[mulogn,sigmalogn],T,1);
y1=y1/1e6;
plot(t1,y1,'r')
hold on

[t2,y2]=simNHPPRPRT(u,theta,1,[a,b],'lognormal',[mulogn,sigmalogn],dantime,danval,T);
y2=y2/1e6;
plot(t2,y2)

[t3,y3]=simNHPPRP(u,theta,1,[a,b],'lognormal',[mulogn,sigmalogn],T,nsim);
[tq,valq]=quantiles(t3,y3,0.2,q);
plot(tq,valq/1e6,'k')
xlim([0,T]);
xlabel('Time (years)','FontSize',16,'FontWeight','Bold')
ylabel('Capital (DKK million)','FontSize',16,'FontWeight','Bold');
set(gca,'LineWidth',1.6,'FontSize',16,'FontWeight','Bold');
box on

% to save the plot in pdf or png please uncomment next 2 lines:
%print -painters -dpdf -r600 STFloss10_01.pdf
%print -painters -dpng -r600 STFloss10_01.png

%Pareto
alphaPareto=1.3127;
lambdaPareto=4.0588e+005;

[t1,y1]=simNHPPRP(u,theta,1,[a,b],'Pareto',[alphaPareto,lambdaPareto],T,1);
 y1=y1/1e6;
 figure
plot(t1,y1,'r')
hold on

[t2,y2]=simNHPPRPRT(u,theta,1,[a,b],'Pareto',[alphaPareto,lambdaPareto],dantime,danval,T);
 y2=y2/1e6;
plot(t2,y2)

[t3,y3]=simNHPPRP(u,theta,1,[a,b],'Pareto',[alphaPareto,lambdaPareto],T,nsim);
[tq,valq]=quantiles(t3,y3,0.2,q);
plot(tq,valq/1e6,'k')
xlim([0,T]);
xlabel('Time (years)','FontSize',16,'FontWeight','Bold')
ylabel('Capital (DKK million)','FontSize',16,'FontWeight','Bold');
set(gca,'LineWidth',1.6,'FontSize',16,'FontWeight','Bold');
box on

% to save the plot in pdf or png please uncomment next 2 lines:
%print -painters -dpdf -r600 STFloss10_02.pdf
%print -painters -dpng -r600 STFloss10_02.png

