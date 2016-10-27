clear all
close all
clc

 RandStream.setGlobalStream(RandStream('mt19937ar','seed',99));

lambda    = 0;                        % default Poisson process NHPP1 
parlambda = [35.32,2.32*2*pi,-0.2];   % NHPP1
distr     = 'lognormal'; % default distribution
params    = [18.3806,1.1052];
T         = 10;          % max time (years)
N         = 100;         % number of trajectories in NHPPALP, default 5000
qu        = [0.05,0.95]; % quantiles
step      = 0.05;        % step in quantile plot

%=========== aggregate process - single realization 
y = simNHPPALP(lambda,parlambda,distr,params,T,1);
y(:,:,2) = y(:,:,2)/1e+9;

%=========== real PCS trajectory 
c    = load('ncl.dat');
t2   = c(ceil((1:2*size(c(:,2),1))/2),2);
t2   = [0;t2;T];
PCS  = cumsum(c(:,3))/1e+9;
PCS2 = PCS(ceil((1:2*size(PCS,1))/2));
PCS2 = [0;0;PCS2];
z    = [t2,PCS2];

%=========== mean of aggregate loss process (only for NHPP1 and lognormal loss size distribution)
t    = [(0:100*T)/100]';
RP   = exp(params(1)+params(2).^2/2).*(parlambda(1)*t-parlambda(2)/2/pi*(cos(2*pi*(t+parlambda(3))) - cos(2*pi*parlambda(3))));
me   = [t,RP/1e+9];

%=========== quantiles
vqu  = quantilelines(simNHPPALP(lambda,parlambda,distr,params,T,N),step,qu');

%=========== plot

plot(y(:,:,1),y(:,:,2),'b-')
hold on
plot(me(:,1),me(:,2),'r-.')
plot(z(:,1),z(:,2),'g-','LineWidth',1.5)
plot(y(:,:,1),60,'k-','LineWidth',1.5)
plot(vqu(:,1),vqu(:,2)/1e+9, 'm--')
plot(vqu(:,1),vqu(:,3)/1e+9, 'm--')
xlim([0,10])
ylim([0,120])
xlabel('Years','FontSize',16,'FontWeight','bold')
ylabel('Aggregate loss process (USD billion)','FontSize',16,'FontWeight','bold')
hold off
box on
set(gca,'FontSize',16,'LineWidth',2,'FontWeight','bold');
% print -painters -dpdf -r600 STFcat02.pdf
% print -painters -dpng -r600 STFcat02.png
