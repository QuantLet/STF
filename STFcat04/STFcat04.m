%----------------------------------------------------------------------------
% Book:         STF
%----------------------------------------------------------------------------
% See also:     BondZeroCoupon, Burrrnd, mixexprnd, simHPP, simNHPP, Paretornd,
%               simNHPPALP
%----------------------------------------------------------------------------
% Quantlet:     STFcat04
%----------------------------------------------------------------------------
% Description:  Plot of the zero-coupon CAT bond price for the Burr claim 
%               amounts and the non-homogeneous Poisson process governing 
%               the flow of losses.
%---------------------------------------------------------------------------- 
% Author:       Awdesch Melzer 20130727
%---------------------------------------------------------------------------- 

clear all
close all
clc

RandStream.setGlobalStream(RandStream('mt19937ar','seed',100));

%parlambda = 34.2;     % HPP
parlambda  = [35.32,2.32*2*pi,-0.2]; % NHPP1
%parlambda = [35.22,0.224,-0.16];    % NHPP2
distr      = 'Burr';   % 'lognormal'
%params    = [18.3806,1.1052];       % lognormal
params     = [0.4801,3.9495*1e16,2.1524]; %Burr

data       = load('ncl.dat');
A          = mean(data(:,3))*(34.2/4);
na         = 41; % default 41
AE         = 12*A;
D          = (A:(12*A-A)/(na-1):AE)';
B          = 0.25;
nb         = 41; % default 41
BE         = 8*B;
T          = (B:(8*B-B)/(nb-1):BE)';
Tmax       = max(T);
lambda     = 0;
N          = 1000; % default 1000
r          = log(1.025);
Z          = 1.06;
d1         = BondZeroCoupon(Z,D,T,r,lambda,parlambda,distr,params,Tmax,N);
y          = d1(:,1);
x          = d1(:,2)/1e9;
z          = d1(:,3);
yy         = reshape(y,na,nb);
xx         = reshape(x,na,nb);
zz         = reshape(z,na,nb);

% Plot
mesh(xx,yy,zz)
xlim([min(x),max(x)])
ylim([min(y),max(y)])
zlim([min(z),max(z)])
 
set(gca,'FontSize',16,'LineWidth',2,'FontWeight','bold');

% print -painters -dpdf -r600 STFcat04.pdf
% print -painters -dpng -r600 STFcat04.png
