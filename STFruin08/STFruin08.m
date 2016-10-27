clear all;
close all,
clc;
format long;

u=[0;10^9;5*10^9;10^10;2*10^10;5*10^10]; % initial capital of insurance company (in USD)
theta=0.3;                               % relative safety loading
dparameters1=[0.0584 , 3.5900e-10 ; 0.9416 , 7.5088e-09]; % weights (first column) and exponential parameters (second column)

m=moments(1,dparameters1);  % 1st raw moment
m2=moments(2,dparameters1); % 2nd raw moment
m3=moments(3,dparameters1); % 3nd raw moment

delta1=(2*m*theta)/(m2+((4*m*m3/3/m2)-m2)*theta) % scale parameter of gamma distribution function
delta2=(1+theta)/(1+((4*m*m3/3/m2^2)-1)*theta)   % shape parameter of gamma distribution function

% the Beekman-Bowers approximation for mixture of 2 exponentials claims with beta1=3.5900e-10, beta2=7.5088e-09, alpha=0.0584 and theta=0.3 (u in USD) 
psi=(1-gamcdf(u,delta2,1/delta1))/(1+theta)