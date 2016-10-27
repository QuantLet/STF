clear all;
close all,
clc;
format long;

u1=[0;10^9;5*10^9;10^10;2*10^10;5*10^10]; % initial capital of insurance company (in USD)
theta1=0.3;                               % relative safety loading
dparameters1=[0.0584 , 3.5900e-10 ; 0.9416 , 7.5088e-09]; % weights (second vector) and exponential parameters (first vector)

m=moments(1,dparameters1);  % 1st raw moment
m2=moments(2,dparameters1); % 2nd raw moment
m3=moments(3,dparameters1); % 3nd raw moment
	
R=adjR(theta1,dparameters1,m,m2,m3);    % adjustment coefficient R for mixture of 2 exponentials distribution claims
mgfprim=mgfs(R,1,dparameters1); % moment generating function for mixture of 2 exponentials distribution claims 

C=(theta1*m)/(mgfprim-m*(1+theta1));

% the Cramer-Lundberg approximation for mixture of 2 exponentials claims with beta1=3.5900e-10, beta2=7.5088e-09, alpha=0.0584 and theta=0.3 (u in USD)
psi=C*exp(-R*u1)
