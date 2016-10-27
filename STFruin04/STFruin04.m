clear all;
close all,
clc;
format long;

% produces the exact ruin probability in infinite time for insurance
% collective risk model with mixture of 2 exponentials distribution claims

u1=[0;10^9;5*10^9;10^10;2*10^10;5*10^10]; % initial capital of insurance company (in USD)
theta1=0.3;                               % relative safety loading
dparameters1=[3.5900e-10 , 0.0584 ; 7.5088e-09 , 0.9416]; % exponential parameters (first column) and weights (second column)

% ruin probability for mixture of 2 exponentials claims with exponetial parameters 
%beta1=3.5900e-10, beta2=7.5088e-09, alpha=0.0584 and theta=0.3 (u in USD)

psi=ruinmix2exps(u1,theta1,dparameters1)
