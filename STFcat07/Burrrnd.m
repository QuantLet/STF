function y=Burrrnd(alpha,lambda,tau,n,m)
%BURRRND Random arrays from Burr distribution.
%   R = BURRRND(ALPHA,LAMBDA,TAU,N,M) returns an M-by-N array of random numbers 
%   chosen from the Burr distribution with parameters ALPHA, LAMBDA, TAU.
%
%   The default values for the parameters ALPHA, LAMBDA, TAU, M, N are
%   1, 1, 2, 1, 1, respectively.
%
%   BURRRND uses the inversion method.


  if nargin<5
    m=1;
  end
  if nargin<4
    n=1;
  end
  if nargin<3
    tau=2;
  end
  if nargin<2
    lambda=1;
  end
  if nargin<1
    alpha=1;
  end
  
  y=(lambda*(rand(n,m).^(-1/alpha)-1)).^(1/tau);

