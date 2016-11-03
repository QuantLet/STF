function y=Paretornd(alpha,lambda,n,m)
%PARETORND Random arrays from Pareto distribution.
%   Y = PARETORND(ALPHA,LAMBDA,N,M) returns an M-by-N array of random numbers 
%   chosen from the Pareto distribution with parameters ALPHA, LAMBDA.
%
%   The default values for ALPHA, LAMBDA, N, M 1, 1, 1, 1, respectively.
%
%   PARETORND uses the inversion method.

  if nargin<4
    m=1;
  end
  if nargin<3
    n=1;
  end
  if nargin<2
    lambda=1;
  end
  if nargin<1
    alpha=1;
  end
 
  y=lambda*(rand(n,m).^(-1/alpha)-1);
