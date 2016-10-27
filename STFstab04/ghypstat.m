function [ex,varx] = ghypstat(lambda,alpha,beta,delta,mu)
%GHYPSTAT Returns expected value and variance of generalized hyperbolic 
%   distribution
%   SYNTAX: [ex,varx] = ghypstat(lambda, alpha, beta, delta, mu)

% Written by Rafal Weron, 20100730

psi = alpha^2 - beta^2;
% chi = delta^2
zeta = delta * sqrt(psi);

% ex = mu + (beta*delta/sqrt(psi))*(mbessel3(lambda+1,zeta)/mbessel3(lambda,zeta));
% varx = mbessel3(lambda+1,zeta)/(zeta*mbessel3(lambda,zeta)); 
% varx = varx + (beta*beta/psi)*( mbessel3(lambda+2,zeta)/mbessel3(lambda,zeta) - (mbessel3(lambda+1,zeta)/(zeta*mbessel3(lambda,zeta)))^2);
% varx = delta^2 * varx;

ex = mu + (beta*delta/sqrt(psi))*(besselk(lambda+1,zeta)/besselk(lambda,zeta));
varx = besselk(lambda+1,zeta)/(zeta*besselk(lambda,zeta));
varx = varx + (beta*beta/psi)*( besselk(lambda+2,zeta)/besselk(lambda,zeta) - (besselk(lambda+1,zeta)/(zeta*besselk(lambda,zeta)))^2);
varx = delta^2 * varx;