function y = pdfHeston(x,theta,kappa,sigma,rho,t,par)
%PDFHESTON Probability density function in the Heston stochastic volatility model.
%   Y=PDFHESTON(X,THETA,KAPPA,SIGMA,RHO,T) returns the Heston pdf 
%   values at points X given long-run variance THETA, level of mean reversion
%   KAPPA, volatility of variance (vol of vol) SIGMA, correlation RHO 
%   and time lag T. 
%   Y=PDFHESTON(X,THETA,KAPPA,SIGMA,RHO,T,1) additionally plots the pdf.
%
%   Sample use:
%     >> pdfHeston((-2:2)/2,.04,2,.3,-.05,1)
%
%   Reference(s):   
% 	[1] A.A.Dragulescu, V.M.Yakovenko (2002) Probability distribution 
%       of returns in the Heston model with stochastic volatility, 
%       Quantitative Finance 2: 443-453.
%   [2] A.Janek, T.Kluge, R.Weron, U.Wystup (2010) FX smile in the
%       Heston model, see http://ideas.repec.org/p/pra/mprapa/25491.html
%       {Chapter prepared for the 2nd edition of "Statistical Tools for 
%       Finance and Insurance", P.Cizek, W.Härdle, R.Weron (eds.), 
%       Springer-Verlag, forthcoming in 2011.} 

% Written by Rafal Weron (2004.07.30)
% Revised by Agnieszka Janek (2010.07.07, 2010.10.20)
% Revised by Rafal Weron (2010.11.16, 2010.12.27)
 
if (nargin < 7)
  par = 0;
end
  
n = length(x); 
% Create an array of pdf values
y = zeros(1,n); 
% Set the limits of integration
pxm=100;
% Compute the marginal pdf
for i=1:n
  % quadgk - evaluate the integral numerically using 
  % the adaptive Gauss-Kronrod quadrature on the interval [-100,100]
  y(i) = 1/(2*pi)*quadgk(@(s) pdfHestonInt(x(i),theta,kappa,sigma,rho,t,s),-100,100,'RelTol',1e-8);
end

% Plot marginal pdf if par==1
if (par == 1)
  plot(x,y)
end
  
%%%%%%%%%%%%% INTERNALLY USED ROUTINE %%%%%%%%%%%%%
  
function y=pdfHestonInt(x,theta,kappa,sigma,rho,t,px)
%PDFHESTONINT Auxiliary function used by PDFHESTON.
%   Y=PDFHESTONINT(X,THETA,KAPPA,SIGMA,RHO,T,PX) returns the values of the 
%   auxiliary function evaluated at points X, integrated  numerically  
%   at points PX, given long-run variance THETA,level of mean reversion 
%   KAPPA, volatility of variance (vol of vol) SIGMA, correlation RHO
%   and time lag T.

Gamma = kappa + 1i*rho*sigma*px;
Gamma2 = Gamma.*Gamma;
Omega = sqrt(Gamma2+sigma^2*(px.^2-1i*px));
Omega2 = Omega.*Omega;
Ft1 = kappa*theta/(sigma^2)*Gamma*t;
  
ncoeff = Omega2 - Gamma2+2*kappa*Gamma;
coeff = ncoeff./(2*kappa*Omega);
Ot = Omega*t/2;
Err = cosh(Ot);

Ft2 = log(Err + coeff.*sinh(Ot));
Ft = Ft1 - 2*kappa*theta/(sigma^2)*Ft2;
Eksp = exp(Ft + 1i*x.*px);
y = real(Eksp);
