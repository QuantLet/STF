function P= HestonVanilla(cp,s,k,v0,vv,rd,rf,tau,kappa,theta,lambda,rho)
%HESTONVANILLA European FX option price in the Heston model.
%   P=HESTONVANILLA(CP,S,K,V0,VV,RD,RF,TAU,KAPPA,THETA,LAMBDA,RHO) returns 
%   the price of a European call (CP=1) or put (CP=-1) option given 
%   spot price S, strike K, initial volatility V0, vol of vol VV, 
%   domestic interest rate RD, foreign interest rate RF, time to maturity 
%   (in years) TAU, level of mean reversion KAPPA, long-run variance THETA, 
%   market price of volatility risk LAMBDA and correlation RHO.
     
%   Sample use:
%       >> HestonVanilla(1,1.03,1,.01,.02,.05,.03,.25,10,.01,0,.5) 
%
%   Reference(s):   
%   [1] H.Albrecher, P.Mayer, W.Schoutens, J.Tistaert (2006) The little
%       Heston trap, Wilmott Magazine, January: 83–92.
%   [2] A.Janek, T.Kluge, R.Weron, U.Wystup (2010) FX smile in the
%       Heston model, see http://ideas.repec.org/p/pra/mprapa/25491.html
%       {Chapter prepared for the 2nd edition of "Statistical Tools for 
%       Finance and Insurance", P.Cizek, W.Härdle, R.Weron (eds.), 
%       Springer-Verlag, forthcoming in 2011.} 

%   Written by Agnieszka Janek and Rafal Weron (2010.07.07)
%   Revised by Rafal Weron (2010.10.08, 2010.12.27)

% Integrate using adaptive Gauss-Kronod quadrature
P1 = 0.5+1/pi.*quadgk(@(phi) HestonVanillaInt(phi,1,s,k,v0,vv,rd,rf,tau,kappa,theta,lambda,rho),0,inf,'RelTol',1e-8);
P2 = 0.5+1/pi.*quadgk(@(phi) HestonVanillaInt(phi,2,s,k,v0,vv,rd,rf,tau,kappa,theta,lambda,rho),0,inf,'RelTol',1e-8);

% Calculate Pplus and Pminus
Pplus = (1-cp)/2 + cp.*P1;
Pminus = (1-cp)/2 + cp.*P2;	

% Calculate option price
P = cp.*( s.*exp(-rf.*tau).*Pplus - k.*exp(-rd.*tau).*Pminus);

%%%%%%%%%%%%% INTERNALLY USED ROUTINE %%%%%%%%%%%%%
  
function F=HestonVanillaInt(phi,m,s,k,v0,vv,rd,rf,tau,kappa,theta,lambda,rho)
%HESTONVANILLAINT Auxiliary function used by HESTONVANILLA.
%   F=HESTONVANILLAINT(phi,m,s,k,v0,vv,rd,rf,tau,kappa,theta,lambda,rho)
%   returns the values of the auxiliary function evaluated at points PHI, 
%   given spot price S, strike K, initial volatility V0, vol of vol VV, 
%   domestic interest rate RD, foreign interest rate RF, time to maturity 
%   (in years) TAU, level of mean reversion KAPPA, long-run variance THETA, 
%   market price of volatility risk LAMBDA and correlation RHO.  

a = kappa.*theta;
u = [0.5 -0.5];
b = [kappa + lambda - rho.*vv kappa + lambda];
x = log(s);
y = log(k);

d = sqrt((1i*rho.*vv.*phi - b(m)).^2 - vv.^2*(2*u(m).*phi.*1i - phi.^2));
g = (b(m) - rho.*vv.*phi.*1i - d)./(b(m) - rho.*vv.*phi.*1i + d);
D = (b(m) - rho.*vv.*phi.*1i - d)./(vv.^2).*((1 - exp(-d.*tau))./(1 - g.*exp(-d.*tau)));
C = (rd - rf).*phi.*1i.*tau + a./(vv.^2).*((b(m) - rho.*vv.*phi.*1i - d).*tau - 2*log((1 - g.*exp(-d.*tau))./(1-g)));
f = exp(C + D.*v0 + 1i.*phi.*x);
F = real(exp(-1i.*phi.*y).*f./(1i.*phi));  