function [v0,vv,kappa,theta,rho,IV,SSE]=HestonVanillaFitSmile(delta,marketvols,spot,rd,rf,tau,cp)
%HESTONVANILLAFITSMILE Fit the Heston model to the FX market implied volatility smile.
%   [V0,VV,KAPPA,THETA,RHO,IV,SSE]=HESTONVANILLAFITSMILE(DELTA,MARKETVOLS,SPOT,RD,RF,TAU,PHI) 
%   returns initial volatility V0, vol of vol VV, mean reversion KAPPA, 
%   long-run mean THETA, correlation RHO, vector of Garman-Kohlhagen  
%   implied volatilities IV, sum of squared errors SSE given a vector of  
%   spot delta values DELTA, vector of the market implied volatilities 
%   MARKETVOLS, spot price SPOT, domestic and foreign interest rates RD
%   and RF, time to maturity (in years) TAU and option type CP:
%     CP = 1 --> call option, CP = -1 --> put option
%        
%   Sample use:
%     >> [v0,vv,kappa,theta,rho,IV,SSE]=HestonVanillaFitSmile(0.2:0.2:0.6,0.1:0.1:0.3,1,0.05,0.05,1,1)
% 
%   Reference(s):
%   [1] A.Janek, T.Kluge, R.Weron, U.Wystup (2010) FX smile in the
%       Heston model, see http://ideas.repec.org/p/pra/mprapa/25491.html
%       {Chapter prepared for the 2nd edition of "Statistical Tools for 
%       Finance and Insurance", P.Cizek, W.Härdle, R.Weron (eds.), 
%       Springer-Verlag, forthcoming in 2011.} 

%   Written by Rafal Weron (2004.07.30) 
%   Revised by Agnieszka Janek (2010.07.07, 2010.10.20)
%   Written by Rafal Weron (2010.12.27) 

% Calculate strikes for market deltas
strikes = GarmanKohlhagen(spot,delta,marketvols,rd,rf,tau,cp,2);
nostrikes = length(strikes);
% Set the initial variance v0 = implied ATM market volatility
v0 = (marketvols(floor(nostrikes/2)+1)).^2;
% Set mean reversion to constant value = 1.5
kappa = 1.5; 		
% Set initial values for vv, theta, rho
initparam = [2*sqrt(v0) 2*v0 0];
% Set default maximum number of iterations, termination
% tolerance on the function value and on params
opt = optimset('MaxIter',1e3,'TolFun', 1e-8, 'TolX', 1e-8);
 
% Find a local minimizer 'param' of the function HestonVanillaSSE with initial parameter
% 'initparam' (Nelder-Mead's method) 
nn = fminsearch(@(param) HestonVanillaSSE(param,cp,spot,strikes,v0,marketvols,rd,rf,tau,kappa),initparam,opt);
  
vv = nn(1);
theta = nn(2);
rho = nn(3);
% Create an array of option values 
P = zeros(1,nostrikes);
% Create an array of implied volatilities
IV = zeros(1,nostrikes);
for i=1:nostrikes
  % Calculate option prices with respective strikes in the Heston
  % stochastic volatility model with lambda = 0;
  P(i) = HestonVanilla(cp,spot,strikes(i),v0,vv,rd,rf,tau,kappa,theta,0,rho); 
  % Calculate implied volatilities for options with respective strikes
  IV(i) = ImplVolFX(P(i),spot,strikes(i),rd,rf,tau,cp);
end

% Calculate SSE of the fit for options with respective strikes
SSE = HestonVanillaSSE([vv theta rho],cp,spot,strikes,v0,marketvols,rd,rf,tau,kappa);

%%%%%%%%%%%%% INTERNALLY USED ROUTINES %%%%%%%%%%%%%
  
function E=HestonVanillaSSE(param,cp,spot,strikes,v0,marketvols,rd,rf,tau,kappa)
%HESTONVANILLASSE Computes the sum of squared errors (SSE) between the
%market and model implied volatilities.
%   E = HESTONVANILLASSE(param,cp,spot,strikes,v0,marketvols,rd,rf,tau,kappa)
%   returns the sum of squared errors of the fit given 3-element vector PARAM:
%     param(1) = vol of vol VV, 
%     param(2) = long-run variance THETA,
%     param(3) = correlation RHO,
%   spot price SPOT, a vector of market implied volatilities MARKETVOLS,
%   a vector of strikes STRIKES, initial volatility V0, domestic interest 
%   rate (annualized) RD, foreign interest rate (annualized) RF, time to 
%   maturity (in years) TAU, level of mean reversion KAPPA and option type CP:
%     CP = 1 --> call option, CP = -1 --> put option

% Set vv, theta, rho
vv = param(1);
theta = param(2);
rho = param(3);

nostrikes = length(strikes);
% Check whether vv, theta and rho are admissible values
if ( vv < 0 || theta < 0|| abs(rho)>1)
  E = Inf;
else
  % Create an array of option values 
  P = zeros(1,nostrikes);
  % Create an array of implied volatilities
  IV = zeros(1,nostrikes);
  for i=1:nostrikes
    % Calculate option prices with respective strikes in Heston's 
    % stochastic volatility model with lambda =0;
    P(i) = HestonVanilla(cp,spot,strikes(i),v0,vv,rd,rf,tau,kappa,theta,0,rho);
    % Calculate implied volatilities for options with respective strikes
    IV(i) = ImplVolFX(P(i),spot,strikes(i),rd,rf,tau,cp);
  end
  % Calculate SSE of the fit 
  E = sum((marketvols - IV).^2 );
end

function vol=ImplVolFX(X,S,K,rd,rf,tau,cp)
%IMPLVOLFX Implied volatility assuming the Garman-Kohlhagen model for 
%a European style FX option.
%   VOL = IMPLVOLFX(X, S, K, rd, rf, tau, cp) determines implied volatility 
%   of a European FX option given its price X, spot price SPOT, strike K,
%   domestic interest rate RD, foreign interest rate RF, time to maturity 
%   (in years) TAU and option type CP:
%     CP = 1 --> call option, CP = -1 --> put option

% Find zero of a function GK(vol) - X 
vol = fzero(@(vol) (GarmanKohlhagen(S, K, vol , rd, rf, tau, cp, 0)-X),0.001,optimset('TolX',1e-8));
