function smile=HestonVanillaSmile(cp,spot,strikes,v0,vv,rd,rf,tau,kappa,theta,lambda,rho)
%HESTONVANILLASMILE The volatility smile implied by the Heston model.
%   SMILE=HESTONVANILLASMILE(CP,SPOT,STRIKES,V0,VV,RD,RF,TAU,KAPPA,THETA,LAMBDA,RHO)
%   returns a vector of volatilities given a vector of strikes STRIKES, 
%   spot price SPOT, initial volatility V0, vol of vol VV, domestic and 
%   foreign interest rates RD and RF, time to maturity (in years) TAU, mean 
%   reversion KAPPA, long-run mean THETA, market price of risk LAMBDA, 
%   correlation RHO and option type CP: 
%     CP=1 -> call option, CP=-1 -> put option.
%     
%   Sample use:
%     >> HestonVanillaSmile(1,0.8,0.2:0.2:1,0,0.2,0.05,0.05,1,2,0.04,0,0.5)
% 
%   Reference(s):
%   [1] A.Janek, T.Kluge, R.Weron, U.Wystup (2010) FX smile in the
%       Heston model, see http://ideas.repec.org/p/pra/mprapa/25491.html
%       {Chapter prepared for the 2nd edition of "Statistical Tools for 
%       Finance and Insurance", P.Cizek, W.Härdle, R.Weron (eds.), 
%       Springer-Verlag, forthcoming in 2011.}  

%   Written by Rafal Weron (2004.07.30)
%   Revised by Agnieszka Janek (2010.07.07, 2010.10.20)
%   Revised by Rafal Weron (2010.12.27)

nostrikes = length(strikes);
% Create an array of option values 
P = zeros(1,nostrikes);
% Create an array of implied volatilities
IV = zeros(1,nostrikes);

for m=1:nostrikes
  lambda = 0;
  % Calculate option prices with respective strikes in the Heston 
  % stochastic volatility model with lambda = 0;
  P(m) = HestonVanilla(cp,spot,strikes(m),v0,vv,rd,rf,tau,kappa,theta,lambda,rho);
  % Calculate implied volatilities for options with respective strikes
  IV(m) = ImplVolFX(P(m),spot,strikes(m),rd,rf,tau,cp);
end
smile = IV;
 
%%%%%%%%%%%%% INTERNALLY USED ROUTINE %%%%%%%%%%%%%

function vol=ImplVolFX(X,S,K,rd,rf,tau,cp)
%IMPLVOLFX Implied volatility assuming the Garman-Kohlhagen model for a
%European style FX option.
%   VOL=IMPLVOLFX(X,S,K,RD,RF,TAU,CP) determines implied volatility 
%   of a European FX option given its price X, spot price SPOT, strike K,
%   domestic interest rate RD, foreign interest rate RF, time to maturity 
%   (in years) TAU and option type CP:
%     CP = 1 --> call option, CP = -1 --> put option

% Find zero of a function GK(vol) - X 
vol = fzero(@(vol) (GarmanKohlhagen(S, K, vol , rd, rf, tau, cp, 0)-X),0.001,optimset('TolX',1e-8));
