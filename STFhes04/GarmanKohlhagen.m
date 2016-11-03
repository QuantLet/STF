function opv = GarmanKohlhagen(S,K,vol,rd,rf,tau,cp,task)
%GARMANKOHLHAGEN European FX option pricing formula.
%   OPV=GARMANKOHLHAGEN(S,K,VOL,RD,RF,TAU,CP,TASK) returns option price, 
%   (spot) delta or strike depending on the value of the TASK parameter:
%     TASK = 0 --> option price (default), 
%     TASK = 1 --> option (spot) delta,
%     TASK = 2 --> option strike (given its delta as parameter K),
%   in the Garman and Kohlhagen (1983) option pricing model.
%   The remaining input parameters are: FX spot S, strike/spot delta K, 
%   volatility VOL, domestic and foreign riskless interest rates RD and RF 
%   (annualized), time to expiry (in years) TAU and option type CP:
%     CP = 1 --> call option (default), CP = -1 --> put option.
%      
%   Sample use:
%     >> C = GarmanKohlhagen(1,1,0.1,0.05,0.03,0.25,1,0)
%  
%   Reference(s):
%   [1] M.B.Garman, S.W.Kohlhagen (1983) Foreign currency option values, 
%       Journal of International Money and Finance 2, 231-237. 
%   [2] A.Janek, T.Kluge, R.Weron, U.Wystup (2010) FX smile in the
%       Heston model, see http://ideas.repec.org/p/pra/mprapa/25491.html
%       {Chapter prepared for the 2nd edition of "Statistical Tools for 
%       Finance and Insurance", P.Cizek, W.Härdle, R.Weron (eds.), 
%       Springer-Verlag, forthcoming in 2011.}  

%   Written by Rafal Weron (2004.06.20)
%   Revised by Agnieszka Janek (2010.07.07)
%   Revised by Rafal Weron (2010.11.16, 2010.12.27)

if (nargin < 6)
  error ('Wrong number of input arguments.')
else
  % Set default values:
  if (nargin < 8)
    task = 0; 	%calculate option price
    if (nargin == 6)
      cp = 1; % call option
    end
  end
end

d1 = ( log(S./K)+(rd-rf+0.5.*vol.^2).*tau )./ ( vol.*sqrt(tau));
d2 = d1 - vol.*sqrt(tau);
switch task
  case 1		% calculate spot delta
    opv = cp.*exp(-rf.*tau).*normcdf(cp.*d1,0,1); 
  case 2  	% calculate strike given spot delta
    opv = S.*exp(-cp.*norminv(K.*exp(rf*tau),0,1).*vol.*sqrt(tau)+(rd-rf+0.5.*vol.^2).*tau);
  otherwise	% calculate option price
    opv = cp.*( S.*exp(-rf.*tau).*normcdf(cp.*d1,0,1)-K.*exp(-rd.*tau).*normcdf(cp.*d2,0,1));
end




