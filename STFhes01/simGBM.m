function x = simGBM(n,x0,mu,sigma,delta,no,method)
%  
%   X = SIMGBM(N,X0,MU,SIGMA,DELTA,NO,METHOD) returns a vector of a sample 
%   trajectory of Geometric Brownian Motion on the time interval [0,N]: 
%     dX(t) = MU*X(t)*dt + SIGMA*X(t)*dW(t), 
%   
%   time step size DELTA, array of normally distributed pseudorandom 
%   numbers NO (array NO is simulated if not provided as an input variable) 
%   and method:
%     METHOD = 0 (default) - direct integration 
%     METHOD = 1 - Euler scheme
%     METHOD = 2 - Milstein scheme 
%     METHOD = 3 - 2nd order Milstein scheme 
%
%   Sample use:
%     >> simGBM(1,.84,.02,sqrt(.1),1/100,normrnd(0,1,100,1),1)
%
%   Reference(s):
%   [1] P.E.Kloeden, E.Platen (1995) Numerical Solution
%       of Stochastic Differential Equations, Springer.
%   [2] R.Korn, E.Korn, G.Kroisandt (2010) Monte Carlo Methods and
%       Models in Finance and Insurance, CRC Press.

%   Written by Rafal Weron (2004.05.21)
%   Revised and added functionality by Rafal Weron (2010.12.27)

% Set default value of 'method'
if (nargin<7)
  method = 0; %direct integration
end
% Generate normal random numbers if not provided in 'no'
if (nargin<6)
  no = normrnd(0,1,ceil(n/delta),1); 
% Check whether length of 'no' is appropriate
else
  no = no(:);
  if (length(no) ~= ceil(n/delta))
    error ('Error: length(no) <> n/delta');
  end
end
    
switch method
  case 1 % Euler scheme
    x = x0*cumprod(1 + mu*delta + sigma.*delta^0.5.*no);
  case 2 % Milstein scheme
    x = x0*cumprod(1 + mu*delta + sigma*delta^0.5.*no ...
      + 0.5*sigma^2*delta*(no.^2-1));
  case 3 % 2nd order Milstein scheme
    x = x0*cumprod(1 + mu*delta + sigma*delta^0.5.*no ...
      + 0.5*sigma^2*delta*(no.^2-1) ...
      + mu*sigma.*no.*(delta^1.5) + 0.5*(mu^2).*(delta^2));
  otherwise % Direct integration
    x = x0*exp(cumsum((mu - 0.5*sigma^2)*delta + sigma.*delta^0.5.*no));
end
% Add starting value
x = [x0; x]; 