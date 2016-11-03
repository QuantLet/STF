function x = simHeston(n,s0,v0,mu,kappa,theta,sigma,rho,delta,no,method)
% SIMHESTON sample trajectories of the spot price and volatility processes 
%   in the Heston model.
%   X = SIMHESTON(n,s0,v0,mu,kappa,theta,sigma,rho,delta,no,method) returns 
%   a 2-column array, containing the simulated trajectories of the spot 
%   price S(t) and volatility v(t) for t=0:DELTA:N, in the model:
% 		dS(t) = mu*S(t)*dt + v^0.5*S(t)*dW1(t) 
% 		dv(t) = kappa*(theta - v(t))*dt + sigma*(v(t)^0.5)*dW2(t)
% 		Cov[dW1(t),dW2(t)] = rho*dt
%   
%   Reference(s):
%   [1] L.Andersen (2007) Efficient Simluation of the Heston Stochastic 
%       Volatility Model, Working Paper.
%   [2] A.Janek, T.Kluge, R.Weron, U.Wystup, "FX Smile in the
%       Heston Model", in Cizek et al. (2011) "Statistical Tools for
%       Finance and Insurance", Springer.


% Set default value of method (Quadratic-Exponential scheme)
if (nargin < 11)
  method = 0; 	
end

% Set default no  
if (nargin < 10)
  no = normrnd(0,1,n/delta,2); 
end
    
% Check whether length(no) is appropriate
if (length(no(:,1))~=(n/delta))
  error ('Size of no is inappropriate. Length of no(:,1) should be equal to n/delta.')
end

% Construct the time vector
t = 0:ceil(n/delta);
sizet = length(t);
x = zeros(sizet,2);
% Construct the correlation matrix
C = reshape([1 rho rho 1],[2 2]);
% Compute correlated pseudo random normal numbers for the spot price (Euler
% scheme) and volatility processes
u = normalcorr(C,no) * sqrt(delta);
if (method == 0) 
  % Use Quadratic-Exponential scheme  
  % Initialize x with starting points log(s0) and v0
  x(1,:) = [log(s0) v0];
  % Initialize phiC (phiC in [1,2])
  phiC = 1.5;
  % Compute x values at each time point
  for i=2:sizet
    m = theta + (x(i-1,2)-theta)*exp(-kappa*delta);
    s2 = x(i-1,2)*sigma^2*exp(-kappa*delta)/kappa*(1-exp(-kappa*delta))+theta*sigma^2/(2*kappa)*(1-exp(-kappa*delta))^2;
    phi = s2/m^2;
    gamma1 = 0.5;
    gamma2 = 0.5;
    K0 = -rho*kappa*theta/sigma*delta;
    K1 = gamma1*delta*(kappa*rho/sigma-0.5)-rho/sigma;
    K2 = gamma2*delta*(kappa*rho/sigma-0.5)+rho/sigma;
    K3 = gamma1*delta*(1-rho^2);
    K4 = gamma2*delta*(1-rho^2);
    if (phi <= phiC)
      b2 = 2*phi^(-1)-1+sqrt(2*phi^(-1))*sqrt(2*phi^(-1)-1);
      a = m/(1+b2);
      % Calculate values of the volatility process 
      x(i,2) = a.*(sqrt(b2)+no(i-1,2)).^2;
    else
      p = (phi-1)/(phi+1);
      beta = (1-p)/m;
      if (0<=cdf('norm',no(i-1,2),0,1) && cdf('norm',no(i-1,2),0,1)<=p)
        x(i,2) = 0;
      else
        if (p<cdf('norm',no(i-1,2),0,1) && cdf('norm',no(i-1,2),0,1)<=1)
          % Calculate values of the volatility process 
          x(i,2) = beta^(-1)*log((1-p)/(1-cdf('norm',no(i-1,2),0,1)));
        end
      end
    end
    % Calculate values of the log(spot) price process
    x(i,1) = x(i-1,1) + mu*delta + K0 + K1*x(i-1,2) + K2*x(i,2) + ...
      sqrt(K3*x(i-1,2) + K4*x(i,2))*no(i-1,1);
  end
  % Calculate values of the spot price process
  x(:,1) = exp(x(:,1));    
else
  % Euler scheme with absorption or reflection for the volatility process
  % Initialize x with starting points s0 (spot price) and v0 (volatility process)
  x(1,:) = [s0 v0];
  for i=2:sizet
    if (method==1)
      % Absorption variant
      x(i-1,2)= max([x(i-1,2) 0]);
    else
      % Reflection variant
      x(i-1,2)= max([x(i-1,2) -x(i-1,2)]);
    end
    x(i,2) = x(i-1,2) + kappa*(theta - x(i-1,2))*delta + sigma*sqrt(x(i-1,2))*u(i-1,2);
    % Use standard Euler scheme for the spot price process
    x(i,1) = x(i-1,1) + x(i-1,1)*(mu*delta + sqrt(x(i-1,2))*u(i-1,1));
  end
end

%%%%%%%%%%%%% INTERNALLY USED ROUTINE %%%%%%%%%%%%%

function x = normalcorr(C,no)
% NORMALCORR Generate correlated pseudo random normal variates 
% using the Cholesky factorization.
%   x = NORMALCORR(n,C,no) returns (size(no,1), sqrt(size(C)))-dimensional 
%   array x containing the simulated trajectories that approximately has
%   a correlation structure defined by matrix C. 
%   NO - matrix of normally distributed pseudorandom numbers.
%     
%   Sample use:
%       >> C = reshape([1 .5 .5 1],[2 2]);
%       >> normalcorr(C,normrnd(0,1,200,2))

M = chol(C)';
x = (M*no')';

