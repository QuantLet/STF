function y = mef(param1, param2, param3, x, ind);
% MEF Mean excess function.
%       y = mef (param1, param2, param3, x, ind) returns the values of the
%       mean excess function for the distribution specified by IND (1 -
%       lognormal, 2 - gamma, 3 - Weibull, 4 - Pareto, 5 - Burr, 6 -
%       mixture of exponential distributions with parameters PARAM1,
%       PARAM2, PARAM3. 

switch ind
    case 1
      % Lognormal
      mu    = param1;
      sigma = param2;
      FP = exp(mu+sigma^2/2) * (1-normcdf((log(x) - mu - sigma^2)/sigma));
      FI = normcdf((log(x)-mu)/sigma);
      y=FP./(1-FI)-x;
    case 2
      % Gamma
      alpha = param1;
      beta  = param2;
      FP = alpha / beta * (1 - gamcdf(x, alpha+1, 1/beta));
      FI = gamcdf (x, alpha, 1/beta);
      y=FP./(1-FI)-x;
    case 3
      % Weibull, c.d.f = 1-exp(-t^alpha * beta)
      alpha = param1;
      beta  = param2;
      y = beta^(-1/alpha)* gamma(1+1/alpha).*(1-gammainc(x.^alpha*beta,1+1/alpha)).*exp(beta*x.^alpha)-x;
    case 4
      % Pareto
      alpha   = param1;
      lambda  = param2; 
      y= (x + lambda) ./(alpha-1);      
    case 5
      % Burr
      alpha   = param1;
      lambda  = param2; 
      tau     = param3;
      FP1 = lambda^(1/tau) * gamma(alpha-1/tau) * gamma(1+1/tau) / gamma(alpha).*(lambda./(lambda+x.^tau)).^(-alpha);
      FP2 = 1-betainc(x.^tau ./ (lambda+x.^tau) , 1+1/tau, alpha - 1/tau);
      y = FP1.* FP2-x
    case 6
      % mixture of exponentials
      beta1   = param1;
      beta2   = param2; 
      c=param3;
      FP1 = c./ beta1.*exp(-beta1*x)+(1-c)./beta2.*exp(-beta2*x);
      FP2 = c.*exp(-beta1*x)+(1-c).*exp(-beta2*x);
      y=FP1./FP2; 
end

