function y=Burrpdf(x,alpha,lambda,tau)
%BURRPDF Burr probability density function (pdf).
%   Y = BURRPDF(X,ALPHA,LAMBDA,TAU) returns the pdf of the Burr
%   distribution with parameters ALPHA,LAMBDA,TAU, evaluated at the values in X.
%   For CONTROL=0 the error message is displayed, if the parmeters are
%   negative.
%   
%   The default values for the parameters ALPHA, LAMBDA, TAU, CONTROL are
%   1, 1, 2, 0, respectively.

  if nargin<4
     tau=2;
  end 
  if nargin<3
    lambda=1;
  end
  if nargin<2
    alpha=1;
  end
  if nargin<1
      error('stats:normpdf:TooFewInputs','Input argument X is undefined.');
  end
  if tau<=0
      error('Non-positive tau!');
  end
  if lambda<=0
      error('Non-positive sigma!');
  end
  if alpha<=0
      error('Non-positive alpha!');
  end

  y=zeros(size(x));
  pos=find(x>0);
  y(pos)=tau*alpha*lambda^alpha*x(pos).^(tau-1).*(lambda+x(pos).^tau).^(-alpha-1);

